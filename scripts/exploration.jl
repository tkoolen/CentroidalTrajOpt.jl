# To do:
#
# * Test that CoPs are in region
# * Test that friction cone constraints are satisfied
# * Move utilities to appropriate places
# * better CoM kinematic constraints
# * Initial condition modification

## Packages
using CentroidalTrajOpt
using LinearAlgebra

using MultilinearOpt
using Polyhedra: hrep, polyhedron
using CentroidalTrajOpt: extrude

using RigidBodyDynamics
using RigidBodyDynamics.PDControl
using RigidBodyDynamics.Contact
using MeshCatMechanisms
using GeometryTypes: Point, GLNormalMesh
using QPControl
using QPWalkingControl
using AtlasRobot
using RigidBodySim
using Rotations
using OSQP
using MathOptInterface

using QPControl.Trajectories: PointTrajectory

using MeshCat: RGBA

const MOI = MathOptInterface

include("optimizers.jl")

function create_atlas()
    urdf = AtlasRobot.urdfpath()
    mechanism = parse_urdf(urdf, floating=true, remove_fixed_tree_joints=false)
    link_colors = Dict(map(body -> string(body) => RGBA(0.7f0, 0.7f0, 0.7f0, 0.5f0), bodies(mechanism)))
    visuals = URDFVisuals(AtlasRobot.urdfpath(); package_path=[AtlasRobot.packagepath()], link_colors=link_colors)
    remove_fixed_tree_joints!(mechanism)
    foot_points = AtlasRobot.foot_contact_points(mechanism)
    sole_frames = AtlasRobot.add_sole_frames!(mechanism)
    nominal_state = MechanismState(mechanism)
    AtlasRobot.setnominal!(nominal_state)
    floating_joint = first(joints(mechanism))
    configuration(nominal_state, floating_joint)[end] += -0.0028061189941; # FIXME
    pelvis = findbody(mechanism, "pelvis");
    mechanism, nominal_state, foot_points, sole_frames, floating_joint, pelvis, visuals
end

function center_of_mass_velocity(state::MechanismState)
    h = momentum(state)
    FreeVector3D(h.frame, linear(h)) / mass(state.mechanism)
end

function create_environment()
    region_data = ContactRegion{Float64}[]
    push!(region_data, ContactRegion(
            AffineMap(one(RotMatrix{3}), zero(SVector{3})),
            0.7,
            0.0,
            Float64[1 0; 0 1; -1 0; 0 -1],
            0.2 * ones(4)
    ))
    push!(region_data, ContactRegion(
            AffineMap(one(RotMatrix{3}) * RotXYZ(0.1, -0.2, 0.3), SVector(1.0, 0.3, 0.2)),
            0.7,
            0.0,
            Float64[1 0; 0 1; -1 0; 0 -1],
            0.4 * ones(4)
    ))
    push!(region_data, ContactRegion(
            AffineMap(one(RotMatrix{3}) * RotXYZ(0.1, -0.2, 0.3), SVector(0.0, 1.0, 0.2)),
            0.7,
            0.0,
            Float64[1 0; 0 1; -1 0; 0 -1],
            0.3 * ones(4)
    ))
    region_data
end

function create_contact_model(
        mechanism::Mechanism,
        foot_points::AbstractDict{BodyID, <:AbstractVector{<:Point3D}},
        region_data::Vector{<:ContactRegion})
    contact_model = ContactModel()
    normal_model = hunt_crossley_hertz(; k=500e3)
    k_tangential = 20e3
    b_tangential = 100.#2 * sqrt(k_tangential * mass(mechanism) / 10)
    world_frame = root_frame(mechanism)
    foot_collision_elements = CollisionElement[]
    for (bodyid, points) in foot_points
        body = findbody(mechanism, bodyid)
        for point in points
            push!(foot_collision_elements, CollisionElement(body, point.frame, Point(point.v)))
        end
    end
    push!(contact_model, foot_collision_elements)
    region_thickness = 0.2
    for (i, region) in enumerate(region_data)
        geometry = GLNormalMesh(polyhedron(extrude(hrep(region.A, region.b), region_thickness)))
        frame = CartesianFrame3D("region_$i")
        add_frame!(root_body(mechanism), Transform3D(frame, world_frame, region.transform))
        element = CollisionElement(root_body(mechanism), frame, geometry)
        group = CollisionElement[element]
        push!(contact_model, group)
        tangential_model = ViscoelasticCoulombModel(region.μ, k_tangential, b_tangential)
        contact_force_model = SplitContactForceModel(normal_model, tangential_model)
        set_contact_force_model!(contact_model, foot_collision_elements, group, contact_force_model)
    end
    return contact_model
end

function create_controller(
        mechanism::Mechanism,
        floating_joint::Joint,
        foot_points::AbstractDict{BodyID, <:AbstractVector{<:Point3D}},
        μ::Number,
        pelvis::RigidBody,
        nominal_state::MechanismState,
        com_trajectory,
    )
    # Low level controller
    optimizer = OSQP.Optimizer(verbose=false, eps_abs=1e-5, eps_rel=1e-5, max_iter=5000, adaptive_rho_interval=25)
    lowlevel = MomentumBasedController{4}(mechanism, optimizer, floatingjoint = floating_joint);
    for (bodyid, points) in foot_points
        body = findbody(mechanism, bodyid)
        for point in points
            normal = FreeVector3D(default_frame(body), 0.0, 0.0, 1.0)
            contact = addcontact!(lowlevel, body, point, normal, μ)
            contact.maxnormalforce[] = 1e6 # TODO
            contact.weight[] = 1e-3
        end
    end

    # Linear momentum controller
    world_frame = root_frame(mechanism)
    com_gains = QPWalkingControl.critically_damped_gains(10.0)
    linear_momentum_controller = PDCoMController(com_gains, PointTrajectory(world_frame, com_trajectory), mass(mechanism))

    # State machine
    contacts = Dict(BodyID(body) => contact for (body, contact) in lowlevel.contacts)
    state_machine = CoMTrackingStateMachine(mechanism, contacts)

    # High level controller
    HumanoidQPController(lowlevel, pelvis, nominal_state,
        state_machine, collect(values(state_machine.end_effector_controllers)), linear_momentum_controller)
end

function MeshCatMechanisms.setelement!(mvis::MechanismVisualizer, contact_model::ContactModel, args...; base_name="collision element ")
    mechanism = mvis.state.mechanism
    i = 1
    for group in contact_model.collision_groups
        for element in group
            geometry = element.geometry
            if geometry isa MeshCatMechanisms.GeometryLike
                frame = element.transform.from
                setelement!(mvis, frame, element.geometry, args...)
                i += 1
            end
        end
    end
end

## Robot setup
mechanism, state0, foot_points, sole_frames, floating_joint, pelvis, visuals = create_atlas()

## Environment
region_data = create_environment()

## Collision setup
contact_model = create_contact_model(mechanism, foot_points, region_data)

## Initial conditions
c0 = center_of_mass(state0).v
ċ0 = center_of_mass_velocity(state0).v
contacts0 = map(values(sole_frames)) do sole_frame # TODO
    p0 = translation(transform_to_root(state0, sole_frame))
    @assert norm(c0 - p0) <= 1.2 # TODO: integrate with kinematic limits
    region_data[1] => p0
end
for i in eachindex(contacts0)
    p_i_0 = last(contacts0[i])
    for j in eachindex(contacts0)
        p_j_0 = last(contacts0[j])
        if i != j
            @assert norm(p_i_0 - p_j_0) >= 0.1 # TODO: integrate with kinematic limits
        end
    end
end

## Additional settings
g = mechanism.gravitational_acceleration.v
max_cop_distance = 0.07

## Optimizer
# optimizer_factory = baron_optimizer_factory()
optimizer_factory = scip_optimizer_factory()

## Problem
problem = CentroidalTrajectoryProblem(optimizer_factory, region_data, c0, ċ0, contacts0;
    g=g, max_cop_distance=max_cop_distance, num_pieces=2, c_degree=3,
    # objective_type=ObjectiveTypes.MIN_EXCURSION);
    objective_type=ObjectiveTypes.FEASIBILITY);

disallow_jumping!(problem)

# fix.(problem.z_vars[:, :, 1], [1.0 1.0; 0.0 0.0; 0.0 0.0])
# fix.(problem.z_vars[:, :, 2], [0.0 0.0; -0.0 -0.0; 1.0 1.0])

if optimizer_factory.constructor == SCIP.Optimizer
    problem.model.optimize_hook = function (model)
        mscip = backend(model).optimizer.model.mscip
        # SCIP.SCIPsetEmphasis(mscip, SCIP.SCIP_PARAMEMPHASIS_FEASIBILITY, true)
        # SCIP.SCIPsetPresolving(mscip, SCIP.SCIP_PARAMSETTING_AGGRESSIVE, true)
        # SCIP.SCIPsetHeuristics(mscip, SCIP.SCIP_PARAMSETTING_AGGRESSIVE, true)
        MOI.optimize!(backend(model))
        return
    end
end

relax = optimizer_factory.constructor == Gurobi.Optimizer || optimizer_factory.constructor == CPLEX.Optimizer
if relax
    @info "Relaxing bilinearities."
    relaxbilinear!(problem.model, method=:Logarithmic1D, disc_level=17, constraints_to_skip=problem.friction_cone_quadratic_constraints)
end

## Create visualizer
using MeshCat
newvis = false
if newvis || (!@isdefined vis) || isempty(vis.core.scope.pool.connections)
    vis = Visualizer()
    open(vis)
    wait(vis)
end
delete!(vis)

## Centroidal trajectory visualization
cvis = CentroidalTrajectoryVisualizer(vis, region_data, norm(g), length(contacts0))

## Robot visualization
mvis = MechanismVisualizer(mechanism, visuals, vis)
copyto!(mvis, state0)

## Environment visualization
setelement!(mvis, contact_model)

## Feasibility
result = CentroidalTrajOpt.solve!(problem);

## Result visualization
set_com_trajectory!(cvis, result)
set_state!(cvis, result, 0.0)

if backend(problem.model).optimizer.model isa SCIP.Optimizer
    mscip = backend(problem.model).optimizer.model.mscip
    # SCIP.print_statistics(mscip)
    # SCIP.print_heuristic_statistics(mscip)
end

## Optimality

# Setting objective function straight away somehow prevents SCIP from finding a feasible point, so:
# SCIP.SCIPsetEmphasis(problem.model.moi_backend.optimizer.model.mscip, SCIP.SCIP_PARAMEMPHASIS_OPTIMALITY, true);
# set_objective(problem.model, MOI.MAX_SENSE, sum(problem.z_vars))
# result = solve!(problem);
# set_objective(problem.model, MOI.MIN_SENSE, sum(problem.Δts))
# result = solve!(problem);

if optimizer_factory.constructor == BARON.Optimizer
    @info "BARON problem file: $(backend(problem.model).optimizer.model.optimizer.inner.problem_file_name)"
    @info "BARON result file: $(backend(problem.model).optimizer.model.optimizer.inner.result_file_name)"
end

## Tests
using Test
c = result.center_of_mass
ċ = map_subfunctions(derivative, c)
c̈ = map_subfunctions(derivative, ċ)
rs = result.centers_of_pressure
fs = result.contact_forces
# τns = result.contact_normal_torques
ps = result.contact_positions
ns = normals(problem)

# Times
@test issorted(result.break_times)
@test all(x -> x >= 0, result.break_times)

T = last(result.break_times)

# Initial and final conditions
@test c(0) ≈ c0 atol=1e-12
@test ċ(0) ≈ ċ0 atol=1e-12
@test ċ(T) ≈ zeros(3) atol=1e-12
@test c̈(T) ≈ zeros(3) atol=1e-12

for t in range(0, T, length=100)
    ftot = sum(fs[contact](t) for contact in problem.contacts.val)

    # Dynamics
    @test c̈(t) ≈ g + ftot atol=1e-5

    # Torque about CoM
    n = ns[problem.regions(1)]
    # τ = sum((rs[j](t) - c(t)) × fs[j](t) + n * τns[j](t) for j in problem.contacts.val)
    τ = sum((rs[j](t) - c(t)) × fs[j](t) for j in problem.contacts.val)

    if !relax
        @test τ ≈ zeros(3) atol=1e-4
    end

    # Friction cones
    for j in problem.contacts.val
        zs = result.contact_indicators[j](t)
        @test sum(zs) <= 1 + 1e-6
        region_idx = findfirst(isequal(1), zs)

        f = fs[j]
        # τn = τns[j]
        p = ps[j]
        r = rs[j]

        f_t = f(t)
        fn_t = f_t ⋅ n
        # τn_t = τn(t)

        @test fn_t >= -1e-7
        @test norm(p(t) - r(t)) < max_cop_distance + 1e-5

        if region_idx == nothing
            @test f_t ≈ zero(f_t) atol=1e-4
            # @test τn_t ≈ zero(τn_t) atol=1e-4
        else
            region = region_data[region_idx]
            μ = region.μ
            @test norm(f_t - n * fn_t) <= μ * fn_t + 1e-6
            # @test -μrot * fn_t <= τn_t
            # @test τn_t <= μrot * fn_s
        end
    end
end

# Continuity around breaks
for t in result.break_times
    @test c(t - 1e-8) ≈ c(t + 1e-8) atol=1e-5
    @test ċ(t - 1e-8) ≈ ċ(t + 1e-8) atol=1e-5
end

## Plan Visualization
plan_animation = Animation(cvis, result)

## Mode sequence
# value.(problem.z_vars)

## Controller
μ_control = 0.7
controller = create_controller(mechanism, floating_joint, foot_points, μ_control, pelvis, state0, result.center_of_mass);

## Simulation
simulate = true
if simulate
    gui = GUI(mvis);
    open(gui)
    state = MechanismState(mechanism)
    copyto!(state, state0)
    Δt = 1 / 500
    pcontroller = PeriodicController(similar(velocity(state)), Δt, controller)
    damping = JointDamping{Float64}(mechanism, AtlasRobot.urdfpath())
    dynamics = Dynamics(
        mechanism,
        SumController(similar(velocity(state)), (pcontroller, damping));
        contact_model=contact_model)
    callback = CallbackSet(RealtimeRateLimiter(poll_interval=pi / 100), CallbackSet(gui; max_fps=60))
    T = last(result.break_times)
    tspan = (0., T)
    contact_state = SoftContactState(contact_model)
    odeproblem = ODEProblem(dynamics, (state, contact_state), tspan; callback=callback)

    # simulate
    @time sol = RigidBodySim.solve(odeproblem, Tsit5(), abs_tol = 1e-8, dt = 1e-6, dtmax=1e-3);
    sim_animation = Animation(mvis, sol)
end

setanimation!(vis, plan_animation);
setanimation!(vis, sim_animation);
setanimation!(vis, merge(plan_animation, sim_animation));
