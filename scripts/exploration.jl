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
using QPWalkingControl: PosePlan, set_pose_plan!

import Polyhedra
import Blink

using MeshCat: RGBA

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

include("optimizers.jl")
include("util.jl")

function create_atlas()
    urdf = AtlasRobot.urdfpath()
    mechanism = parse_urdf(urdf, floating=true, remove_fixed_tree_joints=false)
    link_colors = Dict(map(body -> string(body) => RGBA(0.7f0, 0.7f0, 0.7f0, 0.5f0), bodies(mechanism)))
    visuals = URDFVisuals(AtlasRobot.urdfpath(); package_path=[AtlasRobot.packagepath()], link_colors=link_colors)
    remove_fixed_tree_joints!(mechanism)
    foot_points = AtlasRobot.foot_contact_points(mechanism)
    sole_frames = AtlasRobot.add_sole_frames!(mechanism)
    nominal_state = MechanismState(mechanism)
    AtlasRobot.setnominal!(nominal_state; kneebend=1.6)
    floating_joint = first(joints(mechanism))
    Δz = translation(transform_to_root(nominal_state, first(values(sole_frames))))[3]
    configuration(nominal_state, floating_joint)[end] -= Δz
    setdirty!(nominal_state)
    pelvis = findbody(mechanism, "pelvis");
    mechanism, nominal_state, foot_points, sole_frames, floating_joint, pelvis, visuals
end

function create_environment()
    region_data = ContactRegion{Float64}[]
    push!(region_data, ContactRegion(
            AffineMap(one(RotMatrix{3}), zero(SVector{3})),
            0.7,
            0.0,
            Float64[1 0; 0 1; -1 0; 0 -1],
            0.15 * ones(4)
    ))
    push!(region_data, ContactRegion(
            AffineMap(one(RotMatrix{3}) * RotXYZ(0.1, -0.2, 0.3), SVector(0.7, 0.3, 0.3)),
            0.7,
            0.0,
            Float64[1 0; 0 1; -1 0; 0 -1],
            0.15 * ones(4)
    ))
    push!(region_data, ContactRegion(
            AffineMap(one(RotMatrix{3}) * RotXYZ(-0.1, 0.2, 0.3), SVector(0.5, 1.0, 0.2)),
            0.7,
            0.0,
            Float64[1 0; 0 1; -1 0; 0 -1],
            0.15 * ones(4)
    ))
    push!(region_data, ContactRegion(
        AffineMap(one(RotMatrix{3}) * RotXYZ(0.0, 0.0, 0.0), SVector(1.1, 1.2, 0.1)),
        0.7,
        0.0,
        Float64[1 0; 0 1; -1 0; 0 -1],
        0.1 * ones(4)
    ))
    region_data
end

function create_contact_model(
        mechanism::Mechanism,
        foot_points::AbstractDict{BodyID, <:AbstractVector{<:Point3D}},
        region_data::Vector{<:ContactRegion}; region_offset) # TODO
    contact_model = ContactModel()
    normal_model = hunt_crossley_hertz(; k=500e3)
    k_tangential = 7e3
    b_tangential = 250.#2 * sqrt(k_tangential * mass(mechanism) / 10)
    world_frame = root_frame(mechanism)
    foot_collision_elements = CollisionElement[]
    for (bodyid, points) in foot_points
        body = findbody(mechanism, bodyid)
        for point in points
            push!(foot_collision_elements, CollisionElement(body, point.frame, Point(point.v)))
        end
    end
    push!(contact_model, foot_collision_elements)
    debug_with_flat_ground = false
    if debug_with_flat_ground
        geometry = HalfSpace(SVector(0., 0., 1.), 0.0)
        frame = world_frame
        group = CollisionElement[CollisionElement(root_body(mechanism), frame, geometry)]
        push!(contact_model, group)
        tangential_model = ViscoelasticCoulombModel(0.8, k_tangential, b_tangential)
        contact_force_model = SplitContactForceModel(normal_model, tangential_model)
        set_contact_force_model!(contact_model, foot_collision_elements, group, contact_force_model)
    else
        region_thickness = 0.3
        for (i, region) in enumerate(region_data)
            b̄ = region.b .+ flipsign.(region_offset, region.b) # FIXME
            geometry = Contact.HRep(extrude(hrep(region.A, b̄), region_thickness))
            frame = CartesianFrame3D("region_$i")
            add_frame!(root_body(mechanism), Transform3D(frame, world_frame, region.transform))
            element = CollisionElement(root_body(mechanism), frame, geometry)
            group = CollisionElement[element]
            push!(contact_model, group)
            tangential_model = ViscoelasticCoulombModel(region.μ, k_tangential, b_tangential)
            contact_force_model = SplitContactForceModel(normal_model, tangential_model)
            set_contact_force_model!(contact_model, foot_collision_elements, group, contact_force_model)
        end
    end
    return contact_model
end

function create_controller(
        mechanism::Mechanism,
        contact_body_ids::AbstractVector{BodyID},
        floating_joint::Joint,
        foot_points::AbstractDict{BodyID, <:AbstractVector{<:Point3D}},
        μ::Number,
        pelvis::RigidBody,
        state0::MechanismState,
        nominal_state::MechanismState,
        com_trajectory,
        pose_plans::AbstractDict{BodyID, <:PosePlan}
    )
    # Low level controller
    optimizer = OSQP.Optimizer(verbose=false, eps_abs=1e-5, eps_rel=1e-5, max_iter=5000, adaptive_rho_interval=25)
    lowlevel = MomentumBasedController{4}(mechanism, optimizer, floatingjoint = floating_joint);
    for bodyid in contact_body_ids
        points = foot_points[bodyid]
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
    for (bodyid, plan) in pose_plans
        set_pose_plan!(state_machine, bodyid, plan)
    end

    # High level controller
    HumanoidQPController(lowlevel, pelvis, nominal_state,
        state_machine, collect(values(state_machine.end_effector_controllers)), linear_momentum_controller)
end

## Robot setup
mechanism, state0, foot_points, sole_frames, floating_joint, pelvis, visuals = create_atlas()
contact_body_ids = sort(collect(keys(sole_frames)), by=x -> x.value) # establishes order once and for all

## Environment
region_data = create_environment()

## Parameters
g = mechanism.gravitational_acceleration.v
inside_foot_radius = Inf
outside_foot_radius = 0.0
for (bodyid, points) in foot_points
    body = findbody(mechanism, bodyid)
    sole_frame = sole_frames[bodyid]
    for point in points
        point_sole_frame = fixed_transform(body, point.frame, sole_frame) * point
        dist = norm(point_sole_frame.v)
        global inside_foot_radius = min(inside_foot_radius, dist) # TODO: not technically correct
        global outside_foot_radius = max(outside_foot_radius, dist)
    end
end
max_cop_distance = inside_foot_radius
max_com_to_contact_distance = 1.07
min_inter_contact_distance = 0.2# TODO: 2 * outside_foot_radius
region_offset = outside_foot_radius + 0.04

## Collision setup
contact_model = create_contact_model(mechanism, foot_points, region_data, region_offset=region_offset)

## Load results
load = true
if load
    if @isdefined result
        @error "result already defined."
    else
        import StaticArrays, StaticUnivariatePolynomials, QPControl.Trajectories # for field types of CentroidalTrajectoryResult
        result = load_result();
    end
end

## Initial conditions
c0 = center_of_mass(state0).v
ċ0 = center_of_mass_velocity(state0).v

contacts0 = map(contact_body_ids) do bodyid # TODO
    sole_frame = sole_frames[bodyid]
    p0 = translation(transform_to_root(state0, sole_frame))
    region_data[1] => p0
end

## Final conditions
cf = nothing
# cf = c0# + SVector(0.8, 0.5, 0.05)
# cf = c0 + SVector(0.05, 0.0, -0.05)

## Basic initial state feasibility checks
for i in eachindex(contacts0)
    p_i_0 = last(contacts0[i])
    @assert norm(c0 - p_i_0) <= max_com_to_contact_distance # TODO: integrate with kinematic limits

    for j in eachindex(contacts0)
        p_j_0 = last(contacts0[j])
        if i != j
            @assert norm(p_i_0 - p_j_0) >= min_inter_contact_distance # TODO: integrate with kinematic limits
        end
    end
end

## Create visualizer
using MeshCat
newvis = false
if newvis || (!@isdefined vis) || isempty(vis.core.scope.pool.connections)
    vis = Visualizer()
    # wait(vis)
end
delete!(vis)

## Centroidal trajectory visualization
cvis = CentroidalTrajectoryVisualizer(vis, region_data, norm(g), length(contact_body_ids))

## Robot visualization
mvis = MechanismVisualizer(mechanism, visuals, vis)
copyto!(mvis, state0)
for sole_frame in values(sole_frames)
    setelement!(mvis, sole_frame, 0.1)
end

## Environment visualization
@time setelement!(mvis, contact_model)#, MeshLambertMaterial(color=RGBA(0.9, 0.9, 0.5, 0.95)))

## GUI
gui = GUI(mvis)
if isempty(vis.core.scope.pool.connections)
    open(gui)
end
sleep(1.)

## Optimizer
optimizer_factory = baron_optimizer_factory()
# optimizer_factory = scip_optimizer_factory()

## Problem
problem = CentroidalTrajectoryProblem(optimizer_factory, region_data, c0, ċ0, contacts0;
    cf=cf, g=g,
    max_cop_distance=max_cop_distance, max_com_to_contact_distance=max_com_to_contact_distance, min_inter_contact_distance=min_inter_contact_distance,
    num_pieces=9, c_degree=3,
    # objective_type=ObjectiveTypes.MIN_EXCURSION);
    objective_type=ObjectiveTypes.FEASIBILITY);

disallow_jumping!(problem)

## Final region constraint
fix.(problem.z_vars[problem.pieces(problem.pieces[end]), problem.regions(problem.regions[end])], 1.0)

if optimizer_factory.constructor == SCIP.Optimizer
    problem.model.optimize_hook = scip_optimize_hook
end

relax = optimizer_factory.constructor == Gurobi.Optimizer || optimizer_factory.constructor == CPLEX.Optimizer
if relax
    @info "Relaxing bilinearities."
    relaxbilinear!(problem.model, method=:Logarithmic1D, disc_level=17, constraints_to_skip=problem.friction_cone_quadratic_constraints)
end

## Feasibility
result = CentroidalTrajOpt.solve!(problem);

if backend(problem.model).optimizer.model isa SCIP.Optimizer
    mscip = backend(problem.model).optimizer.model.mscip
    # SCIP.print_statistics(mscip)
    SCIP.print_heuristic_statistics(mscip)
    # SCIP.SCIPprintVersion(mscip, Libc.FILE(0))
end

reoptimize = false
if reoptimize
    ## Final region constraint
    fix.(problem.z_vars[problem.pieces(problem.pieces[end]), problem.regions(problem.regions[2])], 1.0)

    # ## Re-optimize with final CoM constraint
    # cx_desired = c0[1] + 0.5
    # # cy_desired = c0[2] + 0.5
    # let
    #     pieces = problem.pieces
    #     c_coeffs = problem.c_coeffs
    #     cf = problem.c_vars[pieces(last(pieces)), c_coeffs(last(c_coeffs))]
    #     JuMP.fix(cf[1], cx_desired, force=true)
    #     # JuMP.fix(cf[2], cy_desired, force=true)
    #     # JuMP.fix.(cf, cf_desired, force=true)
    # end
    result = CentroidalTrajOpt.solve!(problem)
end

## Result visualization
set_com_trajectory!(cvis, result)
set_state!(cvis, result, 0.0)
plan_animation = Animation(cvis, result)
setanimation!(vis, plan_animation)

## Optimality
# Setting objective function straight away somehow prevents SCIP from finding a feasible point, so:
# SCIP.SCIPsetEmphasis(problem.model.moi_backend.optimizer.model.mscip, SCIP.SCIP_PARAMEMPHASIS_OPTIMALITY, true);
# set_objective!(problem.model, MOI.MIN_SENSE,
# set_objective(problem.model, MOI.MAX_SENSE, sum(problem.z_vars))
# result = solve!(problem);
# set_objective(problem.model, MOI.MIN_SENSE, sum(problem.Δts))
# result = CentroidalTrajOpt.solve!(problem);

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
# @test ċ(T) ≈ zeros(3) atol=1e-12
# @test c̈(T) ≈ zeros(3) atol=1e-12

for t in range(0, T, length=100)
    ftot = sum(fs[contact](t) for contact in problem.contacts.val)

    # Dynamics
    @test c̈(t) ≈ g + ftot atol=1e-5

    # Torque about CoM
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
        @test norm(p(t) - r(t)) < max_cop_distance + 1e-5

        if region_idx == nothing
            @test f_t ≈ zero(f_t) atol=1e-4
            # @test τn_t ≈ zero(τn_t) atol=1e-4
        else
            region = region_data[region_idx]
            n = ns[problem.regions(region_idx)]
            fn_t = f_t ⋅ n
            μ = region.μ
            @test fn_t >= -1e-7
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

## Simulation
simulate = true
if simulate
    # Pose plans
    pose_plans = Dict{BodyID, PosePlan{Float64}}()
    for (i, bodyid) in enumerate(contact_body_ids)
        pose0 = transform_to_root(state0, sole_frames[bodyid])
        pose_plans[bodyid] = PosePlan(pose0, result.contact_positions[i], result.contact_indicators[i], region_data)
    end

    # Controller
    μ_control = 0.7
    nominal_state = state0 # TODO
    controller = create_controller(mechanism, contact_body_ids, floating_joint, foot_points, μ_control, pelvis,
        state0, nominal_state, result.center_of_mass, pose_plans);

    state = MechanismState(mechanism)
    copyto!(state, state0)
    Δt = 1 / 500
    pcontroller = PeriodicController(similar(velocity(state)), Δt, controller)
    damping = JointDamping{Float64}(mechanism, AtlasRobot.urdfpath())
    dynamics = Dynamics(
        mechanism,
        pcontroller,#SumController(similar(velocity(state)), (pcontroller, damping));
        contact_model=contact_model)
    # callback = CallbackSet(RealtimeRateLimiter(poll_interval=pi / 100), )
    callback = CallbackSet(gui; max_fps=30)
    T = last(result.break_times)
    tspan = (0., T + 2)
    contact_state = SoftContactState(contact_model)
    odeproblem = ODEProblem(dynamics, (state, contact_state), tspan)#, callback=callback)

    # simulate
    @time sol = RigidBodySim.solve(odeproblem, Tsit5(), abs_tol = 1e-8, dt = 1e-6, dtmax=1e-3);
    sim_animation = Animation(mvis, sol)

    # setanimation!(vis, plan_animation);
    setanimation!(vis, sim_animation);
    setanimation!(vis, merge(plan_animation, sim_animation));
end

## Set up visualizer for creating videos
video = false
if video
    open(vis, Blink.Window(Dict(:width => 1280, :height => 720, :useContentSize => true)))
    setanimation!(vis, merge(plan_animation, sim_animation));
end

## Save results
CentroidalTrajOpt.Serialization.save_result(result)
