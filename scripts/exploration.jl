# To do:
#
# * Test that CoPs are in region
# * Test that friction cone constraints are satisfied
# * Move utilities to appropriate places
# * better CoM kinematic constraints
# * Initial condition modification

using Pkg
Pkg.activate(@__DIR__)

using CentroidalTrajOpt
using LinearAlgebra
using AmplNLWriter
using Ipopt
using SCIP
using BARON
# using Alpine
# using Juniper
using CPLEX
using Gurobi
using JuMP
using Rotations
using MultilinearOpt

# Environment setup
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

# Initial conditions
c0 = SVector(-0.05, 0.05, 1.0)
ċ0 = SVector(0.5, 0.4, 0.0)
# ċ0 = SVector(0.05, 0.0, 0.0)
contacts0 = [
    region_data[1] => SVector(0.0, 0.15, 0.0),
    region_data[1] => SVector(0.0, -0.15, 0.0),
];

# Additional settings
g = SVector(0.0, 0.0, -9.81);
max_cop_distance = 0.07

# optimizer_factory = with_optimizer(SCIP.Optimizer, limits_gap=0.05, limits_time=10 * 60 * 60, display_verblevel=5,
#     display_width=120, history_valuebased=true, lp_threads=10, branching_preferbinary=true, lp_scaling=false,
#     branching_allfullstrong_priority=536870911, heuristics_multistart_freq=20, heuristics_multistart_onlynlps=false, heuristics_mpec_priority=536870911)#heuristics_subnlp_priority=536870911)#, nlp_solver="ipopt", heuristics_nlpdiving_priority=536870911)#,;
# optimizer_factory = with_optimizer(AmplNLWriter.Optimizer, "/home/twan/code/bonmin/Bonmin-1.8.7/build/bin/bonmin")
# optimizer_factory = with_optimizer(AmplNLWriter.Optimizer, "/home/twan/code/couenne/couenne")
# optimizer_factory = with_optimizer(Ipopt.Optimizer)
# optimizer_factory = with_optimizer(Alpine.Optimizer, nlp_optimizer=Ipopt.Optimizer(print_level=0), mip_optimizer=Gurobi.Optimizer(OutputFlag=0))
# optimizer_factory = let params = Dict{Symbol,Any}()
#     params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level=0)
#     # params[:mip_solver] = with_optimizer(Gurobi.Optimizer, OutputFlag=0)
#     params[:mip_solver] = with_optimizer(CPLEX.Optimizer, CPX_PARAM_SCRIND=0)
#     # Note to self: if you get "Cannot set bounds because variable is of type: BINARY", use LinQuadOptInterface master (needs #91)
#     params[:feasibility_pump] = false
#     with_optimizer(Juniper.Optimizer, params)
# end
# optimizer_factory = with_optimizer(BARON.Optimizer;
#     threads=Sys.CPU_THREADS ÷ 2, MaxTime=10 * 60.0, PrTimeFreq=5., AllowFilterSD=1, AllowFilterSQP=1, AllowIpopt=1#=, NumLoc=20, LocRes=1=#)
optimizer_factory = with_optimizer(Gurobi.Optimizer)

problem = CentroidalTrajectoryProblem(optimizer_factory, region_data, c0, ċ0, contacts0;
    g=g, max_cop_distance=max_cop_distance, num_pieces=5, c_degree=3, objective_type=ObjectiveTypes.FEASIBILITY);

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
    relaxbilinear!(problem.model, method=:Logarithmic1D, disc_level=11, constraints_to_skip=problem.friction_cone_quadratic_constraints)
end

## Visualization
using MeshCat
using LoopThrottle
newvis = false
if newvis || (!@isdefined vis) || isempty(vis.core.scope.pool.connections)
    vis = Visualizer()
    open(vis)
    # wait(vis)
end
delete!(vis)
cvis = CentroidalTrajectoryVisualizer(vis, region_data, norm(g), length(contacts0))
set_objects!(cvis)
sleep(0.1)

## Feasibility
result = solve!(problem);

if backend(problem.model).optimizer.model isa SCIP.Optimizer
    mscip = backend(problem.model).optimizer.model.mscip
    # SCIP.print_statistics(mscip)
    SCIP.print_heuristic_statistics(mscip)
end

## Optimality

# Setting objective function straight away somehow prevents SCIP from finding a feasible point, so:
# SCIP.SCIPsetEmphasis(problem.model.moi_backend.optimizer.model.mscip, SCIP.SCIP_PARAMEMPHASIS_OPTIMALITY, true);
# set_objective(problem.model, MOI.MAX_SENSE, sum(problem.z_vars))
# result = solve!(problem);
# set_objective(problem.model, MOI.MIN_SENSE, sum(problem.Δts))
# result = solve!(problem);

if optimizer_factory.constructor == BARON.Optimizer
    @show backend(problem.model).optimizer.model.optimizer.inner.problem_file_name
    @show backend(problem.model).optimizer.model.optimizer.inner.result_file_name
end

## Initial state
set_com_trajectory!(cvis, result)
set_state!(cvis, result, 0.0)

## Visualization
setanimation!(cvis, result);

## Tests
using Test
c = result.center_of_mass
ċ = map_subfunctions(x -> map_elements(derivative, x), c)
c̈ = map_subfunctions(x -> map_elements(derivative, x), ċ)
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

# let
#     T = last(result.break_times)
#     max_rate = 1 / 2
#     # max_rate = 1
#     @throttle t for t in range(0, T, length = round(Int, 60 * T / max_rate))
#         set_state!(cvis, result, t)
#     end max_rate=max_rate
# end

##
# set_state!(cvis, result, 1.0)
# set_state!(cvis, result, 3 / 3 * (last(result.break_times) - first(result.break_times)))

## Mode sequence
# value.(problem.z_vars)
