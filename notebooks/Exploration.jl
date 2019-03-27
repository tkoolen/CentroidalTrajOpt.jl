#!/usr/bin/env python
# coding: utf-8

# In[3]:


using Pkg
Pkg.activate(@__DIR__)


# To do:
#
# * Test that CoPs are in region
# * Test that friction cone constraints are satisfied
# * Move utilities to appropriate places
# * better CoM kinematic constraints
# * Initial condition modification

# In[4]:


using CentroidalTrajOpt
using LinearAlgebra
using SCIP
using JuMP


# In[5]:


# Environment setup
num_regions = 1#2
region_data = ContactRegion{Float64}[]
for i = 1 : num_regions
    region = ContactRegion(
        AffineMap(one(RotMatrix{3}), zero(SVector{3})),
        0.7,
        0.0,
        Float64[1 0; 0 1; -1 0; 0 -1],
        5.0 * ones(4)
    )
    push!(region_data, region)
end


# In[6]:


# Initial conditions
c0 = SVector(-0.05, 0.05, 0.9)
ċ0 = SVector(0.9, 0.1, 0.1)
contacts0 = [
    region_data[1] => SVector(0.0, 0.15, 0.0),
    region_data[1] => SVector(0.0, -0.15, 0.0),
];


# In[7]:


# Additional settings
g = SVector(0.0, 0.0, -9.81);
max_cop_distance = 0.1


# In[11]:


optimizer_factory = with_optimizer(SCIP.Optimizer, limits_gap=0.05, limits_time=300, display_verblevel=5,
    display_width=120, history_valuebased=true);
# display_lpinfo=true
# optimizer_factory = with_optimizer(Ipopt.Optimizer)


# In[12]:


problem = CentroidalTrajectoryProblem(optimizer_factory, region_data, c0, ċ0, contacts0;
    g=g, max_cop_distance=max_cop_distance, num_pieces=5, c_degree=3);


# In[13]:


result = solve!(problem);


# In[74]:


# set_objective(problem.model, MOI.MAX_SENSE, sum(problem.z_vars))
# result = solve!(problem);


# In[75]:


# Setting objective function straight away somehow prevents SCIP from finding a feasible point, so:
# set_objective(problem.model, MOI.MIN_SENSE, sum(problem.Δts))
# result = solve!(problem);


# In[76]:


using Test
c = result.center_of_mass
ċ = map_subfunctions(x -> map_elements(derivative, x), c)
c̈ = map_subfunctions(x -> map_elements(derivative, x), ċ)
rs = result.centers_of_pressure
fs = result.contact_forces
τns = result.contact_normal_torques
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
    τ = sum((rs[j](t) - c(t)) × fs[j](t) + n * τns[j](t) for j in problem.contacts.val)
    @test τ ≈ zeros(3) atol=1e-4

    # Friction cones
    for j in problem.contacts.val
        zs = result.contact_indicators[j](t)
        @test sum(zs) <= 1
        region = findfirst(isequal(1), zs)

        f = fs[j]
        τn = τns[j]
        p = ps[j]
        r = rs[j]

        f_t = f(t)
        fn_t = f_t ⋅ n
        τn_t = τn(t)

        @test fn_t >= 0
        @test norm(p(t) - r(t)) < max_cop_distance + 1e-5

        # TODO: μ checks (first, get region)
#         @test norm(f_t - n * fn_t) <= μ * fn_t + 1e-6
#         @test -μrot * fn_t <= τn_t
#         @test τn_t <= μrot * fn_s

        if region == nothing
            @test f_t ≈ zero(f_t) atol=1e-4
            @test τn_t ≈ zero(τn_t) atol=1e-4
        end
    end
end

# Continuity around breaks
for t in result.break_times
    @test c(t - 1e-8) ≈ c(t + 1e-8) atol=1e-5
    @test ċ(t - 1e-8) ≈ ċ(t + 1e-8) atol=1e-5
end


# In[77]:


objective_value(problem.model)


# In[78]:


# using Plots


# In[79]:


# gr()
# tvals = range(0, T; length=100)
# cvals = [c(t) for t in tvals]
# rvals = [[rs[j](t) for t in tvals] for j in problem.contacts.val]
# ftotvals = [sum(fs[j](t) for j in problem.contacts.val) for t in tvals]
# plt = plot(getindex.(cvals, 1), getindex.(cvals, 2), getindex.(cvals, 3), label="CoM",
#     xlims=[-0.2, 0.2], ylims=[-0.2, 0.2], zlims=[0, 1.2], aspect_ratio=1)
# for j in problem.contacts.val
#     plot!(plt, getindex.(rvals[j], 1), getindex.(rvals[j], 2), getindex.(rvals[j], 3), label="CoP $j")
# end
# plt


# In[80]:


using MeshCat
using LoopThrottle
newvis = false
if newvis || (!@isdefined vis) || (!@isdefined cvis)
    vis = Visualizer()
    cvis = CentroidalTrajectoryVisualizer(vis, g, length(contacts0))
    set_objects!(cvis)
    open(vis)
end


# In[81]:


let
    T = last(result.break_times)
#     max_rate = 0.3
    max_rate = 1
    @throttle t for t in range(0, T, length = round(Int, 60 * T / max_rate))
        set_state!(cvis, result, t)
    end max_rate=max_rate
end


# In[ ]:
