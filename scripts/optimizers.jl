using JuMP

using SCIP
function scip_optimizer_factory()
    with_optimizer(SCIP.Optimizer)#, limits_time=2 * 60, heuristics_alns_priority=1000000, heuristics_subnlp_priority=1000000)
    # with_optimizer(SCIP.Optimizer,
    #     limits_gap=0.05, limits_time=10 * 60 * 60, display_verblevel=5,
    #     display_width=120, history_valuebased=true, lp_threads=10, branching_preferbinary=true, lp_scaling=false,
    #     branching_allfullstrong_priority=536870911, heuristics_multistart_freq=20, heuristics_multistart_onlynlps=false,
    #     heuristics_mpec_priority=536870911)#heuristics_subnlp_priority=536870911)#, nlp_solver="ipopt", heuristics_nlpdiving_priority=536870911)#,;
end

function scip_optimize_hook(model::JuMP.Model)
    caching_optimizer = backend(model)
    if caching_optimizer.mode == MOIU.AUTOMATIC && caching_optimizer.state == MOIU.EMPTY_OPTIMIZER
        MOIU.attach_optimizer(caching_optimizer)
    end
    mscip = caching_optimizer.optimizer.model.mscip
    SCIP.SCIPsetEmphasis(mscip, SCIP.SCIP_PARAMEMPHASIS_FEASIBILITY, true)#false); println()
    SCIP.SCIPsetPresolving(mscip, SCIP.SCIP_PARAMSETTING_AGGRESSIVE, true)#false); println()
    SCIP.SCIPsetHeuristics(mscip, SCIP.SCIP_PARAMSETTING_AGGRESSIVE, true)#false)
    # SCIP.set_parameter(mscip, "heuristics/alns/freq", 5)
    # SCIP.set_parameter(mscip, "heuristics/subnlp/freq", 5)
    # SCIP.set_parameter(mscip, "heuristics/mpec/freq", 5)
    # SCIP.set_parameter(mscip, "heuristics/feaspump/freq", 5)
    SCIP.set_parameter(mscip, "heuristics/rens/freq", -1)
    SCIP.set_parameter(mscip, "reoptimization/enable", true)
    SCIP.set_parameter(mscip, "limits/time", 2 * 60)
    # SCIP.set_parameter(mscip, "limits/time", 2 * 60 * 60)
    ret = MOI.optimize!(caching_optimizer)
    @assert caching_optimizer.optimizer.model.mscip === mscip
    return ret
end

using AmplNLWriter
AmplNLWriter.setdebug(true)

function bonmin_optimizer_factory()
    with_optimizer(AmplNLWriter.Optimizer, "/home/twan/code/bonmin/Bonmin-1.8.7/build/bin/bonmin")
end

function couenne_optimizer_factory()
    with_optimizer(AmplNLWriter.Optimizer, "/home/twan/code/couenne/couenne")
end


using Ipopt
# ENV["JULIA_IPOPT_LIBRARY_PATH"] = "/home/twan/code/ipopt/CoinIpopt/build/lib"
# ENV["JULIA_IPOPT_EXECUTABLE_PATH"] = "/home/twan/code/ipopt/CoinIpopt/build/bin"
# import Pkg; Pkg.build("Ipopt")

function ipopt_optimizer_factory()
    with_optimizer(Ipopt.Optimizer)
end

function complementarity_reformulate!(model::JuMP.Model)
    for var in JuMP.all_variables(problem.model)
        if JuMP.is_binary(var)
            @constraint model var * (1 - var) == 0
            JuMP.unset_binary(var)
        end
    end
end

# using Alpine
# function alpine_optimizer_factory()
#     with_optimizer(Alpine.Optimizer, nlp_optimizer=Ipopt.Optimizer(print_level=0), mip_optimizer=Gurobi.Optimizer(OutputFlag=0))
# end

# using Juniper
# function juniper_optimizer_factory()
#     params = Dict{Symbol,Any}()
#     params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level=0)
#     # params[:mip_solver] = with_optimizer(Gurobi.Optimizer, OutputFlag=0)
#     params[:mip_solver] = with_optimizer(CPLEX.Optimizer, CPX_PARAM_SCRIND=0)
#     # Note to self: if you get "Cannot set bounds because variable is of type: BINARY", use LinQuadOptInterface master (needs #91)
#     params[:feasibility_pump] = false
#     with_optimizer(Juniper.Optimizer, params)
# end

using BARON
function baron_optimizer_factory()
    with_optimizer(BARON.Optimizer;
         CplexLibName="/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/libcplex1280.so",
         threads=Sys.CPU_THREADS ÷ 2, MaxTime=5 * 60 * 60, PrTimeFreq=5.,
         AllowFilterSD=1, AllowFilterSQP=1, AllowIpopt=1#=, NumLoc=20, LocRes=1=#)
end

using CPLEX
function cplex_optimizer_factory_continuous()
    # https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html

    # Parameters tried:
    # CPX_PARAM_OPTIMALITYTARGET = 2: CPX_OPTIMALITYTARGET_FIRSTORDER - not available for MIP
    # CPX_PARAM_FPHEUR = 1 - not available for MIQCP
    # CPX_PARAM_DEPIND = 3 # dependency checking: turn on at beginning and at end of preprocessing
    # CPX_PARAM_MIQCPSTRAT = 1 - solve QCPs as subproblems
    # CPXPARAM_MIP_Cuts_LiftProj = 3 - Generate lift-and-project cuts very aggressively
    # CPX_PARAM_MIPEMPHASIS = CPLEX.CPX_MIPEMPHASIS_HIDDENFEAS - applies considerable additional effort toward finding high quality feasible solutions that are difficult to locate

    # Parameters to try:
    # CPX_PARAM_MIQCPSTRAT
    # CPX_MIPEMPHASIS_FEASIBILITY
    # barrier-related options
    with_optimizer(CPLEX.Optimizer,
        CPX_PARAM_OPTIMALITYTARGET = CPLEX.CPX_OPTIMALITYTARGET_AUTO,
        CPX_PARAM_DEPIND = 3,
        CPXPARAM_MIP_Cuts_LiftProj = 3)
end

function cplex_optimizer_factory()
    with_optimizer(CPLEX.Optimizer)
end

using Gurobi
function gurobi_optimizer_factory()
    with_optimizer(Gurobi.Optimizer)
end
