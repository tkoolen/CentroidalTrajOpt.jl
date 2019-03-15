module CentroidalTrajOpt

export
    CentroidalTrajectoryProblem,
    solve!

# Re-exports from other packages
export
    SVector,
    with_optimizer

using JuMP
using LinearAlgebra
using StaticArrays
using StaticUnivariatePolynomials
using CoordinateTransformations
using Rotations
using AxisArrays

import QPControl

using QPControl.Trajectories: BezierCurve, derivative, Piecewise

const SUP = StaticUnivariatePolynomials

include("util.jl")
include("problem.jl")

end # module
