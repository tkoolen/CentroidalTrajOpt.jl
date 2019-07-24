module CentroidalTrajOpt

export
    CentroidalTrajectoryProblem,
    ContactRegion,
    solve!,
    normals,
    disallow_jumping!,
    ObjectiveTypes

# Re-exports from other packages
export
    SVector,
    RotMatrix,
    AffineMap,
    with_optimizer,
    derivative,
    breaks

# Should move
export
    map_subfunctions

# Visualization
export
    CentroidalTrajectoryVisualizer,
    set_com_trajectory!,
    set_state!,
    save_result,
    load_result

using JuMP
using LinearAlgebra
using StaticArrays
using StaticUnivariatePolynomials
using CoordinateTransformations
using Rotations
using AxisArrays

import QPControl

using QPControl.Trajectories: BezierCurve, Piecewise, Constant
using QPControl.Trajectories: breaks
using StaticUnivariatePolynomials: derivative

const SUP = StaticUnivariatePolynomials

import Polyhedra
using Polyhedra: hrep, vrep, polyhedron, Mesh

include("util.jl")
include("region.jl")
include("problem.jl")
include("result.jl")

using GeometryTypes: Point, Vec, HyperSphere, HyperRectangle, GLNormalMesh

import MeshCat
using MeshCat: Visualizer
using MeshCat: RGB, RGBA
using MeshCat: ArrowVisualizer, HyperSphere, PointCloud, LineSegments, Cone
using MeshCat: MeshLambertMaterial, MeshBasicMaterial, Texture, PngImage, LineBasicMaterial
using MeshCat: setobject!, settransform!, setprop!, setvisible!
using MeshCat: Animation, atframe, setanimation!

include("visualization.jl")
include("serialization.jl")

using .Serialization

end # module
