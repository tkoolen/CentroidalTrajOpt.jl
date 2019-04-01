module CentroidalTrajOpt

export
    CentroidalTrajectoryProblem,
    ContactRegion,
    solve!,
    center_of_mass,
    normals,
    disallow_jumping!

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
    map_subfunctions,
    map_elements

# Visualization
export
    CentroidalTrajectoryVisualizer,
    set_objects!,
    set_com_trajectory!,
    set_state!

using JuMP
using LinearAlgebra
using StaticArrays
using StaticUnivariatePolynomials
using CoordinateTransformations
using Rotations
using AxisArrays

import QPControl

using QPControl.Trajectories: BezierCurve, Piecewise, Constant, derivative, breaks
using StaticUnivariatePolynomials: derivative

const SUP = StaticUnivariatePolynomials

include("util.jl")
include("region.jl")
include("problem.jl")
include("result.jl")

using GeometryTypes: HyperSphere, Point, Vec
import MeshCat
using MeshCat: AbstractVisualizer
using MeshCat: RGB, RGBA
using MeshCat: ArrowVisualizer, HyperSphere, PointCloud, LineSegments, Cone
using MeshCat: MeshLambertMaterial, Texture, PngImage, LineBasicMaterial
using MeshCat: setobject!, settransform!
using Polyhedra: hrep, polyhedron, Mesh

include("visualization.jl")

end # module
