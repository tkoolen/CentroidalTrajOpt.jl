# Trajectory utilities
# struct Vectorized{T}
#     trajectories::T
# end

# (vectorized::Vectorized)(x, args...) = map(traj -> traj(x, args...), vectorized.trajectories)
# map_elements(f, vectorized::Vectorized) = Vectorized(map(f, vectorized.trajectories))

function map_subfunctions(f, traj::Piecewise)
    Piecewise(map(f, traj.subfunctions), traj.breaks; clamp=traj.clamp)
end

# Polynomial utilities
function scale_argument(x::Polynomial{N}, s) where {N}
    Polynomial(ntuple(i -> s^(i - 1) * x.coeffs[i], Val(N)))
end

# TODO: type piracy; should generalize to any N
function (p::Polynomial)(x, ::Val{2})
    pd = SUP.derivative(p)
    pdd = SUP.derivative(pd)
    p(x), pd(x), pdd(x)
end

# JuMP utilities
function simplify(model::JuMP.Model, x::JuMP.QuadExpr)
    canonical = JuMP.MOIU.canonical # TODO
    QuadExpr(model, canonical(JuMP.moi_function.(x)))
end

val(x) = JuMP.value(x)
val(x::Number) = x

function bezier(model::JuMP.Model, namefunc, degree::Union{Integer, Val})
    BezierCurve(ntuple(j -> @variable(model, base_name=namefunc(j)), degree + 1))
end

function bezier_svec(model::JuMP.Model, namefunc, degree::Union{Integer, Val}, length::Union{Integer, Val})
    SVector(ntuple(i -> bezier(model, j -> namefunc(i, j), degree), length))
end

function constrain_poly_equal(model, x::Polynomial, y::Polynomial)
    @constraint model SVector(x.coeffs) .== SVector(y.coeffs)
end

function constrain_poly_equal(model, x::BezierCurve, y::BezierCurve)
    @constraint model SVector(x.coeffs) .== SVector(y.coeffs)
end

function constrain_poly_equal(model, x::Polynomial{N}, y::Number) where N
    constrain_poly_equal(model, x, Polynomial{N}(Polynomial(y)))
end

function axis_array_vars(model::Model, namefunc, axes...; lower_bound=nothing, upper_bound=nothing)
    vars = map(Iterators.product(map(axis -> axis.val, axes)...)) do I
        @variable model base_name=namefunc(I...)
    end
    for var in vars
        if lower_bound !== nothing
            set_lower_bound(var, lower_bound)
        end
        if upper_bound !== nothing
            set_upper_bound(var, upper_bound)
        end
    end
    AxisArray(vars, axes...)
end

function constrain_l1_norm(model::JuMP.Model, x::AbstractVector{<:JuMP.AbstractJuMPScalar}, value; add_bounds)
    @assert add_bounds
    s = @variable model [1 : length(x)] lower_bound=0 upper_bound=value
    @constraint model s .>=  x
    @constraint model s .>= -x
    @constraint model sum(s) <= value
    model
end

# Geometry utilities
function extrude(polyhedron::Polyhedra.HRepresentation, thickness::Number; zmax=0)
    A = polyhedron.A
    b = polyhedron.b
    n = size(A, 1)
    Ā = [A zeros(n, 1); 0 0 1; 0 0 -1]
    b̄ = [b; zmax; -zmax + thickness]
    hrep(Ā, b̄, polyhedron.linset)
end

# scrap
# function poly(model::JuMP.Model, namefunc, degree::Union{Integer, Val})
#     Polynomial(ntuple(j -> @variable(model, base_name=namefunc(j)), degree + 1))
# end

# function poly_svec(model::JuMP.Model, namefunc, degree::Union{Integer, Val}, length::Union{Integer, Val})
#     SVector(ntuple(i -> poly(model, j -> namefunc(i, j), degree), length))
# end

# JuMP.value(x::Polynomial) = Polynomial(JuMP.value.(x.coeffs))
# JuMP.value(x::BezierCurve) = BezierCurve(JuMP.value.(x.coeffs))
# JuMP.value(x::Piecewise) = Piecewise(JuMP.value.(x.subfunctions), JuMP.value.(x.breaks), x.clamp)

function align_z_axis(rot::Rotation{3}, z_axis::StaticVector{3})
    rot_z_axis = rot * SVector(0, 0, 1)
    rotation_between(rot_z_axis, z_axis) * rot
end

function closest_pose(previous::AffineMap, origin::StaticVector{3}, z_axis::StaticVector{3})
    rot = align_z_axis(previous.linear, z_axis)
    AffineMap(rot, origin)
end
