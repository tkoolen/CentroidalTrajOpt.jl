# Polynomial utilities
function scale_argument(x::Polynomial{N}, s) where {N}
    Polynomial(ntuple(i -> s^(i - 1) * x.coeffs[i], Val(N)))
end

# JuMP utilities
function simplify(model::JuMP.Model, x::JuMP.QuadExpr)
    canonical = JuMP.MOIU.canonical # TODO
    QuadExpr(model, canonical(JuMP.moi_function.(x)))
end

JuMP.value(x::Polynomial) = Polynomial(JuMP.value.(x.coeffs))
JuMP.value(x::BezierCurve) = BezierCurve(JuMP.value.(x.points))
JuMP.value(x::Piecewise) = Piecewise(JuMP.value.(x.subfunctions), JuMP.value.(x.breaks), x.clamp)

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
    @constraint model SVector(x.points) .== SVector(y.points)
end

function constrain_poly_equal(model, x::Polynomial{N}, y::Number) where N
    constrain_poly_equal(model, x, Polynomial{N}(Polynomial(y)))
end

function axis_array_vars(model::Model, namefunc, axes...)
    vars = map(Iterators.product(map(axis -> axis.val, axes)...)) do I
        @variable model base_name=namefunc(I...)
    end
    AxisArray(vars, axes...)
end

# scrap
# function poly(model::JuMP.Model, namefunc, degree::Union{Integer, Val})
#     Polynomial(ntuple(j -> @variable(model, base_name=namefunc(j)), degree + 1))
# end

# function poly_svec(model::JuMP.Model, namefunc, degree::Union{Integer, Val}, length::Union{Integer, Val})
#     SVector(ntuple(i -> poly(model, j -> namefunc(i, j), degree), length))
# end
