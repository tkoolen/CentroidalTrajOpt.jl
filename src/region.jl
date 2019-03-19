struct ContactRegion{T}
    transform::AffineMap{RotMatrix3{T}, SVector{3, T}}
    μ::T
    μrot::T
    A::Matrix{T}
    b::Vector{T}
end
