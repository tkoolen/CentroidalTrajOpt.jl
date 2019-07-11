module CentroidalTrajOptTest

using CentroidalTrajOpt
using Test
using Random
using StaticArrays
using LinearAlgebra
using Rotations

using CentroidalTrajOpt: align_z_axis

@testset "align_z_axis" begin
    rng = MersenneTwister(2)
    for i = 1 : 100
        unaligned = rand(Quat)
        z_axis = normalize(rand(SVector{3}))
        aligned = align_z_axis(unaligned, z_axis)
        @test aligned * SVector(0, 0, 1) ≈ z_axis atol=1e-14
        θ1 = rotation_angle(unaligned \ aligned)
        θ2 = rotation_angle(rotation_between(unaligned * SVector(0, 0, 1), aligned * SVector(0, 0, 1)))
        @test θ1 ≈ θ2 atol=1e-14
    end
end

end
