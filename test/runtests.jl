module CentroidalTrajOptTest

using CentroidalTrajOpt
using Test
using Random
using StaticArrays
using LinearAlgebra
using Rotations
using CoordinateTransformations

using CentroidalTrajOpt: align_z_axis, closest_pose

@testset "align_z_axis" begin
    rng = MersenneTwister(2)
    for i = 1 : 100
        unaligned = rand(rng, Quat)
        z_axis = normalize(rand(rng, SVector{3}))
        aligned = align_z_axis(unaligned, z_axis)
        @test aligned * SVector(0, 0, 1) ≈ z_axis atol=1e-14
        θ1 = rotation_angle(unaligned \ aligned)
        θ2 = rotation_angle(rotation_between(unaligned * SVector(0, 0, 1), aligned * SVector(0, 0, 1)))
        @test θ1 ≈ θ2 atol=1e-14
    end
end

@testset "closest_pose" begin
    rng = MersenneTwister(2)
    for i = 1 : 100
        previous = AffineMap(rand(rng, Quat), rand(rng, SVector{3}))
        z_axis = normalize(rand(SVector{3}))
        origin = rand(rng, SVector{3})
        next = closest_pose(previous, origin, z_axis)
        @test next.translation === origin
        @test next.linear * SVector(0, 0, 1) ≈ z_axis atol=1e-14
    end
end

end
