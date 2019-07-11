## center_of_mass_velocity
using RigidBodyDynamics
function center_of_mass_velocity(state::MechanismState)
    h = momentum(state)
    FreeVector3D(h.frame, linear(h)) / mass(state.mechanism)
end

## setelement! for ContactModel
using MeshCatMechanisms, RigidBodyDynamics
function MeshCatMechanisms.setelement!(mvis::MechanismVisualizer, contact_model::ContactModel, args...; base_name="collision element ")
    mechanism = mvis.state.mechanism
    for group in contact_model.collision_groups
        for element in group
            geometry = element.geometry
            if geometry isa MeshCatMechanisms.GeometryLike
                frame = element.transform.from
                setelement!(mvis, frame, element.geometry, args...)
            end
        end
    end
end

## pose_trajectory
using QPControl.Trajectories: Piecewise, Constant
using CentroidalTrajOpt: closest_pose
function pose_trajectory(
        pose0::Transform3D,
        point_trajectory::Piecewise{<:Constant{<:SVector{3}}},
        contact_indicators::Piecewise{<:Constant{<:AbstractVector{Bool}}},
        region_data::AbstractVector{<:ContactRegion})
    @assert point_trajectory.breaks == contact_indicators.breaks
    n = length(point_trajectory.subfunctions)
    subfunctions = Vector{Union{Constant{typeof(pose0)}, Constant{Nothing}}}(undef, n)
    subfunctions[1] = Constant(pose0)
    previous = pose0
    for i = 2 : n
        time = point_trajectory.breaks[i]
        region_index = findfirst(contact_indicators.subfunctions[i].value)
        subfunctions[i] = if region_index !== nothing
            point = point_trajectory.subfunctions[i].value
            region = region_data[region_index]
            normal = region.transform.linear * SVector(0, 0, 1)
            pose = Transform3D(previous.from, previous.to, closest_pose(AffineMap(previous), point, normal))
            previous = pose
            Constant(pose)
        else
            Constant(nothing)
        end
    end
    return Piecewise(subfunctions, point_trajectory.breaks, clamp=true)
end
