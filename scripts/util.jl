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

## PosePlan
using QPControl.Trajectories: Piecewise, Constant
using CentroidalTrajOpt: closest_pose

function QPWalkingControl.PosePlan(pose0::Transform3D,
        point_trajectory::Piecewise{<:Constant{<:SVector{3}}},
        contact_indicators::Piecewise{<:Constant{<:AbstractVector{Bool}}},
        region_data::AbstractVector{<:ContactRegion})
    @assert point_trajectory.breaks == contact_indicators.breaks
    T = eltype(pose0)
    plan = PosePlan{T}()
    previous_region_index = findfirst(contact_indicators.subfunctions[1].value)
    previous_pose = pose0
    start_time = nothing
    n = length(point_trajectory.subfunctions)
    for i = 2 : n
        region_index = findfirst(contact_indicators.subfunctions[i].value)
        if region_index == previous_region_index
            continue
        elseif region_index === nothing
            start_time = point_trajectory.breaks[i]
            previous_region_index = region_index
            continue
        else
            point = point_trajectory.subfunctions[i].value
            region = region_data[region_index]
            normal = region.transform.linear * SVector(0, 0, 1)
            final_pose = Transform3D(previous_pose.from, previous_pose.to, closest_pose(AffineMap(previous_pose), point, normal))
            duration = point_trajectory.breaks[i] - start_time
            push!(plan, start_time, duration, final_pose)
            previous_pose = final_pose
            previous_region_index = region_index
        end
    end
    return plan
end
