## center_of_mass_velocity
using RigidBodyDynamics
function center_of_mass_velocity(state::MechanismState)
    h = momentum(state)
    FreeVector3D(h.frame, linear(h)) / mass(state.mechanism)
end

## hrep conversion
function RigidBodyDynamics.Contact.HRep(hr::Polyhedra.HRepresentation)
    M = Polyhedra.nhalfspaces(hr)
    N = Polyhedra.fulldim(hr)
    return map(Polyhedra.halfspaces(hr)) do halfspace
        RigidBodyDynamics.Contact.HalfSpace(SVector{N}(halfspace.a), halfspace.β)
    end |> SVector{M} |> RigidBodyDynamics.Contact.HRep
end

function Polyhedra.hrep(hr::RigidBodyDynamics.Contact.HRep)
    # hrep(Vector([Polyhedra.HalfSpace(halfspace.outward_normal, halfspace.offset) for halfspace in hr.halfspaces]))
    hrep([Polyhedra.HalfSpace(halfspace.outward_normal, halfspace.offset) for halfspace in hr.halfspaces])
end

## setelement! for ContactModel
using MeshCatMechanisms, RigidBodyDynamics
function MeshCatMechanisms.setelement!(mvis::MechanismVisualizer, contact_model::ContactModel, args...; base_name="collision element ")
    mechanism = mvis.state.mechanism
    for group in contact_model.collision_groups
        for element in group
            geometry = element.geometry
            frame = element.transform.from
            if geometry isa MeshCatMechanisms.GeometryLike
                setelement!(mvis, frame, element.geometry, args...)
            elseif geometry isa Contact.HRep
                setelement!(mvis, frame, GLNormalMesh(polyhedron(hrep(geometry))))
                setobject!(mvis[frame][:triad], Triad(0.1))
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

function check_kinematic_constraints_satisfied(contacts, c0, max_com_to_contact_distance, min_inter_contact_distance)
    ## Basic initial state feasibility checks
    for i in eachindex(contacts)
        p_i_0 = last(contacts[i])
        norm(c0 - p_i_0) <= max_com_to_contact_distance || error()
        for j in eachindex(contacts)
            p_j_0 = last(contacts[j])
            if i != j
                norm(p_i_0 - p_j_0) >= min_inter_contact_distance || error()
            end
        end
    end
end
