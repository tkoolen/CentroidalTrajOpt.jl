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
