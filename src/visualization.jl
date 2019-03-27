struct CentroidalTrajectoryVisualizer{V<:AbstractVisualizer}
    com_visualizer::V
    force_visualizers::Vector{ArrowVisualizer{V}}
    contact_position_visualizers::Vector{V}
    gravity_mag::Float64
end

function CentroidalTrajectoryVisualizer(vis::AbstractVisualizer, g::AbstractVector, num_contacts::Number)
    com_visualizer = vis[:com]
    force_visualizers = [ArrowVisualizer(vis["f_$i"]) for i = 1 : num_contacts]
    contact_position_visualizers = [vis["p_$i"] for i = 1 : num_contacts]
    gravity_mag = norm(g)
    CentroidalTrajectoryVisualizer(
        com_visualizer, force_visualizers, contact_position_visualizers, gravity_mag)
end

function set_objects!(vis::CentroidalTrajectoryVisualizer)
    # CoM
    # img = PngImage(joinpath(@__DIR__, "..", "assets", "checkerboard.png"))
    # com_material = MeshLambertMaterial(map=Texture(image=img, wrap=(1, 1), repeat=(1, 1)))
    setobject!(vis.com_visualizer, HyperSphere(Point(0., 0, 0), 0.05))#, com_material)

    # Forces
    force_orange = RGB(243 / 255, 118 / 255, 32 / 255)
    for force_visualizer in vis.force_visualizers
        setobject!(force_visualizer, MeshLambertMaterial(; color=force_orange))
    end

    # Contact positions
    for contact_position_visualizer in vis.contact_position_visualizers
        setobject!(contact_position_visualizer, HyperSphere(Point(0., 0, 0), 0.03))
    end
    # TODO: Contact regions
    vis
end

function set_state!(vis::CentroidalTrajectoryVisualizer, result::CentroidalTrajectoryResult, t::Number)
    # CoM
    c = result.center_of_mass(t)
    settransform!(vis.com_visualizer, Translation(c))

    # Forces
    for (force_vis, cop, force) in zip(vis.force_visualizers, result.centers_of_pressure, result.contact_forces)
        r = cop(t)
        f = force(t)
        settransform!(force_vis, Point(r), Vec(f / vis.gravity_mag))
    end

    # Contact positions
    for (position_vis, position, indicator) in zip(vis.contact_position_visualizers, result.contact_positions, result.contact_indicators)
        p = position(t)
        zs = indicator(t)
        if sum(zs) > 0.5
            settransform!(position_vis, Translation(p))
        else
            settransform!(position_vis, LinearMap(Diagonal(SVector(0, 0, 0))))
        end
    end
    vis
end
