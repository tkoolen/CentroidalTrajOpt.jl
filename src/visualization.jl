hide!(vis::AbstractVisualizer) = settransform!(vis, LinearMap(0I))

struct CentroidalTrajectoryVisualizer{V<:AbstractVisualizer}
    vis::V
    com_visualizer::V
    region_data::Vector{ContactRegion{Float64}}
    region_visualizers::Vector{V}
    force_visualizers::Vector{ArrowVisualizer{V}}
    cone_visualizers::Vector{V}
    contact_position_visualizers::Vector{V}
    gravity_mag::Float64
end

function CentroidalTrajectoryVisualizer(vis::AbstractVisualizer,
        region_data::Vector{ContactRegion{Float64}}, g::AbstractVector, num_contacts::Number)
    com_visualizer = vis[:com]
    region_visualizers = [vis["region_$i"] for i in eachindex(region_data)]
    force_visualizers = [ArrowVisualizer(vis["f_$i"]) for i = 1 : num_contacts]
    cone_visualizers = [vis["cone_$i"] for i = 1 : num_contacts]
    contact_position_visualizers = [vis["p_$i"] for i = 1 : num_contacts]
    gravity_mag = norm(g)
    CentroidalTrajectoryVisualizer(vis, com_visualizer, region_data, region_visualizers,
        force_visualizers, cone_visualizers, contact_position_visualizers, gravity_mag)
end

function set_objects!(vis::CentroidalTrajectoryVisualizer)
    # Contact regions
    for (region_vis, region) in zip(vis.region_visualizers, vis.region_data)
        thickness = 0.05
        A, b = region.A, region.b
        global_to_local = inv(region.transform)
        n = size(A, 1)
        Ā = [A zeros(n, 1); 0 0 1; 0 0 -1]
        b̄ = [b; 0; thickness]
        h = hrep(Ā * global_to_local.linear, b̄ - Ā * global_to_local.translation)
        setobject!(region_vis, Mesh(polyhedron(h)))
    end

    # CoM
    # img = PngImage(joinpath(@__DIR__, "..", "assets", "ps-neutral.png"))
    # com_material = MeshLambertMaterial(map=Texture(image=img))#, wrap=(1, 1), repeat=(1, 1)))
    setobject!(vis.com_visualizer, HyperSphere(Point(0., 0, 0), 0.05))#, com_material)
    setvisible!(vis.com_visualizer, false)

    # Forces / contact force cones
    force_orange = RGB(243 / 255, 118 / 255, 32 / 255)
    force_orange_transparent = RGBA(force_orange, 0.2)
    for (force_visualizer, cone_visualizer) in zip(vis.force_visualizers, vis.cone_visualizers)
        setobject!(force_visualizer, MeshLambertMaterial(; color=force_orange))
        cone_height = 0.2
        cone = Cone(Point(0., 0., cone_height), Point(0., 0., 0.), cone_height)
        setobject!(cone_visualizer, cone, MeshLambertMaterial(color=force_orange_transparent))
        setvisible!(cone_visualizer, false)
    end

    # Contact positions
    for contact_position_visualizer in vis.contact_position_visualizers
        setobject!(contact_position_visualizer, HyperSphere(Point(0., 0, 0), 0.02), MeshLambertMaterial(color=RGB(0.1, 0.1, 0.1)))
        setvisible!(contact_position_visualizer, false)
    end
    vis
end

function set_com_trajectory!(vis::CentroidalTrajectoryVisualizer, result::CentroidalTrajectoryResult)
    ts = range(first(result.break_times), last(result.break_times); length=100)
    geometry = PointCloud([Point(result.center_of_mass(t)) for t in ts])
    setobject!(vis.vis[:com_trajectory], LineSegments(geometry, LineBasicMaterial()))
    vis
end

function set_state!(vis::CentroidalTrajectoryVisualizer, result::CentroidalTrajectoryResult, t::Number)
    # CoM
    c = result.center_of_mass(t)
    setvisible!(vis.com_visualizer, true)
    settransform!(vis.com_visualizer, Translation(c))

    # Forces
    for (force_vis, cone_vis, cop, force, indicator) in zip(vis.force_visualizers, vis.cone_visualizers, result.centers_of_pressure, result.contact_forces, result.contact_indicators)
        r = cop(t)
        f = force(t)
        settransform!(force_vis, Point(r), Vec(f / vis.gravity_mag))
        zs = indicator(t)
        if sum(zs) > 0.5
            setvisible!(cone_vis, true)
            region_index = argmax(zs)
            region = vis.region_data[region_index]
            n = region.transform.linear[:, 3]
            μ = region.μ
            R = rotation_between(SVector(0, 0, 1), n)
            settransform!(cone_vis, Translation(r) ∘ LinearMap(R) ∘ LinearMap(Diagonal(SVector(1, 1, 1 / μ))))
        else
            setvisible!(cone_vis, false)
        end
    end

    # Contact positions
    for (position_vis, position, indicator) in zip(vis.contact_position_visualizers, result.contact_positions, result.contact_indicators)
        p = position(t)
        zs = indicator(t)
        if sum(zs) > 0.5
            setvisible!(position_vis, true)
            settransform!(position_vis, Translation(p))
        else
            setvisible!(position_vis, false)
        end
    end
    vis
end
