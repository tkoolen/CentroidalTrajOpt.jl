# hide!(vis::Visualizer) = settransform!(vis, LinearMap(0I))

struct CentroidalTrajectoryVisualizer
    vis::Visualizer
    com_visualizer::Visualizer
    region_data::Vector{ContactRegion{Float64}}
    force_visualizers::Vector{ArrowVisualizer{Visualizer}}
    cone_visualizers::Vector{Visualizer}
    contact_position_visualizers::Vector{Visualizer}
    gravity_mag::Float64
end

function CentroidalTrajectoryVisualizer(vis::Visualizer,
        region_data::Vector{ContactRegion{Float64}}, gravity_mag, num_contacts::Number)

    # CoM
    com_visualizer = vis[:com]
    img = PngImage(joinpath(@__DIR__, "..", "assets", "checkerboard.png"))
    com_material = MeshBasicMaterial(map=Texture(image=img))#, wrap=(1, 1), repeat=(1, 1)))
    setobject!(com_visualizer, HyperSphere(Point(0., 0, 0), 0.05), com_material)
    setvisible!(com_visualizer, false)

    # Forces / contact force cones
    force_visualizers = [ArrowVisualizer(vis["f_$i"]) for i = 1 : num_contacts]
    cone_visualizers = [vis["cone_$i"] for i = 1 : num_contacts]
    force_orange = RGB(243 / 255, 118 / 255, 32 / 255)
    force_orange_transparent = RGBA(force_orange, 0.2)
    for (force_visualizer, cone_visualizer) in zip(force_visualizers, cone_visualizers)
        setobject!(force_visualizer, MeshLambertMaterial(; color=force_orange))
        cone_height = 0.2
        cone = Cone(Point(0., 0., cone_height), Point(0., 0., 0.), cone_height)
        setobject!(cone_visualizer, cone, MeshLambertMaterial(color=force_orange_transparent))
        setvisible!(cone_visualizer, false)
    end

    # Contact positions
    contact_position_visualizers = [vis["p_$i"] for i = 1 : num_contacts]
    for contact_position_visualizer in contact_position_visualizers
        setobject!(contact_position_visualizer, HyperSphere(Point(0., 0, 0), 0.02), MeshLambertMaterial(color=RGB(0.1, 0.1, 0.1)))
        setvisible!(contact_position_visualizer, false)
    end

    # Regions
    for (i, region) in enumerate(region_data)
        extruded = extrude(hrep(region.A, region.b), 1e-4; zmax=1e-4)
        world_to_body = inv(region.transform)
        geometry = GLNormalMesh(polyhedron(hrep(extruded.A * world_to_body.linear, extruded.b - extruded.A * world_to_body.translation)))
        setobject!(vis["region_$i"], geometry, MeshLambertMaterial(color=RGB(0.5, 0.5, 0.5)))
    end

    CentroidalTrajectoryVisualizer(vis, com_visualizer, region_data,
        force_visualizers, cone_visualizers, contact_position_visualizers, gravity_mag)
end

function set_com_trajectory!(vis::CentroidalTrajectoryVisualizer, result::CentroidalTrajectoryResult)
    ts = range(first(result.break_times), last(result.break_times); length=100)
    geometry = PointCloud([Point(result.center_of_mass(t)) for t in ts])
    setobject!(vis.vis[:com_trajectory], LineSegments(geometry, LineBasicMaterial()))
    vis
end

function set_state!(vis::CentroidalTrajectoryVisualizer, result::CentroidalTrajectoryResult, t::Number; show_cones::Bool=false)
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
        if sum(zs) > 0.5 && show_cones
            setvisible!(cone_vis, true)
            region_index = argmax(zs)
            region = vis.region_data[region_index]
            n = region.transform.linear[:, 3]
            μ = region.μ
            R = rotation_between(SVector(0, 0, 1), n)
            settransform!(cone_vis, Translation(r) ∘ LinearMap(R) ∘ LinearMap(Diagonal(SVector(1, 1, 1 / μ))))
        else
            setvisible!(cone_vis, false)
            settransform!(cone_vis, LinearMap(0I))
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
            settransform!(position_vis, LinearMap(0I))
        end
    end
    vis
end

function MeshCat.Animation(
        cvis::CentroidalTrajectoryVisualizer,
        result::CentroidalTrajectoryResult;
        fps::Number=30)
    t0 = first(result.break_times)
    tf = last(result.break_times)
    animation = Animation(fps)
    num_frames = floor(Int, (tf - t0) * fps)
    for frame in 0 : num_frames
        t = t0 + frame / fps
        atframe(animation, frame) do
            set_state!(cvis, result, t)
        end
    end
    return animation
end
