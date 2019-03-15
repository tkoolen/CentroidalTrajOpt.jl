struct CentroidalTrajectoryProblem
    model::JuMP.Model

    # Axes
    c_coeffs
    f_coeffs
    contacts
    regions
    pieces
    coords
    coord2ds

    # Variables
    c_vars
    f_vars
    f̄_vars
    p_vars
    r_vars
    r̄_vars
    τn_vars
    Δts
    z_vars

    # Other data
    normals
end

function CentroidalTrajectoryProblem(optimizer_factory::JuMP.OptimizerFactory;
        c0::AbstractVector,
        ċ0::AbstractVector,
        c_degree = 3,
        num_pieces = 2,
        num_contacts = 2, # TODO
        g = SVector(0.0, 0.0, -9.81),
        max_cop_distance = 0.1
    )

    c_num_coeffs = c_degree + 1
    f_num_coeffs = c_num_coeffs - 2

    num_regions = 1 # TODO

    # Axes
    c_coeffs = Axis{:c_coeff}(1 : c_num_coeffs)
    f_coeffs = Axis{:f_coeff}(1 : f_num_coeffs)
    contacts = Axis{:contact}(1 : num_contacts)
    regions = Axis{:region}(1 : num_regions)
    pieces = Axis{:piece}(1 : num_pieces)
    coords = Axis{:coord}(SVector(:x, :y, :z))
    coord2ds = Axis{:coord2d}(SVector(:x, :y))

    # Environment data. TODO: rework, constructor args.
    transforms = AxisArray(fill(AffineMap(one(RotMatrix{3}), zero(SVector{3})), num_regions), regions)
    μs = AxisArray(fill(0.7, num_regions), regions)
    μrots = AxisArray(fill(0.0, num_regions), regions)
    As = AxisArray(fill([1 0; 0 1; -1 0; 0 -1], num_regions), regions)
    bs = AxisArray(fill([0.5, 0.5, 0.5, 0.5], num_regions), regions);

    # Time stuff
    Δt_tolerance = 0.05 # TODO: constructor arg

    model = Model(optimizer_factory)

    # Indexing convention:
    # i: piece index
    # j: contact index
    # k: coordinate index
    # l: coefficient index
    # m: region index

    Δts = AxisArray(fill(2.0, length(pieces)), pieces)
    # const Δts = nothing
    if Δts === nothing
        Δts = axis_array_vars(model, i -> "Δt[$i]", pieces)
        Δtsqs = axis_array_vars(model, i -> "Δtsq[$i]", pieces)
        for (Δt, Δtsq) in zip(Δts, Δtsqs)
            Δtmin = Δt_tolerance # TODO
            Δtmax = 5.0 # TODO
            set_lower_bound(Δt, Δtmin)
            set_upper_bound(Δt, Δtmax)
            set_lower_bound(Δtsq, Δtmin^2)
            set_upper_bound(Δtsq, Δtmax^2)
            @constraint model Δtsq == Δt^2
        end
        @objective model Min sum(Δts)
    else
        Δtsqs = map(x -> x^2, Δts)
    end

    # Contact normals
    ns = map(transform -> transform.linear[:, 3], transforms)

    c_vars = axis_array_vars(model, (i, k, l) -> "C[p$i, $k, $l]", pieces, coords, c_coeffs)
    f_vars = axis_array_vars(model, (i, j, k, l) -> "F[p$i, c$j, $k, $l]", pieces, contacts, coords, f_coeffs)
    f̄_vars = axis_array_vars(model, (i, j, k, l) -> "F̄[p$i, c$j, $k, $l]", pieces, contacts, coords, f_coeffs)
    p_vars = axis_array_vars(model, (i, j, k) -> "P[p$i, c$j, $k]", pieces, contacts, coords)
    p̄_vars = axis_array_vars(model, (i, j, k) -> "P̄[p$i, c$j, $k]", pieces, contacts, coord2ds)
    r_vars = axis_array_vars(model, (i, j, k, l) -> "R[p$i, c$j, $k, $l]", pieces, contacts, coords, c_coeffs)
    r̄_vars = axis_array_vars(model, (i, j, k, l) -> "R̄[p$i, c$j, $k, $l]", pieces, contacts, coord2ds, c_coeffs)
    τn_vars = axis_array_vars(model, (i, j, l) -> "Tn[p$i, c$j, $l]", pieces, contacts, c_coeffs)
    z_vars = axis_array_vars(model, (i, j, m) -> "z[p$i, c$j, r$m]", pieces, contacts, regions)
    foreach(set_binary, z_vars)

    cprev = nothing
    c′prev = nothing

    for i in 1 : num_pieces
        # Contact/region assignment constraints

        for j in 1 : num_contacts
            # Each contact can be assigned to at most one region
            @constraint model z_vars[pieces(i), contacts(j)] in MOI.SOS1(collect(Float64, 1 : num_regions))
        end

        # CoM and derivatives
        c = [BezierCurve(c_vars[pieces(i), coords(k)]...) for k in 1 : length(coords)]
        c′ = derivative.(c)
        c′′ = derivative.(c′)

        # Initial conditions
        if i == 1
            @constraint model map(x -> x(0), c) .== c0
            @constraint model map(x -> x(0), c′) .== ċ0 .* Δts[i]
        end

        # Final conditions
        if i == num_pieces
            @constraint model map(x -> x(1), c) .== c0 # TODO
            @constraint model map(x -> x(1), c′) .== 0
            @constraint model map(x -> x(1), c′′) .== 0
        end

        # CoM continuity
        if i > 1
            @constraint model map(x -> x(1), cprev) .== map(x -> x(0), c)
            @constraint model map(x -> x(1), c′prev) .* Δts[i] .== map(x -> x(0), c′) .* Δts[i - 1]
        end

        # Forces (global and local)
        fs = [[BezierCurve(f_vars[pieces(i), contacts(j), coords(k)]...) for k in coords.val] for j in contacts.val]
        f̄s = [[BezierCurve(f̄_vars[pieces(i), contacts(j), coords(k)]...) for k in coords.val] for j in contacts.val]

        # Contact positions (global and local)
        ps = [p_vars[pieces(i), contacts(j)] for j in contacts.val]
        p̄s = [p̄_vars[pieces(i), contacts(j)] for j in contacts.val]

        # CoPs (global and local)
        rs = [[BezierCurve(r_vars[pieces(i), contacts(j), coords(k)]...) for k in coords.val] for j in contacts.val]
        r̄s = [[BezierCurve(r̄_vars[pieces(i), contacts(j), coord2ds(k)]...) for k in coord2ds.val] for j in contacts.val]

        # Normal forces
        τns = [BezierCurve(τn_vars[pieces(i), contacts(j)]...) for j in contacts.val]

        # Dynamics
        ftot = sum(fs)
        constrain_poly_equal.(model, c′′, identity((g + ftot) .* Δtsqs[i]))

        # Torque balance
        # TODO: assumes single region
        τ = sum(Polynomial.(rs[j]) × Polynomial.(fs[j]) + ns[regions(1)] .* Polynomial(τns[j]) for j = 1 : num_contacts) -
            Polynomial.(c) × Polynomial.(sum(fs))
        τ = map(x -> Polynomial(simplify.(model, x.coeffs)), τ)
        constrain_poly_equal.(model, τ, 0)

        # CoP and contact position constraints
        for (p, p̄, r, r̄) in zip(ps, p̄s, rs, r̄s)
            @assert length(regions) == 1 # TODO
            A = As[regions(1)]
            b = bs[regions(1)]
            transform = transforms[regions(1)]

            # Contact position constraints
            @constraint model A * p̄ .<= b # TODO: could lose this constraint
            @constraint model p .== transform([p̄; zero(eltype(p̄))])

            # TODO: initial contact position constraints

            # CoP constraints
            for control_point in [map(x -> x.points[l], r̄) for l = 1 : c_num_coeffs]
                @constraint model A * control_point .<= b
                Δ = control_point - p̄
                @constraint model Δ ⋅ Δ <= max_cop_distance^2
            end
            constrain_poly_equal.(model, r, transform([r̄; zero(eltype(r̄))]))
        end

        # Contact force constraints
        for (f, f̄, τn) in zip(fs, f̄s, τns)
            @assert length(regions) == 1 # TODO
            region = regions(1)
            transform = transforms[region]
            n = ns[region]
            μ = μs[region]
            μrot = μrots[region]
            R = transform.linear
            constrain_poly_equal.(model, f, R * f̄)
            for l in 1 : f_num_coeffs
                f̄_control_point = map(x -> x.points[l], f̄)
                fₜ = f̄_control_point[1 : 2]
                fₙ = f̄_control_point[3]
                μfₙ = @variable model lower_bound=0
                @constraint model μfₙ == μ * fₙ
        #         @constraint model [μfₙ; fₜ] in SecondOrderCone()
                @constraint model fₜ ⋅ fₜ <= μfₙ^2
                # TODO: use a SOS condition?
                @constraint model -μrot * fₙ <= τn.points[l]
                @constraint model τn.points[l] <= μrot * fₙ
            end
        end

        # CoM kinematic constraints
        for transform in transforms
            p = transform.translation
            for l in 1 : c_num_coeffs
                c_control_point = map(x -> x.points[l], c)
                Δ = c_control_point - p
                @constraint model Δ ⋅ Δ <= 1.5 # TODO
            end
        end

        cprev = c
        c′prev = c′
    end

    CentroidalTrajectoryProblem(model,
        c_coeffs, f_coeffs, contacts, regions, pieces, coords, coord2ds,
        c_vars, f_vars, f̄_vars, p_vars, r_vars, r̄_vars, τn_vars, Δts, z_vars,
        ns
    )
end

normals(problem) = problem.normals
