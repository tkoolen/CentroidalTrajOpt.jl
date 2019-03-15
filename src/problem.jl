struct CentroidalTrajectoryProblem
    model::JuMP.Model

    # Axes
    c_coeffs
    f_coeffs
    contacts
    pieces
    coords
    coord2ds

    # Variables
    c_vars
    f_vars
    f̄_vars
    r_vars
    r̄_vars
    τn_vars
    Δts

    # Other data
    normals
end

function CentroidalTrajectoryProblem(optimizer_factory::JuMP.OptimizerFactory;
        c0::AbstractVector,
        ċ0::AbstractVector,
        c_degree = 3,
        num_pieces = 2,
        num_contacts = 2, # TODO
        g = SVector(0.0, 0.0, -9.81)
    )

    c_num_coeffs = c_degree + 1
    f_num_coeffs = c_num_coeffs - 2

    # Axes
    c_coeffs = Axis{:c_coeff}(1 : c_num_coeffs)
    f_coeffs = Axis{:f_coeff}(1 : f_num_coeffs)
    contacts = Axis{:contact}(1 : num_contacts)
    pieces = Axis{:piece}(1 : num_pieces)
    coords = Axis{:coord}(SVector(:x, :y, :z))
    coord2ds = Axis{:coord2d}(SVector(:x, :y))

    # Environment data. TODO: rework, constructor args.
    transforms = AxisArray(fill(AffineMap(one(RotMatrix{3}), zero(SVector{3})), num_contacts), contacts)
    μs = AxisArray(fill(0.7, num_contacts), contacts)
    μrots = AxisArray(fill(0.0, num_contacts), contacts)
    As = AxisArray(fill([1 0; 0 1; -1 0; 0 -1], num_contacts), contacts)
    bs = AxisArray(fill([0.5, 0.5, 0.5, 0.5], num_contacts), contacts);

    # Time stuff
    Δt_tolerance = 0.05 # TODO: constructor arg

    model = Model(optimizer_factory)

    # Indexing convention:
    # i: piece index
    # j: contact index
    # k: coordinate index
    # l: coefficient index

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

    c_vars = axis_array_vars(model, (i, k, l) -> "C[$i, $k, $l]", pieces, coords, c_coeffs)
    f_vars = axis_array_vars(model, (i, j, k, l) -> "F[$i, $j, $k, $l]", pieces, contacts, coords, f_coeffs)
    f̄_vars = axis_array_vars(model, (i, j, k, l) -> "F̄[$i, $j, $k, $l]", pieces, contacts, coords, f_coeffs)
    r_vars = axis_array_vars(model, (i, j, k, l) -> "R[$i, $j, $k, $l]", pieces, contacts, coords, c_coeffs)
    r̄_vars = axis_array_vars(model, (i, j, k, l) -> "R̄[$i, $j, $k, $l]", pieces, contacts, coord2ds, c_coeffs)
    τn_vars = axis_array_vars(model, (i, j, l) -> "Tn[$i, $j, $l]", pieces, contacts, c_coeffs)

    cprev = nothing
    c′prev = nothing

    for i in 1 : num_pieces
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

        # CoPs (global and local)
        rs = [[BezierCurve(r_vars[pieces(i), contacts(j), coords(k)]...) for k in coords.val] for j in contacts.val]
        r̄s = [[BezierCurve(r̄_vars[pieces(i), contacts(j), coord2ds(k)]...) for k in coord2ds.val] for j in contacts.val]

        # Normal forces
        τns = [BezierCurve(τn_vars[pieces(i), contacts(j)]...) for j in contacts.val]

        # Dynamics
        ftot = sum(fs)
        constrain_poly_equal.(model, c′′, identity((g + ftot) .* Δtsqs[i]))

        # Torque balance
        τ = sum(Polynomial.(rs[j]) × Polynomial.(fs[j]) + ns[j] .* Polynomial(τns[j]) for j = 1 : num_contacts) -
            Polynomial.(c) × Polynomial.(sum(fs))
        τ = map(x -> Polynomial(simplify.(model, x.coeffs)), τ)
        constrain_poly_equal.(model, τ, 0)

        # CoP constraints
        for (r̄, A, b) in zip(r̄s, As, bs)
            for control_point in [map(x -> x.points[l], r̄) for l = 1 : c_num_coeffs]
                @constraint model A * control_point .<= b
            end
        end
        for (r, r̄, transform) in zip(rs, r̄s, transforms)
            constrain_poly_equal.(model, r, transform([r̄; zero(eltype(r̄))]))
        end

        # Contact force constraints
        for (f, f̄, τn, transform, n, μ, μrot) in zip(fs, f̄s, τns, transforms, ns, μs, μrots)
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
        c_coeffs, f_coeffs, contacts, pieces, coords, coord2ds,
        c_vars, f_vars, f̄_vars, r_vars, r̄_vars, τn_vars, Δts,
        ns
    )
end

normals(problem) = problem.normals
