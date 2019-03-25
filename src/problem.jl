struct CentroidalTrajectoryProblem
    model::JuMP.Model

    # Axes
    c_coeffs
    f_coeffs
    contacts
    regions
    pieces
    coords

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

    # Environment data. TODO: constructor arg.
    num_regions = 2
    region_data = ContactRegion{Float64}[]
    for i = 1 : num_regions
        region = ContactRegion(
            AffineMap(one(RotMatrix{3}), zero(SVector{3})),
            0.7,
            0.0,
            Float64[1 0; 0 1; -1 0; 0 -1],
            0.3 * ones(4)
        )
        push!(region_data, region)
    end

    # Axes
    c_coeffs = Axis{:c_coeff}(1 : c_num_coeffs)
    f_coeffs = Axis{:f_coeff}(1 : f_num_coeffs)
    contacts = Axis{:contact}(1 : num_contacts)
    regions = Axis{:region}(1 : num_regions)
    pieces = Axis{:piece}(1 : num_pieces)
    coords = Axis{:coord}(SVector(:x, :y, :z))

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
    #Δts = nothing
    if Δts === nothing
        Δts = axis_array_vars(model, i -> "Δt[$i]", pieces)
        Δtsqs = axis_array_vars(model, i -> "Δtsq[$i]", pieces)
        for (Δt, Δtsq) in zip(Δts, Δtsqs)
            Δtmin = 0.5
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

    # Contact indicators. z[pieces(i), contacts(j), regions(m)] == 1 implies that
    # during piece i, contact j is in contact with region m.
    z_vars = axis_array_vars(model, (i, j, m) -> "z[p$i, c$j, r$m]", pieces, contacts, regions)
    foreach(set_binary, z_vars)

    # Contact normals
    Rs = AxisArray(map(data -> data.transform.linear, region_data), regions)
    ns = map(R -> R[:, 3], Rs)

    # Continuous variables
    c_vars = axis_array_vars(model, (i, k, l) -> "C[p$i, $k, $l]", pieces, coords, c_coeffs)
    f_vars = axis_array_vars(model, (i, j, k, l) -> "F[p$i, c$j, $k, $l]", pieces, contacts, coords, f_coeffs)
    f̄_vars = axis_array_vars(model, (i, j, k, l, m) -> "F̄[p$i, c$j, $k, $l, r$m]", pieces, contacts, coords, f_coeffs, regions)
    p_vars = axis_array_vars(model, (i, j, k) -> "P[p$i, c$j, $k]", pieces, contacts, coords)
    p̄_vars = axis_array_vars(model, (i, j, k) -> "P̄[p$i, c$j, $k]", pieces, contacts, coords)
    r_vars = axis_array_vars(model, (i, j, k, l) -> "R[p$i, c$j, $k, $l]", pieces, contacts, coords, c_coeffs)
    r̄_vars = axis_array_vars(model, (i, j, k, l) -> "R̄[p$i, c$j, $k, $l]", pieces, contacts, coords, c_coeffs)
    τn_vars = axis_array_vars(model, (i, j, l) -> "Tn[p$i, c$j, $l]", pieces, contacts, f_coeffs)

    #@variable model objval
    #@constraint model objval >= sum((c_vars[pieces(i), c_coeffs(l)] - c0) ⋅ (c_vars[pieces(i), c_coeffs(l)] - c0) for i in 1 : num_pieces, l in 1 : c_num_coeffs)
    #@objective model Min objval

    cprev = nothing
    c′prev = nothing

    for i in 1 : num_pieces
        piece = pieces(i)

        # CoM and derivatives
        c = [BezierCurve(c_vars[piece, coords(k)]...) for k in 1 : length(coords)]
        c′ = derivative.(c)
        c′′ = derivative.(c′)

        # Contact/region assignment constraints
        for j in 1 : num_contacts
            # Each contact can be assigned to at most one region
            contact = contacts(j)
            @constraint model z_vars[piece, contact] in MOI.SOS1(collect(Float64, 1 : num_regions))
            #@constraint model sum(z_vars[piece, contact]) <= 1
            if i > 1
                # Each contact must be unassigned for one piece before it can be reassigned to a region.
                # Let Δzᵢ,ⱼ,ₘ = zᵢ,ⱼ,ₘ - zᵢ₋₁,ⱼ,ₘ. Then there are three cases for ∑ₘ |Δzᵢ,ⱼ,ₘ|:
                # 1) ∑ₘ |Δzᵢ,ⱼ,ₘ| = 0: the assignment during piece i is the same as the assignment during piece i - 1
                # 2) ∑ₘ |Δzᵢ,ⱼ,ₘ| = 1: unassigned during piece i and assigned during piece i - 1 or vice versa
                # 3) ∑ₘ |Δzᵢ,ⱼ,ₘ| = 2: different assignment during piece i and piece i - 1
                # The third case must be prevented, so we want ∑ₘ |Δzᵢ,ⱼ,ₘ| ≤ 1.
                # This is an ℓ₁-norm constraint.
                Δz = z_vars[piece, contact] - z_vars[pieces(i - 1), contact]
                constrain_l1_norm(model, Δz, 1; add_bounds=true)

                # Contact position p may only change when the contact is not assigned to a region,
                # i.e., when ∑ₘ zᵢ,ⱼ,ₘ = 0.
                Δp = p_vars[piece, contact] - p_vars[pieces(i - 1), contact]
                Δpmax = 0.5 # TODO
                @constraint model sum(x -> x^2, Δp) <= (1 - sum(z_vars[piece, contact])) * Δpmax^2
            end
        end

        # Initial conditions
        if i == 1
            @constraint model map(x -> x(0), c) .== c0
            @constraint model map(x -> x(0), c′) .== ċ0 .* Δts[i]
        end

        if i == 1
            for j in 1 : num_contacts
                contact = contacts(j)
                initial_region = regions(1) # TODO: make configurable
                @constraint model z_vars[piece, contact, initial_region] == 1

                # TODO: initial contact position constraint
            end
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

        # Forces (global)
        fs = [[BezierCurve(f_vars[piece, contacts(j), coords(k)]...) for k in coords.val] for j in contacts.val]

        # CoPs (global)
        rs = [[BezierCurve(r_vars[piece, contacts(j), coords(k)]...) for k in coords.val] for j in contacts.val]

        # Normal torques (global)
        τns = [BezierCurve(τn_vars[piece, contacts(j)]...) for j in contacts.val]

        # Dynamics
        ftot = sum(fs)
        constrain_poly_equal.(model, c′′, identity((g + ftot) .* Δtsqs[i]))

        # Torque balance
        τ = sum(Polynomial.(rs[j]) × Polynomial.(fs[j]) for j = 1 : num_contacts) +
            sum(Polynomial.(sum(ns[regions(m)] .* z_vars[piece, regions(m), contacts(j)] for m in 1 : num_regions) .* τns[j]) for j = 1 : num_contacts) -
            Polynomial.(c) × Polynomial.(ftot)
        τ = map(x -> Polynomial(simplify.(model, x.coeffs)), τ)
        constrain_poly_equal.(model, τ, 0)

        # CoP and contact position constraints
        Mr = 5 # TODO
        for j in 1 : num_contacts
            contact = contacts(j)
            p = p_vars[piece, contact]
            p̄ = p̄_vars[piece, contact]
            p̄xy = p̄[1 : 2]
            p̄z = p̄[3]
            for m in 1 : num_regions
                region = regions(m)
                A = region_data[m].A
                b = region_data[m].b
                transform = region_data[m].transform
                z = z_vars[piece, contact, region]

                # Contact position constraints
                @constraint model A * p̄xy .<= b + Mr * (1 - z)
                @constraint model  p̄z <= Mr * (1 - z)
                @constraint model -p̄z <= Mr * (1 - z)
                @constraint model p .== transform(p̄)

                # CoP constraints
                for l = 1 : c_num_coeffs
                    coeff = c_coeffs(l)
                    r = r_vars[piece, contact, coeff]
                    r̄ = r̄_vars[piece, contact, coeff]
                    r̄xy = r̄[1 : 2]
                    r̄z = r̄[3]
                    @constraint model A * r̄xy .<= b + Mr * (1 - z)
                    @constraint model  r̄z <= Mr * (1 - z)
                    @constraint model -r̄z <= Mr * (1 - z)
                    @constraint model r .== transform(r̄)
                    @constraint model sum(x -> x^2, r̄ - p̄) <= max_cop_distance^2
                end
            end
        end

        # Contact force/torque constraints
        Mf = 5 * norm(g) # TODO
        Mτ = norm(g)
        for j in 1 : num_contacts
            contact = contacts(j)
            for m in 1 : num_regions
                region = regions(m)
                μ = region_data[m].μ
                μrot = region_data[m].μrot
                R = Rs[region]
                n = ns[region]
                z = z_vars[piece, contact, region]
                for l in 1 : f_num_coeffs
                    coeff = f_coeffs(l)
                    f̄ = f̄_vars[piece, contact, coeff, region]
                    f̄xy = f̄[1 : 2]
                    f̄z = f̄[3]
                    set_lower_bound(f̄z, 0)
                    @constraint model f̄z <= Mf * z
                    μf̄z = @variable model lower_bound=0 # needed for SecondOrderCone constraint
                    @constraint model μf̄z == μ * f̄z
                    @constraint model [μf̄z; f̄xy] in SecondOrderCone()
                    #@constraint model f̄xy ⋅ f̄xy <= μf̄z^2
                    # TODO: use a SOS condition?
                    # TODO: crude model
                    τn = τn_vars[piece, contact, coeff]
                    @constraint model  τn <= μrot * f̄z + Mτ * z
                    @constraint model -τn <= μrot * f̄z + Mτ * z
                end
            end
            for l in 1 : f_num_coeffs
                coeff = f_coeffs(l)
                @constraint model f_vars[piece, contact, coeff] .== sum(Rs[regions(m)] * f̄_vars[piece, contact, coeff, regions(m)] for m in 1 : num_regions)
            end
        end

        # CoM kinematic constraints
        for j in 1 : num_contacts
            contact = contacts(j)
            p = p_vars[piece, contact]
            for l in 1 : c_num_coeffs
                cpoint = map(x -> x.points[l], c)
                @constraint model sum(x -> x^2, cpoint - p) <= 1.5 # TODO
                @constraint model cpoint[3] - p[3] >= 0.6 # TODO
            end
        end

        cprev = c
        c′prev = c′
    end

    CentroidalTrajectoryProblem(model,
        c_coeffs, f_coeffs, contacts, regions, pieces, coords,
        c_vars, f_vars, f̄_vars, p_vars, r_vars, r̄_vars, τn_vars, Δts, z_vars,
        ns
    )
end

normals(problem) = problem.normals
