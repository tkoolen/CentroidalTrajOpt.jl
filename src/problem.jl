struct CentroidalTrajectoryProblem
    model::JuMP.Model

    # Axes
    c_coeffs
    f_coeffs
    contacts
    regions
    pieces
    coords
    coords2d

    # Variables
    c_vars
    f_vars
    f̄_vars
    r_vars
    p_vars
    p̄_vars
    Δts
    z_vars

    # Other data
    normals
    friction_cone_quadratic_constraints
end

module ObjectiveTypes
export ObjectiveType
@enum ObjectiveType begin
    FEASIBILITY
    MIN_TIME
    MIN_EXCURSION
    MAX_HEIGHT
end
end
using .ObjectiveTypes

function CentroidalTrajectoryProblem(optimizer_factory::JuMP.OptimizerFactory,
        region_data::AbstractVector{<:ContactRegion},
        c0::AbstractVector{<:Number},
        ċ0::AbstractVector{<:Number},
        contacts0::AbstractVector{<:Union{Pair{<:ContactRegion, <:AbstractVector}, Nothing}};
        cf = nothing,
        com_degree = 3,
        cop_degree = com_degree,
        num_pieces = 2,
        g = SVector(0.0, 0.0, -9.81),
        max_cop_distance = 0.1,
        min_com_to_contact_distance,
        max_com_to_contact_distance,
        min_inter_contact_distance,
        objective_type::ObjectiveType = FEASIBILITY,
        min_Δt = 0.6,
        max_Δt = 1.5,
        c_margin_xy = 0.5,
        c_margin_z_min = 0.5,
        c_margin_z_max = 1.2,
        max_force = 3 * norm(g),
        Δrmax = 0.7 # TODO
    )

    c_num_coeffs = com_degree + 1
    f_num_coeffs = c_num_coeffs - 2
    p_num_coeffs = cop_degree + 1
    num_regions = length(region_data)
    num_contacts = length(contacts0)

    # Axes
    c_coeffs = Axis{:c_coeff}(1 : c_num_coeffs)
    f_coeffs = Axis{:f_coeff}(1 : f_num_coeffs)
    p_coeffs = Axis{:r_coeff}(1 : p_num_coeffs)
    contacts = Axis{:contact}(1 : num_contacts)
    regions = Axis{:region}(1 : num_regions)
    pieces = Axis{:piece}(1 : num_pieces)
    coords = Axis{:coord}(SVector(:x, :y, :z))
    coords2d = Axis{:coord}(SVector(:x, :y))

    # Time stuff
    Δt_tolerance = 0.05 # TODO: constructor arg

    model = Model(optimizer_factory)
    optimizer = backend(model)#.optimizer.model

    optimizer_does_sos1 = MOI.supports_constraint(optimizer, MOI.VectorOfVariables, MOI.SOS1{Float64})
    optimizer_does_sos1 || @info "Optimizer does not support SOS1 constraints; using sum <= 1 instead."

    optimizer_does_soc = MOI.supports_constraint(optimizer, MOI.VectorOfVariables, MOI.SecondOrderCone)
    optimizer_does_soc || @info "Optimizer does not support second order cone constraints; using quadratic constraints instead."

    optimizer_does_indicator_constraints = false#optimizer.optimizer.model isa SCIP.Optimizer # TODO: generalize

    # Indexing convention:
    # i: piece index
    # j: contact index
    # k: coordinate index
    # l: coefficient index
    # m: region index

    if min_Δt == max_Δt
        Δts = AxisArray(fill(min_Δt, length(pieces)), pieces)
        Δtsqs = map(x -> x^2, Δts)
    else
        Δts = axis_array_vars(model, i -> "Δt[$i]", pieces)
        Δtsqs = axis_array_vars(model, i -> "Δtsq[$i]", pieces)
        for (Δt, Δtsq) in zip(Δts, Δtsqs)
            set_lower_bound(Δt, min_Δt)
            set_upper_bound(Δt, max_Δt)
            set_lower_bound(Δtsq, min_Δt^2)
            set_upper_bound(Δtsq, max_Δt^2)
            @constraint model Δtsq == Δt^2
        end
        if objective_type == ObjectiveTypes.MIN_TIME
            @objective model Min sum(Δts)
        end
    end
    # Δts = AxisArray([ 1.0910477104335914
    #     0.6545247110231425
    #     0.6014271582441553
    #     0.848420003941298 ], pieces)
    # Δtsqs = map(x -> x^2, Δts)

    # Contact indicators. z[pieces(i), contacts(j), regions(m)] == 1 implies that
    # during piece i, contact j is in contact with region m.
    z_vars = axis_array_vars(model, (i, j, m) -> "z[p$i, c$j, r$m]", pieces, contacts, regions)
    foreach(set_binary, z_vars)

    # Contact normals
    Rs = AxisArray(map(data -> data.transform.linear, region_data), regions)
    ns = map(R -> R[:, 3], Rs)

    # Contact location extrema (local coordinates)
    r̄_extrema = AxisArray(map(region -> polyhedron_extrema(polyhedron(hrep(region.A, region.b))), region_data), regions)
    p̄_extrema = map(((p̄_min, p̄_max),) -> (p̄_min .- max_cop_distance, p̄_max .+ max_cop_distance), r̄_extrema)
    r_vertices = mapreduce(vcat, 1 : num_regions) do m
        let region = region_data[m]
            map(Polyhedra.points(polyhedron(hrep(region.A, region.b)))) do p̄_vertex
                region.transform([p̄_vertex; 0])
            end
        end
    end
    r_min, r_max = polyhedron_extrema(polyhedron(vrep(r_vertices)))
    p_min, p_max = r_min .- max_cop_distance, r_max .+ max_cop_distance
    c_min, c_max = p_min - SVector(c_margin_xy, c_margin_xy, -c_margin_z_min), p_max + SVector(c_margin_xy, c_margin_xy, c_margin_z_max)

    # Continuous variables
    c_vars = axis_array_vars(model, (i, k, l) -> "C[p$i, $k, $l]",pieces, coords, c_coeffs)
    f_vars = axis_array_vars(model, (i, j, k, l) -> "F[p$i, c$j, $k, $l]", pieces, contacts, coords, f_coeffs)
    f̄_vars = axis_array_vars(model, (i, j, k, l, m) -> "F̄[p$i, c$j, $k, $l, r$m]", pieces, contacts, coords, f_coeffs, regions)
    r_vars = axis_array_vars(model, (i, j, k) -> "P[p$i, c$j, $k]", pieces, contacts, coords)
    r̄_vars = axis_array_vars(model, (i, j, k, m) -> "P̄[p$i, c$j, $k, r$m]", pieces, contacts, coords2d, regions)
    p_vars = axis_array_vars(model, (i, j, k, l) -> "R[p$i, c$j, $k, $l]", pieces, contacts, coords, p_coeffs)
    p̄_vars = axis_array_vars(model, (i, j, k, l, m) -> "R̄[p$i, c$j, $k, $l, r$m]", pieces, contacts, coords2d, p_coeffs, regions)

    #@variable model objval
    if objective_type == ObjectiveTypes.MIN_EXCURSION
        # @constraint model objval >= sum((c_vars[pieces(i), c_coeffs(l)] - c0) ⋅ (c_vars[pieces(i), c_coeffs(l)] - c0) for i in 1 : num_pieces, l in 1 : c_num_coeffs)
        @objective model Min sum(x -> x ⋅ x, (c_vars[pieces(i), c_coeffs(l)] - c0) for i in 1 : num_pieces, l in 1 : c_num_coeffs)
    elseif objective_type == ObjectiveTypes.MAX_HEIGHT
        @objective model Max sum(c_vars[c_coeffs(3)])
    else
        @objective model Min 0
    end
    # @objective model Max sum(z_vars)

    cprev = nothing
    c′prev = nothing
    friction_cone_quadratic_constraints = []

    for i in 1 : num_pieces
        piece = pieces(i)

        # CoM and derivatives
        c = [BezierCurve(c_vars[piece, coords(k)]...) for k in 1 : length(coords)]
        c′ = derivative.(c)
        c′′ = derivative.(c′)

        # CoM bounds
        for l in 1 : c_num_coeffs
            cpoint = map(x -> x.coeffs[l], c)
            set_lower_bound.(cpoint, c_min)
            set_upper_bound.(cpoint, c_max)
        end

        # Contact/region assignment constraints
        for j in 1 : num_contacts
            contact = contacts(j)

            # Each contact can be assigned to at most one region
            if optimizer_does_sos1
                @constraint model z_vars[piece, contact] in MOI.SOS1(collect(Float64, 1 : num_regions))
            else
                @constraint model sum(z_vars[piece, contact]) <= 1
            end

            r = r_vars[piece, contact]
            set_lower_bound.(r, r_min)
            set_upper_bound.(r, r_max)
            if i == 1
                # Initial contact assignment.
                if contacts0[j] === nothing
                    @constraint model sum(z_vars[piece, contact]) == 0
                else
                    region_data0, r0 = contacts0[j]
                    m = findfirst(isequal(region_data0), region_data)
                    z_var = z_vars[piece, contact, regions(m)]
                    fix(z_var, 1, force=true)
                    # unset_binary(z_var) # to make Alpine happpy
                    fix.(r, r0, force=true)
                end
            else
                # Each contact must be unassigned for one piece before it can be reassigned to a region.
                # Let Δzᵢ,ⱼ,ₘ = zᵢ,ⱼ,ₘ - zᵢ₋₁,ⱼ,ₘ. Then there are three cases for ∑ₘ |Δzᵢ,ⱼ,ₘ|:
                # 1) ∑ₘ |Δzᵢ,ⱼ,ₘ| = 0: the assignment during piece i is the same as the assignment during piece i - 1
                # 2) ∑ₘ |Δzᵢ,ⱼ,ₘ| = 1: unassigned during piece i and assigned during piece i - 1 or vice versa
                # 3) ∑ₘ |Δzᵢ,ⱼ,ₘ| = 2: different assignment during piece i and piece i - 1
                # The third case must be prevented, so we want ∑ₘ |Δzᵢ,ⱼ,ₘ| ≤ 1.
                # This is an ℓ₁-norm constraint.
                Δz = z_vars[piece, contact] - z_vars[pieces(i - 1), contact]
                constrain_l1_norm(model, Δz, 1; add_bounds=true)

                # Two swing segments in a row is not allowed.
                @constraint model sum(z_vars[piece, contact]) + sum(z_vars[pieces(i - 1), contact]) >= 1

                # Contact reference position r may only change when the contact is not assigned to a region,
                # i.e., when ∑ₘ zᵢ,ⱼ,ₘ = 0.
                Δr = r - r_vars[pieces(i - 1), contact]
                # @constraint model sum(x -> x^2, Δp) <= (1 - sum(z_vars[piece, contact])) * Δpmax^2
                z_no_contact = @variable model
                set_binary(z_no_contact)
                @constraint model z_no_contact == 1 - sum(z_vars[piece, contact])
                @constraint model  Δr .<= z_no_contact * Δrmax
                @constraint model -Δr .<= z_no_contact * Δrmax
            end
        end

        # Initial conditions
        if i == 1
            fix.(map(x -> x.coeffs[1], c), c0, force=true)
            @constraint model map(x -> x(0), c′) .== ċ0 .* Δts[i]
        end

        # Final conditions
        if i == num_pieces
            if cf !== nothing
                fix.(map(x -> x.coeffs[end], c), cf, force=true)
            end
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
        ps = [[BezierCurve(p_vars[piece, contacts(j), coords(k)]...) for k in coords.val] for j in contacts.val]

        # Dynamics
        ftot = sum(fs)
        constrain_poly_equal.(model, c′′, identity((g + ftot) .* Δtsqs[i]))

        # Torque balance
        τ = sum(Polynomial.(ps[j]) × Polynomial.(fs[j]) for j = 1 : num_contacts) -
            Polynomial.(c) × Polynomial.(ftot)
        τ = map(x -> Polynomial(simplify.(model, x.coeffs)), τ)
        constrain_poly_equal.(model, τ, 0)

        # CoP and contact position constraints
        # TODO: could derive tighter big-M from constraints on p and p̄, but it's unclear whether that will actually help performance,
        # because changing Mr to
        Mr = 2.0
        for j in 1 : num_contacts
            contact = contacts(j)
            r = r_vars[piece, contact]
            for m in 1 : num_regions
                region = regions(m)
                r̄ = r̄_vars[piece, contact, region]
                r̄_min, r̄_max = r̄_extrema[region]
                set_lower_bound.(r̄, r̄_min)
                set_upper_bound.(r̄, r̄_max)
                A = region_data[m].A
                b = region_data[m].b
                transform = region_data[m].transform
                z = z_vars[piece, contact, region]

                # Contact position constraints
                @constraint model A * r̄ .<= b
                r_aux = @variable model [1 : 3]
                @constraint model r_aux .== r - transform([r̄; 0])
                if optimizer_does_indicator_constraints
                    for i in eachindex(r_aux)
                        MOI.add_constraint(optimizer, MOI.VectorOfVariables([z; r_aux]), SCIP.IndicatorSet(setindex!(zeros(3), +1, i), 0.0))
                        MOI.add_constraint(optimizer, MOI.VectorOfVariables([z; r_aux]), SCIP.IndicatorSet(setindex!(zeros(3), -1, i), 0.0))
                    end
                else
                    @constraint model  r_aux .<= Mr * (1 - z)
                    @constraint model -r_aux .<= Mr * (1 - z)
                end

                # CoP constraints
                p̄_min, p̄_max = p̄_extrema[region]
                for l = 1 : p_num_coeffs
                    coeff = p_coeffs(l)
                    p = p_vars[piece, contact, coeff]
                    set_lower_bound.(p, p_min)
                    set_upper_bound.(p, p_max)
                    p̄ = p̄_vars[piece, contact, coeff, region]
                    set_lower_bound.(p̄, p̄_min)
                    set_upper_bound.(p̄, p̄_max)
                    # @constraint model A * p̄ .<= b
                    p_aux = @variable model [1 : 3]
                    @constraint model p_aux .== p - transform([p̄; 0])
                    if optimizer_does_indicator_constraints
                        for i in eachindex(p_aux)
                            MOI.add_constraint(optimizer, MOI.VectorOfVariables([z; p_aux]), SCIP.IndicatorSet(setindex!(zeros(3), +1, i), 0.0))
                            MOI.add_constraint(optimizer, MOI.VectorOfVariables([z; p_aux]), SCIP.IndicatorSet(setindex!(zeros(3), -1, i), 0.0))
                        end
                    else
                        @constraint model  p_aux .<= Mr * (1 - z)
                        @constraint model -p_aux .<= Mr * (1 - z)
                    end
                    if max_cop_distance == 0
                        @constraint model p .== r
                    else
                        @constraint model sum(x -> x^2, p - r) <= max_cop_distance^2
                        # constrain_l1_norm(model, r - p, √2 / 2 * max_cop_distance; add_bounds=true)
                    end
                end
            end
        end

        # Contact force/torque constraints
        for j in 1 : num_contacts
            contact = contacts(j)
            μmax = 0.0
            for m in 1 : num_regions
                region = regions(m)
                μ = region_data[m].μ
                μmax = max(μ, μmax)
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
                    set_upper_bound(f̄z, max_force)
                    set_lower_bound.(f̄xy, -μ * max_force)
                    set_upper_bound.(f̄xy, +μ * max_force)
                    # auxiliary variable needed for SecondOrderCone constraint:
                    μf̄z = @variable model lower_bound=0 upper_bound=μ * upper_bound(f̄z)
                    @constraint model μf̄z == μ * f̄z
                    if optimizer_does_indicator_constraints
                        z_opposite = @variable model
                        set_binary(z_opposite)
                        @constraint model z_opposite == 1 - z
                        MOI.add_constraint(optimizer, MOI.VectorOfVariables([z_opposite; f̄z]), SCIP.IndicatorSet([1.0], 0.0))
                    else
                        # @constraint model μf̄z <= max_force * μ * z
                        @constraint model f̄z <= max_force * z
                    end
                    if optimizer_does_soc
                        @constraint model [μf̄z; f̄xy] in SecondOrderCone()
                    else
                        friction_cone_quadratic_constraint = @constraint model f̄xy ⋅ f̄xy <= μf̄z^2
                        push!(friction_cone_quadratic_constraints, friction_cone_quadratic_constraint)
                    end
                end
            end
            for l in 1 : f_num_coeffs
                coeff = f_coeffs(l)
                f = f_vars[piece, contact, coeff]
                set_lower_bound.(f, -max(μmax, 1) * max_force)
                set_upper_bound.(f, +max(μmax, 1) * max_force)
                @constraint model f .== sum(Rs[regions(m)] * f̄_vars[piece, contact, coeff, regions(m)] for m in 1 : num_regions)
            end
        end

        # CoM kinematic constraints
        for j in 1 : num_contacts
            contact = contacts(j)
            r = r_vars[piece, contact]
            for l in 1 : c_num_coeffs
                cpoint = map(x -> x.coeffs[l], c)
                # @constraint model cpoint[3] - r[3] >= 0.75 # TODO
                @constraint model sum(x -> x^2, cpoint - r) >= min_com_to_contact_distance^2
                # set_lower_bound(cpoint[3], 0.7)
                # @constraint model cpoint - r .<= 1.5
                # @constraint model r - cpoint .<= 1.5
                # @constraint model sum(x -> x^2, cpoint - r) >= 0.5^2 # TODO
                @constraint model sum(x -> x^2, cpoint - r) <= max_com_to_contact_distance^2
            end
        end


        # Inter-contact kinematic constraints
        if i > 1
            for j1 in 1 : num_contacts
                contact1 = contacts(j1)
                r1 = r_vars[piece, contact1]
                for j2 in j1 + 1 : num_contacts
                    contact2 = contacts(j2)
                    r2 = r_vars[piece, contact2]
                    @constraint model sum(x -> x^2, (r1 - r2)) >= min_inter_contact_distance^2
                end
            end
            # FIXME: hack:
            @constraint model r_vars[piece, contacts(1), coords(2)] - r_vars[piece, contacts(2), coords(2)] >= -max_cop_distance
        end

        cprev = c
        c′prev = c′
    end

    CentroidalTrajectoryProblem(model,
        c_coeffs, f_coeffs, contacts, regions, pieces, coords, coords2d,
        c_vars, f_vars, f̄_vars, r_vars, p_vars, p̄_vars, Δts, z_vars,
        ns, friction_cone_quadratic_constraints
    )
end

function disallow_jumping!(problem::CentroidalTrajectoryProblem)
    model = problem.model
    pieces = problem.pieces
    for i in 1 : length(pieces)
        @constraint model sum(problem.z_vars[pieces(i)]) >= 1
    end
end

normals(problem) = problem.normals
