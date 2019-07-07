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
    p_vars
    r_vars
    r̄_vars
    # τn_vars
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
end
end
using .ObjectiveTypes

function CentroidalTrajectoryProblem(optimizer_factory::JuMP.OptimizerFactory,
        region_data::AbstractVector{<:ContactRegion},
        c0::AbstractVector{<:Number},
        ċ0::AbstractVector{<:Number},
        contacts0::AbstractVector{<:Union{Pair{<:ContactRegion, <:AbstractVector}, Nothing}};
        c_degree = 3,
        num_pieces = 2,
        g = SVector(0.0, 0.0, -9.81),
        max_cop_distance = 0.1,
        objective_type::ObjectiveType = FEASIBILITY
    )

    c_num_coeffs = c_degree + 1
    f_num_coeffs = c_num_coeffs - 2
    r_num_coeffs = c_num_coeffs
    num_regions = length(region_data)
    num_contacts = length(contacts0)

    # Axes
    c_coeffs = Axis{:c_coeff}(1 : c_num_coeffs)
    f_coeffs = Axis{:f_coeff}(1 : f_num_coeffs)
    r_coeffs = Axis{:r_coeff}(1 : r_num_coeffs)
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

    # Indexing convention:
    # i: piece index
    # j: contact index
    # k: coordinate index
    # l: coefficient index
    # m: region index

    # Δts = AxisArray(fill(0.3, length(pieces)), pieces)
    Δts = nothing
    if Δts === nothing
        Δts = axis_array_vars(model, i -> "Δt[$i]", pieces)
        Δtsqs = axis_array_vars(model, i -> "Δtsq[$i]", pieces)
        for (Δt, Δtsq) in zip(Δts, Δtsqs)
            Δtmin = 0.3
            Δtmax = 1.0 # TODO
            set_lower_bound(Δt, Δtmin)
            set_upper_bound(Δt, Δtmax)
            set_lower_bound(Δtsq, Δtmin^2)
            set_upper_bound(Δtsq, Δtmax^2)
            @constraint model Δtsq == Δt^2
        end
        if objective_type == ObjectiveTypes.MIN_TIME
            @objective model Min sum(Δts)
        end
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
    c_vars = axis_array_vars(model, (i, k, l) -> "C[p$i, $k, $l]",pieces, coords, c_coeffs;
        lower_bound=-2, upper_bound=2)
    f_vars = axis_array_vars(model, (i, j, k, l) -> "F[p$i, c$j, $k, $l]", pieces, contacts, coords, f_coeffs;
        lower_bound=-3 * norm(g), upper_bound=3 * norm(g))
    f̄_vars = axis_array_vars(model, (i, j, k, l, m) -> "F̄[p$i, c$j, $k, $l, r$m]", pieces, contacts, coords, f_coeffs, regions;
        lower_bound=-3 * norm(g), upper_bound=3 * norm(g))
    p_vars = axis_array_vars(model, (i, j, k) -> "P[p$i, c$j, $k]", pieces, contacts, coords;
        lower_bound=-2, upper_bound=2)
    p̄_vars = axis_array_vars(model, (i, j, k, m) -> "P̄[p$i, c$j, $k, r$m]", pieces, contacts, coords2d, regions;
        lower_bound=-1, upper_bound=1)
    r_vars = axis_array_vars(model, (i, j, k, l) -> "R[p$i, c$j, $k, $l]", pieces, contacts, coords, r_coeffs;
        lower_bound=-2, upper_bound=2)
    r̄_vars = axis_array_vars(model, (i, j, k, l, m) -> "R̄[p$i, c$j, $k, $l, r$m]", pieces, contacts, coords2d, r_coeffs, regions;
        lower_bound=-1, upper_bound=1)
    # τn_vars = axis_array_vars(model, (i, j, l) -> "Tn[p$i, c$j, $l]", pieces, contacts, f_coeffs;
    #     lower_bound=-1 * norm(g), upper_bound=1 * norm(g))

    #@variable model objval
    if objective_type == ObjectiveTypes.MIN_EXCURSION
        # @constraint model objval >= sum((c_vars[pieces(i), c_coeffs(l)] - c0) ⋅ (c_vars[pieces(i), c_coeffs(l)] - c0) for i in 1 : num_pieces, l in 1 : c_num_coeffs)
        @objective model Min sum(x -> x ⋅ x, (c_vars[pieces(i), c_coeffs(l)] - c0) for i in 1 : num_pieces, l in 1 : c_num_coeffs)
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

        # Contact/region assignment constraints
        for j in 1 : num_contacts
            # Each contact can be assigned to at most one region
            contact = contacts(j)

            if optimizer_does_sos1
                @constraint model z_vars[piece, contact] in MOI.SOS1(collect(Float64, 1 : num_regions))
            else
                @constraint model sum(z_vars[piece, contact]) <= 1
            end
            if i == 1
                # Initial contact assignment.
                if contacts0[j] === nothing
                    @constraint model sum(z_vars[piece, contact]) == 0
                else
                    region_data0, p0 = contacts0[j]
                    m = findfirst(isequal(region_data0), region_data)
                    z_var = z_vars[piece, contact, regions(m)]
                    fix(z_var, 1, force=true)
                    unset_binary(z_var) # to make Alpine happpy
                    fix.(p_vars[piece, contact], p0, force=true)
                end
            end
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
                Δpmax = 0.8 # TODO
                # @constraint model sum(x -> x^2, Δp) <= (1 - sum(z_vars[piece, contact])) * Δpmax^2
                @constraint model  Δp .<= (1 - sum(z_vars[piece, contact])) * Δpmax
                @constraint model -Δp .<= (1 - sum(z_vars[piece, contact])) * Δpmax
            end
        end

        # Initial conditions
        if i == 1
            @constraint model map(x -> x(0), c) .== c0
            @constraint model map(x -> x(0), c′) .== ċ0 .* Δts[i]
        end

        # Final conditions
        if i == num_pieces
            # @constraint model map(x -> x(1), c) .== c0 # TODO
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
        # τns = [BezierCurve(τn_vars[piece, contacts(j)]...) for j in contacts.val]

        # Dynamics
        ftot = sum(fs)
        constrain_poly_equal.(model, c′′, identity((g + ftot) .* Δtsqs[i]))

        # Torque balance
        # τ = sum(Polynomial.(rs[j]) × Polynomial.(fs[j]) for j = 1 : num_contacts) +
        #     sum(Polynomial.(sum(ns[regions(m)] .* z_vars[piece, regions(m), contacts(j)] for m in 1 : num_regions) .* τns[j]) for j = 1 : num_contacts) -
        #     Polynomial.(c) × Polynomial.(ftot)
        τ = sum(Polynomial.(rs[j]) × Polynomial.(fs[j]) for j = 1 : num_contacts) -
            Polynomial.(c) × Polynomial.(ftot)
        τ = map(x -> Polynomial(simplify.(model, x.coeffs)), τ)
        constrain_poly_equal.(model, τ, 0)

        # CoP and contact position constraints
        Mr = 2 # TODO
        for j in 1 : num_contacts
            contact = contacts(j)
            p = p_vars[piece, contact]
            for m in 1 : num_regions
                region = regions(m)
                p̄ = p̄_vars[piece, contact, region]
                A = region_data[m].A
                b = region_data[m].b
                transform = region_data[m].transform
                z = z_vars[piece, contact, region]

                # Contact position constraints
                @constraint model A * p̄ .<= b
                @constraint model  (p - transform([p̄; 0])) .<= Mr * (1 - z)
                @constraint model -(p - transform([p̄; 0])) .<= Mr * (1 - z)

                # CoP constraints
                for l = 1 : r_num_coeffs
                    coeff = r_coeffs(l)
                    r = r_vars[piece, contact, coeff]
                    r̄ = r̄_vars[piece, contact, coeff, region]
                    @constraint model A * r̄ .<= b
                    @constraint model  (r - transform([r̄; 0])) .<= Mr * (1 - z)
                    @constraint model -(r - transform([r̄; 0])) .<= Mr * (1 - z)
                    if max_cop_distance == 0
                        @constraint model r .== p
                    else
                        @constraint model sum(x -> x^2, r - p) <= max_cop_distance^2
                        # @constraint model  (r - p) .<= max_cop_distance
                        # @constraint model -(r - p) .<= max_cop_distance
                    end
                end
            end
        end

        # Contact force/torque constraints
        Mf = 3 * norm(g) # TODO
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
                    # auxiliary variable needed for SecondOrderCone constraint:
                    μf̄z = @variable model lower_bound=0 upper_bound=μ * upper_bound(f̄z)
                    @constraint model μf̄z == μ * f̄z
                    @constraint model μf̄z <= Mf * μ * z
                    if optimizer_does_soc
                        @constraint model [μf̄z; f̄xy] in SecondOrderCone()
                    else
                        friction_cone_quadratic_constraint = @constraint model f̄xy ⋅ f̄xy <= μf̄z^2
                        push!(friction_cone_quadratic_constraints, friction_cone_quadratic_constraint)
                    end
                    # TODO: use a SOS condition?
                    # TODO: crude model
                    # τn = τn_vars[piece, contact, coeff]
                    # if μrot == 0
                    #     fix(τn, 0, force=true)
                    # else
                    #     @constraint model  τn <= μrot * f̄z + Mτ * z
                    #     @constraint model -τn <= μrot * f̄z + Mτ * z
                    # end
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
                cpoint = map(x -> x.coeffs[l], c)
                @constraint model cpoint[3] - p[3] >= 0.6 # TODO
                # set_lower_bound(cpoint[3], 0.7)
                # @constraint model cpoint - p .<= 1.5
                # @constraint model p - cpoint .<= 1.5
                # @constraint model sum(x -> x^2, cpoint - p) >= 0.5^2 # TODO
                @constraint model sum(x -> x^2, cpoint - p) <= 1.2^2 # TODO
            end
        end

        # Inter-contact kinematic constraints
        for j1 in 1 : num_contacts
            contact1 = contacts(j1)
            p1 = p_vars[piece, contact1]
            for j2 in j1 + 1 : num_contacts
                contact2 = contacts(j2)
                p2 = p_vars[piece, contact2]
                @constraint model sum(x -> x^2, (p1 - p2)) >= 0.1^2 # TODO: parameterize
            end
        end

        cprev = c
        c′prev = c′
    end

    CentroidalTrajectoryProblem(model,
        c_coeffs, f_coeffs, contacts, regions, pieces, coords, coords2d,
        c_vars, f_vars, f̄_vars, p_vars, r_vars, r̄_vars, #=τn_vars,=# Δts, z_vars,
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
