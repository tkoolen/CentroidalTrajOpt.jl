struct CentroidalTrajectoryResult
    break_times
    center_of_mass
    contact_forces
    contact_positions
    centers_of_pressure
    # contact_normal_torques
    contact_indicators
end

function CentroidalTrajectoryResult(problem::CentroidalTrajectoryProblem)
    # Axes
    pieces = problem.pieces
    coords = problem.coords
    contacts = problem.contacts
    regions = problem.regions

    # Times
    Δt_vals = val.(problem.Δts)
    break_times = pushfirst!(cumsum(Δt_vals), 0)

    # Center of mass
    c_vals = map(val, problem.c_vars)
    c_subfunctions = [poly_cat(map(k -> poly_piece_val(c_vals, Δt_vals[i], pieces(i), coords(k)), coords.val)) for i in pieces.val]
    c = Piecewise(c_subfunctions, break_times)

    # Contact forces
    f_vals = map(val, problem.f_vars)
    fs = map(contacts.val) do contact
        subfunctions = [poly_cat(map(k -> poly_piece_val(f_vals, Δt_vals[i], pieces(i), coords(k), contacts(contact)), coords.val)) for i in pieces.val]
        Piecewise(subfunctions, break_times)
    end

    # Contact positions
    r_vals = map(val, problem.r_vars)
    rs = map(contacts.val) do contact
        subfunctions = [Constant(map(k -> r_vals[pieces(i), contacts(contact), coords(k)], coords.val)) for i in pieces.val]
        Piecewise(subfunctions, break_times)
    end

    # CoPs
    p_vals = map(val, problem.p_vars)
    ps = map(contacts.val) do contact
        subfunctions = [poly_cat(map(k -> poly_piece_val(p_vals, Δt_vals[i], pieces(i), coords(k), contacts(contact)), coords.val)) for i in pieces.val]
        Piecewise(subfunctions, break_times)
    end

    # # Normal torques
    # τn_vals = map(val, problem.τn_vars)
    # τns = map(contacts.val) do contact
    #     subfunctions = [poly_piece_val(τn_vals, Δt_vals[i], pieces(i), contacts(contact)) for i in pieces.val]
    #     Piecewise(subfunctions, break_times)
    # end

    # Contact indicators
    z_vals = map(val, problem.z_vars)
    zs = map(contacts.val) do contact
        subfunctions = [Constant(z_vals[pieces(i), contacts(contact)] .> 0.5) for i in pieces.val]
        Piecewise(subfunctions, break_times)
    end

    CentroidalTrajectoryResult(break_times, c, fs, rs, ps, #=τns,=# zs)
end

function solve!(problem::CentroidalTrajectoryProblem; ignore_optimize_hook::Bool=problem.model.optimize_hook === nothing)
    optimize!(problem.model, ignore_optimize_hook=ignore_optimize_hook)
    CentroidalTrajectoryResult(problem)
end

function poly_piece_val(bezier_coeffs, Δt, indices...)
    scale_argument(Polynomial(BezierCurve(bezier_coeffs[indices...]...)), 1 / Δt)
end

function poly_cat(polys::StaticVector{N, P}) where {N, M, P<:Polynomial{M}}
    coeffs = ntuple(i -> SVector(ntuple(j -> polys[j].coeffs[i], Val(N))), Val(M))
    Polynomial(coeffs)
end
