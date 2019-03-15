struct CentroidalTrajectoryResult
    break_times
    center_of_mass
    contact_forces
    contact_positions
    centers_of_pressure
    contact_normal_torques
end

function CentroidalTrajectoryResult(problem::CentroidalTrajectoryProblem)
    # Axes
    pieces = problem.pieces
    coords = problem.coords
    contacts = problem.contacts

    # Times
    Δt_vals = val.(problem.Δts)
    break_times = pushfirst!(cumsum(Δt_vals), 0)

    # Center of mass
    c_vals = map(val, problem.c_vars)
    c_subfunctions = [Vectorized([poly_piece_val(c_vals, Δt_vals[i], pieces(i), coords(k)) for k in coords.val]) for i in pieces.val]
    c = Piecewise(c_subfunctions, break_times)

    # Contact forces
    f_vals = map(val, problem.f_vars)
    fs = map(contacts.val) do contact
        subfunctions = [Vectorized([poly_piece_val(f_vals, Δt_vals[i], pieces(i), coords(k), contacts(contact)) for k in coords.val]) for i in pieces.val]
        Piecewise(subfunctions, break_times)
    end

    # Contact positions
    p_vals = map(val, problem.p_vars)
    ps = map(contacts.val) do contact
        subfunctions = [Constant([p_vals[pieces(i), contacts(contact), coords(k)] for k in coords.val]) for i in pieces.val]
        Piecewise(subfunctions, break_times)
    end

    # CoPs
    r_vals = map(val, problem.r_vars)
    rs = map(contacts.val) do contact
        subfunctions = [Vectorized([poly_piece_val(r_vals, Δt_vals[i], pieces(i), coords(k), contacts(contact)) for k in coords.val]) for i in pieces.val]
        Piecewise(subfunctions, break_times)
    end

    # Normal torques
    τn_vals = map(val, problem.τn_vars)
    τns = map(contacts.val) do contact
        subfunctions = [poly_piece_val(τn_vals, Δt_vals[i], pieces(i), contacts(contact)) for i in pieces.val]
        Piecewise(subfunctions, break_times)
    end

    CentroidalTrajectoryResult(break_times, c, fs, ps, rs, τns)
end

function solve!(problem::CentroidalTrajectoryProblem)
    optimize!(problem.model)
    CentroidalTrajectoryResult(problem)
end

function poly_piece_val(bezier_coeffs, Δt, indices...)
    scale_argument(Polynomial(BezierCurve(bezier_coeffs[indices...]...)), 1 / Δt)
end




# function center_of_mass(problem::CentroidalTrajectoryProblem)
#     pieces = problem.pieces
#     coords = problem.coords
#     Δt_vals = val.(problem.Δts)
#     break_times = pushfirst!(cumsum(Δt_vals), 0)

#     Piecewise(subfunctions, break_times)
# end

# function contact_trajectory(problem::CentroidalTrajectoryProblem, bezier_coeffs, contact_index)
#     pieces = problem.pieces
#     coords = problem.coords
#     contacts = problem.contacts
#     Δt_vals = val.(problem.Δts)
#     break_times = pushfirst!(cumsum(Δt_vals), 0)
#     vals = map(val, bezier_coeffs)
#     subfunctions = [Vectorized([poly_piece_val(vals, Δt_vals[i], pieces(i), coords(k), contacts(contact_index)) for k in coords.val]) for i in pieces.val]
#     Piecewise(subfunctions, break_times)
# end

# function contact_force(problem::CentroidalTrajectoryProblem, contact_index)

# end

# function contact_forces(problem::CentroidalTrajectoryProblem)
#     AxisArray([contact_force(problem, i) for i in problem.contacts.val], problem.contacts)
# end
