## Save/load results
module Serialization

export save_result, load_result

import ..CentroidalTrajOpt: CentroidalTrajectoryResult
import JLD2
import FileIO

function save_result(result::CentroidalTrajectoryResult, filename::Union{AbstractString, Nothing}=nothing)
    if filename === nothing
        default_dir = get(pwd, ENV, "CENTROIDAL_TRAJ_OPT_RESULT_DIR")
        println("Directory [$default_dir]?")
        dir = readline()
        isempty(dir) && (dir = default_dir)
        println("File name?")
        basename = readline()
        filename = joinpath(dir, basename)
    end
    FileIO.save(filename, Dict("result" => result))
    return filename
end

function load_result(filename::AbstractString)
    return FileIO.load(filename)["result"]
end

end
