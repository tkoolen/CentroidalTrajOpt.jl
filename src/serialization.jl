## Save/load results
module Serialization

export save_result, load_result

import ..CentroidalTrajOpt: CentroidalTrajectoryResult
import Gtk
import JLD2
import FileIO

function save_result(result::CentroidalTrajectoryResult, filename::Union{AbstractString, Nothing}=nothing)
    if filename === nothing
        filename = Gtk.save_dialog_native("Save as...", Gtk.Null(), (Gtk.GtkFileFilter("*.jld2", name="All supported formats"), "*.jld2"))
    end
    if !isempty(filename)
        FileIO.save(filename, Dict("result" => result))
    end
    return filename
end

function load_result(filename::Union{AbstractString, Nothing}=nothing)
    if filename === nothing
        filename = Gtk.open_dialog_native("Open result file", Gtk.GtkNullContainer(), ("*.jld2", Gtk.GtkFileFilter("*.jld2", name="All supported formats")))
    end
    if !isempty(filename)
        return FileIO.load(filename)["result"]
    end
    return nothing
end

end
