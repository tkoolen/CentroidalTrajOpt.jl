using PackageCompiler

scripts_dir = joinpath(@__DIR__, "..", "scripts")
toml_file = joinpath(scripts_dir, "Project.toml")
script = joinpath(scripts_dir, "exploration.jl")
precompiles_file = joinpath(scripts_dir, "precompiles.jl")
blacklist = ["CentroidalTrajOpt", "QPWalkingControl", "QPControl"]#, "AtlasRobot"]
sysimg_dest = joinpath(scripts_dir, "sysimg.so")

@info "Generating precompile statements"
PackageCompiler.snoop(nothing, toml_file, script, precompiles_file, false, blacklist)
@info "Done generating precompile statements"

sleep(1)
@info "Compiling system image"
new_sysimg, old_sysimg = compile_incremental(toml_file, precompiles_file)
@info "Done compiling system image"
cp(new_sysimg, sysimg_dest, force=true)
