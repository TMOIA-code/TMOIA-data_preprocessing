# Simple Model Interpretability
module InterpretModel

using BSON, Flux, Random, Tables

if !isdefined(@__MODULE__, :MyUtils)
    include("utils.jl")
    using .MyUtils
end

export interpreter_inDir, interpreter

function interpreter(path_model::String, bias::Bool=false)
    BSON.@load path_model model
    path_wt = string(dirname(path_model), "/", "interp_weight_", basename(path_model)[1:(end-5)], ".txt")
    my_write_table(model.layers[1].weight, path_wt, toTable=true)
    if bias
        path_bias = string(dirname(path_model), "/", "interp_bias_", basename(path_model)[1:(end-5)], ".txt")
        my_write_table(model.layers[1].bias, path_bias, toTable=true)
    end
    return nothing
end

function interpreter_inDir(dir_models::String, bias::Bool=true)
    grepMd = ".bson" |> Regex
    fList = readdir(dir_models)
    mds = fList[occursin.(grepMd, fList)]
    println(mds)
    Threads.@threads for mdn in eachindex(mds)
        interpreter(string(dir_models, "/", mds[mdn]), bias)
    end
    return nothing
end

end # of InterpretModel
