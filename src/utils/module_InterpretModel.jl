# Simple Model Interpretability
module InterpretModel

using BSON, Flux, Random, Tables, CSV

export interpreter_inDir, interpreter

function interpreter(path_model::String, bias::Bool=false)
    BSON.@load path_model model
    path_wt = string(dirname(path_model), "/", "interp_weight_", basename(path_model)[1:(end-5)], ".txt")
    CSV.write(path_wt, Tables.table(model.layers[1].weight), delim = "\t", header = false)
    if bias
        path_bias = string(dirname(path_model), "/", "interp_bias_", basename(path_model)[1:(end-5)], ".txt")
        CSV.write(path_bias, Tables.table(model.layers[1].bias), delim = "\t", header = false)
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
