# Dimension reduction
using Flux, CUDA, BSON, Dates
__precompile__()
include("utils.jl")
include("utils_dl.jl")
include("utils_dr.jl")
include("models.jl")


function DimReduction(whichSp::String = "01", whX::String = "SNP", whY::String = "_pheno_zscore",
                        batch_size::Int64 = 32, epoch_max::Int64 = 550, lr::Float64=0.0001, es_delay::Int64 = 70,
                        dirDataSplit::String = "path to r10", calc_whLayer::Int64 = 1;
                        path_Gene2SNP::String = "none", isSimpleMd::Bool=true, subDir::String="", md_dimL2::Int64=2000, isZipModel::Bool=true,
                        bias::Bool=false, selectY::Bool=false, ySelected::Vector{Int64}=[0, 0], stdX::Bool=false, stdY::Bool=false, stdFunc=std_zscore!)
    #
    if subDir == ""; subDir = string("dr_", md_dimL2, "/"); end
    w_pref = string(dirDataSplit, "/", subDir, whichSp, "_")
    mkpath(dirname(w_pref))
    path_wrt_r = string(w_pref, "drResult_", whX, ".txt")
    path_save_model = string(w_pref, "drModel_", whX, ".bson")
    write_csvRows(path_wrt_r, [Dates.now(Dates.UTC), dirDataSplit, whichSp, batch_size, lr], false)
    ##
    mx_trn, mx_val, mx_tst = get_loader(whichSp, dirDataSplit, true, batch_size, whX, whY, selectY, y_selected=ySelected, std_x=stdX , std_y=stdY, std_func=stdFunc)
    num_feat = size(mx_tst[1])[1]
    num_trait = size(mx_tst[2])[1]
    ##
    isFreezeMdParam = false
    where2Freeze = zeros(Float32, 2, 2)
    if whX == "SNP"
        isFreezeMdParam = true
        model = md_snpS(path_Gene2SNP, num_trait, isSimpleMd)
        where2Freeze = whGradsToDel(path_Gene2SNP)
    else
        model = md_fc(num_feat, num_trait, md_dimL2, isSimpleMd, bias)
    end
    println(model)
    println("\n")
    #
    trn_loader, val_loader, tst_loader = get_loader(whichSp, dirDataSplit, false, batch_size, whX, whY, selectY, y_selected=ySelected, std_x=stdX , std_y=stdY, std_func=stdFunc)
    model = myTrain(model, trn_loader, val_loader, tst_loader, epoch_max, lr, es_delay, isFreezeMdParam, true, true,
                    path_w_rec=path_wrt_r, path_save_model=path_save_model, isZipModel=isZipModel, wh2freeze=where2Freeze)
    trn_loader, val_loader, tst_loader, where2Freeze = nothing, nothing, nothing, nothing, nothing
    #
    io = open(path_wrt_r, "a")
    write(io, string("\n", Dates.now(Dates.UTC), "\n\n", model, "\n"))
    close(io)
    @views begin
        g_trn = calc_RD(model, mx_trn[1], calc_whLayer, bias) |> transpose |> Tables.table ; CSV.write(string(w_pref, "trn_dr_", whX, ".txt"), g_trn, delim = "\t", header = false) ; g_trn = nothing;
        g_val = calc_RD(model, mx_val[1], calc_whLayer, bias) |> transpose |> Tables.table ; CSV.write(string(w_pref, "val_dr_", whX, ".txt"), g_val, delim = "\t", header = false) ; g_val = nothing;
        g_tst = calc_RD(model, mx_tst[1], calc_whLayer, bias) |> transpose |> Tables.table ; CSV.write(string(w_pref, "tst_dr_", whX, ".txt"), g_tst, delim = "\t", header = false) ; g_tst = nothing;
    end
    return nothing
end


function run_DR(whichX::Vector{String}, whichY::String, traitsIn::Vector{String}, nameSplits::Vector{String}, dirDataset::String,
                batch_size::Int64, epoch_max::Int64, lr::Float64, patienceX::Int64, whLayerRD::Int64, MdDimLayer2::Int64=2000;
                pathG2S::String="none", suffixTrait::String="_AllOmics", subDir::String="/r10", isSimpleMd::Bool=true, isZipModel::Bool=true,
                bias::Bool=true, selectY::Bool=false, ySelected::Vector{Int64}=[0, 0], stdX::Bool=false, stdY::Bool=false, stdFunc=std_zscore!)
    for wx in eachindex(whichX)
        for tnx in eachindex(traitsIn)
            for wp in nameSplits
                println("\n", "---- Omic: ", whichX[wx], "\n", "---- Trait: ", traitsIn[tnx], "\n", "---- Random split: ", wp, "\n")
                DimReduction(wp, whichX[wx], whichY, batch_size, epoch_max, lr, patienceX,
                             string(dirDataset, traitsIn[tnx], suffixTrait, subDir),
                             whLayerRD, path_Gene2SNP=pathG2S, md_dimL2=MdDimLayer2, isSimpleMd=isSimpleMd, isZipModel=isZipModel, bias=bias,
                             selectY=selectY, ySelected=ySelected, stdX=stdX, stdY=stdY, stdFunc=stdFunc)
            end
        end
    end
    return nothing
end
