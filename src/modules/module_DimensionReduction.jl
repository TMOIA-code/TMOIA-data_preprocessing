module DimensionReduction

using Flux, CUDA, BSON, Dates, Tables
include("utils.jl")
include("utils_dl.jl")
include("models.jl")
using .MyUtils

export pick_SNPsInGene
export DimReduction
export run_DR

## Pick SNPs with gene found
function pick_SNPsInGene_ut(L_sta::Int64, L_end::Int64,
                            FileSNP::String, snp_remain::Vector{Int64}, IdPos::Vector{Int64}, w_name::String)
    #write 1st row
    @views tmpLine = my_readline(FileSNP, snp_remain[L_sta], delim=',')[Not([1,2,3,4])][IdPos]
    my_write_table(tmpLine, w_name, toTable=true)
    @views for lnN in (L_sta + 1):L_end
        tmpLine = my_readline(FileSNP, snp_remain[lnN], delim=',')[Not([1,2,3,4])][IdPos]
        my_write_table(tmpLine, w_name, isAppend=true, toTable=true)
    end
    return nothing
end

function pick_SNPsInGene(FileGeneToSNP::String, FileSNP::String, FileIdWhLine::String)
    IdPos = parse.(Int64, my_readline(FileIdWhLine, 1, delim='\t'))
    g2s = my_read_table(FileGeneToSNP, String, '\t')
    @views snp_remain = parse.(Int64, unique(g2s[:, 1])) |> sort
    println("\n--- SNPs with gene found: ", length(snp_remain), "\n")
    ## Multi-thread
    threadNum = Threads.nthreads()
    startPts, endPts = multiThreadSplit(length(snp_remain), threadNum)
    ## MT begins
    Threads.@threads for nT in 1:threadNum
        pick_SNPsInGene_ut(startPts[nT], endPts[nT], FileSNP, snp_remain, IdPos,
                           string(dirname(FileSNP), "/", ".tmp_nt", format_numLen(nT, 2), "_", basename(FileSNP), ".txt"))
    end
    ## Conc parts
    conc_fileParts(string(dirname(FileSNP), "/.tmp_nt*"),
                   string(dirname(FileSNP), "/", "_lastest_", basename(FileSNP), ".txt"))
    return nothing
end


## Count GeneToSNP
function count_gene2snp(fpath::String, returnWhLine::Bool = false)
    fin = my_read_table(fpath, String, '\t')
    ## Must sort SNP by which line because dataset is sorted. !!!!!
    @views orderL = sortperm(parse.(Int64, fin[:,1]))
    @views fin = fin[orderL, :]
    ##
    @views gene_uniq = unique(fin[:,2])
    numsSNP = [] |> Vector{Int64}
    positions = []
    @views for gn in eachindex(gene_uniq)
        pos = findall(x -> x == gene_uniq[gn], fin[:,2])
        push!(numsSNP, length(pos))
        push!(positions, pos)
    end
    if returnWhLine; return positions, size(fin)[1]; end;
    return numsSNP
end

function whGradsToDel(pathGeneToSNP::String)
    geneL, num_snp = count_gene2snp(pathGeneToSNP, true)
    wh2Keep = zeros(Float32, length(geneL), num_snp)
    #seq_1toN = collect(StepRange(1, Int64(1), num_snp))
    Threads.@threads for gn in eachindex(geneL)
        wh2Keep[gn, geneL[gn]] .= 1.0
    end
    return wh2Keep
end

# =======================================

function DimReduction(whichSp::String = "01", whX::String = "SNP", whY::String = "_pheno_zscore",
                        batch_size::Int64 = 32, paramTrain::myStructParamTrain = myStructParamTrain(550, 1e-4, 70),
                        dirDataSplit::String = "path to r10", calc_whLayer::Int64 = 1;
                        path_Gene2SNP::String = "none", isSimpleMd::Bool=true, subDir::String="", md_dimL2::Int64=2000, isZipModel::Bool=true,
                        bias::Bool=false, selectY::Bool=false, ySelected::Vector{Int64}=[0, 0], stdX::Bool=false, stdY::Bool=false, stdFunc=std_zscore!)
    #
    if subDir == ""; subDir = string("dr_", md_dimL2, "/"); end
    w_pref = string(dirDataSplit, "/", subDir, whichSp, "_")
    mkpath(dirname(w_pref))
    path_wrt_r = string(w_pref, "drResult_", whX, ".txt")
    path_save_model = string(w_pref, "drModel_", whX, ".bson")
    marks = [Dates.now(Dates.UTC), dirDataSplit, whichSp, batch_size, lr]
    my_write_table(marks, path_wrt_r, toTable=true)
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
    dataIn = myStructDataIn(trn_loader, tst_loader, val_loader)
    trn_loader, val_loader, tst_loader = nothing, nothing, nothing
    model = myTrain(model, dataIn, paramTrain, isFreezeMdParam, true, true,
                    path_w_rec=path_wrt_r, path_save_model=path_save_model, isZipModel=isZipModel, wh2freeze=where2Freeze)
    where2Freeze = nothing
    #
    io = open(path_wrt_r, "a")
    write(io, string("\n", Dates.now(Dates.UTC), "\n\n", model, "\n"))
    close(io)
    @views begin
        g_trn = calc_RD(model, mx_trn[1], calc_whLayer, bias) |> transpose |> Tables.table ; my_write_table(g_trn, string(w_pref, "trn_dr_", whX, ".txt")) ; g_trn = nothing;
        g_val = calc_RD(model, mx_val[1], calc_whLayer, bias) |> transpose |> Tables.table ; my_write_table(g_val, string(w_pref, "val_dr_", whX, ".txt")) ; g_val = nothing;
        g_tst = calc_RD(model, mx_tst[1], calc_whLayer, bias) |> transpose |> Tables.table ; my_write_table(g_tst, string(w_pref, "tst_dr_", whX, ".txt")) ; g_tst = nothing;
    end
    return nothing
end


function run_DR(whichX::Vector{String}, whichY::String, traitsIn::Vector{String}, nameSplits::Vector{String}, dirDataset::String,
                batch_size::Int64, paramTrain::myStructParamTrain, whLayerRD::Int64, MdDimLayer2::Int64=2000;
                pathG2S::String="none", suffixTrait::String="_AllOmics", subDir::String="/r10", isSimpleMd::Bool=true, isZipModel::Bool=true,
                bias::Bool=true, selectY::Bool=false, ySelected::Vector{Int64}=[0, 0], stdX::Bool=false, stdY::Bool=false, stdFunc=std_zscore!)
    for wx in eachindex(whichX)
        for tnx in eachindex(traitsIn)
            for wp in nameSplits
                println("\n", "---- Omic: ", whichX[wx], "\n", "---- Trait: ", traitsIn[tnx], "\n", "---- Random split: ", wp, "\n")
                DimReduction(wp, whichX[wx], whichY, batch_size, paramTrain,
                             string(dirDataset, traitsIn[tnx], suffixTrait, subDir),
                             whLayerRD, path_Gene2SNP=pathG2S, md_dimL2=MdDimLayer2, isSimpleMd=isSimpleMd, isZipModel=isZipModel, bias=bias,
                             selectY=selectY, ySelected=ySelected, stdX=stdX, stdY=stdY, stdFunc=stdFunc)
            end
        end
    end
    return nothing
end

end # of DimensionReduction
