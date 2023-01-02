## Filter
__precompile__()
include("utils.jl")
include("utils_searchSNP.jl")


## Filter samples: with 3 omics? in tags? (line2 removed) (sorted by id_Xomics)
function filter_cols(origFilePaths::Vector{String}, tags::Vector{String}, grepLen::Int64=16, delim::String="\t")
    id_Xomics = sampl_own_Xomics(origFilePaths, grepLen, delim, false)
    Threads.@threads for nf in eachindex(origFilePaths)
        filter_cols_byfile(origFilePaths[nf], tags, id_Xomics, delim)
    end
    return nothing
end


## Filter features by threshold NA and threshold Var, paralleled.
function filter_rows(filePath::String, threshold_var::Float64=0.001, threshold_NA::Float64=0.25; delim::String="\t", headersN::Int64=2)## remove 2-row headers
    ##
    fileDir = string(dirname(filePath), "/")
    fileName = basename(filePath)
    # Slice for parallel
    num_row = my_count_lines(filePath) - headersN
    @show num_row
    threadNum = Threads.nthreads()
    startPts, endPts = multiThreadSplit(num_row, threadNum)
    ## Calc vars, multi-thread
    Threads.@threads for nT in 1:threadNum
        calc_var_A2Blines(filePath,
                          string(fileDir, ".tmp_running_THvar-", threshold_var, "_THna-", threshold_NA, "_p-", format_numLen(nT, 2), "_", fileName),
                          (headersN + startPts[nT]), (headersN + endPts[nT]),  threshold_var, threshold_NA, delim)
    end
    ### Write header
    write_csvRows(string(fileDir, ".tmp_0000_header_", fileName, "_THvar-", threshold_var, "_THna-", threshold_NA, ".txt"), my_readline(filePath, 1, delim=delim), false)
    ## Conc parts
    path_w = string(fileDir, "_f-X_", fileName, "_THvar-", threshold_var, "_THna-", threshold_NA, ".txt")
    grep_w = string(fileDir, ".tmp_*")
    conc_fileParts(grep_w, path_w)
    return nothing
end


##==================================##

## Filter SNPs by threshold Var, paralleled.
function filter_rows_SNP(filePath::String, threshold_var::Float64=0.01, delim_r::String=",")
    ##
    fileDir = dirname(filePath) * "/"
    fileName = basename(filePath)
    # Slice for parallel
    num_row = my_count_lines(filePath)
    @show num_row
    threadNum = Threads.nthreads()
    startPts, endPts = multiThreadSplit(num_row, threadNum)
    # 
    @show threadNum
    ## Calc vars, multi-thread
    Threads.@threads for nT in 1:threadNum
        calc_var_A2Blines_s(filePath, *(fileDir, ".tmp_running_THvar-", threshold_var, "_p-", format_numLen(nT, 2), "_", fileName[1:end-4], ".txt"),
                            startPts[nT], endPts[nT],  threshold_var, delim_r=delim_r)
    end
    conc_fileParts(*(fileDir, ".tmp_*"), *(fileDir, "_", fileName[1:end-4], "_THvar-", threshold_var, ".txt"))
    return nothing
end


function search_gene_of_SNP(pathSNPmx::String, dirGTF::String, intergenic::Bool=true; skipRow::Int64=0, slimResult::Bool=true)
    ##sh TAIR10.1/split_chr.sh
    ## Count lines for GTF of each Chr
    GTFs = readdir(dirGTF)
    GTFs = GTFs[occursin.(Regex(".gtf"), GTFs)]
    nRowGTFs = Vector{Int64}(undef, length(GTFs))
    Threads.@threads for gn in eachindex(GTFs)
        nRowGTFs[gn] = my_count_lines(dirGTF * "/" * GTFs[gn])
    end
    ## Run
    num_thread = Threads.nthreads()
    num_SNP = my_count_lines(pathSNPmx) - skipRow
    staL, endL = multiThreadSplit(num_SNP, num_thread)
    if slimResult
        Threads.@threads for nT in 1:num_thread
            search_gene_by_snp(pathSNPmx, dirGTF, nRowGTFs, (skipRow + staL[nT]), (skipRow + endL[nT]), nT, intergenic, skipRow)
        end
    else
        Threads.@threads for nT in 1:num_thread
            search_region_by_snp(pathSNPmx, dirGTF, nRowGTFs, (skipRow + staL[nT]), (skipRow + endL[nT]), nT, intergenic, skipRow)
        end
    end
    ## Combine file parts
    dirW = dirname(pathSNPmx) * "/"
    if intergenic
        pathGTFout = *(dirW, "_gtf-", "withIntergenic_", basename(pathSNPmx)[1:end-4], "_", basename(dirGTF), ".gtf")
        pathGeneOut = *(dirW, "_gene-", "withIntergenic_", basename(pathSNPmx)[1:end-4], "_", basename(dirGTF), ".txt")
    else
        pathGTFout = *(dirW, "_gtf-", "onlyGene_", basename(pathSNPmx)[1:end-4], "_", basename(dirGTF), ".gtf")
        pathGeneOut = *(dirW, "_gene-", "onlyGene_", basename(pathSNPmx)[1:end-4], "_", basename(dirGTF), ".txt")
    end
    #@show pathGeneOut pathGTFout
    grep_gtf = dirW * ".p-gtf-*"
    grep_gene = dirW * ".p-gene-*"
    if !slimResult
        conc_fileParts(grep_gtf, pathGTFout)
    else
        conc_fileParts(grep_gene, pathGeneOut)
    end
    return nothing
end


################################# Filter other omics #############################

function filter_ath(fPath::String, thresholdVar::Float64, getNumFeat::Bool=false, nheader::Int64=1)
    colToDel = Int64[]
    mx_df = CSV.File(fPath, header=nheader) |> DataFrame
    mxIn = mx_df |> Matrix
    dimx = size(mxIn)
    for nc in 1:dimx[2]
        if var(mxIn[:, nc]) < thresholdVar
            push!(colToDel, nc)
        end
    end
    num_feat = dimx[2] - length(colToDel)
    println("---- num feature: ", num_feat)
    if getNumFeat
        return num_feat
    else
        return mx_df[:, Not(colToDel)]
    end
end

# Get proper number of features by Binary search
function control_num_feat(fPath::String, numFeature::Int64=1000, isReturn::Bool=false)
    guess0 = 1.0
    guess1 = 0.5
    dynam = filter_ath(fPath, guess0, true)
    while dynam != numFeature
        if dynam < numFeature
            guess0 = guess0 - guess1
        elseif dynam > numFeature
            guess0 = guess0 + guess1
        end
        dynam = filter_ath(fPath, guess0, true)
        @show guess0
        guess1 = guess1 / 2
    end
    mxOut = filter_ath(fPath, guess0)
    newF = string(dirname(fPath), "/flt/", "_flt-", numFeature, "_", basename(fPath))
    CSV.write(newF, mxOut)
    if isReturn
        return mxOut, guess0
    end
    return nothing
end
