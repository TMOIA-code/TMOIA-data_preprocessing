module MyFilters

include("utils.jl")

export filter_cols, filter_rows
export filter_rows_SNP, filter_ath
export control_num_feat


## Q: Whose var?
function Var_minmax_normalize(vecin::Vector, vec_noNA::Vector{Float64}, containNA::Bool=false, returnVar::Bool=true)
    vecin = string.(vecin)
    vecin = parse.(Float64, vecin)
    if containNA
        meann = mean(vec_noNA)
        minn = minimum(vec_noNA, init=Inf)
        rangein = maximum(vec_noNA, init=-1000) - minn
        patchn = (meann-minn)/rangein
    else
        minn = minimum(vecin)
        rangein = maximum(vecin) - minn
    end
    @inbounds @simd for num_g in eachindex(vecin)
        if containNA && vecin[num_g] == -1
            vecin[num_g] = patchn
        else
            vecin[num_g] = (vecin[num_g] - minn) / rangein
        end
    end
    if returnVar
        return var(vecin), vecin
    else
        return vecin
    end
end

function procc_RowsElem(rowswh::String, tmp_row::Vector{Float64}, tmp_noNA::Vector{Float64})
    pp_spl = split(rowswh, "e")
    if rowswh == "NA"
        pp = -1 |> Float64
    ## If num is too little, e.g. 1e-16
    elseif length(pp_spl) > 1 && parse(Int64, pp_spl[2]) < -15
        pp = 0 |> Float64
        append!(tmp_noNA, pp)
    else
        pp = parse(Float64, rowswh)
        append!(tmp_noNA, pp)
    end
    append!(tmp_row, pp)
    return tmp_row, tmp_noNA
end
function procc_byRow(oneRow::Vector{String})
    tmp_row = [] |> Vector{Float64}
    tmp_noNA = [] |> Vector{Float64}
    for elm in eachindex(oneRow)
        tmp_row, tmp_noNA = procc_RowsElem(string(oneRow[elm]), tmp_row, tmp_noNA)
    end
    cNA = false
    if length(tmp_noNA) < length(tmp_row); cNA = true; end;
    return tmp_row, tmp_noNA, cNA
end

function calc_var_A2Blines(filePath::String, writePath::String,
                            lineBegin::Int64, lineEnd::Int64,
                            threshold_var::Float64=0.001, threshold_NA::Float64=0.25,
                            delim::String="\t", skipCol::Vector{Int64}=[1])
    mark_new = true
    for rn in lineBegin:lineEnd
        row = my_readline(filePath, rn, delim=delim, maxLen=8192*32)
        tmp_row, tmp_noNA, cNA = procc_byRow(row[Not(skipCol)])
        if (1 - threshold_NA) * (length(row) - length(skipCol)) > length(tmp_noNA)
            continue
        end
        varN, vecFilled = Var_minmax_normalize(tmp_row, tmp_noNA, cNA)
        if threshold_var < varN
            out = [string(row[skipCol])]
            append!(out, string.(vecFilled))
            if !mark_new
                write_csvRows(writePath, out, true, delim=delim)
            else
                write_csvRows(writePath, out, false, delim=delim)
                mark_new = false
            end
        end
    end
    return nothing
end

function calc_var_A2Blines_s(filePath::String, writePath::String,
                               lineBegin::Int64, lineEnd::Int64,
                               threshold_var::Float64=0.001;
                               skipCol::Vector{Int64}=[1,2,3,4],
                               delim_r::String=",", delim_w::String="\t")
    mark_new = true
    for rn in lineBegin:lineEnd
        row = my_readline(filePath, rn, delim=delim_r, maxLen=8192*32)
        linePure = parse.(Float32, row[Not(skipCol)]) |> Vector{Float32}
        if threshold_var < var(linePure)
            if !mark_new
                write_csvRows(writePath, row, true, delim=delim_w)
            else
                write_csvRows(writePath, row, false, delim=delim_w)
                mark_new = false
            end
        end
    end
    return nothing
end

function filter_cols_byfile(origFilePath::String, tags::Vector{String}, id_Xomics::Vector{String}, delim::String="\t")
    grepLen = length(id_Xomics[1])
    line1 = my_readline(origFilePath, 1, delim=delim)
    cols_inTags = my_readline(origFilePath, 2, delim=delim, maxLen=8192*32) .∈ (tags,)
    cols_Xomics = cut_vec_str(line1, grepLen) .∈ (id_Xomics,)
    avail_cols  = (cols_inTags + cols_Xomics) .> 1
    line1_avail = line1[avail_cols]
    println("Available col num: ", sum(avail_cols))
    ## sort cols
    col_sort = Int64[]
    for idom in eachindex(id_Xomics)
        push!(col_sort, findfirst(x -> x == id_Xomics[idom], cut_vec_str(line1_avail, grepLen)))
    end
    ## header
    o1 = [string(line1[1])]
    append!(o1, line1_avail[col_sort])
    write_csvRows(string(dirname(origFilePath), "/", ".fltCol_", basename(origFilePath)), o1, false, delim=delim)
    ##
    for rown in CSV.Rows(origFilePath, delim=delim, header=2)
        rown = string.(rown)
        out = [string(rown[1])]
        append!(out, rown[avail_cols][col_sort])
        write_csvRows(string(dirname(origFilePath), "/", ".fltCol_", basename(origFilePath)), out, true, delim=delim)
    end
    return nothing
end

##==================================##

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

end # of MyFilters
