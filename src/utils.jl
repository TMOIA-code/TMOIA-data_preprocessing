## Utils
__precompile__()
using CSV, Statistics, StatsBase, Tables, DataFrames, DelimitedFiles, Random
#using StatsBase: ZScoreTransform, UnitRangeTransform, transform, fit
#using PProf, Profile
#using BenchmarkTools



function my_count_lines(filePath::String)::Int64## The result may differ to `wc -l filePath`
    num_l::Int64 = 0
    ioin = open(filePath, "r")
    num_l = countlines(ioin)
    close(io)
    return num_l
end


function util_my_readline(filePath::String, nth::Int64, maxLen::Int64=8192*16, eol::AbstractChar='\n', windowSize::Int64=8192)::String
    aeol = UInt8(eol)
    a = Vector{UInt8}(undef, windowSize)
    nl = nb = 0
    nthd = nth - 1
    numElem = maxLen + 1
    outStrU = Vector{UInt8}(undef, numElem)
    UFilled = 0
    io = open(filePath, "r", lock=false)
    while !eof(io)
        nb = readbytes!(io, a)
        @views for i=1:nb#@simd
            @inbounds nl += a[i] == aeol
            if nl == nthd
                UFilled += 1
                outStrU[UFilled] = a[i]
            end
            if nl == nth || UFilled == numElem
                break
            end
        end
        if nl == nth || UFilled == numElem
            break
        end
    end
    close(io)
    outStr = ""
    if UFilled > 1
        if nth > 1
            outStr = join(Char.(outStrU[2:UFilled]))
        else
            outStr = join(Char.(outStrU[1:UFilled]))
        end
    else
        if nth < 2
            outStr = Char(outStrU[1]) * ""
        end
    end
    return outStr
end
##
function my_readline(filePath::String, nth::Int64; maxLen::Int64=8192*16, delim::Char='\t', isSplit::Bool=true, nElem::Int64=0, eol::AbstractChar='\n', windowSize::Int64=8192)::Union{String, Vector{String}}
    outStr = util_my_readline(filePath, nth, maxLen, eol, windowSize)
    if isSplit && length(outStr) > 1
        o_split::Vector{String} = split(outStr, delim)
        if nElem > 0
            return o_split[1:nElem]
        else
            return o_split
        end
    else
        return outStr
    end
end


#= Depend on `sed`
function read_Xth_line(filePath::String, delim::String="\t", whichLine::Int64=1, whichCol::Vector{Int64}=[0]; NoSplit::Bool=false, cut::Bool=false, str_start::Int64=0, str_len::Int64=0)
    if cut
        mySh = string("mywow=`sed -n '", whichLine, "p;", whichLine, "q' ", filePath, "`\n",
                        "echo \${mywow: ", str_start, ": ", str_len, "}\n")
        pathSh = tempname()
        io = open(pathSh, "w"); write(io, mySh); close(io);
        vout_str = read(`sh $pathSh`, String)
    else
        n2, n3 = "p;", "q"
        vout_str = read(`sed -n $whichLine$n2$whichLine$n3 $filePath`, String)
    end
    if NoSplit
        if vout_str[end] == '\n'; vout_str = vout_str[1:end-1]; end;
        return vout_str
    else
        vout = split(vout_str, delim)
        if vout[end][end] == '\n'; vout[end] = vout[end][1:end-1]; end;
        if minimum(whichCol) > 0; return vout[whichCol] |> Vector{String}; end;
        return vout |> Vector{String}
    end
    
end
=#


function my_read_table(Xdir::String, type::DataType=Float32, delim::Char=',', isTranspose::Bool=false)
    txt = open(Xdir) do file; read(file, String); end;
    out = readdlm(IOBuffer(txt), delim, type)
    if isTranspose; out = transpose(out) |> Matrix; end;
    return out
end


function std_zscore!(mx::Matrix) ## By row
    dt = StatsBase.fit(StatsBase.ZScoreTransform, mx)
    StatsBase.transform!(dt, mx)
end
function std_normalize!(mx::Matrix)
    dt = StatsBase.fit(StatsBase.UnitRangeTransform, mx)
    StatsBase.transform!(dt, mx)
end
#
function std_tvt(p1::Matrix, p2::Matrix, p3::Matrix, stdFunc=std_zscore!)
    conc = hcat(p1, p2, p3)
    ncol1, ncol2, ncol3 = size(p1)[2], size(p2)[2], size(p3)[2]
    stdFunc(conc)
    return conc[:, 1:ncol1], conc[:, (ncol1+1):(ncol1+ncol2)], conc[:, (ncol1+ncol2+1):end]
end

function conc_fileParts(grepF::String, path_write::String)
    path_sh = tempname(dirname(path_write))
    path_tmp = tempname(dirname(path_write))
    tmp_sh = string("cat ", grepF, " > ", path_tmp, "\n",
                    "rm ", grepF, "\n",
                    "mv ", path_tmp, " ", path_write, "\n")
    io = open(path_sh, "w"); write(io, tmp_sh); close(io);
    run(`sh $path_sh`)
    rm(path_sh, force=true)
    return nothing
end


function format_numLen(InNum::Int64, lenO::Int64)::String
    ns = string(InNum)
    while length(ns) < lenO
        ns = string("0", ns)
    end
    return ns
end


function fill_NA(vecin::Vector{String})
    vecNums = Float64[]
    vecNew = zeros(Float64, length(vecin))
    @inbounds for nx in eachindex(vecin)
        if vecin[nx] == "NA"
            continue
        else
            push!(vecNums, parse(Float64, vecin[nx]))
        end
    end
    meann = mean(vecNums)
    minn = minimum(vecNums, init=Inf)
    rangen = maximum(vecNums, init=-1000) - minn
    patchn = (meann - minn) / rangen
    ## Fill
    Threads.@threads for nx in eachindex(vecin)
        if vecin[nx] == "NA"
            vecNew[nx] = patchn
        else
            vecNew[nx] = parse(Float64, vecin[nx])
        end
    end
    return reshape(vecNew, 1, size(vecNew)[1])
end


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


function write_csvRows(filename::String, strIn::Vector, isAppend::Bool=false; lineBreakHead::Bool=false, lineBreakTail::Bool=true, delim::String="\t")
    if lineBreakHead; tmp_str = "\n"; else tmp_str = ""; end;
    for tt in 1:length(strIn)
        if tt==1
            tmp_str = string(tmp_str, strIn[tt])
        else
            tmp_str = string(tmp_str, delim, strIn[tt])
        end
    end
    if lineBreakTail; tmp_str = string(tmp_str, "\n"); end;
    if isAppend
        io = open(filename, "a")
    else
        io = open(filename, "w")
    end
    write(io, tmp_str)
    close(io)
    return nothing
end


function multiThreadSplit(numRC::Int64, numThread::Int64)
    pLen = numRC / numThread |> floor
    startPts = [] |> Vector{Int64}
    endPts = [] |> Vector{Int64}
    for nT in 1:numThread
        append!(startPts, Int64(1 + (nT - 1) * pLen))
        append!(endPts, Int64(nT * pLen))
    end
    endPts[numThread] = numRC
    return startPts, endPts
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


function cut_vec_str(vec_str::Vector, len::Int64=15)
    out = [] |> Vector{String}
    for strN in eachindex(vec_str)
        append!(out, [first(vec_str[strN], len)])
    end
    return out
end

#=
function my_split3(num::Int64=500, seed::Int64=1234, ratio_A::Float64=0.7, ratio_B::Float64=0.15)#ratio_C=0.15
    shuffled_indices = Random.shuffle(MersenneTwister(seed), 1:num)
    n_divA = round(Int, ratio_A * num)
    n_divB = round(Int, ratio_B * num)
    indices_pA = shuffled_indices[1:n_divA]
    indices_pB = shuffled_indices[(n_divA + 1):(n_divA + n_divB)]
    indices_pC = shuffled_indices[(n_divA + n_divB + 1):end]
    return indices_pA, indices_pB, indices_pC
end
=#
function my_split_bootstrap(num::Int64=600, seed::Int64=1234)### Size of test&validation proportion is dynamic!
    trnO = sample(MersenneTwister(seed), collect(StepRange(1, Int64(1), num)), num, replace=true, ordered=false)
    num_trn_uniq = length(unique(trnO))
    num_tstval = num - num_trn_uniq
    num_tst = ceil(num_tstval / 2) |> Int64
    tstval = setdiff(collect(StepRange(1, Int64(1), num)), unique(trnO))
    tstO = tstval[1:num_tst]
    valO = tstval[(num_tst + 1):end]
    return trnO, valO, tstO
end
function write_splits(fPaths::Vector{String}, splitN::Int64=1, seed::Int64=1234, dlm::Char='\t')
    num_sample = my_count_lines(fPaths[1])
    idx_trn, idx_val, idx_tst = my_split_bootstrap(num_sample, seed)
    MarkSplit = format_numLen(splitN, 2)
    @views for pn in eachindex(fPaths)
        fin = my_read_table(fPaths[pn], String, dlm)
        CSV.write(string(dirname(fPaths[pn]), "/r10/", MarkSplit, "_trn_", basename(fPaths[pn])), Tables.table(fin[idx_trn, :]), header=false, delim="\t")
        CSV.write(string(dirname(fPaths[pn]), "/r10/", MarkSplit, "_val_", basename(fPaths[pn])), Tables.table(fin[idx_val, :]), header=false, delim="\t")
        CSV.write(string(dirname(fPaths[pn]), "/r10/", MarkSplit, "_tst_", basename(fPaths[pn])), Tables.table(fin[idx_tst, :]), header=false, delim="\t")
    end
    return nothing
end
function write_randomRepeats(fPaths::Vector{String}; repN::Int64=10, seedInit::Int64=2345, isNThreads::Bool=false, dlm::Char='\t')
    mkpath(string(dirname(fPaths[1]), "/r10/"))
    if isNThreads
        Threads.@threads for rn in 1:repN
            write_splits(fPaths, rn, (seedInit + rn), dlm)
        end
    else
        for rn in 1:repN
            write_splits(fPaths, rn, (seedInit + rn), dlm)
        end
    end
    return nothing
end

function pipe_randomSplits(dirIn::String, traits::Vector{String}, fOmics::Vector{String}, FolderSuffix::String)
    for trts in eachindex(traits)
        for oms in eachindex(fOmics)
            write_randomRepeats([string(dirIn, traits[trts], FolderSuffix, "/", fOmics[oms])], isNThreads=true)
        end
    end
    return nothing
end


## Pick samples' id owning meth, mirna and mrna
function sampl_own_Xomics(origFilePaths::Vector{String}, len_id::Int64=16, delim::String="\t", isWrite::Bool=false, wrt_pref::String=".ownAll_")
    out = [] |> Vector{String}
    for fpath in eachindex(origFilePaths)
        if fpath < 2
            out = my_readline(origFilePaths[fpath], 1, delim=delim)[Not(1)]
            if len_id > 1; out = cut_vec_str(out, len_id); end;
            continue
        end
        tcgaB = my_readline(origFilePaths[fpath], 1, delim=delim)[Not(1)]
        if len_id > 1; tcgaB = cut_vec_str(tcgaB, len_id); end;
        out = intersect(unique(out), unique(tcgaB))
    end
    if isWrite; CSV.write(string(wrt_pref, "_len_", len_id, ".csv"), Tables.table(out), header=false); end;
    return out
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


function x_in_vecInt(vecIn::Vector{Int64})::Int64
    mtx = 1
    for nElem in eachindex(vecIn)
        mtx = mtx * vecIn[nElem]
    end
    return mtx
end


function rename_files(dirX::String, keyW::String, changeTo::String, excludeW::String="")
    if dirX[end] == '/'
        dirX = dirX[1:end-1]
    end
    fs = readdir(dirX)
    grepW = keyW |> Regex
    if excludeW != ""
        exclW = excludeW |> Regex
        fs = fs[Not(occursin.(exclW, fs))]
    end
    fs = fs[occursin.(grepW, fs)]
    println(fs)
    for fsn in eachindex(fs)
        tmp_new = replace(fs[fsn], keyW => changeTo)
        tmp_path = string(dirX, "/", tmp_new)
        mv(string(dirX, "/", fs[fsn]), tmp_path)
    end
    return nothing
end
#e.g.
#rename_files("r10/", "val", "x-val", "phen")
#rename_files("r10/", "val", "y-val", "x-")
