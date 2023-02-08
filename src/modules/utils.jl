module MyUtils

export my_count_lines, my_read_table, my_readline, my_write_table
export std_normalize!, std_zscore!
export conc_fileParts, format_numLen, cut_vec_str, multiThreadSplit, rename_files

using StatsBase, Tables, DelimitedFiles
#using StatsBase: ZScoreTransform, UnitRangeTransform, transform, fit
#using PProf, Profile
#using BenchmarkTools
#using Statistics, DataFrames

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


function my_read_table(XPath::String, type::DataType=Float32, delim::Char=',', isTranspose::Bool=false)
    txt = open(XPath) do file; read(file, String); end;
    out = readdlm(IOBuffer(txt), delim, type)
    if isTranspose; out = transpose(out) |> Matrix; end;
    return out
end

function my_write_table(tableO, XPath::String; delim::Char='\t', isAppend::Bool=false, toTable::Bool=false)
    if isAppend
        isapd = "a"
    else
        isapd = "w"
    end
    if toTable
        open(XPath, isapd) do io
            writedlm(io, Tables.table(tableO), delim)
        end
    else
        open(XPath, isapd) do io
            writedlm(io, tableO, delim)
        end
    end
    return nothing
end

function std_zscore!(mx::Matrix) ## By row
    dt = StatsBase.fit(StatsBase.ZScoreTransform, mx)
    StatsBase.transform!(dt, mx)
end
function std_normalize!(mx::Matrix)
    dt = StatsBase.fit(StatsBase.UnitRangeTransform, mx)
    StatsBase.transform!(dt, mx)
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

#=
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
=#

#=
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
=#

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

function cut_vec_str(vec_str::Vector, len::Int64=15)
    out = [] |> Vector{String}
    for strN in eachindex(vec_str)
        append!(out, [first(vec_str[strN], len)])
    end
    return out
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
# e.g.
# rename_files("xxx/r10/", "val", "x-val", "phen")
# rename_files("xxx/r10/", "val", "y-val", "x-")

end
