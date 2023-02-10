module SearchSNP

if !isdefined(@__MODULE__, :MyUtils)
    include("utils.jl")
    using .MyUtils
end

export search_gene_of_SNP

function x_in_vecInt(vecIn::Vector{Int64})::Int64
    mtx = 1
    for nElem in eachindex(vecIn)
        mtx = mtx * vecIn[nElem]
    end
    return mtx
end

## Warn: Trashing this!!!
function BS_approx_num(goal::Int64, pathF::String, nrowF::Int64=my_count_lines(pathF), inWhichCol::Vector{Int64}=[4,5])::Int64
    guessO::Float64 = ceil(nrowF / 2)#floor??
    guessΔ::Float64 = ceil(nrowF / 4)
    xTry = Vector{Int64}(undef, length(inWhichCol))
    xTryStr = my_readline(pathF, Int64(guessO), maxLen=50)
    #if only(xTryStr) == ""; xTryStr = my_readline(pathF, Int64(guessO)-1, maxLen=50); end;
    xTry .= parse.(Int64, xTryStr[inWhichCol])
    mark_cantfind::Int64 = 0
    while x_in_vecInt(sign.(goal .- xTry)) > 0
        if sum(goal .- xTry) > 0
            guessO += guessΔ
        elseif sum(goal .- xTry) < 0
            guessO -= guessΔ
        end
        xTryStr = my_readline(pathF, Int64(guessO), maxLen=50)
        if length(xTryStr) == 0
            guessO -= 1
            xTryStr = my_readline(pathF, Int64(guessO), maxLen=50)
        end
        xTry .= parse.(Int64, xTryStr[inWhichCol])
        guessΔ = ceil(guessΔ / 2)
        if guessΔ < 2
            mark_cantfind += 1
            if mark_cantfind > 4 # Can't find goal
                return 0
            end
        end
    end
    return Int64(guessO)
end

function util_search_LbyL(lineApprox::Int64, pathGTF::String, nChr::String, nPos::Int64, markL::Int64)
    ## grep lines by gene id
    tmp_gene = split(my_readline(pathGTF, lineApprox)[9], "\"")[2]
    tmp_lines_gtf = split(read(`grep $tmp_gene $pathGTF`, String), "\n") |> Vector{String}
    if tmp_lines_gtf[end] == ""; pop!(tmp_lines_gtf); end;
    ## Count ranges
    tmp_lengths = zeros(Int64, length(tmp_lines_gtf))
    tmp_lines = Vector{Vector{String}}(undef, length(tmp_lines_gtf))
    for ns in eachindex(tmp_lines_gtf)
        tmp_l = split(tmp_lines_gtf[ns], "\t") |> Vector{String}
        staEnd = parse.(Int64, tmp_l[4:5])
        DnUp = nPos .- staEnd
        if (DnUp[1] * DnUp[2]) ≤ 0
            tmp_lines[ns] = tmp_l
            tmp_lengths[ns] = staEnd[2] - staEnd[1] + 1
        else
            tmp_lines[ns] = ["0", "0"]
        end
    end
    tmp_lines_gtf, tmp_l = nothing, nothing
    ## Find SNP's position accurately. ???????????? findfirst or findall ?????????
    tmp_gl = tmp_lines[findfirst(x -> length(x)>2 && x[3]=="gene", tmp_lines)]# |> only
    best_line = Vector{String}(undef, length(tmp_gl))
    gene_line = tmp_gl
    tmp_gl = nothing
    ### trick
    #if minimum(tmp_lengths) < 0; error("\n!!!!!!!!!!!!!!!!!!!! Error !!!!!!!!!!!!!!!!!\n"); end;
    tmp_lengths[findall(x -> x==0, tmp_lengths)] .= typemax(Int64)
    best_line = tmp_lines[findmin(tmp_lengths)[2]]
    tmp_lines = nothing
    gLo, gUp = parse(Int64, gene_line[4]), parse(Int64, gene_line[5])
    xLo, xUp = parse(Int64, best_line[4]), parse(Int64, best_line[5])
    gene_line = nothing
    ## Write
    outGTF = Vector{String}(undef, (1 + length(best_line)))
    outGTF[1] = string(markL)
    outGTF[2:end] .= best_line
    outGene = Vector{String}(undef, 13)
    #### | which line in SNPmx | gene | seq id | type | [type]'s length | chr | SNP's position in chr | SNP's dist to xLo | SNP's dist to xUp | xPerc | SNP's dist to gLo | SNP's dist to gUp | gPerc |
    outGene .= string.([markL, tmp_gene, best_line[1], best_line[3], minimum(tmp_lengths), nChr, nPos, (nPos-xLo), (xUp-nPos), ((nPos-xLo)/minimum(tmp_lengths)), (nPos-gLo), (gUp-nPos), ((nPos-gLo)/maximum(tmp_lengths))]) |> Vector{String}
    return outGene, outGTF
end

function util_search_gene_by_snp(nChr::String, nPos::Int64, markL::Int64, pathGTF_chrX::String, nRowGTFx::Int64, intergenic::Bool)::String
    outVec = Vector{String}(undef, 5)
    ## BS
    lineApprox = BS_approx_num(nPos, pathGTF_chrX, nRowGTFx, Int64[4,5])
    if lineApprox == 0
        if intergenic
            geneSelfDef = string("interG_", nChr, "_", nPos)
            outVec = string.([markL, geneSelfDef, nChr, nPos, ""])
        end
    else
        tmpO = my_readline(pathGTF_chrX, lineApprox)[[1,9]]
        gene! = split(tmpO[2], "\"")[2]
        #### | which line in SNPmx | gene | nChr | nPos | seq id |
        outVec = string.([markL, gene!, nChr, nPos, tmpO[1]])
    end
    if isassigned(outVec, 1)
        out = join(outVec, "\t")
    else
        out = ""
    end
    return out
end

function util_search_region_by_snp(nChr::String, nPos::Int64, markL::Int64, pathGTF_chrX::String, nRowGTF::Vector{Int64}, intergenic::Bool)
    outGene = [""]
    outGTF = [""]
    ## BS
    lineApprox = BS_approx_num(nPos, pathGTF_chrX, nRowGTF[parse(Int64, nChr)], Int64[4,5])
    if lineApprox == 0
        if intergenic
            outGene = string.([markL, string("interG_", nChr, "_", nPos), nChr, nPos])
        end
    else
        outGene, outGTF = util_search_LbyL(lineApprox, pathGTF_chrX, nChr, nPos, markL)
    end
    if length(outGene) < 2
        o_gene = ""
    else
        o_gene = join(outGene, "\t")
    end
    if length(outGTF) < 2
        o_gtf = ""
    else
        o_gtf = join(outGTF, "\t")
    end
    return o_gene, o_gtf
end

##==============================================- Main utils -==============================================##

# Give simple results: | which line in SNPmx | gene | Chr | SNP Position | seq id |
function search_gene_by_snp(pathSNPmx::String, dirGTF::String, nRowGTF::Vector{Int64},
                                LnBegin::Int64, LnEnd::Int64, nT::Int64,
                                intergenic::Bool=true, skipRow::Int64=0)
    ##
    dirW = dirname(pathSNPmx) * "/"
    pathGeneOut = string(dirW, ".p-gene-", format_numLen(nT, 2), "-", basename(pathSNPmx)[1:end-4], "_", basename(dirGTF), ".txt")
    rm(pathGeneOut, force=true)
    tmpG::String = "NA"
    ##
    io = open(pathGeneOut, "a")
    for Lx in LnBegin:LnEnd
        markL = Lx - skipRow
        posPart = my_readline(pathSNPmx, Lx, delim=',', maxLen=36, nElem=2)
        nChr = posPart[1][4:end] |> String
        nPos = parse(Int64, posPart[2])
        pathGTF_chrX = string(dirGTF, "/", "AT", nChr, "G.gtf") ## The GTF file to read
        ## Search
        tmpG = util_search_gene_by_snp(nChr, nPos, markL, pathGTF_chrX, nRowGTF[parse(Int64, nChr)], intergenic)
        if length(tmpG) > 0; write(io, tmpG, "\n"); end;
    end
    close(io)
    return nothing
end

# Give full information
function search_region_by_snp(pathSNPmx::String, dirGTF::String, nRowGTF::Vector{Int64},
                                LnBegin::Int64, LnEnd::Int64, nT::Int64,
                                intergenic::Bool=true, skipRow::Int64=0)
    ##
    dirW = dirname(pathSNPmx) * "/"
    pathGTFout = string(dirW, ".p-gtf-", format_numLen(nT, 2), "-", basename(pathSNPmx)[1:end-4], "_", basename(dirGTF), ".gtf")
    pathGeneOut = string(dirW, ".p-gene-", format_numLen(nT, 2), "-", basename(pathSNPmx)[1:end-4], "_", basename(dirGTF), ".txt")
    ##
    oGene, oGTF = Vector{String}(undef, LnEnd-LnBegin+1), Vector{String}(undef, LnEnd-LnBegin+1)
    for Lx in LnBegin:LnEnd
        markL = Lx - skipRow
        posPart = my_readline(pathSNPmx, Lx, delim=',', maxLen=36, nElem=2)
        nChr = posPart[1][4:end] |> String
        nPos = parse(Int64, posPart[2])
        pathGTF_chrX = string(dirGTF, "/", "AT", nChr, "G.gtf") ## The GTF file to read
        ##
        oGene[Lx-LnBegin+1], oGTF[Lx-LnBegin+1] = util_search_region_by_snp(nChr, nPos, markL, pathGTF_chrX, nRowGTF, intergenic)
    end
    rm(pathGeneOut, force=true)
    rm(pathGTFout, force=true)
    open(pathGeneOut, "a") do io
        for gx in eachindex(oGene)
            if length(oGene[gx]) > 1; write(io, oGene[gx], "\n"); end;
        end
    end
    open(pathGTFout, "a") do io
        for gx in eachindex(oGTF)
            if length(oGTF[gx]) > 1; write(io, oGTF[gx], "\n"); end;
        end
    end
    return nothing
end

function search_gene_of_SNP(pathSNPmx::String, dirGTF::String, intergenic::Bool=true; skipRow::Int64=0, slimResult::Bool=true)
    ##sh split_chr.sh
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

end # of SearchSNP
