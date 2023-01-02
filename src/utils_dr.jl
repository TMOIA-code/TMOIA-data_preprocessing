# Utils for dimension reduction
__precompile__()
include("utils.jl")

#################################### reduce SNPs ####################################

## Pick SNPs with gene found
function pick_SNPsInGene_ut(L_sta::Int64, L_end::Int64,
                            FileSNP::String, snp_remain::Vector{Int64}, IdPos::Vector{Int64}, w_name::String)
    #write 1st row
    @views tmpLine = my_readline(FileSNP, snp_remain[L_sta], delim=',')[Not([1,2,3,4])][IdPos]
    write_csvRows(w_name, tmpLine, false)
    @views for lnN in (L_sta + 1):L_end
        tmpLine = my_readline(FileSNP, snp_remain[lnN], delim=',')[Not([1,2,3,4])][IdPos]
        write_csvRows(w_name, tmpLine, true)
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
