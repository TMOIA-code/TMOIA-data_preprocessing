__precompile__(true)

module TMOIA

include("modules/utils.jl")
using .MyUtils
export rename_files, my_count_lines, my_read_table, my_readline, my_write_table

include("modules/module_MyFilters.jl")
using .MyFilters
export filter_cols, filter_rows, filter_rows_SNP, filter_ath

include("modules/module_SplitData.jl")
using .SplitData
export write_randomRepeats, my_split_bootstrap, pipe_randomSplits

include("modules/module_SearchSNP.jl")
using .SearchSNP
export search_gene_of_SNP

include("modules/module_DimensionReduction.jl")
using .DimensionReduction
export pick_SNPsInGene, DimReduction, run_DR
export myStructDataIn, myStructParamTrain

include("modules/module_InterpretModel.jl")
using .InterpretModel
export interpreter_inDir, interpreter

end # module TMOIA
