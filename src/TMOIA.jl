__precompile__(true)

module TMOIA

#include("utils/utils.jl")

include("utils/module_MyFilters.jl")
using .MyFilters
export filter_cols, filter_rows, filter_rows_SNP, filter_ath

include("utils/module_SplitData.jl")
using .SplitData
export write_randomRepeats, my_split_bootstrap, pipe_randomSplits

include("utils/module_SearchSNP.jl")
using .SearchSNP
export search_gene_of_SNP

include("utils/module_DimensionReduction.jl")
using .DimensionReduction
export pick_SNPsInGene, DimReduction, run_DR
export myStructDataIn, myStructParamTrain

include("utils/module_InterpretModel.jl")
using .InterpretModel
export interpreter_inDir, interpreter

end # module TMOIA
