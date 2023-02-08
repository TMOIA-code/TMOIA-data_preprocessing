#
argv <- commandArgs(trailingOnly = TRUE)
#path_ids_1135 = "ath_ids_1135.txt"
#setwd("ath")
library(foreach)

pick_traits = function(traits=c("FT16", "DTF1"), pathToPheno="ath/pheno") {
  traits = sort(traits)
  files = dir(pathToPheno)[which(as.vector(gregexec("CSV", dir(pathToPheno)), mode="numeric") > -1)] |> sort()
  name_tra = strsplit(files, ".CSV") |> as.vector(mode="character") |> sort()
  fileWh = which(name_tra %in% traits)
  fs = list()
  for (fn in fileWh) {
    fs[[length(fs)+1]] = read.csv(paste0(pathToPheno, "/", files[fn]))
  }
  names(fs) = traits
  phMerge = fs[[1]]
  if(length(traits) > 1){
    for (fnn in 2:length(fs)) {
      phMerge = merge(phMerge, fs[[fnn]], by.x = "ID", by.y = "ID")
    }
  }
  colnames(phMerge) = c("ID", traits)
  phMerge = phMerge[!duplicated(phMerge$ID), ]
  rownames(phMerge) = NULL
  phMerge = phMerge[order(phMerge[,1]), ]
  rownames(phMerge) = NULL
  ## z-score std for phenotypes
  phMerge_std = cbind(phMerge[,1], scale(phMerge[, 2:ncol(phMerge)]))
  colnames(phMerge_std) = colnames(phMerge)
  return(list(pheno=phMerge, pheno_std=phMerge_std))
}
#write.table(, "ath_phenos.txt", quote = FALSE, sep = "\t", row.names = FALSE)

prepr_g2s = function(path_g2s="_gene_onlyGene__f-X_snp_matrix.csv_THvar-0.25_ATXG.txt",
                     path_w="_snp2gene_THvar025_onlyGene.txt"){
  g2ss = readLines(path_g2s)[-1]
  g2ss = foreach(sl = 1:length(g2ss), .combine="rbind") %do% {
    strsplit(g2ss[sl], "\t")[[1]][1:2]
  }
  colnames(g2ss) = NULL
  write.table(g2ss, path_w, quote=F, sep="\t", row.names=F, col.names=F)
  return(g2ss)
}

pick_snp_g2s = function(path_g2s="THvar025/snp/_snp2gene_THvar025.txt"){
  g2s = read.table(path_g2s, sep = "\t")
  return(sort(g2s[,1]))
}


pick_IDsSNP_main = function(traits = c("FT16"),
                            dir_write = "ath/dataset/.test/",
                            dirPheno = "ath/pheno",
                            pathSNP = "ath/snp/_f-X_snp_matrix.csv_THvar-0.25.csv",
                            pathG2SNP = "_snp2gene_THvar025_onlyGene.txt",
                            pathMeths = c("ath/meth_exp/627_mCG.csv",
                                          "ath/meth_exp/627_mCHG.csv",
                                          "ath/meth_exp/627_mCHH.csv"),
                            name_meth = c("mCG", "mCHG", "mCHH"),
                            pathExp = "ath/meth_exp/627_exp.csv",
                            pathIDs1135 = "ath/ath_ids_1135.txt",
                            isReturnID = FALSE){
  #
  meths = list()
  for (mCXX in 1:length(pathMeths)) {
    meths[[mCXX]] = read.csv(pathMeths[mCXX])
    if (mCXX == 1) { ids_meth = meths[[1]][,1] }
    if (mCXX > 1) { ids_meth = intersect(ids_meth, meths[[mCXX]][,1]) }
  }
  ids_meth = sort(ids_meth)
  names(meths) = name_meth
  #
  expF = read.csv(pathExp)
  ids_exp = expF[,1] |> sort()
  #
  phenotypesList = pick_traits(traits, dirPheno)
  ids_pheno = phenotypesList[[1]][,1] |> sort()
  #
  ids_1135 = read.table(pathIDs1135, sep = "\t")[, 1]# |> sort()
  ids_inters = intersect(ids_pheno, ids_1135) |> intersect(ids_meth) |> intersect(ids_exp) |> sort()
  print(paste0("Num of samples: ", length(ids_inters)))
  if (isReturnID) { return(ids_inters) }
  whRowMeth = which(ids_meth %in% ids_inters)
  #
  SNPsLine = pick_snp_g2s(pathG2SNP)
  SNP = read.table(pathSNP, sep = ",")
  whColSnp = which(ids_1135 %in% ids_inters) |> sort()
  SNP = SNP[SNPsLine, 5:ncol(SNP)][, whColSnp] |> t()
  #
  whRowPhe = which(ids_pheno %in% ids_inters)
  phenotypes = phenotypesList[[1]][whRowPhe, ]
  phenotypes_std = phenotypesList[[2]][whRowPhe, ]
  #
  if(!dir.exists(dir_write)){
    dir.create(dir_write)
  }
  # Write ID & feature
  write.table(ids_inters, paste0(dir_write, "id.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(colnames(meths[[1]])[-1], paste0(dir_write, "gene_meth.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(colnames(expF)[-1], paste0(dir_write, "gene_exp.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  # Write X
  write.table(SNP, paste0(dir_write, "SNP.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(expF[which(ids_exp %in% ids_inters), -c(1)], paste0(dir_write, "exp.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  for (mcx in 1:length(pathMeths)) {
    write.table(meths[[mcx]][whRowMeth, -c(1)], paste0(dir_write, name_meth[mcx], ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  # Write Y
  write.table(colnames(phenotypes)[-1], paste0(dir_write, "trait.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(phenotypes[, -c(1)], paste0(dir_write, "pheno.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(phenotypes_std[, -c(1)], paste0(dir_write, "pheno_zscore.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


if(TRUE){
  trait = strsplit(argv[1], "[+]")[[1]]
  pick_IDsSNP_main(trait, argv[2])#"FT16_AllOmics/"
}


#================================== Save RData =================================
if(FALSE){
  source("utils.R")
  save.image_multiThreads("ath/dataset/var_0.25", 3)
  #write.table(matrix(which(ids_1135 %in% pick_ids), 1), "snp_id_pos.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  #write.table(matrix(ids_1135, 1), "snp_THvar0.25_withID.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  #write.table(snp_var0.25, "snp_THvar0.25_withID.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
#================================ write genes ==================================
if(FALSE){
  path_g2s="geneToSNP_THvar0.25.txt"
  g2s = read.table(path_g2s, sep = "\t")
  g2s = g2s[order(g2s[,1]), ]
  genes = g2s[!duplicated(g2s[,2]), 2:3]
  write.table(g2s,
              paste0("ath/dataset/THvar0.25_snpdata/", "gene2SNP_sorted_THvar0.25.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(genes,
              paste0("ath/dataset/THvar0.25_snpdata/", "genes_THvar0.25.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

