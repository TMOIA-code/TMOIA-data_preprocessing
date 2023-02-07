#library(XML)
library(foreach)
setwd('~/_prj/')

## Download subtypes
if (TRUE) {
  library(TCGAbiolinks)
  subtypes = PanCancerAtlas_subtypes()
  write.csv(subtypes, '_subtype.csv', row.names = F)
}

pick_subtype = function(id="TCGA-E9-A1RI-11A-41R-A169-07",
                        coln=3,
                        subtypes="subtype.csv") {
  sbs = read.csv(subtypes)
  out = c()
  for (tps in 1:nrow(sbs)) {
    len_id = nchar(id)
    len_sd = nchar(sbs[tps, 1])
    if (len_id == len_sd) {
      if (id == sbs[tps, 1]) {
        return(sbs[tps, coln])
      } else { next }
    } else if (len_id < len_sd) {
      next
    } else {
      if (substring(id, 1, len_sd) == sbs[tps, 1]) {
        out = c(out, sbs[tps, coln])
      }
    }
  }
  if (length(out) == 0) { return(NA) }
  return(out)
}

pick_type_TCGA = function(file.TCGA, subtypeN=3, subtypeCSV="subtype.csv",
                          save_csv=TRUE, saveNum=FALSE) {
  type_id = read.table(file.TCGA, nrows = 1, sep = "\t")[-1] |> as.vector(mode='character') |> unique()
  types = foreach(id = type_id, .combine = 'c') %do% {
    pick_subtype(id, subtypeN, subtypeCSV)
  }
  if (save_csv) {
    dirF = strsplit(file.TCGA, '/')[[1]]
    dirF = dirF[-length(dirF)] |> paste(collapse = "/")
    ## original
    save_nam0 = paste0(dirF, '/labels_sbt', subtypeN, '_asN', 0, '_',
                      tail(strsplit(file.TCGA, '/')[[1]], 1), '.csv')
    write.table(as.matrix(types), save_nam0, 
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ',')
    ## coded
    if (saveNum) {
      save_nam1 = paste0(dirF, '/labels_sbt', subtypeN, '_asN', 1, '_',
                         tail(strsplit(file.TCGA, '/')[[1]], 1), '.csv')
      write.table(as.matrix(as.numeric(as.factor(types))), save_nam1, 
                  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ',')
    }
  }
  #return(type_id)
}


if(FALSE){
  pick_type_TCGA = function(file.TCGA, isNum=TRUE, save_csv=TRUE){
    if (isNum) {
      type_id = read.table(file.TCGA, nrows = 1, sep = "\t")[-1] |> as.vector(mode='character') |> unique() |> substring(14, 15) |> as.numeric()
    } else {
      type_id = read.table(file.TCGA, nrows = 1, sep = "\t")[-1] |> as.vector(mode='character') |> unique() |> substring(14, 16)
    }
    if (save_csv) {
      dirF = strsplit(file.TCGA, '/')[[1]]
      dirF = dirF[-length(dirF)] |> paste(collapse = "/")
      save_nam = paste0(dirF, '/labels_', tail(strsplit(file.TCGA, '/')[[1]], 1), '.csv')
      write.table(as.matrix(type_id), save_nam, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ',')
    }
    #return(type_id)
  }
  
  download_TCGA_xml = function(file.TCGA, prj="TCGA-KIRC") {
    TCGAs = read.table(file.TCGA, nrows = 1, sep = "\t")[-1] |> as.vector(mode='character') |> unique() |> substring(1, 12)
    ## Download xml
    queryXX = GDCquery(project = prj, legacy = TRUE,
                       data.category = 'Biospecimen', barcode = TCGAs)
    #data.type = '', #platform = 'IlluminaGA_miRNASeq')
    GDCdownload(queryXX)
    queryXX = queryXX[[1]][[1]]
    queryXX = data.frame(id=queryXX$id, cases=queryXX$cases,
                         project=queryXX$project, file_name=queryXX$file_name)
    ## Save index:id-tcga
    file.TCGA.name = strsplit(file.TCGA, '/')[[1]] |> tail(1)
    print(file.TCGA.name)
    print(prj)
    if(!dir.exists('idx')) { dir.create('idx') }
    grep_exi = grep(file.TCGA.name, dir('idx'))
    if(length(grep_exi) > 0) {
      queryXX = rbind(read.csv(paste0('idx/', dir('idx')[grep_exi[1]])), queryXX)
      write.csv(queryXX, paste0('idx/', 'index_', file.TCGA.name, '.csv'), row.names=F, quote=F)
    } else {
      write.csv(queryXX, paste0('idx/', 'index_', file.TCGA.name, '.csv'), row.names=F, quote=F)
    }
  }
  
  ## unused
  find_TCGA_type = function(idx_csv, prj='TCGA-KIRC', TCGA="TCGA-A3-3307",
                            dire='GDCdata/') {
    idx = read.csv(idx_csv)
    idx = idx$id[which(idx$cases == TCGA)]
    dire = paste0(dire, prj, '/legacy/Biospecimen/Biospecimen_Supplement/', idx)
    xmlF = dir(dire)[1]
    tmp_sample_type = xmlToList(paste0(dire, '/', xmlF))[["patient"]][["samples"]][["sample"]]
    sample_type = tmp_sample_type[["sample_type"]][["text"]]
    sample_type_id = tmp_sample_type[["sample_type_id"]][["text"]] |> as.numeric()
    return(list(type_id = sample_type_id, type_ = sample_type))
  }
}
