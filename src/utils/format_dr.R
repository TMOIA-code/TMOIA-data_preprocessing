# Proccess dr

ids_split = function(dir_id_r10) {
  files = dir(dir_id_r10)
  ids = c()
  for (nId in files) {if (as.numeric(regexec('id', nId)[[1]][1]) != -1) {ids = c(ids, nId)}}
  return(ids)
}
dr_split_snp = function(dir_id_r10, dir_snp_r10) {
  paths_ids = ids_split(dir_id_r10)
  files_snp = dir(dir_snp_r10)
  fs_snp = c()
  for (ns in files_snp) {if (as.numeric(regexec('snp', ns)[[1]][1]) != -1) {fs_snp = c(fs_snp, ns)}}
  print(fs_snp)
  for (f_snp in 1:length(fs_snp)) {
    path_snp = paste0(dir_snp_r10, "/", fs_snp[f_snp])
    tmp_snp = read.table(path_snp, TRUE, ",")
    tmp_fs_ids = paths_ids[(3*(f_snp-1)+1):(3*f_snp)]
    for (tvt in 1:3) {
      tmp_path = paste0(dir_id_r10, "/", tmp_fs_ids[tvt])
      tmp_ids = read.table(tmp_path, sep="\t")[,1]
      whL = c()
      for (nid in tmp_ids) {whL = c(whL, which(tmp_snp[,1] == nid))}
      write.table(tmp_snp[whL, -1],
                  paste0(dir_snp_r10, "/", substring(tmp_fs_ids[tvt], 1, 6), "_dr_snp.txt"),
                  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }
}


dr_split_other = function(dir_r10_id, dir_r10_dr) {
  paths_ids = ids_split(dir_r10_id)
  files_om = dir(dir_r10_dr)
  fs_om = c()
  for (ns in files_om) {if (as.numeric(regexec('m', ns)[[1]][1]) != -1) {fs_om = c(fs_om, ns)}}
  print(fs_om)
  for (f_n in 1:length(fs_om)) {
    om = gsub("2000.csv", "", fs_om[f_n])
    path_om = paste0(dir_r10_dr, "/", fs_om[f_n])
    tmp_om = read.table(path_om, TRUE, ",")
    #drs = substr(paths_ids, 1,2) |> unique()
    for (drs in 1:10) {
      tmp_fs_ids = paths_ids[(3*(drs-1)+1):(3*drs)]
      for (tvt in 1:3) {
        tmp_path = paste0(dir_r10_id, "/", tmp_fs_ids[tvt])
        tmp_ids = read.table(tmp_path, sep="\t")[,1]
        whL = c()
        for (nid in tmp_ids) {whL = c(whL, which(tmp_om[,1] == nid))}
        write.table(tmp_om[whL, -1],
                    paste0(dir_r10_dr, "/", substring(tmp_fs_ids[tvt], 1, 6), "_dr_", om, ".txt"),
                    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      }
    }
    
  }
}
