
ABC_all_calc <- function(scores,
                         pred_file,
                         key_file = NA,
                         output_cell,
                         output_bed = "temp.bed"){
 # library(data.table)
 # cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
 # tabb = data.frame(fread(paste0(cell, "/", "AllPredictions.AvgHiC.ABC0.015.500.encodeFormat.txt.gz")))
  tabb = data.frame(paste0(pred_file))
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")

  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ABC_blood_calc <- function(scores,
                           pred_file,
                           key_file,
                           output_cell,
                           output_bed = "temp.bed"){
 # library(data.table)
 #  cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
 # tabb = data.frame(fread(paste0(cell, "/", "AllPredictions.AvgHiC.ABC0.015.500.encodeFormat.txt.gz")))
 #  tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC.txt", header=F)[,1])
  tabb = data.frame(paste0(pred_file))
  tissuenames2 = as.character(read.table(paste0(key_file)))
  tissuenames2 = intersect(unique(tabb$CellType), tissuenames2)
  tabb=tabb[which(tabb$CellType %in% tissuenames2==T), ]
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")

  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ABC_brain_calc <- function(scores,
                           pred_file,
                           key_file=NA,
                           output_cell,
                           output_bed = "temp.bed"){
  #library(data.table)
  #cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
  #tabb = data.frame(fread(paste0(cell, "/", "AllPredictions.AvgHiC.ABC0.015.500.encodeFormat.txt.gz")))
  tabb = data.frame(fread(paste0(pred_file)))
  tissuenames2 = c("neur", "Neur", "astro", "spinal", "Brain", "brain")
  tissue_ids = as.numeric(unlist(sapply(tissuenames2, function(x) return(grep(x, tabb$CellType)))))
  tabb = tabb[tissue_ids, ]
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")

  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ChromHMM_all_calc <- function(scores,
                              pred_file,
                              key_file,
                              output_cell,
                              output_bed = "temp.bed"){
  #chromhmm_biosamples = read.delim("/n/groups/price/kushal/ENCODE/data/EG_predictions/ChromHMM_biosample_key.tsv")
  #library(data.table)
  #cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
  #tabb_pre = data.frame(fread(paste0(cell, "/", "AllPredictions_ChromHMM2017_encodeFormat.tsv.gz")))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  tabb = tabb_pre
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  chromhmm_biosamples = read.delim(paste0(key_file))
  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


ChromHMM_blood_calc <- function(scores,
                                pred_file,
                                key_file,
                                output_cell,
                                output_bed = "temp.bed"){
  #chromhmm_biosamples = read.delim("/n/groups/price/kushal/ENCODE/data/EG_predictions/ChromHMM_biosample_key.tsv")
  chromhmm_biosamples = read.delim(paste0(key_file))
  blood_eids = chromhmm_biosamples[grep("BLD", chromhmm_biosamples[,3]), 1]
  #cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
  #tabb_pre = data.frame(fread(paste0(cell, "/", "AllPredictions_ChromHMM2017_encodeFormat.tsv.gz")))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  tabb = tabb_pre[which(tabb_pre$CellType %in% blood_eids == T), ]
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")

  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ChromHMM_brain_calc <- function(scores,
                                pred_file,
                                key_file,
                                output_cell,
                                output_bed = "temp.bed"){
  #chromhmm_biosamples = read.delim("/n/groups/price/kushal/ENCODE/data/EG_predictions/ChromHMM_biosample_key.tsv")
  chromhmm_biosamples = read.delim(paste0(key_file))
  brain_eids = chromhmm_biosamples[grep("BRN", chromhmm_biosamples[,3]), 1]
  #library(data.table)
  #cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
  #tabb_pre = data.frame(fread(paste0(cell, "/", "AllPredictions_ChromHMM2017_encodeFormat.tsv.gz")))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  tabb = tabb_pre[which(tabb_pre$CellType %in% brain_eids == T), ]
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")

  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

EpiMap_all_calc <- function(scores,
                            pred_file,
                            key_file,
                            output_cell,
                            output_bed = "temp.bed"){
  #epimap_biosamples = read.delim("/n/groups/price/kushal/ENCODE/data/EG_predictions/EpiMap_biosample_key.tsv")
  epimap_biosamples = read.delim(paste0(key_file))
  #library(data.table)
  #cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
  #tabb_pre = data.frame(fread(paste0(cell, "/", "EpiMap_Appended.txt.gz")))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  tabb = tabb_pre
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

EpiMap_blood_calc <- function(scores,
                              pred_file,
                              key_file,
                              output_cell,
                              output_bed = "temp.bed"){
  #epimap_biosamples = read.delim("/n/groups/price/kushal/ENCODE/data/EG_predictions/EpiMap_biosample_key.tsv")
  epimap_biosamples = read.delim(paste0(key_file))
  blood_eids = epimap_biosamples$biosample[c(grep("Blood & T-cell", epimap_biosamples$sampleCategory),
                                             grep("HSC & B-cell", epimap_biosamples$sampleCategory),
                                             grep("Lymphoblastoid", epimap_biosamples$sampleCategory))]
  #library(data.table)
  #cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
  #tabb_pre = data.frame(fread(paste0(cell, "/", "EpiMap_Appended.txt.gz")))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  tabb = tabb_pre[which(tabb_pre$CellType %in% blood_eids == T), ]
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

EpiMap_brain_calc <- function(scores,
                              pred_file,
                              key_file,
                              output_cell,
                              output_bed = "temp.bed"){
  #epimap_biosamples = read.delim("/n/groups/price/kushal/ENCODE/data/EG_predictions/EpiMap_biosample_key.tsv")
  epimap_biosamples = read.delim(paste0(key_file))
  brain_eids = epimap_biosamples$biosample[c(grep("Brain", epimap_biosamples$sampleCategory),
                                             grep("Neurosph", epimap_biosamples$sampleCategory))]
  #library(data.table)
  #cell="/n/groups/price/kushal/ENCODE/data/EG_predictions"
  #tabb_pre = data.frame(fread(paste0(cell, "/", "EpiMap_Appended.txt.gz")))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  tabb = tabb_pre[which(tabb_pre$CellType %in% brain_eids == T), ]
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")

  matched_ids = match(tabb2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(tabb2[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
