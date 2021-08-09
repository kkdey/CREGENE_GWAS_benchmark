## Generate bedgraph files from a gene score file and CRE-gene links files in blood

# open log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(data.table)
library(R.utils)

merge_functions_file <- file.path(snakemake@scriptdir, "ALL_CRE_gene_linkfiles.R")
suppressPackageStartupMessages(source(merge_functions_file))

message("Reading in the gene program directory:", snakemake@input$programdir)

ll = list.files(snakemake@input$programdir, pattern = ".txt")
annot_names = as.character(sapply(ll, function(x) return(strsplit(x, ".txt")[[1]][1])))

for(numl in 1:length(ll)){
  score_file = paste0(snakemake@input$programdir, "/", ll[numl])
  gene_scores = read.delim(score_file, header=F)

  if(!dir.exists(paste0(snakemake@output$outdir, "/", annot_names[numl]))){
    dir.create(paste0(snakemake@output$outdir, "/", annot_names[numl]))
  }
  scores = gene_scores[,2]
  names(scores) = gene_scores[,1]

  out1 = ABC_blood_calc(scores,
                      pred_file = snakemake@params$pred_file_ABC,
                      key_file = snakemake@params$key_file_ABC,
                      output_cell = paste0(snakemake@output$outdir, "/", annot_names[numl]),
                      output_bed = paste0("ABC_blood.bed"))

  out2 = ChromHMM_blood_calc(scores,
                           pred_file = snakemake@params$pred_file_ChromHMM,
                           key_file = snakemake@params$key_file_ChromHMM,
                           output_cell = paste0(snakemake@output$outdir, "/", annot_name),
                           output_bed = paste0("ChromHMM_blood.bed"))

  out3 = EpiMap_blood_calc(scores,
                         pred_file = snakemake@params$pred_file_EpiMap,
                         key_file = snakemake@params$key_file_EpiMap,
                         output_cell = paste0(snakemake@output$outdir, "/", annot_name),
                         output_bed = paste0("EpiMap_blood.bed"))
}
