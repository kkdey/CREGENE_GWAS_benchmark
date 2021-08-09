# rules to perform comparisons of CRE predictions to CRISPR data

def get_pred_config(wildcards):
  pred_config = config["comparisons"][wildcards.comparison]["pred_config"]
  if pred_config is None:
    comparison = wildcards.comparison
    pred_config = "results/" + comparison + "/pred_config.txt"
  return pred_config

## RULES -------------------------------------------------------------------------------------------

rule build_bedgraphs_from_programs_all:
  input:
    programs = lambda wildcards: config["comparisons"][wildcards.comparison]["programs"]["all"],
  output:
    outdir = "results/{comparison}/bedfiles/"
  log: "results/{comparison}/logs/build_bedgraphs_from_programs_all.log"
  params:
    pred_file_ABC = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"]["ABC"],
    pred_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"]["ChromHMM"],
    pred_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["pred]["EpiMap"],
    key_file_ABC = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["ABC"],
    key_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["ChromHMM"],
    key_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["EpiMap"]
  conda: "../envs/r_packages.yml"
  script:
   "../../workflow/scripts/build_bedgraphs_from_programs_all.R"


rule build_bedgraphs_from_programs_blood:
  input:
    programs = lambda wildcards: config["comparisons"][wildcards.comparison]["programs"]["all"],
  output:
    outdir = "results/{comparison}/bedfiles/"
  log: "results/{comparison}/logs/build_bedgraphs_from_programs_blood.log"
  params:
    pred_file_ABC = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"]["ABC"],
    pred_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"]["ChromHMM"],
    pred_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["pred]["EpiMap"],
    key_file_ABC = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["ABC"],
    key_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["ChromHMM"],
    key_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["EpiMap"]
  conda: "../envs/r_packages.yml"
  script:
   "../../workflow/scripts/build_bedgraphs_from_programs_blood.R"


rule build_bedgraphs_from_programs_brain:
  input:
    programs = lambda wildcards: config["comparisons"][wildcards.comparison]["programs"]["all"],
  output:
    outdir = "results/{comparison}/bedfiles/"
  log: "results/{comparison}/logs/build_bedgraphs_from_programs_brain.log"
  params:
    pred_file_ABC = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"]["ABC"],
    pred_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"]["ChromHMM"],
    pred_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["pred]["EpiMap"],
    key_file_ABC = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["ABC"],
    key_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["ChromHMM"],
    key_file_ChromHMM = lambda wildcards: config["comparisons"][wildcards.comparison]["key"]["EpiMap"]
  conda: "../envs/r_packages.yml"
  script:
   "../../workflow/scripts/build_bedgraphs_from_programs_brain.R"


rule process_bedfiles:
  input:
    bed_dir = "results/{comparison}/bedfiles/ALL/"
  output:
    outdir = "results/{comparison}/bedfiles/ALL"
  log: "results/{comparison}/logs/process_bedfiles.log"
  conda: "../envs/r_packages.yml"
  script:
   "../../workflow/scripts/process_bedfiles.sh"


rule create_annot_from_bedgraph:
  input:
    bed_dir = "results/{comparison}/bedfiles/ALL/"
  output:
    annot_dir = "results/{comparison}/annotations/ALL"
  log: "results/{comparison}/logs/create_annot_from_bedgraph.log"
  conda: "../envs/ldsc_python.yml"
  params:
    bimpath= lambda wildcards: config["comparisons"][wildcards.comparison]["bimpath"],
  script:
   "../../workflow/scripts/create_annot_from_bedgraph.sh"


rule ldsc_compute_ldscores:
  input:
    annot_dir = "results/{comparison}/annotations/ALL"
  output:
    annot_dir = "results/{comparison}/annotations/ALL"
  log: "results/{comparison}/logs/ldsc_compute_ldscores.log"
  conda: "../envs/ldsc_python.yml"
  params:
    bfile_path= lambda wildcards: config["comparisons"][wildcards.comparison]["bfile_path"],
    hapmap_path= lambda wildcards: config["comparisons"][wildcards.comparison]["hapmap_path"],
    ld_wind_cm= 1
  script:
   "../../workflow/scripts/ldsc_compute_ldscores.sh"


rule ldsc_h2_analysis:
  input:
    annot_dir = "results/{comparison}/annotations/ALL",
    trait = lambda wildcards: config["comparisons"][wildcards.comparison]["trait"]
  output:
    results_dir = "results/{comparison}/ldsc_results/ALL"
  log: "results/{comparison}/logs/ldsc_h2_analysis.log"
  conda: "../envs/ldsc_python.yml"
  params:
    baseline_cell= lambda wildcards: config["comparisons"][wildcards.comparison]["baseline_cell"],
    baseline_version= "baseline_Epi",
    weights_path= lambda wildcards: config["comparisons"][wildcards.comparison]["weights_path"],
    freq_path= lambda wildcards: config["comparisons"][wildcards.comparison]["freq_path"]
  script:
   "../../workflow/scripts/ldsc_h2_analysis.sh"
