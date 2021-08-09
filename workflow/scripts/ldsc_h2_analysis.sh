annot_cell=snakemake@input$annot_dir
baseline_cell=snakemake@params$baseline_cell
baseline_version=snakemake@params$baseline_version
weights_path=snakemake@params$weights_path
freq_path=snakemake@params$freq_path
sumstat=snakemake@input$trait
output_cell_pre=snakemake@output$results_dir

IFS="
"

module load conda2
source activate ldsc

if [ ! -d $output_cell_pre ]
then
    mkdir $output_cell_pre
fi

output_cell=$output_cell_pre/$baseline_version

if [ ! -d $output_cell ]
then
    mkdir $output_cell
fi

echo $output_cell
for line in `ls $annot_cell | awk '{print $1}' | sort | uniq`;
do
    annot_module=`echo $line | awk '{print $1}'`
    echo $annot_cell $annot_module
    if [ ! -d $annot_cell/$annot_module ]
    then
        echo "Error: annotation module directory not found" > ldsc_logfile.log
        exit 100
    fi
    if [ ! -d $output_cell/$annot_module ]
    then
        mkdir $output_cell/$annot_module
    fi
    for ll in `ls $annot_cell/$annot_module | awk '{print $1}' | sort | uniq`;
    do
      	annot_dir=`echo $ll | awk '{print $1}'`
      	echo $annot_dir
      	if [ ! -d $annot_cell/$annot_module/$annot_dir ]
      	then
                  echo "Error: annotation module directory not found" > ldsc_logfile.log
                  exit 101
      	fi
      	if [ ! -d $output_cell/$annot_module/$annot_dir ]
      	then
                  mkdir $output_cell/$annot_module/$annot_dir
      	fi
      	if [ ! -f $output_cell/$annot_module/$annot_dir/$sumstat.results ]
        then
        python ldsc.py --h2 $sumstat --ref-ld-chr $annot_cell/$annot_module/$annot_dir/$annot_dir.,$baseline_cell/$baseline_version/baselineLD.  --frqfile-chr $freq_path/1000G.EUR.QC. --w-ld-chr $weights_path/weights.hm3_noMHC. --overlap-annot --print-coefficients --print-delete-vals --out $output_cell/$annot_module/$annot_dir/$sumstat
        fi
     done
done



