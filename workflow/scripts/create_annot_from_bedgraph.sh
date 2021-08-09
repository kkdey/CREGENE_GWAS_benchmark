module load gcc/6.2.0
module load conda2
source activate ldsc
bed_cell=snakemake@input$bed_dir
bim_cell=snakemake@params$bimpath
annot_cell=snakemake@output$annot_dir

IFS="
"

TASKFILE=/n/groups/price/kushal/ENCODE/data/pops.txt
for line in `ls $bed_cell | awk '{print $1}' | sort | uniq`;
do
    name=`echo $line | awk '{print $1}'`
    if [ ! -d $annot_cell/$name ]
    then
	    mkdir $annot_cell/$name
    fi
    for bedline in `ls $bed_cell/$name/ | cat | sort | uniq | cut -f 1 -d '.'`;
    do
    	bedname=`echo $bedline | awk '{print $1}'`
    	if [ ! -d $annot_cell/$name/$bedname ]
    	then
    	    mkdir $annot_cell/$name/$bedname
    	fi
    	if [ ! -f $annot_cell/$name/$bedname/$bedname.22.annot.gz ]
    	then
    	    python  make_annot_combine_from_bedgraph.py --bedname $bedname --bedfile_path $bed_cell/$name --bimfile_path $bim_cell --annot_path $annot_cell/$name/$bedname
    	fi
    done
done
