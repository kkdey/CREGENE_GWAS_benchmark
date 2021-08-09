annot_cell=snakemake@output$annot_dir
bfile_path=snakemake@params$bfile_path
hapmap_path=snakemake@params$hapmap_path
ld_window=snakemake@params$ld_wind_cm

IFS="
"

TASKFILE=/n/groups/price/kushal/ENCODE/data/pops2.txt

module load conda2
source activate ldsc

for line in `ls $annot_cell | awk '{print $1}' | sort | uniq`;
do
   annot_module=`echo $line | awk '{print $1}'`
   echo $annot_cell $annot_module
   for ll in `ls $annot_cell/$annot_module | awk '{print $1}' | sort | uniq`;
   do
       annot_dir=`echo $ll | awk '{print $1}'`
       echo $annot_dir
       if [ ! -d $annot_cell/$annot_module/$annot_dir ]
       then
	   mkdir $annot_cell/$annot_module/$annot_dir
       fi
       for chrom in {1..22}
       do
       if [ ! -f $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom.l2.ldscore.gz ]
       then
           ~/.conda/envs/ldsc/bin/python $ldsc_path/ldsc.py --bfile $bfile_path/1000G.EUR.QC.$chrom --l2 --ld-wind-cm $ld_window --yes-really --annot $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom.annot.gz --print-snps $hapmap_path/hm.$chrom.snp --out $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom"
       fi
    done
  done
done

