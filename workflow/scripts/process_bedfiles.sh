bed_cell=snakemake@input$bed_dir
TASKFILE=/n/groups/price/kushal/ENCODE/data/bonemarrow.txt

for line in `ls $bed_cell | awk '{print $1}' | sort | uniq`;
do
   annot_name=`echo $line | awk '{print $1}'`
   input_cell=$bed_cell/$annot_name
   echo  $input_cell
   names=`ls $input_cell | cut -f 1 -d '.'`
   for name in $names
   do
       bedtools sort -i $input_cell/$name.bed > $input_cell/$name.2.bed
       bedtools merge -i $input_cell/$name.2.bed -c 4 -o max > $input_cell/$name.3.bed
       mv $input_cell/$name.3.bed $input_cell/$name.bed
       rm $input_cell/$name.2.bed
   done
done



