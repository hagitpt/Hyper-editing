
genome_fasta=$1
source_dir=$2
pre_out=$3
statistic=$4

SamToolsView_soft=$5 # full path + tool name


analyse_mm_script="analyse_mm.pl" # if needed insert the proper path before the script name
sort_R_script="sort_R_read.pl" # if needed insert the proper path before the script name

> temp_file # clear existing file or open new if not exist


for bamFILE in $source_dir/*FilterReads*.bam; do

$analyse_mm_script $genome_fasta $bamFILE $SamToolsView_soft>> temp_file

done


$sort_R_script temp_file temp_statistic > $pre_out.analyseMM
cat temp_statistic >> $statistic


rm temp_file temp_statistic
