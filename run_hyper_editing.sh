

#command line: $run_he_script $genome_bwa_ind $genome_trans_bwa_ind $genome_fasta $Q $PE $GAP $dir_pre $bwa_run $ue_detect_args $source_dir $bwa_aln_soft $bwa_mem_soft $SamToFastq_soft $SamTools_soft
##########################################################################################

bwa_aln_soft=${17} # full path + tool name
bwa_mem_soft=${18} # full path + tool name
SamToFastq_soft=${19} # full path + tool name
SamTools_soft=${20} # full path + tool name

source_dir=${16}  # line per each src_fastq_file: path/full_name[tab]out_pre #### for PE files each path+full_name will appear in different line with same out_pre for both files!!! 
genome_bwa_ind=$1 # path+prefix of the index genome expected 5 files like: prefix.amb, prefix.ann, prefix.bwt, prefix.pac, prefix.sa
genome_trans_bwa_ind=$2 # path+prefix of the index transformed genome: for each of the 6 transformed (a2c a2g a2t g2c t2c t2g)=>tt: 5 files: prefix.tt.amb, prefix.tt.ann, prefix.tt.bwt, prefix.tt.pac, prefix.tt.sa + 1 fasta file prefix.tt.fa => tot 36 files
genome_fasta=$3 # path+full_name of the fasta file of the original genome
Q=$4 # PHRED score of base quality (33 or 64)
PE=$5 # 0 => single end, 1 => paired end
GAP=$6 # gap max size between the pairs
dir_pre=$7 # path+full_name of the output directory
bwa_run=$8 # 0 => pipeline will not run BWA mapping if analyseMM file exists, 1 => first run of the source files, pipeline will begin with BWA mapping 
ue_detect_args="$9 ${10} ${11} ${12} ${13} ${14} ${15}"  # args meaning: -Min of edit sites at Ultra-Edit read -Min fraction of edit sites/mm sites -Min sequence quality for counting editing event -Max fraction of same letter in cluster -Min of cluster length -Max initiate index of cluster -Min ending index of cluster


##########################################################################################

#scripts
unmapped_script="pre_unmapped.sh" # if needed insert the proper path before the script name
TransRun_script="TransRun.sh" # if needed insert the proper path before the script name
analyse_mm_script="analyse_mm.sh" # if needed insert the proper path before the script name
detectUE_script="detect_ue.pl" # if needed insert the proper path before the script name

#files and directories
unmap_dir="$dir_pre/unMap"
Trans_run_dir="$dir_pre/TransRun" 
analyseMM_dir="$dir_pre/AnalyseMM" 
stat_files="$dir_pre/statistic"
general_stat="$stat_files/general"

##########################################################################################

args=`echo $ue_detect_args | tr " " "_"`

if [ $PE == 1 ]; then
	args="PE_$args" 
elif [ $PE == 0 ]; then
	args="SE_$args" 
fi
arg_stat_det="$stat_files/detect.$args"
log_file="$dir_pre/log.$args"
UE_detect_dir_pre="$dir_pre/UEdetect.$args"

if [ $bwa_run == 1 ]; then
mkdir $dir_pre
mkdir $stat_files
> $general_stat # clear existing file or open new if not exist
fi

mkdir $UE_detect_dir_pre 
mkdir $Trans_run_dir 
mkdir $unmap_dir  
mkdir $analyseMM_dir 

echo -e "\n----------\n" >> $log_file 
echo -e "\n----------\n" >> $arg_stat_det 


PE_1=1

while read line; do   

file_path=`echo $line | awk '{print $1}'`   
out_name=`echo $line | awk '{print $2}'`  

if [ $PE == 1 ] && [ $PE_1 == 1 ]; then
	out_pre="$out_name-1"     
elif [ $PE == 1 ] && [ $PE_1 == 2 ]; then
	out_pre="$out_name-2"     
elif [ $PE == 0 ]; then
	out_pre="$out_name"
fi

if [ $bwa_run == 1 ] || ([ $bwa_run == 0 ] && [ ! -f "$analyseMM_dir/$out_pre.analyseMM" ]); then 
echo "Start run $file_path in `date`" >> $log_file

$unmapped_script $file_path $genome_bwa_ind $Q $PE $unmap_dir $dir_pre/$out_pre $general_stat $bwa_aln_soft $bwa_mem_soft $SamToFastq_soft $SamTools_soft   
$TransRun_script $genome_trans_bwa_ind $unmap_dir/$out_pre.mem.um.fastq $Trans_run_dir/$out_pre $out_pre $general_stat $bwa_aln_soft $SamTools_soft
$analyse_mm_script $genome_fasta $Trans_run_dir/$out_pre/bamFiles $analyseMM_dir/$out_pre $general_stat $SamTools_soft
echo "End mapping and mm analysis of $file_path in `date`" >> $log_file
fi

if [ $PE == 1 ] && [ $PE_1 == 1 ]; then
	first_out_pre=$out_pre
	PE_1=2
elif [ $PE == 1 ] && [ $PE_1 == 2 ]; then
	echo "Start detection UE of $file_path pair in `date`" >> $log_file
	$detectUE_script $analyseMM_dir/$first_out_pre.analyseMM $UE_detect_dir_pre/$first_out_pre $arg_stat_det $ue_detect_args $SamTools_soft $PE $GAP $unmap_dir/$out_pre.aln.bam $unmap_dir/$out_pre.mem.bam 
	$detectUE_script $analyseMM_dir/$out_pre.analyseMM $UE_detect_dir_pre/$out_pre $arg_stat_det $ue_detect_args $SamTools_soft $PE $GAP $unmap_dir/$first_out_pre.aln.bam $unmap_dir/$first_out_pre.mem.bam 
	echo "End detection UE of $file_path pair in `date`" >> $log_file
	PE_1=1
elif [ $PE == 0 ]; then
	echo "Start detection UE of $file_path in `date`" >> $log_file
	$detectUE_script $analyseMM_dir/$out_pre.analyseMM $UE_detect_dir_pre/$out_pre $arg_stat_det $ue_detect_args $SamTools_soft $PE  
	echo "End detection UE of $file_path in `date`" >> $log_file
fi

done<$source_dir

##### removing the intermediates files and folders
rm -r $unmap_dir $Trans_run_dir $analyseMM_dir



