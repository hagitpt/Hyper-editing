bwa_aln_soft=$8 # full path + tool name
bwa_mem_soft=$9 # full path + tool name
SamToFastq_soft="${10} --fastq" # full path + tool name
SamToolsView_soft="${11} view" # full path + tool name

source_file=$1
genome_ind=$2
Q=$3
PE=$4
pre_dir=$5
path_out_pre=$6
stat_file=$7

out_pre=`basename $path_out_pre` 

# for PHRED-64 => Q=64, for PHERED-33 => Q=33

read_count=`cat $source_file | wc -l`; read_count=$(($read_count/4))
echo -e "\n\nOriginal $source_file , out prefix: $path_out_pre\nreads in original $source_file: $read_count" >> $stat_file

if [ $Q == 64 ]; then 
$bwa_aln_soft aln -I -t 5 $genome_ind $source_file -f $pre_dir/$out_pre.sai
else if [ $Q == 33 ]; then
$bwa_aln_soft aln -t 5 $genome_ind $source_file -f $pre_dir/$out_pre.sai
fi
fi

# n=50, output of 50 suboptimal alignments. for PE matching 
$bwa_aln_soft samse -n 50 $genome_ind $pre_dir/$out_pre.sai $source_file -f $pre_dir/$out_pre.aln.sam

$SamToolsView_soft -bS $pre_dir/$out_pre.aln.sam -o $pre_dir/$out_pre.aln.bam
$SamToolsView_soft -b -f4 $pre_dir/$out_pre.aln.bam -o $pre_dir/$out_pre.aln.um.bam

$SamToFastq_soft -o $pre_dir/$out_pre.aln.um.fastq $pre_dir/$out_pre.aln.um.bam

read_count=`cat $pre_dir/$out_pre.aln.um.fastq | wc -l`; read_count=$(($read_count/4))
echo "bwa aln unmapped reads: $read_count" >> $stat_file

$bwa_mem_soft mem -M -t 5 -k 50 $genome_ind $pre_dir/$out_pre.aln.um.fastq > $pre_dir/$out_pre.mem.sam 

$SamToolsView_soft -bS $pre_dir/$out_pre.mem.sam -o $pre_dir/$out_pre.mem.bam
$SamToolsView_soft -b -f4 $pre_dir/$out_pre.mem.bam -o $pre_dir/$out_pre.mem.um.bam

$SamToFastq_soft -o $pre_dir/$out_pre.mem.um.fastq $pre_dir/$out_pre.mem.um.bam

read_count=`cat $pre_dir/$out_pre.mem.um.fastq | wc -l`; read_count=$(($read_count/4))
echo "bwa mem unmapped reads: $read_count" >> $stat_file

#if [ $PE == 1 ]; then
#$SamToolsView_soft -h -F 4 $pre_dir/$out_pre.aln.bam > $pre_dir/$out_pre.aln.map.sam 
#$SamToolsView_soft -h -F 4 $pre_dir/$out_pre.mem.bam > $pre_dir/$out_pre.mem.map.sam  
#fi

rm $pre_dir/$out_pre.sai $pre_dir/$out_pre.aln.sam $pre_dir/$out_pre.mem.sam $pre_dir/$out_pre.aln.um.bam $pre_dir/$out_pre.mem.um.bam
