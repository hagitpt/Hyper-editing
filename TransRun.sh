
genome_ind=$1
source_file=$2
Trans_run_dir=$3
pre_out=$4
statistic=$5

bwa_aln_soft=$6 # full path + tool name
SamToolsView_soft="$7 view" # full path + tool name

transform_script="fastq_transform.pl" # if needed insert the proper path before the script name
get_orig_script="Get_orig_read.pl" # if needed insert the proper path before the script name
filter_sam_script="filter_sam.pl" # if needed insert the proper path before the script name


Trans_dir="$Trans_run_dir/TransFiles" 
bam_files="$Trans_run_dir/bamFiles" 
reject_files="$Trans_run_dir/Reject" 
filter_statistic="$reject_files/statistic" 


mkdir $Trans_run_dir $Trans_dir $bam_files $reject_files 

#names of alignment files:

align_full="full.$pre_out"
align_origReads="origReads.$pre_out"
align_FilterReads="FilterReads.$pre_out"

#transform orig fastq file

$transform_script a g $source_file > $Trans_dir/a2g.$pre_out.fastq
$transform_script t c $source_file > $Trans_dir/t2c.$pre_out.fastq
$transform_script a c $source_file > $Trans_dir/a2c.$pre_out.fastq
$transform_script t g $source_file > $Trans_dir/t2g.$pre_out.fastq
$transform_script g c $source_file > $Trans_dir/g2c.$pre_out.fastq
$transform_script c g $source_file > $Trans_dir/c2g.$pre_out.fastq
$transform_script a t $source_file > $Trans_dir/a2t.$pre_out.fastq
$transform_script t a $source_file > $Trans_dir/t2a.$pre_out.fastq

echo -e "\nTrans Run" >> $statistic

#run bwa on transform files

#a2g++


echo "statistic of a2g++" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.a2g $Trans_dir/a2g.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.a2g $bam_files/temp.sai $Trans_dir/a2g.$pre_out.fastq > $bam_files/temp.sam 

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2g++.$align_full.bam  ###


#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2g++.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2g++.$align_origReads.sam -o $bam_files/a2g++.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/a2g++.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2g++.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2g++.$align_FilterReads.sam -o $bam_files/a2g++.$align_FilterReads.bam 
cat $filter_statistic.a2g++ >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?




#a2g+-

echo "statistic of a2g+-" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.a2g $Trans_dir/t2c.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.a2g $bam_files/temp.sai $Trans_dir/t2c.$pre_out.fastq > $bam_files/temp.sam 

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2g+-.$align_full.bam ###

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2g+-.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2g+-.$align_origReads.sam -o $bam_files/a2g+-.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/a2g+-.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2g+-.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2g+-.$align_FilterReads.sam -o $bam_files/a2g+-.$align_FilterReads.bam 
cat $filter_statistic.a2g+- >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?


#a2g-+

echo "statistic of a2g-+" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.t2c $Trans_dir/a2g.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.t2c $bam_files/temp.sai $Trans_dir/a2g.$pre_out.fastq > $bam_files/temp.sam 

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2g-+.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2g-+.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2g-+.$align_origReads.sam -o $bam_files/a2g-+.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/a2g-+.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2g-+.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2g-+.$align_FilterReads.sam -o $bam_files/a2g-+.$align_FilterReads.bam 
cat $filter_statistic.a2g-+ >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?



#a2g--

echo "statistic of a2g--" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.t2c $Trans_dir/t2c.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.t2c $bam_files/temp.sai $Trans_dir/t2c.$pre_out.fastq > $bam_files/temp.sam 

#$$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2g--.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2g--.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2g--.$align_origReads.sam -o $bam_files/a2g--.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/a2g--.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2g--.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2g--.$align_FilterReads.sam -o $bam_files/a2g--.$align_FilterReads.bam 
cat $filter_statistic.a2g-- >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?



#a2c++

echo "statistic of a2c++" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.a2c $Trans_dir/a2c.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.a2c $bam_files/temp.sai $Trans_dir/a2c.$pre_out.fastq > $bam_files/temp.sam 

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2c++.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2c++.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2c++.$align_origReads.sam -o $bam_files/a2c++.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/a2c++.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2c++.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2c++.$align_FilterReads.sam -o $bam_files/a2c++.$align_FilterReads.bam 
cat $filter_statistic.a2c++ >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?


#a2c+-

echo "statistic of a2c+-" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.a2c $Trans_dir/t2g.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.a2c $bam_files/temp.sai $Trans_dir/t2g.$pre_out.fastq > $bam_files/temp.sam 

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2c+-.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2c+-.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2c+-.$align_origReads.sam -o $bam_files/a2c+-.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/a2c+-.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2c+-.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2c+-.$align_FilterReads.sam -o $bam_files/a2c+-.$align_FilterReads.bam 
cat $filter_statistic.a2c+- >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?


#a2c-+

echo "statistic of a2c-+" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.t2g $Trans_dir/a2c.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.t2g $bam_files/temp.sai $Trans_dir/a2c.$pre_out.fastq > $bam_files/temp.sam 

#$$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2c-+.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2c-+.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2c-+.$align_origReads.sam -o $bam_files/a2c-+.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/a2c-+.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2c-+.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2c-+.$align_FilterReads.sam -o $bam_files/a2c-+.$align_FilterReads.bam 
cat $filter_statistic.a2c-+ >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?



#a2c--

echo "statistic of a2c--" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.t2g $Trans_dir/t2g.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.t2g $bam_files/temp.sai $Trans_dir/t2g.$pre_out.fastq > $bam_files/temp.sam 

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2c--.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2c--.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2c--.$align_origReads.sam -o $bam_files/a2c--.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/a2c--.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2c--.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2c--.$align_FilterReads.sam -o $bam_files/a2c--.$align_FilterReads.bam 
cat $filter_statistic.a2c-- >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?



#g2c++

echo "statistic of g2c++" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.g2c $Trans_dir/g2c.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.g2c $bam_files/temp.sai $Trans_dir/g2c.$pre_out.fastq > $bam_files/temp.sam 

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/g2c++.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 


#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/g2c++.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/g2c++.$align_origReads.sam -o $bam_files/g2c++.$align_origReads.bam 


#filtering the reads
$filter_sam_script $bam_files/g2c++.$align_origReads.sam $reject_files $filter_statistic > $bam_files/g2c++.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/g2c++.$align_FilterReads.sam -o $bam_files/g2c++.$align_FilterReads.bam 
cat $filter_statistic.g2c++ >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?


#g2c+-

echo "statistic of g2c+-" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.g2c $Trans_dir/c2g.$pre_out.fastq -f $bam_files/temp.sai

$bwa_aln_soft samse -n 50 $genome_ind.g2c $bam_files/temp.sai $Trans_dir/c2g.$pre_out.fastq > $bam_files/temp.sam

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/g2c+-.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 

#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/g2c+-.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/g2c+-.$align_origReads.sam -o $bam_files/g2c+-.$align_origReads.bam 

#filtering the reads
$filter_sam_script $bam_files/g2c+-.$align_origReads.sam $reject_files $filter_statistic > $bam_files/g2c+-.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/g2c+-.$align_FilterReads.sam -o $bam_files/g2c+-.$align_FilterReads.bam 
cat $filter_statistic.g2c+- >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?


#a2t++

echo "statistic of a2t++" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.a2t $Trans_dir/a2t.$pre_out.fastq -f $bam_files/temp.sai

$bwa_aln_soft samse -n 50 $genome_ind.a2t $bam_files/temp.sai $Trans_dir/a2t.$pre_out.fastq > $bam_files/temp.sam

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2t++.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 

#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2t++.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2t++.$align_origReads.sam -o $bam_files/a2t++.$align_origReads.bam 

#filtering the reads
$filter_sam_script $bam_files/a2t++.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2t++.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2t++.$align_FilterReads.sam -o $bam_files/a2t++.$align_FilterReads.bam
cat $filter_statistic.a2t++ >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?


#a2t+-

echo "statistic of a2t+-" >> $statistic

$bwa_aln_soft aln -t 5 -n 2 -o 0 -N $genome_ind.a2t $Trans_dir/t2a.$pre_out.fastq -f $bam_files/temp.sai 

$bwa_aln_soft samse -n 50 $genome_ind.a2t $bam_files/temp.sai $Trans_dir/t2a.$pre_out.fastq > $bam_files/temp.sam

#$SamToolsView_soft -bS $bam_files/temp.sam -o $bam_files/a2t+-.$align_full.bam 

#remove unmapped and map to chrUn and random reads
grep -v "chrUn" $bam_files/temp.sam | grep -v "_random" | grep -v "chrUextra" > $bam_files/temp2.sam
$SamToolsView_soft -F 4 -S -h -o $bam_files/temp3.sam $bam_files/temp2.sam 

#retrieve original reads
$get_orig_script $bam_files/temp3.sam $source_file > $bam_files/a2t+-.$align_origReads.sam
$SamToolsView_soft -bS $bam_files/a2t+-.$align_origReads.sam -o $bam_files/a2t+-.$align_origReads.bam 

#filtering the reads
$filter_sam_script $bam_files/a2t+-.$align_origReads.sam $reject_files $filter_statistic > $bam_files/a2t+-.$align_FilterReads.sam
$SamToolsView_soft -bS $bam_files/a2t+-.$align_FilterReads.sam -o $bam_files/a2t+-.$align_FilterReads.bam
cat $filter_statistic.a2t+- >> $statistic
echo "" >> $statistic

rm $bam_files/*.sa?