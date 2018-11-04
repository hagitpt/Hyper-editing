# Hyper-editing
A genome-wide map of hyper-edited RNA reveals numerous new sites Hagit T. Porath1, Shai Carmi2, Erez Y. Levanon1 1 The Mina and Everard Goodman Faculty of Life Sciences, Bar-Ilan University, Ramat-Gan, 5290002, Israel 2 Department of Computer Science, Columbia University, New York, NY 10027, USA Correspondence should be addressed to E.Y.L. (Erez.Levanon@biu.ac.il)

Running of hyper_editing scripts for detection of hyper editing in RNA-Seq datasets (fastq files):

Input: fastq files to be tested for hyper editing and afasta file of genome reference.

Output: bed files of the hyper edited clusters and some statistics. For running the pipeline do as instructed bellow:

The pipeline uses the following tools: "bwa", "bam2fastx" and "samtools". Make sure that these tools are installed in your system and specify the paths to each tool, as follows: "bwa" tool at TransformIndexBWA_genome.sh, and "bwa", "bam2fastx" and "samtools" tools at config_file.sh.
If needed (i.e. running pipeline from outside folder) insert pipeline's scripts paths to the following files (in the specified lines): TransformIndexBWA_genome.sh, config_file.sh, run_hyper_editing.sh, analyse_mm.sh, TransRun.sh. and also insert the path to the "unique_simple_repeats.txt" file at filter_sam.pl.
Run TransformIndexBWA_genome.sh for bwa indexing of the reference genome and transforming of the genome (and indexing of the transformed genomes). Reference genome must be in fasta format and its extension should be ".fa", the command for running TransformIndexBWA_genome.sh script: TransformIndexBWA_genome.sh prefix_of_the_genome_fasta_file (without ".fa"). Transforming and indexing (this step) should be run once for each analysed genome assembly.
Prepare the "file_list" file which specified the input fastq files and the output prefix for each of them: Each line in the "file_list" file is: full_path_and_fastq_file_name/TAB/output_prefix (If the input files are of paired-end reads, each of the paired files should appear in a separate line. Paired files have to be listed one after another).
Edit the config_file.sh according to your specifications (in each line replace the string in the "XXX" with the proper information):
source_file="file_list" :insert the full path and name of the prepared "file_list".
genome_bwa_ind="genome_prefix": insert the full path+prefix of the index genome (expected 5 files like: prefix.amb, prefix.ann, prefix.bwt, prefix.pac, prefix.sa)
genome_trans_bwa_ind="transformed_genome_prefix": insert the full path+prefix of the index of the transformed genomes (for each of the 6 transformed [a2c a2g a2t g2c t2c t2g]=>tt: 5 files: prefix.tt.amb, prefix.tt.ann, prefix.tt.bwt, prefix.tt.pac, prefix.tt.sa + 1 fasta file prefix.tt.fa => tot 36 files)
genome_fasta="genome_prefix.fa": insert the path+full_name of the fasta file of the original genome
Q="33": insert the base quality PHRED score offset (33 or 64) of the source fastq files
PE="0": 0 => for single end files, 1 => for paired end
GAP=500000: for paired end files: the maximum gap size (bp) between the pairs
dir_pre="output_folder_prefix" : full output folder path
bwa_run="1": 1 => first running of the pipeline on those source files, 0 => rerunning of the pipeline on same source files
ue_detect_args="0.05 0.6 30 0.6 0.1 0.8 0.2": detection of hyper-editing arguments (for more details see the article Methods section) Argument definitions:
Minimum of edit sites at the hyper-edited read (float- fraction of read length)
Minimum of fraction between edit sites/mm sites at the hyper-edited read (float)
Minimum sequence quality for counting site as editing event (int- range between 1-40)
Maximum fraction of same letter in editing-cluster (float- fraction of read length)
Minimum of editing-cluster length (float- fraction of read length)
Maximum of initiate index of editing cluster in the read (float- fraction of read length)
Minimum of ending index of editing cluster in the read (float- fraction of read length)
Output files and information:

The folder UEdetect./detection_arguments/ contains the bed files of the detected clusters and sites, and detailed lists of all possible editing types for each analyzed fastq file.
The folder statistic contains two files: one for general statistic on the mapping and filtering and the other contains the detection statistic.
Log file of the running. 
Good Luck!
