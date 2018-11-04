
#script genome_prefix

bwa_index_soft="bwa" # if needed insert the proper path before the tool name
Trans_script_path="" # if needed insert the proper path to the transform scripts ("a2c.pl" "a2g.pl" "a2t.pl" "g2c.pl" "t2c.pl" "t2g.pl")

genome_prefix=$1

#index the source genome:
$bwa_index_soft index -p $genome_prefix -a bwtsw $genome_prefix.fa

TransTypes="a2c a2g a2t g2c t2c t2g"

for tt in $TransTypes; do
#transform
$Trans_script_path/$tt.pl $genome_prefix.fa > $genome_prefix.$tt.fa
#indexing
$bwa_index_soft index -p $genome_prefix.$tt -a bwtsw $genome_prefix.$tt.fa
done
