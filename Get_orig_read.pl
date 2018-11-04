#!/private/apps/bin/perl516 -w

########################################################################################
#
#     in: sam file of alignments of transformed reads, fastq file of the original reads
#     out: sam file contain: header, just the mapped reads with the original sequence
#     Pay attention that readID (optional begins with @) in fastq file must be the same as readID at sam file
#
######################################################################################## 

use strict;

my $sam_file = $ARGV[0];
my $fastq_file = $ARGV[1];
my $orig_seq = "";

my %reads = ();

open (my $InFastFile,$fastq_file);

while (<$InFastFile>)
{
    # retrieve the 4 lines of each read
    my $IDLine = $_;
    my $SeqLine = <$InFastFile>;
    my $BlankLine = <$InFastFile>;
    my $BQLine = <$InFastFile>;
    
    my $ID = "";
    chomp ($IDLine);
    chomp ($SeqLine);

    # delete the prefix "@" of the ID line
   if ($IDLine =~ /^\@/) 
   {
	$ID = substr ($IDLine,1,length($IDLine)-1);
   }
   $reads{$ID} = $SeqLine;

}

close ($InFastFile);

open (my $InSamFile,$sam_file);
 
while (<$InSamFile>) 
{
    chomp;
    
    if (m/^@/) # header line
    {
	print "$_\n";
	next ;
    }
    
    #my ($qry_name, $flag, $map_ref, $map_b, $mapq, $cigar, $null1, $null2, $null3, $qry_seq, $qry_qual) = split;
   
    my @Fields = split;

    if (!@Fields)
    {
	next;
    }

    
    if ((($Fields[1]>>2) & 1) == 1) # 0x4 is set=> read unmapped
    {
	#print "Unmapped: $_\n";
	next;
    }

    if ((($Fields[1]>>4) & 1) == 1) # 0x10 is set=> read aligned to "-" strand (reverse complemented)
    {
	$orig_seq = revcomp($reads{$Fields[0]});
    }
    else # 0x10 is unset => read aligned to "+" strand
    {
	$orig_seq = $reads{$Fields[0]};
    }
    
    $Fields[9]=$orig_seq;
    my $new_line = join "\t",@Fields;
    print "$new_line\n";

}

close ($InSamFile);


sub revcomp {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}