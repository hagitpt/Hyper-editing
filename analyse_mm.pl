#!/private/apps/bin/perl516 -w

########################################################################################
#
#     in: reference genome (fasta format) ARGV[0]
#         bam file to be analysed for mm (base quality is encode by PHRED-33) ARGV[1]
#         [and path for the samtools tool ARGV[2]]
#     out: list of mm analysis of the reads
#
########################################################################################

use strict;

my $SamToolsView_soft="$ARGV[2] view"; # full path + tool name

sub analyse_mm 
{
    my ($comb,$ref_seq,$read) = @_;

    my %reverse_mm_types = ('A2G' => 'T2C', 'A2C' => 'T2G');

    my $mm_str = "";

    my @rna_seq = split(//, uc $read);
    my @dna_seq = split(//, uc $ref_seq);

    my $mm_count = 0;
    my $expected_edit = "";
    my $edit_type = "";
    my $edit_count=0;
    my $edit_count1=0; 
    my $edit_count2=0;
    my $edit_loc = "";
    my @edit_loc1 = ();
    my @edit_loc2 = ();
    my @mm_loc = ();
    my $mm_locs = "";


    if ($comb =~ /(\S{3})-\S{1}/) # comb options are: A2G++,A2G+-,A2G-+,A2G--,A2C++,A2C+-,A2C-+,A2C--,G2C++,G2C+-,A2T++,A2T+-
    {
	$expected_edit = $reverse_mm_types{$1};
	#print "comb is $comb, expcted is $expected_edit\n"; 
    }
    else
    {
	$expected_edit = substr ($comb,0,3);
	#print "comb is $comb, expcted is $expected_edit\n"; 
    }

    
    for my $i (0..(length($read)-1))
    {				
	if ($rna_seq[$i] ne $dna_seq[$i])
	{
	    #print "mm at $i $dna_seq[$i] to $rna_seq[$i]\n";
	    $mm_count ++;
	    push @mm_loc,$i;
	}
	elsif (($rna_seq[$i] eq "N") and ($dna_seq[$i] eq "N"))
	{
	    #print "NN at $i $dna_seq[$i] to $rna_seq[$i]\n";
	    $mm_count ++;
	    push @mm_loc,$i;
	    next;
	}

	if ($dna_seq[$i] eq (substr ($expected_edit,0,1))) 
	{
	    if ($rna_seq[$i] eq (substr ($expected_edit,2,1)))
	    {
		$edit_count1++;
		push @edit_loc1,$i;
	    }

	}
	elsif  ($dna_seq[$i] eq (substr ($expected_edit,2,1)))
	{
	    if ($rna_seq[$i] eq (substr ($expected_edit,0,1)))
	    {
		$edit_count2++;
		push @edit_loc2,$i;
	    }

	}
	
    }

    if ($edit_count2 > $edit_count1)
    {
	$edit_count = $edit_count2;
	$edit_type = reverse ($expected_edit);
	$edit_loc = join(';',@edit_loc2);

	foreach my $mm (@mm_loc)
	{
		if (!grep(/^$mm$/, @edit_loc2))
		{
			$mm_locs .= "$mm;";
		}
	}
    }
    else
    {
	$edit_count = $edit_count1;
	$edit_type = $expected_edit;
	$edit_loc = join(';',@edit_loc1);
	foreach my $mm (@mm_loc)
	{
		if (!grep(/^$mm$/, @edit_loc1))
		{
			$mm_locs .= "$mm;";
		}
	}

    }
	if ($mm_locs)
	{
		chop $mm_locs;
	}
	else
	{
		$mm_locs = ";";
	}
    $mm_str = "$mm_count\t$edit_type\t$edit_count\t$edit_loc\t$mm_locs";
    return $mm_str;
}


my $hits_amount = 50; # samse -n=50

my $ref_gen = $ARGV[0];
my $sam_path = $ARGV[1];
my @path = split ("/",$ARGV[1]);
my $sam_name = $path[scalar(@path)-1];
my $comb = uc (substr($sam_name,0,5));

my $count_orig_map = 0;
my $count_reject_hit_amount = 0;
my $count_analysed = 0;
my $count_sub_optimals = 0;
my $count_w_suboptimals = 0;
my $count_U =0;

my %ref_hash = (); # map of sequence names to sequence strings
my $name = "";
my $seq = "";

my $ref_seq;
my $mm_str;


# retrieve reference fasta file

open (my $InRef,$ref_gen);


while (<$InRef>) 
{
    chomp;
    
    if (m/>(\S+)/) # header line
    {
	$ref_hash{$name} = $seq if ($name ne "");
	$name = $1;
	$seq = "";
    }
    else 
    {
	$seq .= $_;
    }
}
$ref_hash{$name} = $seq if ($name ne "");

close ($InRef);




open (my $InSamFile,"$SamToolsView_soft $sam_path |");


while (<$InSamFile>)
{
    chomp;
    my $line = $_;

    if (m/^@/) # header line
    {
	next ;
    }
    
    my @Fields = split;

    if (!@Fields)
    {
	next;
    }

    
    if ((($Fields[1]>>2) & 1) == 1) # 0x4 is set=> read unmapped
    {
	next;
    }

    $count_orig_map ++;

    my $hits = 0;
    foreach my $field (@Fields)
    {
	if ($field =~ m/X0:i:(\d+)/)
    	{
		$hits += $1;
	}
	elsif ($field =~ m/X1:i:(\d+)/)
	{
    		$hits += $1;
	}
    }

    if ($hits > ($hits_amount+1)) # reject alignments with more then ($hits_amount+1) hits (optimal and suboptimal)
    {
	$count_reject_hit_amount += 1;
	#print "reject: $Fields[13],$Fields[14]\t\t$hits > 51,$count_reject_hit_amount\n";
	next;
    }

    $count_analysed ++;

    $ref_seq = uc (substr($ref_hash{$Fields[2]},$Fields[3]-1,length($Fields[9])));
    #print "\nlocation: $Fields[2],$Fields[3]\n\nDNA: $ref_seq\nRNA: $Fields[9]\n";
    $mm_str = analyse_mm($comb,$ref_seq,$Fields[9]);

    #addition of 1 base at each end of the ref seq, in case editing sites are located at the ends
    if($Fields[3]-1 == 0) #ref_seq first base is the first base of chr
    {
	$ref_seq = "N".$ref_seq;
    }
    else
    {
	$ref_seq = (uc(substr($ref_hash{$Fields[2]},$Fields[3]-2,1))).$ref_seq;
    }
    if (length($ref_hash{$Fields[2]}) < ($Fields[3]+length($Fields[9])))#ref_seq last base is the last base of chr
    {
	$ref_seq = $ref_seq."N";
    } 
    else
    {
	$ref_seq = $ref_seq.(uc(substr($ref_hash{$Fields[2]},$Fields[3]+length($Fields[9])-1,1)));
    }

    ##$ref_seq= uc (substr($ref_hash{$Fields[2]},$Fields[3]-2,length($Fields[9])+2));
    my $hits2print = "$Fields[0]\t$Fields[2]\t$Fields[3]\t$ref_seq\t$mm_str\n";
     

    if ($hits > 1)
    {	
	my $XA = "";
	foreach my $field (@Fields)
    	{
		if ($field =~ m/^XA:Z:/)
    		{
			$XA = $field;	
		}
	}

	$count_w_suboptimals ++;
	
	$XA =~ s/^XA:Z://;
	my @sub_optimals = split (";",$XA);
	#print "sub optimal @sub_optimals\n\n";

	foreach my $s (@sub_optimals)
	{ 
	    $count_sub_optimals++;
	    my @sub_optimal = split (",",$s);

	    $sub_optimal[1] =~ s/^[+-]//;
	    $ref_seq = uc (substr($ref_hash{$sub_optimal[0]},$sub_optimal[1]-1,length($Fields[9])));
		if (length($ref_seq) < length($Fields[9]))
		{
			$hits --;		
			next;
		}
	    $mm_str = analyse_mm($comb,$ref_seq,$Fields[9]);

	    #addition of 1 base at each end of the ref seq, in case editing sites are located at the ends
	    if($sub_optimal[1]-1 == 0) #ref_seq first base is the first base of chr
	    {
		$ref_seq = "N".$ref_seq;
	    }
	    else
	    {
				$ref_seq = (uc(substr($ref_hash{$sub_optimal[0]},$sub_optimal[1]-2,1))).$ref_seq; 
	    }
		if (length($ref_hash{$sub_optimal[0]}) < ($sub_optimal[1]+length($Fields[9])))#ref_seq last base is the last base of chr
	    {
		$ref_seq = $ref_seq."N";
	    } 
	    else
	    {
				$ref_seq = $ref_seq.(uc(substr($ref_hash{$sub_optimal[0]},$sub_optimal[1]+length($Fields[9])-1,1)));
	    }
	    #$ref_seq= uc (substr($ref_hash{$sub_optimal[0]},$sub_optimal[1]-2,length($Fields[9])+2));
		$hits2print = "$hits2print$Fields[0]\t$sub_optimal[0]\t$sub_optimal[1]\t$ref_seq\t$mm_str\n";
	}
    }
	
   else
	{
		$count_U ++;
	}	

	print "$line\thits:$hits\n$hits2print";

}
