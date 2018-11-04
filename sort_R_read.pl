#!/private/apps/bin/perl516 -w

########################################################################################
#
#     in: a file with list of mm analysis of the reads ARGV[0]
#         [and file to write the statistic ARGV[1]]
#     out: a sorted file of reads (repeat reads are cluster together) 
#
########################################################################################

use strict;

my $source_file = $ARGV[0];
my $max_hits_amount = 50; # samse -n=50
my $statistic = $ARGV[1];

open (my $f_statistic,">$statistic");
select((select($f_statistic), $|=1)[0]);

open (my $InFile,$source_file);
select((select($InFile), $|=1)[0]);

my %reads_count = ();
my %reads_sign = ();
my %reads_hit_amount = ();
my %reads_align_lines = ();

my $sign="";

my $count_orig_reads =0;
my $count_orig_R_reads =0;
my $count_orig_U =0;
my $count_reject_hit_amount =0;
my $count_reject_record=0;
my $count_new_R=0;
my $count_new_R_hits=0;
my $count_not_changed_R=0;
my $count_not_changed_R_hits=0;	
my $count_not_changed_U=0;

while (<$InFile>)
{
    chomp;
    my $align_line = $_;

    my @Fields = split;
    
    $sign="";
    $count_orig_reads ++;
    
    my $hits_amount = 0;

    foreach my $field (@Fields)
    {
		if ($field =~ m/hits:(\d+)/) 
    	{
			$hits_amount = $1;
		}
		elsif ($field =~ m/X1:i:(\d+)/)
		{
    		$hits_amount += $1;
		}
    }

    if($hits_amount > 1)
    {
	$count_orig_R_reads ++;
    }
    else
    {
	$count_orig_U ++;
    }

    if( exists $reads_count{$Fields[0]})
    {
	$reads_count{$Fields[0]} ++;
	$reads_hit_amount{$Fields[0]}+=$hits_amount;
	if($reads_sign{$Fields[0]}!= $Fields[1])
	{
		$sign="-";
	}
	#print "R read: $Fields[0], tot hit amount: $reads_hit_amount{$Fields[0]}, sign: $reads_sign{$Fields[0]}, new sign: $sign\n";	 
	    
    }
    
    else
    {
	$reads_count{$Fields[0]}=1;
	$reads_sign{$Fields[0]}=$Fields[1]; 
	$reads_hit_amount{$Fields[0]}=$hits_amount;
	$reads_align_lines{$Fields[0]}="$align_line\n";
    }

    for (0.. ($hits_amount-1)) 
    {
	my $hit_line = <$InFile>;
	chomp $hit_line;
	my $new_hit_line="$hit_line\t$sign\n";
	$reads_align_lines{$Fields[0]}.=$new_hit_line;
	#print "line: $hit_line\nnew line: $new_hit_line hash rec: $reads_align_lines{$Fields[0]}";
    }
}
    
foreach my $read_id (keys %reads_count)
{
	if($reads_hit_amount{$read_id}>($max_hits_amount+1)) # tot hit >51
	{
		$count_reject_hit_amount ++;
		$count_reject_record+=$reads_count{$read_id};
		#print "reject read: $read_id, tot hit amount: $reads_hit_amount{$read_id} read count: $reads_count{$read_id}\n";
		next;
	}

	if($reads_count{$read_id}>1)
	{
		$count_new_R++;
		$count_new_R_hits+=$reads_hit_amount{$read_id};

		
		#$reads_align_lines{$read_id} =~ s/X0:i:\d+\tX1:i:\d+/X0:i:$reads_hit_amount{$read_id}\tX1:i:0/;
		$reads_align_lines{$read_id} =~ s/hits:\d+\n/hits:$reads_hit_amount{$read_id}\n/;
		#print "R read: $read_id, tot hit amount: $reads_hit_amount{$read_id} read count: $reads_count{$read_id}\n";
		#print "$reads_align_lines{$read_id}";


	}
	else
	{
		if($reads_hit_amount{$read_id}>1)
		{
			$count_not_changed_R++;
			$count_not_changed_R_hits+=$reads_hit_amount{$read_id};

		}
		else
		{		
			$count_not_changed_U ++;
		}
	}

	print "$reads_align_lines{$read_id}";

}

=pod
print $f_statistic "\nStatistic of sort repeat reads:\n";
print $f_statistic "count_orig_reads\tcount_orig_R_reads\tcount_orig_U\tcount_reject_hit_amount\tcount_reject_record\tcount_new_R\tcount_new_R_hits\tcount_not_changed_R\tcount_not_changed_R_hits\tcount_not_changed_U\n";
print $f_statistic "$count_orig_reads\t$count_orig_R_reads\t$count_orig_U\t$count_reject_hit_amount\t$count_reject_record\t$count_new_R\t$count_new_R_hits\t$count_not_changed_R\t$count_not_changed_R_hits\t$count_not_changed_U\n";
=cut


my $tot_R = $count_new_R+$count_not_changed_R;
my $tot_U = $count_not_changed_U;
print $f_statistic "\nStatistic of Analyse MM and sort repeat reads:\n";
print $f_statistic "tot Repeat reads\ttot Uniq reads\n";
print $f_statistic "$tot_R\t$tot_U\n";
