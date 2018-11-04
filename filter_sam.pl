#!/usr/bin/perl -w

########################################################################################
#
#     in: [0]sam file to be filtered (base quality is encode by PHRED-33) [1]folder of the reject_files [2] file for filtering statistic  ###
#     out: filtered sam file
#
######################################################################################## 

use strict;


my $sim_rep_file = "unique_simple_repeats.txt"; # if needed insert the proper path to the file

my $max_letter = 0.6;
my $min_letter = 0.1;
my $max_others = 0.1;
my $repeat_len = 10; # for single repeat_len?repeat_len, and for other repeats repeat_len

my $qual_offset = 33;
my $baseQ_cutoff = 0.1; # for calculating the amount of base qualities to ignore 
my $min_baseQ_avg = 25;

my $sam_file = $ARGV[0];
my @path = split ("/",$ARGV[0]);
my $full_align = $path[scalar(@path)-1];
my @path2 = split (/\./,$full_align);
my $align = $path2[0];

my $reject_files = $ARGV[1];
my $statistic_file = $ARGV[2];
my $count_orig = 0;
my $count_invalid = 0;
my $count_invalid_rep = 0;
my $count_invalid_con = 0;
my $count_invalid_qual = 0;

open (my $f_rep, ">$reject_files/inv_rep.$align");
select((select($f_rep), $|=1)[0]);

open (my $f_con, ">$reject_files/inv_con.$align");
select((select($f_con), $|=1)[0]);

open (my $f_qual, ">$reject_files/inv_qual.$align");
select((select($f_qual), $|=1)[0]);

open (my $f_log, ">$statistic_file.$align");
select((select($f_log), $|=1)[0]);

open (my $f, $sim_rep_file);
my $repeats_str = '';

while (<$f>)
{
    my $rep_line = uc $_;
    chomp $rep_line;
    if (length($rep_line) > 1)
    {
	$repeats_str .= "($rep_line)\{$repeat_len,\}|";
    }
    if (length($rep_line) == 1)
    {
	$repeats_str .= "($rep_line)\{$repeat_len,\}\\w?($rep_line)\{$repeat_len,\}|";
    }

}

chop $repeats_str;

close($f);

open (my $InSamFile,$sam_file);


while (<$InSamFile>)
{
    chomp;
    
    my $line = $_;

    if (m/^@/) # header line
    {
	print "$line\n";
	next ;
    }

    $count_orig += 1;
   
    #my @Fields = split ("\t",$line);
    
    my @Fields = split;
    
    
    if (!@Fields)
    {
	next;
    }
    
# Remove reads with repeats
    if ($Fields[9] =~ /$repeats_str/o)	
    {
       ##print "$SeqLine is repeat $&\n";
       $count_invalid += 1;
       $count_invalid_rep += 1;
       print $f_rep "$line\n";
       next;
    }

# Remove reads due to bases content 
# Remove reads due to low base qualities

    my %nuc_map = ('A' => 0, 'C' => 1, 'G' => 2, 'T' =>3);
    my @bases = split(//,$Fields[9]);
    my @hist = (0,0,0,0);
    my $num_others = 0;
    my $tot = scalar(@bases);

    my @qualities = split(//,$Fields[10]);
    my @known_qual = ();

    for my $i (0 .. ($tot-1))
    {
		if ($bases[$i] eq 'N')		
		{
			$num_others++;
		}
		else
		{
			$hist[$nuc_map{$bases[$i]}]++;
			push (@known_qual , (ord($qualities[$i])-$qual_offset));
		}
	}		
		
	my $sum_known = $tot - $num_others;
	my @hist_s = sort {$a <=> $b} @hist;
	my $min = $hist_s[0] / $sum_known;
	my $max = $hist_s[3] / $sum_known;
	
	if (($min < $min_letter) or ($max > $max_letter) or (($num_others/$tot) > $max_others))
	{
	    #print "$Fields[9] is content invalid\n";
	    #print "A= $hist[0],C= $hist[1],G= $hist[2],T= $hist[3], min- $min, max- $max, other- $num_others\n";
	    #$count_invalid += 1;
	    #$count_invalid_con += 1;
	    $count_invalid += 1;
	    $count_invalid_con += 1;
	    print $f_con "$line\n";
	    next;
	}
    
    my $num2cutoff = int($baseQ_cutoff*$sum_known);
    my @known_qual_s = sort {$a <=> $b} @known_qual; 
    my $sum_qual=0;

    my $qual_str ="";
    for my $j ($num2cutoff .. ($sum_known-1))
    {
	$sum_qual += $known_qual_s[$j];
	$qual_str .= " $known_qual_s[$j]";
    }
    
    my $avg_qual = $sum_qual/($sum_known-$num2cutoff);
    
    if ($avg_qual < $min_baseQ_avg)
    {
	
	#print "$SeqLine\n$BQLine invalid quality\n num2cutoff= $num2cutoff, avg_qual= $avg_qual, sum_qual= $sum_qual\n";
	#print "qulities: @known_qual_s\n";
	#print "after cutoff: $qual_str\n";

	$count_invalid += 1;
	$count_invalid_qual += 1;
	print $f_qual "$line\n";
	next;
    }
	    
    print "$line\n"; 
}
#print "Rep invalid = $count_invalid_rep, Con invalid = $count_invalid_con, Total invalid = $count_invalid\n";

print $f_log "map to chr\t$count_orig\nFilter: Tot invalid\t$count_invalid\nInvalid repeat\t$count_invalid_rep\nInvalid content\t$count_invalid_con\nInvalid qualities\t$count_invalid_qual\n";

close ($InSamFile);
close ($f_rep);
close ($f_con);
close ($f_qual);
close ($f_log);

