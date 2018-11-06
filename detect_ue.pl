#!/private/apps/bin/perl516 -w

########################################################################################
#
#     in: a file with list of mm analysis of the reads ARGV[0]
#	   output prefix ARGV[1]
#         a statistic file ARGV[2]
#		ARGV[10]=> samtools tool path
#         set of arguments for determining Ultra Edit read
#     out: a set of Ultra Edit reads
#
########################################################################################

use strict;
use File::Path qw(mkpath);

my $qual_offset = 33;

my $SamToolsView_soft="$ARGV[10] view"; # full path + tool name

my $source_file = $ARGV[0];
my $out_pre = $ARGV[1];
my $stat_file = $ARGV[2];

my $frac_min_edit_sites = $ARGV[3]; # 4 ##Shai 12
my $min_edit_of_mm = $ARGV[4]; # Shai 0.9
#my $min_edit_of_pot = $ARGV[]; # Shai 0.2
my $min_qual = $ARGV[5];

#my $min_letter = $ARGV[]; # 0.1
my $max_letter = $ARGV[6]; # 0.6
my $frac_min_len = $ARGV[7]; # 8
my $frac_min_int = $ARGV[8]; # 15
my $frac_min_end = $ARGV[9]; # 60


my %mm_types = ('A2G' => 'A2G', 'T2C' => 'A2G', 'G2A' => 'G2A', 'C2T' => 'G2A', 'A2C' => 'A2C', 'T2G' => 'A2C', 'C2A' => 'C2A', 'G2T' => 'C2A', 'G2C' => 'G2C', 'C2G' => 'G2C', 'A2T' => 'A2T', 'T2A' => 'A2T');
my %comp_base = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'N' => 'N');

my $count_analysed_U = 0;
my $count_selected_from_R = 0;
my $count_orig_R = 0;  

my %count_orig_wo_reg_R = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0); ###
my %count_reject_qual_sites = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_min_edit_sites = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_min_edit_of_mm = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
#my %count_reject_min_edit_of_pot = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_len = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_int = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_end = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_content = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_PE = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);

my %count_ue = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_edit_sites = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);

my %count_up_tot_bases = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_up = ('A2G' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'G2A' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'A2C' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'C2A' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'G2C' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'A2T' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},);

my %count_down_tot_bases = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_down = ('A2G' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'G2A' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'A2C' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'C2A' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'G2C' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'A2T' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},);

my %frac_up = ('A' => 0,'C' => 0,'G' => 0,'T' => 0);
my %frac_down = ('A' => 0,'C' => 0,'G' => 0,'T' => 0);

my %lens = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %align_mm = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %clust_mm = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %edit_ind = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);  ##### updated 14/04/14

my $UE_dir = "$out_pre.UE.bed_files";
my $ES_dir = "$out_pre.ES.bed_files";

mkpath ($UE_dir);
mkpath ($ES_dir);

open (my $f_statistic,">>$stat_file");
select((select($f_statistic), $|=1)[0]);

open (my $ue_list,">$out_pre.UE.list");
select((select($ue_list), $|=1)[0]);

open (my $ue_det,">$out_pre.UE.Details");
select((select($ue_det), $|=1)[0]);

open (my $InFile,$source_file);
select((select($InFile), $|=1)[0]);

my %bed_files = ();
my %es_bed_files = ();

foreach my $edit_type (keys %count_ue)
{
 
    my $bed_file =">$UE_dir/$edit_type.bed";

    my $es_bed_file =">$ES_dir/$edit_type.bed";

    open ($bed_files{$edit_type},$bed_file);

    open ($es_bed_files{$edit_type},$es_bed_file);

}

# Hash for the paired reads mapping
my %map_reads = ();

#PE=1 paired end files, PE=0 single end file
my $PE = $ARGV[11];
my $max_GAP;
my $map_file1;
my $map_file2;

if ($PE == 1) 
{ 

    $max_GAP = $ARGV[12];
    $map_file1 = $ARGV[13];
    $map_file2 = $ARGV[14];


    open (my $InMapFile1,"$SamToolsView_soft $map_file1 |");
    select((select($InMapFile1), $|=1)[0]);
    
    open (my $InMapFile2,"$SamToolsView_soft $map_file2 |");
    select((select($InMapFile2), $|=1)[0]);

    while (<$InMapFile1>) 
    {
	chomp;
	
	if (m/^@/) # header line
	{
	    next ;
	}

	my $align_line = $_;
	my @F = split;
	
	if ((($F[1]>>2) & 1) == 1) # 0x4 is set=> read unmapped 
    {
		next;
    }
	  
	if (exists $map_reads{$F[0]})
	{
	    
	    $map_reads{$F[0]}.=$align_line;
	}
	
	else
	{
	    $map_reads{$F[0]} = $align_line;
	}
    }
    

    while (<$InMapFile2>) 
    {
	chomp;
	
	if (m/^@/) # header line
	{
	    next ;
	}

	my $align_line = $_;
	my @F = split;
  
	if ((($F[1]>>2) & 1) == 1) # 0x4 is set=> read unmapped 
    {
		next;
    }
	if (exists $map_reads{$F[0]})
	{
	    
	    $map_reads{$F[0]}.=$align_line;
	}
	
	else
	{
	    $map_reads{$F[0]} = $align_line;
	}
    }

}


while (<$InFile>)
{
    chomp;
    my $align_line = $_;
	
    my @Fields = split;

    my $hits_amount = 0;

     foreach my $field (@Fields)
    {
		if ($field =~ m/hits:(\d+)/)

    	{
			$hits_amount = $1;
		}
	} 
    
    my %hits_ed_frac = ();
    my $hit;
    my ($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,$edit_count,$edit_loc,$mm_loc,$read_sign);
    
    for (0.. ($hits_amount-1)) 
    {
	$hit=<$InFile>;
	chomp $hit;
	($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,$edit_count,$edit_loc,$mm_loc,$read_sign) = split ("\t",$hit);
	
	$hits_ed_frac{$hit}=$edit_count/$mm_count;
    }
    
    if ($hits_amount == 1)
    {
	$count_analysed_U ++;
    }
    
    else
    {
	$count_orig_R ++;
	
	my @sorted_id = sort {$hits_ed_frac{$b} <=> $hits_ed_frac{$a}} keys %hits_ed_frac;
	
	if ($hits_ed_frac{$sorted_id[0]} >= $hits_ed_frac{$sorted_id[1]}+0.1)
	{
	    $hit=$sorted_id[0];
	    ($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,$edit_count,$edit_loc,$mm_loc,$read_sign) = split ("\t",$hit);
	    $count_selected_from_R ++;
	}
	else
	{
	    next;
	}
	
    }
    
	$count_orig_wo_reg_R{$mm_types{$edit_type}}++;

    my @edit_sites = split (";",$edit_loc);
	my @qual = split (//,$Fields[10]);
	my @qualities;
	if($read_sign eq "-") 
    {
		@qualities = reverse @qual;
    }
	else
	{
		@qualities = @qual;
	} 
    my @sites2splice = ();
    
	# record sites with low qual
    for my $i (0..(scalar(@edit_sites)-1))
    {
	if (ord($qualities[$edit_sites[$i]])-$qual_offset < $min_qual)
	{
	    $edit_count--;
	    $mm_count--;
	    $count_reject_qual_sites{$mm_types{$edit_type}}++;
	    push (@sites2splice,$i);
	}
    }
    
	# splice the low qual sites from @edit_sites array
    if (@sites2splice)
    {
	my $spliced_sites = 0;
	foreach my $site2splice (@sites2splice)
	{
	    splice (@edit_sites ,$site2splice-$spliced_sites,1);
	    $spliced_sites ++;
	}
	
    }
    
    my $read_len = length($Fields[9]);
    my $min_edit_sites = $frac_min_edit_sites*$read_len;
    my $min_len = $frac_min_len*$read_len;
    my $min_int = $frac_min_int*$read_len;
    my $min_end = $frac_min_end*$read_len;
    
    
    if ($edit_count < $min_edit_sites)
    {
	$count_reject_min_edit_sites{$mm_types{$edit_type}}++;
	next;
    }
    
    if (($edit_count/$mm_count) < $min_edit_of_mm)
    {
	$count_reject_min_edit_of_mm{$mm_types{$edit_type}}++;
	next;
    }
    
    
    my $edit_clust = substr($Fields[9],$edit_sites[0],$edit_sites[scalar(@edit_sites)-1]-$edit_sites[0]+1);
    my $edit_len = length ($edit_clust);
    
    my @prim_sites = split (";",$edit_loc);

    my $dna_clust = substr($ref_seq,$edit_sites[0]+1,$edit_sites[scalar(@edit_sites)-1]-$edit_sites[0]+1); #ref_seq included 2 bases at the ends
    my $edited_nuc = substr ($edit_type,0,1);
    
	#############################################################
	#### computing the pot editing sites (no. of As) but for now do not use it!
	
	my $potential_editing =  0;
    
    while ($dna_clust =~ /$edited_nuc/g)
    {
		$potential_editing ++; 
    }
    
    
    if (@sites2splice)
    {
		foreach my $site2splice (@sites2splice)
		{
			if (($prim_sites[$site2splice] < $edit_sites[0]) or ($prim_sites[$site2splice] > $edit_sites[scalar(@edit_sites)-1]))
			{
				next;
			}
			$potential_editing --;
		}
    }
	#############################################################
	
    
# Remove reads due to bases content 

    my %nuc_map = ('A' => 0, 'C' => 1, 'G' => 2, 'T' =>3);
    my @bases = split(//,$edit_clust);
    my @hist = (0,0,0,0);
    my $num_others = 0;
    my $tot = scalar(@bases);
    
    for my $i (0 .. ($tot-1))
    {
	if ($bases[$i] eq 'N')		
	{
	    $num_others++;
	}
	else
	{
	    $hist[$nuc_map{$bases[$i]}]++;
	}
    }		
    
    my $sum_known = $tot - $num_others;
    my @hist_s = sort {$a <=> $b} @hist;
    my $min = $hist_s[0] / $sum_known;
    my $max = $hist_s[3] / $sum_known;
    
    if ($max > $max_letter)  #($min < $min_letter) or ($max > $max_letter))
    {
	$count_reject_content{$mm_types{$edit_type}}++;
	next;
    }
    
   
    if ($edit_len < $min_len)
    {
	$count_reject_len{$mm_types{$edit_type}}++;
	next;
    }
    
    if ($edit_sites[0] > $min_int)
    {
	$count_reject_int{$mm_types{$edit_type}}++;
	next;
    }
    
    if ($edit_sites[(scalar(@edit_sites)-1)] < $min_end)
    {
	$count_reject_end{$mm_types{$edit_type}}++;
	next;
    }


    my $read_seq;
    my $dna_sign;
    my $rna_sign;
    my $sign_dt1;
    my $sign_dt2;
    
    if($edit_type eq $mm_types{$edit_type})
    {
	$dna_sign = "+";
    }
    else
    {
	$dna_sign = "-";
    }
    
    if($read_sign eq "-")
    {
	$read_seq = revcomp($Fields[9]);
	$sign_dt1 = -1;
	
    }
    else
    {
	$sign_dt1 = 1;
	$read_seq = $Fields[9];
    }
    
    if ((($Fields[1]>>4) & 1) == 1) # 0x10 is set=> read aligned to "-" strand (reverse complemented)
    {
	$sign_dt2 = -1;
    }
    else
    {
	$sign_dt2 = 1;
    }
    
    if(($sign_dt1*$sign_dt2)==1)
    {
	$rna_sign = "+";
    }
    else
    {
	$rna_sign = "-";	
    }
    
#check PE alignment

    if ($PE == 1) 
    {
	
	 if (not exists $map_reads{$read_id}) # the paired read is unmapped
	 {
	     $count_reject_PE{$mm_types{$edit_type}}++;
	     next;
	 }
	 my @hits = split ("\n",$map_reads{$read_id});
	 
	 if (scalar (@hits) == 1)
	 {
	     chomp $map_reads{$read_id};
	     if ($map_reads{$read_id} =~ m/XA:Z:(\S+$)/)
	     {
		 @hits = split (";",$1);
		 foreach my $i (0..(scalar(@hits)-1))
		 {     
		     my @f = split (",",$hits[$i]);
		     if ($f[1]=~  m/(^[-+])(\d+$)/)
		     {			 
			 $hits[$i] = "$f[0]\t$2\t$1";
		     }
		 }
	     }
	     else
	     {
		 shift(@hits); 
	     }
	    
	     my @F = split ("\t",$map_reads{$read_id});
	     if ((($F[1]>>4) & 1) == 1) # 0x10 is set=> read aligned to "-" strand (reverse complemented)
	     {
		 $F[1]="-";
	     }
	     else
	     {
		 $F[1]="+";
	     }
	     unshift(@hits,"$F[2]\t$F[3]\t$F[1]");
	 }
	 else
	 { 	 
	     foreach my $i (0..(scalar(@hits)-1))
	     {
		 chomp $hits[$i];
		 my @F = split ("\t",$hits[$i]);
		 if ((($F[1]>>4) & 1) == 1) # 0x10 is set=> read aligned to "-" strand (reverse complemented)
		 {
		     $F[1]="-";
		 }
		 else
		 {
		     $F[1]="+";
		 }
		 $hits[$i]="$F[2]\t$F[3]\t$F[1]";
	     }
	 }
	 
	 my $hits_len = scalar (@hits);

	 my $PEmatch = 0;
	 
	 foreach my $hit (@hits)
	 {
	     my @Fhit = split ("\t",$hit);

	     if (($Fhit[0] eq $ref_chr) && (abs($Fhit[1]-$ref_base) < $max_GAP) && ($Fhit[2] ne $rna_sign))
	     {
		 $PEmatch = 1;
		 last;
	     }
	 }
	 
	 if ($PEmatch == 0)
	 {
	     $count_reject_PE{$mm_types{$edit_type}}++;
	     next;
	 }
	     
    }

#the hit is UE
#record the UE


    my $dna_seq = substr($ref_seq,1,$read_len); #ref_seq included 2 bases at the ends

    my $prim_edit_loc = $edit_loc;
    $edit_loc = join (';',@edit_sites);
    $hit = "$read_id\t$ref_chr\t$ref_base\t$ref_seq\t$mm_count\t$edit_type\t$edit_count\t$edit_loc\t$mm_loc";
    
    my $bed_ref_base = $ref_base-1+$edit_sites[0];
    my $end_pos=$ref_base+$edit_sites[scalar(@edit_sites)-1];

    $count_ue{$mm_types{$edit_type}}++;
    $count_edit_sites{$mm_types{$edit_type}}+=$edit_count;

    print {$bed_files{$mm_types{$edit_type}}} "$ref_chr\t$bed_ref_base\t$end_pos\t($rna_sign)$read_id\t$edit_type\t$dna_sign\n";
	
  
     foreach my $ed_site (@edit_sites)
    {
	my $es_bed_ref_base = $ref_base-1+$ed_site;
	my $es_bed_end_base = $es_bed_ref_base+1;
	my $triplet = substr($ref_seq,$ed_site,3); 
	my $real_ed_ind = $ed_site; 
	if($rna_sign eq "-")
    { 
	$real_ed_ind = $read_len-$ed_site-1;
    }

	print {$es_bed_files{$mm_types{$edit_type}}} "$ref_chr\t$es_bed_ref_base\t$es_bed_end_base\t($rna_sign)$read_id;$real_ed_ind\t$triplet\t$dna_sign\n"; 
	$edit_ind{$mm_types{$edit_type}} += $real_ed_ind;
		    
	if($edit_type eq $mm_types{$edit_type})
	{
	    my $up_base = substr($ref_seq,$ed_site,1);
	    if ($up_base ne "N")
	    {
		$count_up_tot_bases{$mm_types{$edit_type}}++;
		$count_up{$mm_types{$edit_type}}{$up_base}++;
	    }
	    
	    my $down_base = substr($ref_seq,$ed_site+2,1);
	    if ($down_base ne "N")
	    {
		$count_down_tot_bases{$mm_types{$edit_type}}++;
		$count_down{$mm_types{$edit_type}}{$down_base}++;
	    }
	}
	    
	else
	{
	    my $orig_up_base = substr($ref_seq,$ed_site+2,1);
	    my $up_base = $comp_base{$orig_up_base};
	    if ($up_base ne "N")
	    {
		$count_up_tot_bases{$mm_types{$edit_type}}++;
		$count_up{$mm_types{$edit_type}}{$up_base}++;
	    }

	    my $orig_down_base = substr($ref_seq,$ed_site,1);
	    my $down_base = $comp_base{$orig_down_base};
	    if ($down_base ne "N")
	    {
		$count_down_tot_bases{$mm_types{$edit_type}}++;
		$count_down{$mm_types{$edit_type}}{$down_base}++;
	    }
	    
	}
    }

    my @prim_edit_sites = split (";",$prim_edit_loc);
    my @mm_sites = split (";",$mm_loc);
    my $dna_end_pos=$ref_base+$read_len-1;

    my $mm_not_edit = $mm_count-$edit_count;    
    my $mm_in_clust = 0;

    if (@mm_sites)
    {
	foreach my $mm (@mm_sites)
	{
	    if ($mm > $edit_sites[0] && $mm < $edit_sites[scalar(@edit_sites)-1])
	    {
		$mm_in_clust++;
	    }	
	}
    }
    
    $lens{$mm_types{$edit_type}} += $edit_len;
    $align_mm{$mm_types{$edit_type}} += $mm_not_edit;
    $clust_mm{$mm_types{$edit_type}} += $mm_in_clust;
    
###print "Hit: $hit\n"; #####################

########### $mm_types{$edit_type} report the real 1 of 12 edit_type instead 1 of 6
    print $ue_list "$edit_type\t($rna_sign)$Fields[0]\t$ref_chr\t$ref_base\t$dna_seq\t$read_seq\t$read_len\t$edit_len\t$edit_count\t$mm_not_edit\n";
    print $ue_det "Edited read: ($rna_sign)$Fields[0], Edit type: $edit_type, Aligns to: $ref_chr:$ref_base-$dna_end_pos\nEdit sites: $edit_count, Mismatches sites in alignment: $mm_not_edit, Mismatches sites in cluster: $mm_in_clust, Alignment length: $read_len, Cluster length: $edit_len\nEdit indexes: $edit_loc, Mismatches indexes: $mm_loc\nDNA: $dna_seq\n     ";
    
    foreach $b (0..$read_len-1)
    {
	if (grep(/^$b$/, @edit_sites))
	{
	    print $ue_det "*";
	    next;
	}
	if (grep(/^$b$/, @prim_edit_sites))
	{
	    print $ue_det "-";
	    next;
	}
	if (grep(/^$b$/, @mm_sites))
	{
	    print $ue_det "X";
	    next;
	}
	print $ue_det "|";

    }
    print $ue_det "\nRNA: $read_seq\n|) Match site,  X) Mismatch site, *) Edit site, -) Edit site with quality < $min_qual (not counted)\n\n";
}


print $f_statistic "\nStatistic of detection UE of $out_pre\n";
print $f_statistic "(Arguments: $ARGV[3]\t$ARGV[4]\t$ARGV[5]\t$ARGV[6]\t$ARGV[7]\t$ARGV[8]\t$ARGV[9])\n";

#print $f_statistic "\ncount_orig_reads\tcount_uniq_loc\n$count_orig_reads\t$count_uniq_loc\n";
#print $f_statistic "\nedit_type\tcount_analysed\tcount_pot_edit_sites\tcount_reject_qual_sites\tcount_reject_ue\tcount_reject_min_edit_sites\tcount_reject_min_edit_of_mm\tcount_edit_sites\tcount_ue\n";

#print header

print $f_statistic

"edit_type".
"\torig_(WO_reg_R)".
"\treject_qual_sites".
"\treject_min_edit_sites".
"\treject_min_edit_of_mm".
"\treject_content".
"\treject_len".
"\treject_int".
"\treject_end";
if ($PE == 1)
{
    print $f_statistic "\treject_PE";
}

print $f_statistic "\t".
"\tUE".
"\tavg_lens".
"\tedit_sites".
"\tavg_es".
"\talign_mm".
"\tavg_align_mm".
"\tclust_mm".
"\tavg_clust_mm".
"\tavg_edit_ind".

"\t".
"\tA_up\tC_up\tG_up\tT_up".
"\tA_down\tC_down\tG_down\tT_down".

"\n";

foreach my $edit_type (keys %count_ue)
{

#Analyse upstream/downstream bases

    foreach my $base (keys %{$count_up{$edit_type}})
    { 
	$frac_up{$base}=0;
	$frac_down{$base}=0;
	$frac_up{$base}=$count_up{$edit_type}{$base}/$count_up_tot_bases{$edit_type} unless ($count_up_tot_bases{$edit_type} ==0); 
	$frac_down{$base}=$count_down{$edit_type}{$base}/$count_down_tot_bases{$edit_type} unless ($count_down_tot_bases{$edit_type} ==0); 
    }

#Analyse Edit clusters

    my $avg_lens = 0;
    $avg_lens = $lens{$edit_type}/$count_ue{$edit_type} unless ($count_ue{$edit_type}==0);

    
#Analyse Edit sites

    my $avg_ed_site = 0;
    $avg_ed_site = $count_edit_sites{$edit_type}/$count_ue{$edit_type} unless ($count_ue{$edit_type}==0);

    my $avg_align_mm = 0;
    $avg_align_mm = $align_mm{$edit_type}/$count_ue{$edit_type} unless ($count_ue{$edit_type}==0);

    my $avg_clust_mm = 0;
    $avg_clust_mm = $clust_mm{$edit_type}/$count_ue{$edit_type} unless ($count_ue{$edit_type}==0);

	my $avg_edit_ind = 0;
    $avg_edit_ind = $edit_ind{$edit_type}/$count_edit_sites{$edit_type} unless ($count_edit_sites{$edit_type}==0);


	
	
#print results

    print $f_statistic

"$edit_type".
"\t$count_orig_wo_reg_R{$edit_type}".
"\t$count_reject_qual_sites{$edit_type}".
"\t$count_reject_min_edit_sites{$edit_type}".
"\t$count_reject_min_edit_of_mm{$edit_type}".
"\t$count_reject_content{$mm_types{$edit_type}}".
"\t$count_reject_len{$edit_type}".
"\t$count_reject_int{$edit_type}".
"\t$count_reject_end{$edit_type}";
if ($PE == 1)
{
    print $f_statistic "\t$count_reject_PE{$edit_type}";
}

print $f_statistic "\t".
"\t$count_ue{$edit_type}".
"\t$avg_lens".
"\t$count_edit_sites{$edit_type}".
"\t$avg_ed_site".
"\t$align_mm{$edit_type}".
"\t$avg_align_mm".
"\t$clust_mm{$edit_type}".
"\t$avg_clust_mm".
"\t$avg_edit_ind".

"\t".
"\t$frac_up{A}\t$frac_up{C}\t$frac_up{G}\t$frac_up{T}".
"\t$frac_down{A}\t$frac_down{C}\t$frac_down{G}\t$frac_down{T}".

"\n";
}


sub revcomp {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}
