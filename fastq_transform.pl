#!/usr/bin/perl -w

use strict;

########################################################################################
#
#     in: A letter to be transformed [0], a letter to transformed to [1] 
#         and fastq file of the reads to be transformed [2]
#     out: fastq file with transformed reads
#
######################################################################################## 


my $fastq_file = $ARGV[2];
my $s_str="";
my $tr_str="";

#print "inputfile $fastq_file\n";

if (($ARGV[0] =~ "[Tt]") and ($ARGV[1] =~ "[Cc]"))
{
    $s_str = 'T';
    $tr_str = 'C';
}
elsif (($ARGV[0] =~ "[Aa]") and ($ARGV[1] =~ "[Gg]"))
{
    $s_str = 'A';
    $tr_str = 'G';
}
elsif (($ARGV[0] =~ "[Tt]") and ($ARGV[1] =~ "[Gg]"))
{
    $s_str = 'T';
    $tr_str = 'G';
}
elsif (($ARGV[0] =~ "[Aa]") and ($ARGV[1] =~ "[Tt]"))
{
    $s_str = 'A';
    $tr_str = 'T';
}
elsif (($ARGV[0] =~ "[Tt]") and ($ARGV[1] =~ "[Aa]"))
{
    $s_str = 'T';
    $tr_str = 'A';
}
elsif (($ARGV[0] =~ "[Aa]") and ($ARGV[1] =~ "[Cc]"))
{
    $s_str = 'A';
    $tr_str = 'C';
}
elsif (($ARGV[0] =~ "[Gg]") and ($ARGV[1] =~ "[Cc]"))
{
    $s_str = 'G';
    $tr_str = 'C';
}
elsif (($ARGV[0] =~ "[Cc]") and ($ARGV[1] =~ "[Gg]"))
{
    $s_str = 'C';
    $tr_str = 'G';
}
else
{   
    print "Illegal letters to transform!\n";
    die;
}

open (my $InFile,$fastq_file);

while (<$InFile>)
{
    my $IDLine = $_;
    my $SeqLine = <$InFile>;
    my $BlankLine = <$InFile>;
    my $BQLine = <$InFile>;

    $SeqLine =~ s/$s_str/$tr_str/g;

    print "$IDLine$SeqLine$BlankLine$BQLine";
}
