#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(&basename &dirname);

####################################################################################
# Realign.pl
#
# Authors: Yi Liao (06/26/2020)
# 
# Copyright (C) Not for commercial use
#
# Realign the .synnet.axt file using EMboss stretcher
# 
# Prerequisite : EMboss stretcher
#
# Usage:  perl 4_PairRealign.pl -axt *.synnet.axt -output *.synnet.religned.axt -stretcher /path/to/stretcher
#
# Options: -axt        [FILE] *.synnet.axt file        [REQUIRED]
#          -output     [FILE] output file              [REQUIRED]
#          -stretcher  [FILE] path to stretcher        [REQUIRED]
#          -help       print this information
#
####################################################################################

my ($axt,$output,$stretcher,$Help);
GetOptions( 'axt=s' => \$axt,
            'output=s' => \$output,
            'stretcher=s' => \$stretcher,
            'help' => \$Help
            );

if ($Help){
print <<"END.";
Usage:  perl 4_PairRealign.pl -axt *.synnet.axt -output *.synnet.religned.axt -stretcher /path/to/stretcher

Options: -axt        [FILE] *.synnet.axt file        [REQUIRED]
         -output     [FILE] output file              [REQUIRED]
         -stretcher  [FILE] path to stretcher        [REQUIRED]
         -help       print this information
END.
exit;
}

##Default
$stretcher ||='';

### Main 
$/ = "\n\n";
open In, "$axt" or die "Provide the .axt file!";
open Out, ">$output" or die "Can't open the output file!";
while (<In>) {
chomp;
my @temp = split (/\n/,$_) unless ($_=~/^$/);
my @unit = split (/\s/,$temp[-3]);
my $target = join ("_",("TARGET",$unit[1],$unit[2]));
my $query = join ("_",("QUERY",$unit[4],$unit[5]));
if (!$temp[-1]=~/-/ and !$temp[-2]=~/-/) {
print Out "$temp[-3]\n$temp[-2]\n$temp[-1]\n\n";
} else {
$temp[-1] =~s/-//g;
$temp[-2] =~s/-//g;
my ($Seq1,$Seq2) = &Realign($target,$query,$temp[-2],$temp[-1]);
print Out "$temp[-3]\n$Seq1\n$Seq2\n\n";
  }
}

##################################
####### Sub Functions ############
##################################
sub Realign {
my ($t,$q,$seq1,$seq2)=@_;
open SEQ1, ">$t.fa" or die "$!";
print SEQ1 ">$t\n$seq1\n";
open SEQ2, ">$q.fa" or die "$!";
print SEQ2 ">$q\n$seq2\n";
`${stretcher}stretcher -asequence $t.fa -bsequence $q.fa -outfile $t.$q.stretcher`;
my ($s1, $s2) = &ParseStretcherOutPut("$t.$q.stretcher");
`rm $t.fa $q.fa $t.$q.stretcher`;
return ($s1,$s2);
}
##########################
sub ParseStretcherOutPut {
my $stretcher = shift;
local $/="\n";
open EMBOSS, "$stretcher" or die "$!";
my @ref;
my @query;
while (<EMBOSS>) {
chomp;
if ($_=~/\#/) {
}elsif ($_=~/^TARGET/) {
my @tmp1 = split (/\s/,$_);
push (@ref,$tmp1[1]);
} elsif ($_=~/^QUERY/) {
my @tmp2 = split (/\s/,$_);
push (@query,$tmp2[1]);
 }
}
my $seq1 = join ("",@ref);
my $seq2 = join ("",@query);
return ($seq1,$seq2);
}
