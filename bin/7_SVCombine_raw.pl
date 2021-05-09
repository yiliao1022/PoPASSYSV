#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

####################################################################################
# SVCombine_raw.pl
#
# Authors: Yi Liao (06/30/2020)
#
# Copyright (C) Not for commercial use
#
# Merge raw SVs from pairwise comparison files
#
# Prerequisite :
#
#  Usage: perl $0 -sv SV.txt -svtype INS
#  -sv         a file combined with all paired files     [REQUIRED]
#  -svtype     SV types [INS|DEL|CNV|INV]                [REQUIRED]
#  -ins_tol    coordinate tolerance either up or downstream to be considered or tested if they are the same insertion event [default : 5]
#  -help       print this
#
####################################################################################
my ($svfile,$svtype,$ins_tol,$del_ov,$Help);

GetOptions(
  "sv=s"=>\$svfile,
  "svtype=s"=>\$svtype,
  "ins_tol=f"=>\$ins_tol,
  "del_ov=f"=>\$del_ov,
  "help"     =>\$Help
);

$ins_tol ||= 5;
$del_ov ||= 0.8;

##################################################################################################
if ($Help){
print <<"END.";
  Usage: perl $0 -sv SV.txt -svtype INS
  -sv         a file combined with all paired files     [REQUIRED]
  -svtype     SV types [INS|DEL|CNV|INV]                [REQUIRED]
  -ins_tol    coordinate tolerance either up or downstream to be considered or tested if they are the same insertion event [default : 5]
  -del_ov     coordinate overlap for deletions that are considered or tested if thet are the same deletion event [default : 0.8] 
  -help       print this
END.
exit;
}

if ($svtype eq "INS") {
&Insertion ($svfile,$ins_tol);
} elsif ($svtype eq "DEL") {
&Deletions ($svfile,$del_ov);
} elsif ($svtype eq "CNV") {
 
} elsif ($svtype eq "INV") {
 
} else {
  print "You don't provide the SV type or the type you provided is not currently supported!";
}


############  Combine Deletions

sub Deletions {
my ($del,$overlap_cutoff) = @_;
`sort -k1,1 -k2,2n -k3,3n -k8,8n $del > $del.sort.txt`;
my %hash;
my $log="NA_0_20000000000000_0_0_0";
open DEL, "$del.sort.txt" or die "$!";
open DELOUT, ">$del.sort.merge.txt" or die "$!";

while (<DEL>) {
  chomp;
  my @tmp = split (/_/,$log);
  my ($ref_chr,$ref_s,$ref_e,$alt_chr,$alt_s,$alt_e,$len,$five,$three,$seq) = (split(/\t/,$_))[0,1,2,4,5,6,7,9,10,11];
  my $DEL = join ("_",($ref_chr,$ref_s,$ref_e,$len,$five,$three));
  my $overlap = &Overlap($ref_s,$ref_e,$tmp[1],$tmp[2]);
  if (($ref_chr eq $tmp[0]) and $overlap>$overlap_cutoff and $len/$tmp[3]>0.95 and $len/$tmp[3]<1.05) {
  my $species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
  push @{$hash{$log}},$species;
  } else {
  $log = $DEL;
  my $species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
  push @{$hash{$log}},$species;
  }
}

my $head = "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\tMERGE_SAMPLES_NUM\tMERGE_SAMPLES\tQUERY\tQ_POS\tQ_END\tFIVE\tTHREE\n";

foreach my $ind (keys %hash) {
      my @tmp = split (/\_/,$ind);
      my $id = ${$hash{$ind}}[0];
      my ($alt,$alt_s,$alt_e) = (split(/_/,$id))[0,1,2];
      my $num = $#{$hash{$ind}} + 1;
      print DELOUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$id\tDEL\t$tmp[3]\t$num\t";

if ($#{$hash{$ind}} == 0) {
      print DELOUT "${$hash{$ind}}[0]\t";
} else {
  my $last_one = pop @{$hash{$ind}};
  foreach my $svind (@{$hash{$ind}}) {
      print DELOUT "$svind,";
  }
      print DELOUT "$last_one\t";
}
      print DELOUT "$alt\t$alt_s\t$alt_e\t$tmp[4]\t$tmp[5]\n";
}

`sort -k1,1 -k2,2n $del.sort.merge.txt > $del.combined.sorted.txt`;
`sed -i '1i $head' $del.combined.sorted.txt`;
`rm $del.sort.txt $del.sort.merge.txt`;
}

###### Combine insertions
sub Insertion {
my ($ins,$tol) = @_;
`sort -k5,5 -k6,6n -k7,7n -k8,8n $ins > $ins.sorted.txt`;
my %hash;
my $log="NA_0_0_0_NNNNNN";

open INS, "$ins.sorted.txt" or die "No insertion file is provided!";
open OUT, ">$ins.sorted.combined.txt" or die "Unable to write  output";
while (<INS>) {
chomp;
my @tmp = split (/_/,$log);
my ($ref_chr,$ref_s,$ref_e,$alt_chr,$alt_s,$alt_e,$len,$five,$three,$seq) = (split(/\t/,$_))[4,5,6,0,1,2,7,9,10,11];
my $sequence = uc ($seq);
my $seqIdentity = &SeqIdentity($tmp[4],$sequence);
my $SV = join ("_",($ref_chr,$ref_s,$ref_e,$len,$sequence,$five,$three));
if ($ref_chr eq $tmp[0] and $ref_s<($tmp[1]+$tol) and $seqIdentity==1) {
  my $Species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
  push @{$hash{$log}},$Species;
} else {
  $log = $SV;
  my $Species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
  push @{$hash{$log}},$Species;
  }
}

my $head = "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\tMERGE_SAMPLES_NUM\tMERGE_SAMPLES\tQUERY\tQ_POS\tQ_END\tFIVE\tTHREE\n";

foreach my $ind (keys %hash) {
      my @tmp = split (/\_/,$ind);
      my $id = ${$hash{$ind}}[0];
      my ($alt,$alt_s,$alt_e) = (split(/_/,$id))[0,1,2];
      my $num = $#{$hash{$ind}} + 1;
      print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$id\tINS\t$tmp[3]\t$num\t";

if ($#{$hash{$ind}} == 0) {
      print OUT "${$hash{$ind}}[0]\t";
} else {
  my $last_one = pop @{$hash{$ind}};
  foreach my $svind (@{$hash{$ind}}) {
      print OUT "$svind,";
  }
      print OUT "$last_one\t";
}
      print OUT "$alt\t$alt_s\t$alt_e\t$tmp[5]\t$tmp[6]\n";
}

`sort -k1,1 -k2,2n $ins.sorted.combined.txt > $ins.combined.sorted.txt`;
`sed -i '1i $head' $ins.combined.sorted.txt`;
`rm $ins.sorted.txt $ins.sorted.combined.txt`;
}

###### Compare identity of the insertion sequences from difference genomes
sub SeqIdentity {
my ($seq1,$seq2) = @_;
my $len1=length($seq1)+0.1;
my $len2=length($seq2)+0.1;
my $ident=0;
if ($seq1 eq $seq2) {
  $ident=1;
} elsif ($len1 > 5 and $len1 < 21 and $len1/$len2 > 0.8 and $len1/$len2<1.2) {
  my @seeds;
  my $num = $len1-4+1;
  my $i;
for ($i=0;$i<$num;$i++) {
my $seed = substr ($seq1,$i,4);
push (@seeds,$seed);
}
my $m=0;
foreach my $sed (@seeds) {
if ($seq2=~/$sed/) {
  $m++;
 }
}
 my $perc = $m/$num;
   if ($perc > 0.5) {
     $ident =1;
  }
} elsif ($len1 > 20 and $len2 > 20 and $len1/$len2 > 0.8 and $len1/$len2<1.2) {
  my @seeds;
  my $num = $len1-10+1;
  my $i;
  for ($i=0;$i<$num;$i++) {
  my $seed = substr ($seq1,$i,10);
  push (@seeds,$seed);
  }
  my $m=0;
  foreach my $sed (@seeds) {
  if ($seq2=~/$sed/) {
    $m++;
   }
  }
  my $perc = $m/$num;
  if ($perc > 0.5) {
    $ident =1;
  }
}
return $ident;
}


############ Calculate the overlap between two coordinate intervals

sub Overlap {
my ($c1,$c2,$c3,$c4) = @_;
my $overlap;

if ($c1>$c4) {
  $overlap =0;
} elsif ($c1>=$c3 and $c1<=$c4 and $c2 <=$c4) {
  $overlap = ($c2-$c1)/($c4-$c3);
} elsif ($c1>=$c3 and $c1<=$c4 and $c2 >$c4) {
  $overlap = ($c4-$c1)/($c4-$c3);
} 
return $overlap;  
}
