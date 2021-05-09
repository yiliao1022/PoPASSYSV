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
my ($svfile,$svtype,$ins_tol,$del_ov,$candidateNum,$Help);

GetOptions(
  "sv=s"=>\$svfile,
  "svtype=s"=>\$svtype,
  "ins_tol=f"=>\$ins_tol,
  "del_ov=f"=>\$del_ov,
  "cant=i"=>\$candidateNum,
  "help"     =>\$Help
);

$ins_tol ||= 13;
$del_ov ||= 0.95;
$candidateNum ||= 20;
##################################################################################################
if ($Help){
print <<"END.";
  Usage: perl $0 -sv SV.txt -svtype INS
  -sv         a file combined with all paired files     [REQUIRED]
  -svtype     SV types [INS|DEL|CNV|INV]                [REQUIRED]
  -ins_tol    coordinate tolerance either up or downstream to be considered or tested if they are the same insertion event [default : 13]
  -del_ov     coordinate overlap for deletions that are considered or tested if thet are the same deletion event [default : 0.95] 
  -cant       number of the nearest SVs are used to check whether the current one belong to [default : 20]
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
my @content;
my $log="NA_0_0_0_NNNNNN_0_0";
push (@content,$log);
open INS, "$ins.sorted.txt" or die "No insertion file is provided!";
open OUT, ">$ins.sorted.combined.txt" or die "Unable to write  output";

while (<INS>) {
chomp;
my ($ref_chr,$ref_s,$ref_e,$alt_chr,$alt_s,$alt_e,$len,$five,$three,$seq) = (split(/\t/,$_))[4,5,6,0,1,2,7,9,10,11];
my $sequence = uc ($seq);
my $SV = join ("_",($ref_chr,$ref_s,$ref_e,$len,$sequence,$five,$three));

my $i;

for ($i=0;$i<=$#content;$i++) {
my $Comparion = &SeqIdentity($content[$i],$SV,$tol);
if ($Comparion ==1) {
  my $Species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
  push @{$hash{$content[$i]}},$Species;
  last;
}

if ($i == $#content) {
 
  if ($#content < $candidateNum) {
    unshift (@content,$SV);
  } else {
    pop (@content);
    unshift (@content,$SV);
  }
  my $Species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
  push @{$hash{$SV}},$Species;
}

  


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
my ($SV1,$SV2,$dis) = @_;
my ($chr1,$beg1,$end1,$len1,$sequence1,$five1,$three1) = (split(/\_/,$SV1))[0,1,2,3,4,5,6];
my ($chr2,$beg2,$end2,$len2,$sequence2,$five2,$three2) = (split(/\_/,$SV2))[0,1,2,3,4,5,6];
my $ident=0;

if (($chr1 eq $chr2) && ($sequence1 eq $sequence2) && (abs($beg1-$beg2) < $dis)) {
  $ident=1;
} elsif (($chr1 eq $chr2) && (abs($beg1-$beg2) < $dis) && ($five1 == $five2 || $three1==$three2)) {
  $ident=1;
} elsif (($chr1 eq $chr2) &&  $len1 > 5 && $len1 < 50 && $len1/$len2 > 0.9 && $len1/$len2<1.1) {

  my @seeds;
  my $num = $len1-4+1;
  my $i;
  for ($i=0;$i<$num;$i++) {
   my $seed = substr ($sequence1,$i,4);
   push (@seeds,$seed);
  }

  my $m=0;
  foreach my $sed (@seeds) {
  if ($sequence2=~/$sed/) {
    $m++;
   }
  }
  
  my $perc = $m/$num;
   if ($perc > 0.5) {
     $ident =1;
   }
  
} elsif (($chr1 eq $chr2) && $len1 > 50 && (abs($beg1-$beg2) < $dis)) {
  my @seeds;
  my $num;
  my $i;
  my ($SEQ1,$SEQ2);
  if ($len1 > $len2) {
    $SEQ1 = $sequence2;
    $SEQ2 = $sequence1;
    $num = $len2-10+1;
  } else {
    $SEQ1 = $sequence1;
    $SEQ2 = $sequence2;
    $num = $len1-10+1;
  }

  
  for ($i=0;$i<$num;$i++) {
  my $seed = substr ($SEQ1,$i,10);
  push (@seeds,$seed);
  }
  
  my $m=0;
  foreach my $sed (@seeds) {
  if ($SEQ2=~/$sed/) {
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













