#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use POSIX qw(strftime);
use File::Basename qw(basename dirname);
use List::Util 'first';

####################################################################################
#
# Authors: Yi Liao (08/06/2020)
#
# Copyright (C) Not for commercial use
#
# Further computationally validated the merge raw SVs
#
# Prerequisite : TBA package
#
####################################################################################
my $time= strftime "%F", localtime;
my ($svraw,$svtype,$glst,$singmaf,$tba,$refsize,$refname,$speciesname,$OUTPUTdir,$extend,$head,$help);

GetOptions(
    "indel=s"=>\$svraw,
    "indeltype=s"=>\$svtype,
    "glst=s"=>\$glst,
    "singmaf=s"=>\$singmaf,
    "tba=s"=>\$tba,
    "refsize=s"=>\$refsize,
    "refname=s"=>\$refname,
    "speciesname=s"=>\$speciesname,
    "extend=s"=>\$extend,
    'outdir=s'=>\$OUTPUTdir,
    'head'=>\$head,
    'help'=>\$help,
);

############# Help information ##########
if ($help){
print <<"END.";
  Usage: perl $0 -indel SV.txt -indeltype INS -glst genome.lst -refsize B73.sizes -singmaf /path/to/sing/ -tba /path/to/tba -refname B73 -speciesname maize -outdir B73_tmp
  -indel        raw merged InDels file                                  [REQUIRED]
  -indeltype    InDel types [INS|DEL]                                   [REQUIRED]
  -singmaf      path to pairwise genome algnment maf files              [REQUIRED]
  -glst         list of sample names file (one sample per line)         [REQUIRED]
  -refsize      reference chromosome sizes file                         [REQUIRED]
  -refname      reference name, e.g. hg38, IRGSP1                       [REQUIRED]
  -speciesname  species name, e.g. human,Drosophila,rice,maize          [REQUIRED]
  -tba          path to tba program                                     [REQUIRED]
  -outdir       output Dir                                              [REQUIRED]
  -head         whether to print VCF header                             [OPTIONAL] 
  -help         print this information
END.
exit;
}

##Define default parameters
$OUTPUTdir ||= ".";
$OUTPUTdir =~ s/\/$//;
$extend ||= 10000;
mkdir($OUTPUTdir) unless(-d $OUTPUTdir);

#######################################################
###################### Main ###########################
#######################################################

open (GENOME,"$glst") or die "Genome/species names do not exist, $!"; ### read genoms/assemblies ID list from a input file (How many genomes are you processed)
my @Genome=<GENOME>;
close GENOME;
chomp @Genome;
my $reference=shift @Genome;
my $ref_glst = \@Genome;

open VCF, ">$OUTPUTdir/$svraw.vcf" or die "$!";
open ERR, ">$OUTPUTdir/$svraw.err" or die "$!";

############ Printing VCF header ###########################
if ($head) {
print VCF "\#\#fileformat=VCFv4.2\n\#\#fileDate=$time\n\#\#source=YiLiao_PopASSYSV_version20200901\n\#\#reference=$refname\n";
open Contig, "$refsize" or die "$!";
while (<Contig>) {
chomp;
my @tmp = split(/\s/,$_);
print VCF "\#\#contig=<ID=$tmp[0],length=$tmp[1],assembly=$refname,md5=XXXXXXXX,species=\"$speciesname\",taxonomy=x>\n";
}
close Contig;
print VCF "##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">
##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">
##ALT=<ID=INS,Description=\"Insertion of novel sequence\">
##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">
##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">
##ALT=<ID=INV,Description=\"Inversion\">
##ALT=<ID=CNV,Description=\"Copy number variable region\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">
";
print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
my $num;
for ($num=0;$num<$#Genome;$num++) {
print VCF "$Genome[$num]\t";
}
print VCF "$Genome[-1]\n";
}

#######################################################
my $first=`head -1 $svraw`;
my $last=`tail -1 $svraw`;
my @five = split (/\t/,$first);
my @three = split (/\t/,$last);
my $beg = $five[1]-$extend;
if ($beg<0 or ($beg !~ /^\d+$/)) {
$beg=0;
}
my $end = $three[2]+$extend;
my %liftover;

foreach my $g (@Genome) {
    unless (-e "$OUTPUTdir/$svraw.$reference.$g.sing.maf") {
    `${tba}mafFind $singmaf/$reference.$g.sing.maf $beg $end > $OUTPUTdir/$svraw.$reference.$g.sing.maf`;
     }
    my $liftover_coord = &liftOver("$OUTPUTdir/$svraw.$reference.$g.sing.maf");
    $liftover{$g} = $liftover_coord;
}

open IN, "$svraw" or die "$!";
while (<IN>) {
chomp;
my ($chr,$beg,$end,$len,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,8,9,10,11,12];
my $ref = (split(/\./,$chr))[0];
my $alt = (split(/\./,$qchr))[0];

eval {
if ($svtype eq "DEL") {
print VCF "$chr\t$beg\tSV$beg\tN\t\<$svtype\>\t22\tPASS\tSVTYPE=INS;END=$end;SVLEN=$len\tGT:GQ\t";
   my $beg0 = $beg + 1;
   my $end0 = $end + 2;
   foreach my $ind (@Genome) {
     if (exists ${$liftover{$ind}}{$beg} && exists ${$liftover{$ind}}{$end0}) {
          my @tmp1 = split (/_/, ${$liftover{$ind}}{$beg});
          my @tmp2 = split (/_/, ${$liftover{$ind}}{$end0}); 
          my $check = abs ($tmp2[1]-$tmp1[1]);
          print "$chr\t$beg\t$ind\t$tmp2[1]\t$tmp1[1]\t$check\n";
          if ($check == $len + 2) {
             if ($ind eq $Genome[-1]) {
              print VCF "0|0:22";
              } else {
              print VCF "0|0:22\t";
              }     
           } elsif ($check <= 2) {
             if ($ind eq $Genome[-1]) {
              print VCF "1|1:22";
              } else {
              print VCF "1|1:22\t";
              }             
           } else {
             if ($ind eq $Genome[-1]) {
              print VCF "\.|\.:22";
              } else {
              print VCF "\.|\.:22\t";
              }
           }    
     } else {
       if ($ind eq $Genome[-1]) {
        print VCF "\.|\.:22";
        } else {
        print VCF "\.|\.:22\t";
        } 
     }
  }    
} elsif ($svtype eq "INS") {
print VCF "$chr\t$beg\tSV$beg\tN\t\<$svtype\>\t22\tPASS\tSVTYPE=INS;END=$end;SVLEN=\-$len\tGT:GQ\t"; 
my $beg0 = $beg - 2;
my $end0 = $end + 2;
foreach my $ind (@Genome) {
   if (exists ${$liftover{$ind}}{$beg0} && exists ${$liftover{$ind}}{$end0}) {
       my @tmp1 = split (/_/, ${$liftover{$ind}}{$beg0});
       my @tmp2 = split (/_/, ${$liftover{$ind}}{$end0}); 
       my $check = abs ($tmp2[1]-$tmp1[1]);
       print "$chr\t$beg\t$ind\t$tmp2[1]\t$tmp1[1]\t$check\n";
       if ($check == $len + 4) {
          if ($ind eq $Genome[-1]) {
           print VCF "1|1:22";
           } else {
           print VCF "1|1:22\t";
           }     
        } elsif ($check <= 4) {
          if ($ind eq $Genome[-1]) {
           print VCF "0|0:22";
           } else {
           print VCF "0|0:22\t";
           }             
        } else {
          if ($ind eq $Genome[-1]) {
           print VCF "\.|\.:22";
           } else {
           print VCF "\.|\.:22\t";
           }
        }    
  } else {
    if ($ind eq $Genome[-1]) {
     print VCF "\.|\.:22";
     } else {
     print VCF "\.|\.:22\t";
     } 
  }
 } 
}
print VCF "\n";
};

if ($@) {
  print ERR "Failed\t$_\n";
}

}

#########################################################
############################ Subroutines ################
#########################################################
sub liftOver {
my $maf = shift;
my %liftover;
local $/="\n\n";
open MAF, "<$maf" or die "$!";
my @first = split (/\n/,<MAF>);
my @array;
foreach my $ele (@first) {
  if ($ele=~/^s/) {
   push (@array,$ele);
  }
}

my ($target,$Ts,$Textend,$Tstrand,$Tsize,$Tseq) = (split(/\s+/,$array[0]))[1,2,3,4,5,6];
my ($query,$Qs,$Qextend,$Qstrand,$Qsize,$Qseq) = (split(/\s+/,$array[1]))[1,2,3,4,5,6];
my @targetSEQ = split(//,$Tseq);
my @querySEQ = split (//,$Qseq);
my $i;
my $n=0;
my $m=0;

for ($i=0;$i<=$#targetSEQ;$i++) {
   if ($targetSEQ[$i]!~/\-/) {
     $n++;
   }
   if ($querySEQ[$i]!~/\-/) {
     $m++;
   }

   if ($targetSEQ[$i]!~/\-/ && $querySEQ[$i]!~/\-/) {
     if ($Qstrand =~ /\+/) {
       my $tcoord =  $Ts + $n + 1;
       my $qcoord =  $Qs + $m + 1;
       my $qchrCoord = join  ("_",($query,$qcoord));
       $liftover{$tcoord} = $qchrCoord;
     } elsif ($Qstrand =~/\-/) {
       my $tcoord =  $Ts + $n + 1;
       my $qcoord =  $Qsize - ($Qs + $m);
       my $qchrCoord = join  ("_",($query,$qcoord));
       $liftover{$tcoord} = $qchrCoord;
     }
   }
}

while (<MAF>) {
next if ($_=~/eof/);
my @tmp = split (/\n/,$_);
my ($target1,$Ts1,$Textend1,$Tstrand1,$Tsize1,$Tseq1) = (split(/\s+/,$tmp[1]))[1,2,3,4,5,6];
my ($query1,$Qs1,$Qextend1,$Qstrand1,$Qsize1,$Qseq1) = (split(/\s+/,$tmp[2]))[1,2,3,4,5,6];
my @targetSEQ1 = split(//,$Tseq1);
my @querySEQ1 = split (//,$Qseq1);
my $i1;
my $n1=0;
my $m1=0;
for ($i1=0;$i1<=$#targetSEQ1;$i1++) {
   if ($targetSEQ1[$i1]!~/\-/) {
     $n1++;
   }
   if ($querySEQ1[$i1]!~/\-/) {
     $m1++;
   }
   if ($targetSEQ1[$i1]!~/\-/ && $querySEQ1[$i1]!~/\-/) {
     if ($Qstrand1 =~ /\+/) {
       my $tcoord1 =  $Ts1 + $n1 + 1;
       my $qcoord1 =  $Qs1 + $m1 + 1;
       my $qchrCoord1 = join  ("_",($query1,$qcoord1));
       $liftover{$tcoord1} = $qchrCoord1;;
     } elsif ($Qstrand1 =~/\-/) {
       my $tcoord1 =  $Ts1 + $n1 + 1;
       my $qcoord1 =  $Qsize1 - ($Qs1 + $m1);
       my $qchrCoord1 = join  ("_",($query1,$qcoord1));
       $liftover{$tcoord1} = $qchrCoord1;
     }
   }
 }
}
my $refhash = \%liftover;
return $refhash;
}
#################################################################
