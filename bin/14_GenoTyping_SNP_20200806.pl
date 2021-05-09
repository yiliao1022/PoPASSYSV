#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use POSIX qw(strftime);
use File::Basename qw(basename dirname);
use List::Util 'first';
use List::MoreUtils qw(uniq);

####################################################################################
#
# Authors: Yi Liao (08/05/2020)
#
# Copyright (C) Not for commercial use
#
# Further computationally validated the merge raw SVs
#
# Prerequisite: TBA package
#
####################################################################################
my $time= strftime "%F", localtime;
my ($glst,$singmaf,$refseq,$tba,$refchr,$refsize,$splitn,$refname,$speciesname,$OUTPUTdir,$head,$help);

GetOptions(
    "glst=s"=>\$glst,
    "singmaf=s"=>\$singmaf,
    "tba=s"=>\$tba,
    "refsize=s"=>\$refsize,
    "splitn=i"=>\$splitn,
    "refname=s"=>\$refname,
    "speciesname=s"=>\$speciesname,
    'refseq=s'=>\$refseq,
    "chr=s"=>\$refchr,
    'outdir=s'=>\$OUTPUTdir,
    'head'=>\$head,
    'help'=>\$help,
);

############# Help information ##########
if ($help){
print <<"END.";
  Usage: perl $0 -glst genome.lst -refsize B73.sizes -singmaf /path/to/sing/ -refseq /path/to/ref -tba /path/to/tba -refname B73 -speciesname maize -chr B73.chr1 -outdir B73_tmp

       -singmaf      path to single coverage maf file                                [REQUIRED]
       -glst         a file with the list of genomes (one line one genome)           [REQUIRED]
       -refseq       reference fasta;only contain the current processed chromosome   [REQUIRED]
       -refsize      reference sizes                                                 [REQUIRED]
       -refname      reference name, e.g. hg38, IRGSP1                               [REQUIRED]
       -speciesname  species name, e.g. human, Drosophila, rice , maize              [REQUIRED]
       -tba          path for tba program                                            [REQUIRED]
       -outdir       define a output dir                                             [REQUIRED]
       -chr          which chromosome to process                                     [REQUIRED]
       -splitn       split the chromosome into n segments [default: 1]               [Optional]
       -head         whether to print the VCF header or not                          [Optional]
       -help         print this information
END.
exit;
}

###########################
$OUTPUTdir ||= ".";
$OUTPUTdir =~ s/\/$//;
$splitn ||= 1;
mkdir($OUTPUTdir) unless(-d $OUTPUTdir);
###########################

open (GENOME,"$glst") or die "Genome/species names do not exist, $!"; ### read genoms/assemblies ID list from a input file (How many genomes are you processed)
my @Genome=<GENOME>;
close GENOME;
chomp @Genome;
my $reference=shift @Genome;
my $ref_glst = \@Genome;


open (SIZE,"$refsize") or die "Genome/species names do not exist, $!"; 
my %fasizes;
while (<SIZE>) {
chomp;
my @temp = split (/\s/,$_);
$fasizes{$temp[0]}=$temp[1];
}

################################

open VCF,">$OUTPUTdir/$refchr.SNP.vcf" or die "$!";

#################   determine whether to print the vcf header  #####################################################
if ($head) {
print OUT "\#\#fileformat=VCFv4.2\n\#\#fileDate=$time\n\#\#source=YiLiao_PopASSYSV_version20200901\n\#\#reference=$refname\n";
open Contig, "$refsize" or die "$!";
while (<Contig>) {
chomp;
my @tmp = split(/\s/,$_);
print OUT "\#\#contig=<ID=$tmp[0],length=$tmp[1],assembly=$refname,md5=XXXXXXXX,species=\"$speciesname\",taxonomy=x>\n";
}
close Contig;
print OUT "##phasing=phased
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">
##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
my $num;
for ($num=0;$num<$#Genome;$num++) {
print VCF "$Genome[$num]\t";
}
print VCF "$Genome[-1]\n";
}


#####################
####### Main ########
#####################

my $length = $fasizes{$refchr};
my $window = int($length/$splitn);

my $i;
for ($i=0;$i<$splitn;$i++) {
my $beg = $i*$window+1;
my $end;
if ($i==$splitn-1) {
  $end = $length;
} else {
  $end = ($i+1)*$window;
}

my %liftover;
my $refbase = &refBase($refseq,$beg,$end);

foreach my $g (@Genome) {
    `${tba}mafFind $singmaf/$reference.$g.sing.maf $beg $end > $OUTPUTdir/$reference.$g.sing.$beg.$end.maf`;
     my $liftover_coord = &liftOver("$OUTPUTdir/$reference.$g.sing.$beg.$end.maf");
     $liftover{$g} = $liftover_coord;
}

my $j;

for ($j=$beg;$j<=$end;$j++) {
my @tmpBase;
    
   foreach my $ind (@Genome) {
     my $key = join ("_",($refchr,$j));
     if (exists ${$liftover{$ind}}{$key}) {
            push (@tmpBase,${$liftover{$ind}}{$key});
     } else {
            push (@tmpBase, "\-");
      }
    }
    
    my @alt = uniq @tmpBase;
    my $element_omitted0 = "-";
    my $element_omitted1 = ${$refbase}{$j};
    my @alt_uniq0 = grep { $_ ne $element_omitted0 } @alt;
    my @alt_uniq1 = grep { $_ ne $element_omitted1 } @alt_uniq0;
   
if ($#alt_uniq1 > -1) {

     my $Qual = 50;
     my $Filter = "PASS";
     my $NS = 0;

     foreach my $ind (@tmpBase) {
       if ($ind ne "\-") {
           $NS++;
       }
     }

     my $DP = 20;
     my @AF;
     foreach my $ele (@alt_uniq1) {
     my $num = 0;
     foreach my $ele1 (@tmpBase) {
     if ($ele1 eq $ele) {
        $num++;
      }
     }
     my $fre = $num/$NS;
     my @tmp = split (//,$fre);
     if ($#tmp>3) {
     $fre = substr($fre,0,4);
     }
     push (@AF,$fre);
     }
     
     my $alt_array = Print_array (@alt_uniq1);
     my $fre_array = Print_array (@AF);
     print VCF "$refchr\t$j\tSNP$j\t${$refbase}{$j}\t$alt_array\t$Qual\t$Filter\tNS=$NS;DP=$DP;AF=$fre_array\tGT:GQ:DP:HQ\t";
     
     my $n;
     for ($n=0;$n<=$#tmpBase;$n++) {
      
     if ($tmpBase[$n] eq ${$refbase}{$j}) {
        if ($n==$#tmpBase) {
        print VCF "0|0:48:10:51,51\n";
        } else {
        print VCF "0|0:48:10:51,51\t";
        }
        
     } elsif ($tmpBase[$n] eq "-") {
        if ($n==$#tmpBase) {
        print VCF "\.|\.:48:10:51,51\n";
        } else {
        print VCF "\.|\.:48:10:51,51\t";
        }
        
     } else {

     my $m;
     for ($m=0;$m<=$#alt_uniq1;$m++){
        my $t = $m + 1;
        if ($tmpBase[$n] eq $alt_uniq1[$m]) {
        if ($n==$#tmpBase) {
        print VCF "$t|$t:48:10:51,51\n";
        } else {
        print VCF "$t|$t:48:10:51,51\t";
        }
      }   
     }
    }
   }
  }    
 }
`rm $OUTPUTdir/*.maf`;
}

##############################################
############    Subroutines    ###############
##############################################
sub Print_array {
my (@array) = @_;
my $newarray =  join ',',@array;
return $newarray;
}
#################################################################
sub refBase {
my ($reffasta,$beg,$end) = @_;
my %hash;
my $seq = &SubSeq($reffasta,$beg,$end);
my @tmp = split(//,$seq);
my $i;
for ($i=0;$i<=$#tmp;$i++) {
my $coord = $beg + $i;
$hash{$coord} = $tmp[$i];
}
my $refhash = \%hash;
return $refhash;
}
#################################################################
sub SubSeq {
my ($seq,$beg,$end) = @_;
local $/="\n>";
open FA, "$seq" or die "$!";
my $seqrange;
   while (my $sequence = <FA>) {
         chomp $sequence;
         my ($id) = $sequence =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
         $sequence =~ s/^>*.+\n//;  # remove FASTA header
         $sequence =~ s/\n//g;  # remove endlines
         my $len = $end - $beg + 1;
         $seqrange = substr $sequence, $beg -1, $len;
    }
return $seqrange;
}
#################################################################
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

my ($target,$Ts,$Textend,$Tstrand,$Tsize,$Tseq0) = (split(/\s+/,$array[0]))[1,2,3,4,5,6];
my ($query,$Qs,$Qextend,$Qstrand,$Qsize,$Qseq0) = (split(/\s+/,$array[1]))[1,2,3,4,5,6];
my $Tseq =uc($Tseq0);
my $Qseq =uc($Qseq0);
my @targetSEQ = split(//,$Tseq);
my @querySEQ = split (//,$Qseq);
my $i;
my $n=0;

for ($i=0;$i<=$#targetSEQ;$i++) {
   if ($targetSEQ[$i]!~/\-/) {
     $n++;
  my $tcoord =  $Ts + $n;
  my $qbase = $querySEQ[$i];
  my $tinf = join ("_",($target,$tcoord));
  $liftover{$tinf} = $qbase;
   }
}

while (<MAF>) {
next if ($_=~/eof/);
my @tmp = split (/\n/,$_);
my ($target1,$Ts1,$Textend1,$Tstrand1,$Tsize1,$Tseq1) = (split(/\s+/,$tmp[1]))[1,2,3,4,5,6];
my ($query1,$Qs1,$Qextend1,$Qstrand1,$Qsize1,$Qseq1) = (split(/\s+/,$tmp[2]))[1,2,3,4,5,6];
$Tseq1 = uc ($Tseq1);
$Qseq1 = uc ($Qseq1);
my @targetSEQ1 = split(//,$Tseq1);
my @querySEQ1 = split (//,$Qseq1);

my $j;
my $m=0;
for ($j=0;$j<=$#targetSEQ1;$j++) {
  if ($targetSEQ1[$j]!~/\-/) {
       $m++;
       my $tcoord1 =  $Ts1 + $m;
       my $qbase1 = $querySEQ1[$j];
       my $tinf1 = join ("_",($target1,$tcoord1));
       $liftover{$tinf1} = $qbase1;
  }
 }
}
my $refhash = \%liftover;
return $refhash;
}
#################################################################
