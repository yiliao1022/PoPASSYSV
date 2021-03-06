#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use List::MoreUtils qw(uniq);

####################################################################################
# MAF2VCF_SNP.pl
#
# Authors: Yi Liao (07/11/2020)
#
# Copyright (C) Not for commercial use
#
# Report SNP from .maf file to .vcf file
#
# Prerequisite : TBA tools; MUSCLE
#
# Usage: perl $0 -sv SV.txt -svtype INS -muscle /path/to/muscle -tba /path/to/tba
# -sv         raw merged SV file                        [REQUIRED]
# -svtype     SV types [INS|DEL|CNV|INV]                [REQUIRED]
# -gseq       path to all genome sequences              [REQUIRED]
# -muscle     path for muscle program
# -tba        path for tba program
# -help       print this information
#
####################################################################################

my $time= strftime "%F", localtime;
my ($input,$fasize,$ref,$query,$species,$glst,$output,$header,$Help);

GetOptions (
"input:s" =>\$input,
"fasize:s" =>\$fasize,
"ref:s" => \$ref,
"query:s" => \$query,
"spec:s" => \$species,
"glst:s"=>\$glst,
"output:s" =>\$output,
"head"=>\$header,
"help" =>\$Help
);

if ($Help) {
print << "END.";
Useage: perl $0 -input maf2fasta.fa -spec rice/human/fruitfly/maize -ref IRGSP1 -fasize chr.sizes -glst glst.txt -head

-input      : maf2fasta files
-fasize     : reference sizes by chromosome
-spec       : species name, e.g. rice, human,maize
-ref        : reference name, e.g. hg38, ISO1, IRGSP1
-glst       : individuals in the population, the first one should be the reference
-head       : if defined then print the vcf header eitherwise not 

END.
exit;
}

open OUT, ">$input.vcf" or die "$!";



#################   determine whether to print the vcf header  #####################################################
if ($header) {
print OUT "\#\#fileformat=VCFv4.2\n\#\#fileDate=$time\n\#\#source=YiLiao_PopASSYSV_version20200901\n\#\#reference=$ref\n";
open Contig, "$fasize" or die "$!";
while (<Contig>) {
chomp;
my @tmp = split(/\s/,$_);
print OUT "\#\#contig=<ID=$tmp[0],length=$tmp[1],assembly=$ref,md5=XXXXXXXX,species=\"$species\",taxonomy=x>\n";
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
}
####################################################################################################################


open GLST, "$glst" or die "$!";
chomp (my @genome=<GLST>);
close GLST;

#my @genome = ("Nipp", "Aaaa", "Bbbb", "Cccc", "Dddd", "Eeee", "Gggg", "IR64", "Jjjj", "Kkkk", "Llll", "MH63", "Mmmm", "NN22", "Nnnn", "Oooo", "R498", "Ssss", "ZS97", "Basm", "N311", "DomS");


my ($sequence,$coordinate,$index)=&Parse_MAF2FastaOutput($input);

my $ref1 = shift @genome;
my $g;
for ($g=0;$g<$#genome;$g++) {
print OUT "$genome[$g]\t";
}
print OUT "$genome[-1]\n";

my $a;
my $size = @$coordinate;
for ($a=0;$a<$size;$a++) {
     my $flag =0;
     my @alt;
     foreach my $ind (@genome) {
     $$sequence{$ref1}->[$$coordinate[$a]] = uc ($$sequence{$ref1}->[$$coordinate[$a]]);
     $$sequence{$ind}->[$$coordinate[$a]] = uc ($$sequence{$ind}->[$$coordinate[$a]]);
     if ($$sequence{$ref1}->[$$coordinate[$a]] ne $$sequence{$ind}->[$$coordinate[$a]]) {
     $flag = 1;
     push (@alt,$$sequence{$ind}->[$$coordinate[$a]]);
     }
     }
     my @alt_uniq = uniq @alt;
     my $element_omitted = "-";
     @alt_uniq = grep { $_ ne $element_omitted } @alt_uniq;
     if ($flag==1 and $#alt_uniq > -1) {
     my $Qual = 50;
     my $Filter = "PASS";
     my $NS = @genome;
     my $DP = 20;
     my @AF;
     foreach my $ele (@alt_uniq) {
     my $num = 0;
     foreach my $ele1 (@alt) {
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
     my $alt_array = Print_array (@alt_uniq);
     my $fre_array = Print_array (@AF);
     print OUT "Chr01\t$$index{$$coordinate[$a]}\tSNP$$index{$$coordinate[$a]}\t$$sequence{$ref1}->[$$coordinate[$a]]\t$alt_array\t$Qual\t$Filter\tNS=$NS;DP=$DP;AF=$fre_array\tGT:GQ:DP:HQ\t";

     my $m;
     for ($m=0;$m<$#genome;$m++) {

     if ($$sequence{$genome[$m]}->[$$coordinate[$a]] eq $$sequence{$ref1}->[$$coordinate[$a]]) {
     print OUT "0|0:48:10:51,51\t";
     } elsif ($$sequence{$genome[$m]}->[$$coordinate[$a]] eq $alt_uniq[0]) {
     print OUT "1|1:48:10:51,51\t";
     } elsif ( $#alt_uniq >0 and ($$sequence{$genome[$m]}->[$$coordinate[$a]] eq $alt_uniq[1])) {
     print OUT "2|2:48:10:51,51\t";
     } else {
     print OUT ".\t";
     }
    }
     if ($$sequence{$genome[-1]}->[$$coordinate[$a]] eq $$sequence{$ref1}->[$$coordinate[$a]]) {
     print OUT "0|0:48:10:51,51\n";
     } elsif ($$sequence{$genome[-1]}->[$$coordinate[$a]] eq $alt_uniq[0]) {
     print OUT "1|1:48:10:51,51\n";
     } elsif ( $#alt_uniq >0 and ($$sequence{$genome[-1]}->[$$coordinate[$a]] eq $alt_uniq[1])) {
     print OUT "2|2:48:10:51,51\n";
     } else {
     print OUT ".\n";
     }



   }

}



##################################
############ Subroutine ##########
##################################

sub Print_array {
my (@array) = @_;
my $newarray =  join ',',@array;
return $newarray;
}

sub Parse_MAF2FastaOutput {
   my ($file) = @_;
   my @data;
   my %hash;
   $/="\n";
   open(DATA, "$file") or die "Unable to open file maf2fasta file, $!";
   while (<DATA>) {
         chomp;
       push (@data, $_);
   }
   close DATA;
## load reference sequence into hash
   my @seq_num = split (/\s/,$data[0]);
   my @ref_name = split (/:/,$data[1]);
   my $ref = $ref_name[0];
   my @ref_seq = split(//,$data[$seq_num[0]+1]);
   $hash{$ref}=\@ref_seq;
 ##
## load query sequence into hash
   my $i;
   for ($i=2;$i<$seq_num[0]+1;$i++) {
   my @query_seq = split(//,$data[$seq_num[0]+$i]);
      $hash{$data[$i]} = \@query_seq;
   }
 ##
## relate the sequence coordinates to Maf2fasta reference coordinates
   my $len_n=1;
   my %index_proxy;
   my $t;
   my @coordinate_proxy;
   for ($t=0;$t<$seq_num[1];$t++) {
       if ($hash{$ref}->[$t] =~/[ATGCatgc]/) {
       $index_proxy{$t}=$len_n;
       $len_n++;
       push (@coordinate_proxy,$t);
       }
   }
##
my $seq=\%hash;
my $cor_proxy=\@coordinate_proxy;
my $index_convert=\%index_proxy;
   return ($seq,$cor_proxy,$index_convert);
}
