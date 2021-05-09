#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Long;

####################################################################################
# SVCombine_evaluated.pl
#
# Authors: Yi Liao (07/01/2020)
#
# Copyright (C) Not for commercial use
#
# Further computationally validated the merge raw SVs
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
my ($SVmerge_raw,$svtype,$muscle,$glst,$tba,$genome,$len,$tbamaf,$fasize,$species,$ref,$header,$Help);

GetOptions(
  "sv=s"=>\$SVmerge_raw,
  "gseq=s"=>\$genome,
  "svtype=s"=>\$svtype,
  "fasize=s"=>\$fasize,
  "spec=s"=>\$species,
  "ref=s"=>\$ref,
  "muscle=s"=>\$muscle,
  "glst=s"=>\$glst,
  "tba=s"=>\$tba,
  "len=i"=>\$len,
  "tbamaf=s"=>\$tbamaf,
  "head"=>\$header,
  "help"=>\$Help
);

$muscle ||= '';
$tba ||= '';
$len ||= 100;
##################################################################################################
if ($Help){
print <<"END.";
  Usage: perl $0 -sv SV.txt -svtype INS -muscle /path/to/muscle -tba /path/to/tba -tbamaf tba.maf -gseq /path/to/genomes -len 200
  -sv         raw merged SV file                                      [REQUIRED]
  -svtype     SV types [INS|DEL|CNV|INV]                              [REQUIRED]
  -gseq       path to all genome sequences                            [REQUIRED]
  -tbamaf     tba maf file (multiple genome alignment)                [REQUIRED]
  -glst       a file with the list of genomes (one line one genome)   [REQUIRED]
  -fasize     reference sizes
  -ref        reference name, e.g. hg38, IRGSP1
  -spec       species name, e.g. human, Drosophila, rice , maize
  -muscle     path for muscle program
  -tba        path for tba program
  -len        sequence length extend from the breakpoint  [default:100]
  -help       print this information
END.
exit;
}
################
$/="\n";
open (GLST,'<:encoding(UTF-8)',$glst) or die "$!\n";
chomp(my @genome_lst = <GLST>);
close GLST;
my $ref_glst=\@genome_lst;

#################   determine whether to print the vcf header  #####################################################
open SV, "$SVmerge_raw" or die "$!";
open VCF, ">$SVmerge_raw.vcf" or die "$!";
open ERR, ">sv.err" or die "$!";

if ($header) {
print VCF "\#\#fileformat=VCFv4.2\n\#\#fileDate=$time\n\#\#source=PopASSYSV_version20200901\n\#\#reference=$ref\n";
open Contig, "$fasize" or die "$!";
while (<Contig>) {
chomp;
my @tmp = split(/\s/,$_);
print VCF "\#\#contig=<ID=$tmp[0],length=$tmp[1],assembly=$ref,md5=XXXXXXXX,species=\"$species\",taxonomy=x>\n";
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
}
#####################################################################################################################
print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
my $ref1 = shift @genome_lst;
my $g;
for ($g=0;$g<$#genome_lst;$g++) {
print VCF "$genome_lst[$g]\t";
}
print VCF "$genome_lst[-1]\n";


################ Genotying start ......
if ($svtype eq "INS") {
####################################################
  while (<SV>) {
     next if $_=~/\#/;
     chomp;

eval {
     my ($chr,$s,$e,$length,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,8,9,10,11,12];
     my $beg = $s - $len;
     my $end = $e + $len;
     `${tba}mafFind $tbamaf $beg $end $chr > $chr.$s.sub.maf`;
     my ($coord,$beg_index,$end_index) = &Maf2fasta("$chr.$s.sub.maf",$qchr,$qs,$qe);
     my $end0 = $s + $length;
     print VCF "$chr\t$s\tSV$s\tN\t\<$svtype\>\t22\tPASS\tSVTYPE=$svtype;END=$end0;SVLEN=$length\tGT:GQ\t";       
     my $m;

     for ($m=0;$m<$#genome_lst;$m++) {
       if (exists ${$coord}{$genome_lst[$m]}{$beg_index} and exists ${$coord}{$genome_lst[$m]}{$end_index}) {
          my $check = abs (${$coord}{$genome_lst[$m]}{$beg_index}-${$coord}{$genome_lst[$m]}{$end_index});
          if ($check < 5) {
            print VCF "0|0:22\t";
          } elsif ($check>0.8*$length and $check<1.2*$length) {      
            print VCF "1|1:22\t";
          } else {
            print VCF "./.:22\t"
          }
       } else {
            print VCF "./.:22\t"       
       }       
     }

     if (exists ${$coord}{$genome_lst[$#genome_lst]}{$beg_index} and exists ${$coord}{$genome_lst[$#genome_lst]}{$end_index}) {
        my $check = abs (${$coord}{$genome_lst[$#genome_lst]}{$beg_index}-${$coord}{$genome_lst[$#genome_lst]}{$end_index});
        if ($check < 5) {
          print VCF "0|0:22";
        } elsif ($check>0.8*$length and $check<1.2*$length) {      
          print VCF "1|1:22";
        } else {
          print VCF "./.:22"
        }
     } else {
          print VCF "./.:22"       
     } 
     print VCF "\n";
     `rm $chr.$s.sub.maf`;           
};

if ($@) {
  print ERR "Failed\t$_\n";
}
    
}
####################################################
}

##############################################
############## Subroutines ###################
##############################################

sub Maf2fasta {
my($maf,$ref,$qs,$qe)=@_;
############################################ section1
open MAF, "$maf" or die "Missing maf file!";
my %coord_proxy;
my $s=1;
my $n=0;
my $ref_end=0;
my $gap_ref = &GetRefCoord("$maf",$ref);
$ref = (split(/\./,$ref))[0];
local $/="\n\n"; # record separator as a blank line when reading
my $ref_strand;
while (<MAF>) {
chomp;
next if ($_=~/eof/);
my  @temp = split(/\n/,$_);
$s = $s + $ref_end + ${$gap_ref}[$n];
foreach my $line (@temp) {
  if ($line =~/^s/) {
    my ($spec,$start,$extend,$strand,$qsize,$seq) = (split(/\s+/,$line))[1,2,3,4,5,6];
    $spec = (split(/\./,$spec))[0];
    if ($spec =~/$ref/) {
    $ref_strand = $strand;
    print "reference strand: $ref_strand \n";
    }   
    my $link =&links($s,$start,$strand,$qsize,$seq); 
    if (!keys %{$coord_proxy{$spec}}) {
    $coord_proxy{$spec}=$link;
    } else {
    my %linkcombine=(%{$coord_proxy{$spec}},%{$link});
    $coord_proxy{$spec}=\%linkcombine;
     }
    }
   }
  $n++;
  my @array = sort {$a<=>$b} (keys %{$coord_proxy{$ref}});
  $ref_end = $array[-1];  
  }
  my $beg_index;
  my $end_index; 
  if ($ref_strand =~/\+/) {
  my $begin = $qs-1;
  my @all_matches = grep {${$coord_proxy{$ref}}{$_} eq $begin} keys %{$coord_proxy{$ref}};
  $beg_index = $all_matches[0];
  my @all_matches1 = grep {${$coord_proxy{$ref}}{$_} eq $qe} keys %{$coord_proxy{$ref}};
  $end_index = $all_matches1[0];
  } elsif ($ref_strand =~/\-/) {
  my $end = $qe+1;
  my @all_matches = grep {${$coord_proxy{$ref}}{$_} eq $qs} keys %{$coord_proxy{$ref}};
  $beg_index = $all_matches[0];
  my @all_matches1 = grep {${$coord_proxy{$ref}}{$_} eq $end} keys %{$coord_proxy{$ref}};
  $end_index = $all_matches1[0];
  }
  print "$beg_index\t$end_index\n"; 
  my $ref_hash = \%coord_proxy;
  return ($ref_hash,$beg_index,$end_index);
}

sub GetRefCoord {
my ($maf,$ref) = @_;
my @gaps;
my $end;
local $/="\n";
open In, "$maf" or die "$!";
  while (<In>) {
  chomp;
  if ($_=~/$ref/ and $_=~/^s/) {
    my ($start,$extend) = (split(/\s+/,$_))[2,3];
    if ($end) {
     my $gap = $start - $end;
     push (@gaps,$gap);
   } else {
     push (@gaps,"0");
   }
    $end = $start + $extend;
   }
  } 
  my $gapref = \@gaps;
  return $gapref; 
}

sub links {
my ($s,$beg,$strand,$qsize,$seq) = @_;
my %hash;
my @tmp = split(//,$seq);
my $i;
my $j=0;
 for ($i=0;$i<=$#tmp;$i++) {
    if ($strand =~/\+/) {
        if ($tmp[$i]=~/[ATGCNatgcn]/) {
            $hash{$i+$s} = $j + $beg;
            $j++;
        }  else {
            $hash{$i+$s} = $j + $beg;
        }
    } elsif ($strand =~/\-/) {
        if ($tmp[$i]=~/[ATGCNatgcn]/) {
          $hash{$i+$s} = $qsize - $beg - $j;
          $j++;
        }  else {
          $hash{$i+$s} = $qsize - $beg - $j;
        }
    }
 }
my $reflinks = \%hash; 
return $reflinks;
}
