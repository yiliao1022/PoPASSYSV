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

open SV, "$SVmerge_raw" or die "$!";
open VCF, ">$SVmerge_raw.vcf" or die "$!";
open ERR, ">sv.err" or die "$!";
#################   determine whether to print the vcf header  #####################################################
if ($header) {
print VCF "\#\#fileformat=VCFv4.2\n\#\#fileDate=$time\n\#\#source=YiLiao_PopASSYSV_version20200901\n\#\#reference=$ref\n";
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
print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
}
#####################################################################################################################

my $ref1 = shift @genome_lst;
my $g;
for ($g=0;$g<$#genome_lst;$g++) {
print VCF "$genome_lst[$g]\t";
}
print VCF "$genome_lst[-1]\n";

#####################################################################################################################

if ($svtype eq "INS") {
####################################################
  while (<SV>) {
     next if $_=~/\#/;
     chomp;
     my ($chr,$s,$e,$length,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,8,9,10,11,12];
     my $beg = $s - $len;
     my $end = $e + $len;
     `${tba}mafFind $tbamaf $beg $end $chr > $chr.$s.sub.maf`;
 
      my $genome_array = &Maf2fasta("$chr.$s.sub.maf",$qchr,$genome);
     `${muscle}muscle3.8.31_i86linux64 -in $chr.$s.sub.maf.fa -out $chr.$s.mulcle.fa -maxiters 2`;
      &Muscle2fasta ($genome_array,"$chr.$s.mulcle.fa");
     `mv $chr.$s.mulcle.fa.reformed.fa $chr.$s.tmp.fa`;
      my ($tag,$hashref) = &ParseMaf2fastaINS("$chr.$s.tmp.fa",$qs,$qe,$ref_glst,$length);

  if ($tag eq "Sucess") {
     my $end0 = $s+$length;
     print VCF "$chr\t$s\tSV$s\tA\t\<$svtype\>\t22\tPASS\tSVTYPE=INS;END=$end0;SVLEN=$length\tGT:GQ\t";       
     my $m;

     for ($m=0;$m<$#genome_lst;$m++) {
       print VCF "${$hashref}{$genome_lst[$m]}|${$hashref}{$genome_lst[$m]}:22\t";
     }
       print VCF "${$hashref}{$genome_lst[$#genome_lst]}|${$hashref}{$genome_lst[$#genome_lst]}:22\t";
              
       } else {
           print ERR "$chr\t$s\t$e\tFailed\n";
       }     
           print VCF "\n";
  }
####################################################

} elsif ($svtype eq "DEL") {

####################################################
    while (<SV>) {
       next if $_=~/\#/;
       chomp;
       my ($chr,$s,$e,$length,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,8,9,10,11,12];
       my $beg = $s-$len;
       my $end = $e+$len;
       `${tba}maf2fasta $genome/$chr $tbamaf $beg $end > $chr.TBA.$beg.$end.fasta`;

       if (-z "$chr.TBA.$beg.$end.fasta") { ## -z The file exists and has zero size
           system "rm $chr.TBA.$beg.$end.fasta";
           print ERR "$_\tFailed!\n";
           next;
        }

       my $ref = &ParseMaf2fastaDEL("$chr.TBA.$beg.$end.fasta",$len,$ref_glst,$length);
       my $end0 = $s+$length;
       print VCF "$chr\t$s\tSV$s\tA\t\<$svtype\>\t22\tPASS\tSVTYPE=INS;END=$end0;SVLEN=\-$length\tGT:GQ\t";       
       my $m;

       for ($m=0;$m<$#genome_lst+1;$m++) {
         print VCF "${$ref}{$genome_lst[$m]}|${$ref}{$genome_lst[$m]}:22\t";
       }

        print VCF "\n";
     #  `rm $chr.TBA.$beg.$end.fasta`;
    }
####################################################

} elsif ($svtype eq "CNV") {


} elsif ($svtype eq "INV") {


} else {
  print "Please set the sv type!";
}



#########################################
############# Sub modules ###############
#########################################
sub Maf2fasta {
my($maf,$ref,$g)=@_;
local $/="\n\n"; # record separator as a blank line when reading
open MAF, "$maf" or die "Missing maf file!";
open OUT, ">$maf.fa" or die "Cann't write the output file!";
my %hash;
while (<MAF>) {
chomp;
my @temp = split (/\n/,$_);
foreach my $line (@temp) {
  if ($line =~/^s/) {
    my ($spec,$start,$extend,$strand,$qsize) = (split(/\s+/,$line))[1,2,3,4,5];
   # $spec=(split(/\./,$spec))[0];
    my $segment = join ("_",($start,$extend,$strand,$qsize));
    push @{$hash{$spec}},$segment;
    }
  }
}

my %interval;
foreach my $species (keys %hash) {
  my $five;
  my $three;
  my $strand;
  my @plus;
  my @minus;
   foreach my $seg (@{$hash{$species}}) {
     my ($s,$ext,$or,$size)=(split(/_/,$seg))[0,1,2,3];
     my $e = $s + $ext;
     if ($or=~/\+/) {
      push (@plus,$s);
      push (@plus,$e);
     } elsif($or=~/-/){
      my $srev = $size - $s;
      my $erev = $size - $e;
      push (@minus,$srev);
      push (@minus,$erev);
     }
  }

  if ($#plus>0 and $#minus<0) {
    my @plus_array = sort { $a <=> $b } @plus;
    $five = $plus_array[0];
    $three = $plus_array[-1];
    $strand = "plus";
  } elsif ($#plus<0 and $#minus>0) {
    my @minus_array = sort { $a <=> $b } @minus;
    $five = $minus_array[0];
    $three = $minus_array[-1];
    $strand = "minus";
  }  elsif ($#minus>0 and $#plus>0) {
  my @plus_array = sort { $a <=> $b } @plus;
  my @minus_array = sort { $a <=> $b } @minus;
  my $plus_len = $plus_array[-1]-$plus_array[0];
  my $minus_len = $minus_array[-1]-$minus_array[0];
    if ($plus_len >= $minus_len) {
      $five = $plus_array[0];
      $three = $plus_array[-1];
      $strand = "plus";
    } else {
      $five = $minus_array[0];
      $three = $minus_array[-1];
      $strand = "minus";
    }
  } else {
    print "There is no sequences available for this species!";
  }
   my $coordinates = join ("_",($five,$three,$strand));
   $interval{$species} = $coordinates;
}

my @genomes;
#foreach my $keys (keys %interval) {
 #   push(@genomes,$keys);
#}

if (exists $interval{$ref}) {
  my ($s,$e,$strand) = (split (/\_/,$interval{$ref}))[0,1,2];
  print "$ref\t$s\t$e\t$strand\n";
  my $genome = (split (/\./,$ref))[0];
  my $seq = &SubSeq ("$g/$genome",$ref,$s,$e,$strand);
  print OUT ">$ref\_$s\_$e\n$seq\n";
    push (@genomes,"$ref\_$s\_$e");
    foreach my $spec_seq (keys %interval) {
            next if ($spec_seq eq $ref);
            my ($s,$e,$strand) = (split (/\_/,$interval{$spec_seq}))[0,1,2];
            print "$spec_seq\t$s\t$e\t$strand\n";
            my $genome = (split (/\./,$spec_seq))[0];
            my $seq = &SubSeq ("$g/$genome",$spec_seq,$s,$e,$strand);
            print OUT ">$spec_seq\_$s\_$e\n$seq\n";
            push (@genomes,"$spec_seq\_$s\_$e");
      }
    }

my $refer = \@genomes;
return $refer;
}

### Extract sequence based on coordinate
sub SubSeq {
my ($seq,$chr,$beg,$end,$strand) = @_;
local $/="\n>";
open FA, "$seq" or die "$!";
my $seqrange;
#open Out, ">$chr\_$beg\_$end";
   while (my $sequence = <FA>) {
         chomp $sequence;
         my ($id) = $sequence =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
#         print "$id\n";
if ($id eq $chr) {
         $sequence =~ s/^>*.+\n//;  # remove FASTA header
         $sequence =~ s/\n//g;  # remove endlines
         my $len = $end - $beg + 1;
         $seqrange = substr $sequence, $beg -1, $len;
#         print Out ">$chr\_$beg\_$end\n$seqrange\n";
         last;
        }
    }

my $sequence;
if ($strand eq "plus") {
$sequence = $seqrange;
} elsif ($strand eq "minus") {
$sequence = &RevCom ($seqrange);
}
return $sequence;
}

### Reverse and complement the sequence
sub RevCom {
my $seq = shift;
my $revcomp = reverse $seq;
$revcomp =~ tr/ATGCNatgcn/TACGNtacgn/;
return $revcomp;
}

### convert the TAB fasta format to normal fasta file for MUSCLE
sub ReadMaf2fasta {
my $seq = shift;
my @sequences;
my @genomes;
local $/="\n";
open In,"$seq" or die "$!";
open Out,">$seq.realign.fa" or die "$!";
  while (<In>) {
       chomp;
       push (@sequences,$_)
  }
my $i;
my @unit = split (/:/,$sequences[1]);
$sequences[1] = $unit[0];
for ($i=1;$i<=$#sequences/2;$i++) {
push (@genomes,$sequences[$i]);
$sequences[$i+$#sequences/2] =~s/\-//g;
my $seq_tmp=uc ($sequences[$i+$#sequences/2]);
print Out ">$sequences[$i]\n$seq_tmp\n";
 }
my $ref = \@genomes;
return $ref;
}

##### convert the MUSCLE output format to TAB fasta format
sub Muscle2fasta {
my ($genomes,$muscle) =@_;
local $/="\n>";
my %hash;
open In, "$muscle" or die "$!";
open Out, ">$muscle.reformed.fa" or die "$!";
my $n = $#{$genomes}+1;
print Out "$n 1000\n";
while (my $seq = <In>) {
chomp $seq;
my ($id) = $seq =~ /^>*(\S+)/;
print "$id\n";
$seq =~ s/^>*.+\n//;
$seq =~ s/\n//g;
$hash{$id} = $seq;
}
foreach my $g (@{$genomes}) {
    print Out "$g\n";
  }

foreach my $gid (@{$genomes}) {
    print Out "$hash{$gid}\n";
  }
}


########### Parse Insertions from maf2fasta output

sub ParseMaf2fastaINS {
my ($file,$qs,$qe,$genome,$len) = @_;
my %hash;
my %hash_identity;
my %hash_deletion_overlap;

local $/="\n";
open (IN,'<:encoding(UTF-8)',$file) or die "$!\n";
chomp(my @data = <IN>);
close IN;

my @seq_num = split (/\s/,$data[0]);
my @ref_name = split (/\./,$data[1]);
my @ref_se = split (/\_/,$data[1]);
my $ref_s = $ref_se[1];
my $ref_e = $ref_se[2];
my $extend5 = $qs - $ref_s;
my $extend3 = $ref_e - $qe;
my $tag;
my %hash_deletion_value=();

if ($extend5 < 0 || $extend3 <0) {
$tag = "Failed";
} else {
$tag = "Sucess";
my $ref = $ref_name[0];
my @ref_seq = split(//,$data[$seq_num[0]+1]);
my $ref_len = length ($data[$seq_num[0]+1]);

$hash{$ref}=\@ref_seq;
my @coordinate_proxy;
my $len_n=0;
my $t;
for ($t=0;$t<$ref_len;$t++) {
 if ($hash{$ref}->[$t] =~/[ATGCatgc]/) {
  $len_n++;
  my $s;
  for ($s=$extend5+1;$s<$extend5+$len+1;$s++) {
  if ($len_n == $s) {
  push (@coordinate_proxy,$t);
   }
  }
 }
}

#### tag0 Read queries sequence to hash
my $i;
for ($i=2;$i<$seq_num[0]+1;$i++) {
my @query_seq = split(//,$data[$seq_num[0]+$i]);
my @unit = split (/\./,$data[$i]);
   $hash{$unit[0]} = \@query_seq;
}
#### tag0

###### tag1  calculate the rate of flanking sequence that is aligned.
foreach my $ind (@{$genome}) {
     if (exists $hash{$ind}) {
       my $j;
       my $n;
       my $m;
       my $start = $coordinate_proxy[0];
       my $end = $coordinate_proxy[$#coordinate_proxy];

for ($j=0;$j<$start;$j++) {
     if ($hash{$ref}->[$j] =~/[ATGCatgc]/) {
         $n++;
         $hash{$ref}->[$j] = uc ($hash{$ref}->[$j]);
         if ($hash{$ind}->[$j] =~/[ATGCatgc]/) {
              $hash{$ind}->[$j] = uc ($hash{$ind}->[$j]);
                  if ($hash{$ref}->[$j] eq $hash{$ind}->[$j]) {
                        $m++;
                  }
           }
      }
 }

for ($j=$end+1;$j<$ref_len;$j++) {
     if ($hash{$ref}->[$j] =~/[ATGCatgc]/) {
          $n++;
          $hash{$ref}->[$j] = uc ($hash{$ref}->[$j]);
        if ($hash{$ind}->[$j] =~/[ATGCatgc]/) {
                $hash{$ind}->[$j] = uc ($hash{$ind}->[$j]);
            if ($hash{$ref}->[$j] eq $hash{$ind}->[$j]) {
                         $m++;
           }
        }
    }
 }

$hash_identity{$ind} = $m/$n+0.1;
} else {
$hash_identity{$ind} = 0;
}
}
####### tag1

####### tag2 calculate the rate of deletion sequence present between reference and the query, 0 indicates a deletion, while 1 indicate no deletion
   foreach my $ind1 (@{$genome}) {
           my $n1=0;
           my $m1=0;
     foreach my $coordin (@coordinate_proxy) {
              $n1++;
         if (exists $hash{$ind1} and ${$hash{$ind1}}[$coordin]=~/[ATGCatgc]/) {
              $m1++;
         }
    }
      $hash_deletion_overlap{$ind1} = $m1/$n1+0.1;
    }
####### tag2

###### tag3 determine whether the deletion occur in each sample based on both the flanking region alignment and the deleted region.
  foreach my $ind2 (@{$genome}) {
        if ($hash_identity{$ind2} < 0.5) {
           $hash_deletion_value{$ind2}=2;
      } elsif ($hash_identity{$ind2}>=0.5 && $hash_deletion_overlap{$ind2} >= 0.5) {
           $hash_deletion_value{$ind2}=1;
      } elsif ($hash_identity{$ind2}>=0.5 && $hash_deletion_overlap{$ind2} < 0.5) {
           $hash_deletion_value{$ind2}=0;
      }
   }
####### tag3

}
my $reference = \%hash_deletion_value;
return ($tag,$reference);
}

###############################################################
################## Parse Deletions for maf2fasta output #######
###############################################################

sub ParseMaf2fastaDEL {
my ($file,$extend,$genome,$len) = @_;
my %hash;
my %hash_identity;
my %hash_deletion_overlap;

local $/="\n";
open (IN,'<:encoding(UTF-8)',$file) or die "$!\n";
chomp(my @data = <IN>);
close IN;

my @seq_num = split (/\s/,$data[0]);
my @ref_name = split (/:/,$data[1]);
my $ref = $ref_name[0];
my @ref_seq = split(//,$data[$seq_num[0]+1]);
$hash{$ref}=\@ref_seq;

my @coordinate_proxy;
my $len_n=0;
my $t;
for ($t=0;$t<$seq_num[1];$t++) {
 if ($hash{$ref}->[$t] =~/[ATGCatgc]/) {
  $len_n++;
  my $s;
  for ($s=$extend+1;$s<$extend+$len+1;$s++) {
  if ($len_n == $s) {
  push (@coordinate_proxy,$t);
   }
  }
 }
}

my $i;
for ($i=2;$i<$seq_num[0]+1;$i++) {
my @query_seq = split(//,$data[$seq_num[0]+$i]);
$hash{$data[$i]} = \@query_seq;
}

foreach my $ind (@{$genome}) {
 if (exists $hash{$ind}) {
###
my $j;
my $n;
my $m;
my $start = $coordinate_proxy[0];
my $end = $coordinate_proxy[$#coordinate_proxy];

for ($j=0;$j<$start;$j++) {
     if ($hash{$ref}->[$j] =~/[ATGCatgc]/) {
         $n++;
         $hash{$ref}->[$j] = uc ($hash{$ref}->[$j]);
       if ($hash{$ind}->[$j] =~/[ATGCatgc]/) {
           $hash{$ind}->[$j] = uc ($hash{$ind}->[$j]);
              if ($hash{$ref}->[$j] eq $hash{$ind}->[$j]) {
                 $m++;
              }
        }
     }
   }

for ($j=$end+1;$j<$seq_num[1];$j++) {
     if ($hash{$ref}->[$j] =~/[ATGCatgc]/) {
         $n++;
         $hash{$ref}->[$j] = uc ($hash{$ref}->[$j]);
       if ($hash{$ind}->[$j] =~/[ATGCatgc]/) {
           $hash{$ind}->[$j] = uc ($hash{$ind}->[$j]);
              if ($hash{$ref}->[$j] eq $hash{$ind}->[$j]) {
                 $m++;
              }
        }
     }
   }

$hash_identity{$ind} = $m/$n+0.1;
###
} else {
    $hash_identity{$ind} = 0;
 }
}

foreach my $ind1 (@{$genome}) {
   my $n1=0;
   my $m1=0;
      foreach my $coordin (@coordinate_proxy) {
         $n1++;
         if (exists $hash{$ind1} and ${$hash{$ind1}}[$coordin]=~/[ATGCatgc]/) {
         $m1++;
         }
       }
$hash_deletion_overlap{$ind1} = $m1/$n1+0.1;
}

my %hash_deletion_value;
foreach my $ind2 (@{$genome}) {
   if ($hash_identity{$ind2} < 0.5) {
       $hash_deletion_value{$ind2}=2; ### missing data
    } elsif ($hash_identity{$ind2}>=0.5 && $hash_deletion_overlap{$ind2} > 0.5) {
       $hash_deletion_value{$ind2}=0; ### similar to reference
    } elsif ($hash_identity{$ind2}>=0.5 && $hash_deletion_overlap{$ind2} < 0.5) {
       $hash_deletion_value{$ind2}=1; ### support a deletion
    }
}
my $reference = \%hash_deletion_value;
return $reference ;
}
