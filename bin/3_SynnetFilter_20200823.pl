#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(&basename &dirname);

####################################################################################
# SynNetFilter.pl
#
# Authors: Yi Liao (06/19/2020)
#
# Copyright (C) Not for commercial use
#
# Prerequisite : Perl, Kent's Utilities, TBA program.
#
# Remove non-syntenic fills from *.chain.filter.tnet.synnet, only keep reliable syntenic and conserved orthologous alignments.
#
####################################################################################

my ($input,$output,$one2mul,$toplen,$seclen,$nonsynlen,$synlen,$invlen,$alignrate,$chain,$tba,$kent,$target,$rawsv,$query,$qFarCutoff,$Help);

# get input parameters
GetOptions( 'input=s'=>\$input,
            'output=s'=>\$output,
            'one2mul=s'=>\$one2mul,
            'alrate=f'=>\$alignrate,
            'toplen=i'=>\$toplen,
            'nonsynlen=i'=>\$nonsynlen,
            'synlen=i'=>\$synlen,
            'invlen=i'=>\$invlen,
            'qFarCutoff=i'=>\$qFarCutoff,
	        'kentUtil=s'=>\$kent,
	        'rawsv=s'=>\$rawsv,
            'tba=s'=>\$tba,
	        'chain=s'=>\$chain,
	        'target=s'=>\$target,
	        'query=s'=>\$query,
            'help'=>\$Help
           );

#### Printing help information
if (defined ($Help) || !defined($input) || !defined($output)){
print <<"END.";
Usage: perl $0 --input <SynNet file> --output SybNet.filter --alrate 0.5 --topCutoff 1000000 --secCutoff 500000 --qFarCutoff 5000000 --target ~/genome/B73 --query ~/genome/NC350 --kentUtil /path/to/kent --tba /path/to/tba --chain /path/to/all.chain.filter -rawsv rawsv.txt

Options
--input         <s> : Chain/Net/Synnet file
--output        <s> : output filtered Chain/Net/Synnet file
--one2mul       <s> : List of sequences that are aligned to the reference in the query
--rawsv         <s> : output the syntenic gaps information
--alrate        <i> : rate cutoff for the aligned sequences to the total range of the top chain or other filled chain that to be considered as syntenic region.
--toplen        <i> : the length cutoff for top chain that is to be considered as transolation.
--nonsynlen     <i> : the length cutoff for secondary or below chain that is to be considered as transolation.
--synlen        <i> : the length cutoff for secondary or below chain that is to be considered as transolation.
--invlen        <i> : the length cutoff for secondary or below chain that is to be considered as transolation.
--qFarCutoff    <i> : far cutoff of the breakpoint of a inversion far from the 
--help              : show this help information
END.
exit;
}

### set default values
$toplen ||=500000;
$nonsynlen ||=500000;
$synlen ||=500000;
$invlen ||=10000;
$alignrate ||=0.4;
$qFarCutoff ||=5000000;
$kent ||='';
$tba ||='';

### Main ###
open In, "$input" or die "$!";
open Out, ">$output" or die "$!";
open Out1, ">$one2mul" or die "$!";

my %hash;
my $tag=0;
my $tag1=0;
my $tag2=0;
my $chr;

while (<In>) {
chomp;
if ($_=~/^net|\#/) {
  print Out "$_\n";
     if ($_=~/^net/) {
         my @tmp = split (/\s/,$_);
         my @unit = split (/\./,$tmp[1]);
            $chr = $unit[1];
            print "Processing $chr\n";
      }
} elsif ($_=~/^\s{1}fill/) {
my @temp = split (/\s/, $_);
my @tmp1 = split (/\./,$temp[4]);
my $query_chr = $tmp1[1];
my $ref = \@temp;
my $ali_idx = &array_search ($ref,"ali");
my $topalrate = $temp[$ali_idx+1]/$temp[3];

 if ($temp[3]>$toplen and $topalrate > $alignrate) { # and ($chr eq $query_chr)) {
 print Out "$_\n";
 push (@{$hash{$chr}},$query_chr);
 $tag2=0;
 } else {
 $tag2=1;
 } 
} elsif ($_=~/^\s{2}gap/ and $tag2==0) {
 print Out "$_\n"; 
} elsif ($_=~/^\s{3}fill/ and $_=~/nonSyn/) {

  my @temp = split (/\s+/, $_);
  my @tmp1 = split (/\./,$temp[4]);
  my $query_chr = $tmp1[1];
  my $ref = \@temp;
  my $ali_idx = &array_search ($ref,"ali");
  my $secalrate = $temp[$ali_idx+1]/$temp[3];

   if ($tag2==0 and $temp[3]>$nonsynlen and $secalrate > $alignrate) { # and ($chr eq $query_chr)) {
   print Out "$_\n";
   push (@{$hash{$chr}},$query_chr);
   $tag=0;
   } else {
   $tag=1;
   }

} elsif ($_=~/^\s{3}fill/ and $_=~/syn|inv/) {

     my @temp = split (/\s+/, $_);
     my @tmp1 = split (/\./,$temp[4]);
     my $query_chr = $tmp1[1];
     my $ref = \@temp;
     my $ali_idx = &array_search ($ref,"ali");
     my $qFar_idx = &array_search ($ref,"qFar");
     my $qOver_idx = &array_search ($ref,"qOver");
     my $secalrate = $temp[$ali_idx+1]/$temp[3];   

     if ($_=~/syn/) {
         if ( $tag2==0 and $temp[3]>$synlen and $secalrate > $alignrate) { # or (4*$temp[3]>$temp[$qFar_idx+1]))) {
              print Out "$_\n"; 
              $tag=0;
         } else {
              $tag=1;
         } 
     } elsif ($_=~/inv/) {
       if ( $tag2 ==0 and $temp[3]>$invlen and $secalrate > $alignrate and $temp[$qOver_idx+1] > 0 and $temp[$qFar_idx+1]<$qFarCutoff) {
            print Out "$_\n"; 
            $tag=0;
       } elsif ( $tag2 ==0 && $temp[3]>5*$invlen && $temp[7] > 5*$invlen && $temp[$qFar_idx+1]<$temp[3]) { 
            print Out "$_\n";
            $tag=0;
       } else {
            $tag=1;
       }     
     }
        
} elsif ($_=~/^\s{4}gap/ and $tag==0 and $tag2==0) {
        print Out "$_\n";
} elsif ($_=~/^\s{5}fill/ and $_=~/nonSyn/) {

  my @temp = split (/\s+/, $_);
  my @tmp1 = split (/\./,$temp[4]);
  my $query_chr = $tmp1[1];
  my $ref = \@temp;
  my $ali_idx = &array_search ($ref,"ali");
  my $secalrate = $temp[$ali_idx+1]/$temp[3];
   if ($tag==0 and $tag2==0 and $temp[3]>$nonsynlen and $secalrate > $alignrate) { # and ($chr eq $query_chr)) {
   print Out "$_\n";
   push (@{$hash{$chr}},$query_chr);
   $tag1=0;
   } else {
   $tag1=1;
   }
        
} elsif ($_=~/^\s{5}fill/ and $_=~/syn|inv/) {

  my @temp = split (/\s+/, $_);
  my @tmp1 = split (/\./,$temp[4]);
  my $query_chr = $tmp1[1];
  my $ref = \@temp;
  my $ali_idx = &array_search ($ref,"ali");
  my $qFar_idx = &array_search ($ref,"qFar"); 
  my $qOver_idx = &array_search ($ref,"qOver");
  my $secalrate = $temp[$ali_idx+1]/$temp[3];   

  if ($_=~/syn/) {
      if ( $tag2==0 and $tag==0 and $temp[3]>$synlen and $secalrate > $alignrate) { #or (4*$temp[3]>$temp[$qFar_idx+1]))) {
           print Out "$_\n"; 
           $tag1=0;
      } else {
           $tag1=1;
      } 
  } elsif ($tag2==0 and $tag==0 and $_=~/inv/) {
    if ($temp[3]>$invlen and $secalrate > $alignrate and $temp[$qOver_idx+1] > 0 and $temp[$qFar_idx+1]<$qFarCutoff) {
         print Out "$_\n"; 
         $tag1=0;
    } elsif ($temp[3]>5*$invlen && $temp[7] > 5*$invlen && $temp[$qFar_idx+1]<$temp[3]) {
         print Out "$_\n";
         $tag1=0;
    } else {
         $tag1=1;
    }     
  }

} elsif ($_=~/^\s{6}gap/ and $tag1==0 and $tag2==0 and $tag==0) {
     print Out "$_\n";
} else {
  next;
  }
}


foreach my $chrom (sort{$hash{$a}<=>$hash{$b}} keys %hash)  {
  print Out1 "$chrom\t@{$hash{$chrom}}\n";
}

close In;
close Out;
close Out1;

### To sing.maf 
my $target_file = basename($target);
my $query_file = basename($query);
my $target_path = dirname($target);
my $query_path = dirname($query);

unless (-s "$target_path/$target_file.sizes") {
        `${kent}faSize -detailed $target > $target_path/$target_file.sizes`;
}
unless (-s "$target_path/$target_file.2bit") {
        `${kent}faToTwoBit $target $target_path/$target_file.2bit`;
}
unless (-s "$query_path/$query_file.sizes") {
        `${kent}faSize -detailed $query > $query_path/$query_file.sizes`;
}
unless (-s "$query_path/$query_file.2bit") {
        `${kent}faToTwoBit $query $query_path/$query_file.2bit`;
}


`${kent}netToAxt $output $chain $target_path/$target_file.2bit $query_path/$query_file.2bit $output.axt`;
`${kent}axtToMaf $output.axt $target_path/$target_file.sizes $query_path/$query_file.sizes $output.axt.maf`;
`${tba}single_cov2 $output.axt.maf > $target_file.$query_file.sing.maf`;
`rm $output.axt.maf`;

### quickly to synnic gaps

open Gap, "$output" or die "$!";
open SYN, ">$rawsv" or die "$!";
my $Chr;
while (<Gap>) {
chomp;
$_ =~ s/^\s+//;
my @splits = split (/\s+/,$_);
if ($_=~/^net/) {
$Chr = $splits[1];
} elsif ($_=~/^gap/) {
my $end = $splits[1]+$splits[2];
print SYN "$Chr\t$splits[1]\t$end\n";
 }
}

close Gap;
close SYN;

##########################
#######sub Function ######
##########################

sub array_search {
    my ($arr,$elem) = @_;
    my $idx;
    for my $i (0..$#$arr) {
        if ($arr->[$i] eq $elem) {
            $idx = $i;
            last;
         }
     }
  return $idx;            
}
