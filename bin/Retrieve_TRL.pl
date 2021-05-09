#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(&basename &dirname);
use List::Util qw/ min max /;
use Data::Dumper;

####################################################################################
# Retrieve_TRL.pl
#
# Authors: Yi Liao (06/19/2020)
#
# Copyright (C) Not for commercial use
#
# Prerequisite : Perl
#
# Retrieve potential translocation events.
#
####################################################################################

my ($input1,$input2,$output,$Help);

# get input parameters
GetOptions( 'in1=s'=>\$input1,
            'in2=s'=>\$input2,
            'output=s'=>\$output,
            'help'=>\$Help
           );
         
my $data = &read2Hash ($input2);

open TRL, "$input1" or die "$!";
open OUT, ">$output" or die "$!";

while (<TRL>) {
   chomp;
   next if ($_=~/\#/);
my @temp = split (/\t/,$_);
my @unit = split (/\./,$temp[3]);

foreach my $key (keys %{${$data}{$unit[1]}}) {
        my @num = split (/\_/,$key);
        my $rate = &intersection ($temp[4],$temp[5],$num[0],$num[1]);
        if ($rate > 0.8 ) {
        my @tmp = split ("_", ${${$data}{$unit[1]}}{$key});
           if ($tmp[0] eq $temp[0]) {
           my $rate0 = &intersection ($temp[1],$temp[2],$tmp[1],$tmp[2]);
              if ($rate0 > 0.8) {
                 print OUT "$_\n";
               }
            }
         last;       
        }
    }
}


#################################
sub read2Hash {         
my $input = shift;
my %hash;
open IN, "$input" or die "$!";
 while (<IN>) {
   chomp;
   next if ($_=~/\#/);
   my @temp = split (/\t/,$_);
   my @unit = split (/\./,$temp[0]);
   my $key = join ("_",($temp[1],$temp[2]));
   my $value = join ("_",($temp[3],$temp[4],$temp[5]));
      $hash{$unit[1]}{$key} = $value;  
 }
   my $refhash = \%hash;
   return $refhash;
}

############ intersection
sub intersection {
my ($n1,$n2,$n3,$n4)=@_;
my $overlap;
my $min = max($n1,$n3);
my $max = min($n2,$n4);
if ($min > $max) {
$overlap = 0;
} else {
$overlap = ($max-$min)/($n2-$n1);
}
return $overlap;
}

