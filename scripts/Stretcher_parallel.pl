#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

####################################################################################
# Stretcher_parallel.pl
#
# Authors: Yi Liao (07/10/2020)
#
# Copyright (C) Not for commercial use
#
# run realign the axt file parallelly
#
# Prerequisite : Stretcher
#
# Usage: perl $0 -axt sample.axt -stretcher /path/to/stretcher -wk ./ -split 10
# -axt           raw axt file                                      [REQUIRED]
# -stretcher     path to stretcher program                         [REQUIRED]
# -wk            working folder                                    [REQUIRED]
# -split         how many subfiles are splited to                  [REQUIRED]
# -help          print this information
#
####################################################################################
my ($axt,$stretcher,$kent,$wk,$split,$sh,$popassysv,$Help);
GetOptions(
  "axt=s"=>\$axt,
  "stretcher=s"=>\$stretcher,
  "wk=s"=>\$wk,
  "split=i"=>\$split,
  "sh=s"=>\$sh,
  "kent=s"=>\$kent,
  "popassysv=s"=>\$popassysv,
  "help"     =>\$Help
);

$stretcher ||= '';
$wk ||= '.';
$split ||= 10000;
##################################################################################################
if ($Help){
print <<"END.";
  Usage: perl $0 -axt sample.axt -stretcher /path/to/stretcher -wk ./ -split 10
  -axt           raw axt file                                      [REQUIRED]
  -stretcher     path to stretcher program                         [REQUIRED]
  -wk            working folder                                    [REQUIRED]
  -split         how many subfiles are splited to                  [REQUIRED]
  -sh            bash output file
  -help       print this information

END.
exit;
}
################

`cat $axt | grep -v "#" > $axt.out`;
#`rm $axt`;
my $folder = "SUBFORDER";
`${kent}splitFile $axt.out $split Sub`;
if (-d "$wk/$folder") {
`rm -r $wk/$folder`;
`mkdir $wk/$folder`;
} else {
`mkdir $wk/$folder`;
}
`cp Sub* $wk/$folder`;
`rm Sub*`;
open SH,">$sh" or die "Error in writint to $sh";
opendir (DIR, "$wk/$folder") || die "Error in opening dir $folder";
while (my $filename = readdir(DIR)) {
  if ($filename =~/Sub/) {
    print SH "${popassysv}4_PairRealign.pl -axt $filename -output $filename.realigned.axt -stretcher ${stretcher}\n";
  }
}

