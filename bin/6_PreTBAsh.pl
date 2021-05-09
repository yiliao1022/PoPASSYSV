#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(&basename &dirname);

####################################################################################
# RunTBA.pl
#
# Authors: Yi Liao (06/19/2020)
# 
# Copyright (C) Not for commercial use
#
# Multiple genome alignments using TBA program  
# 
# 
# Prerequisite : Kent's utilities; TBA
#
# Usage:  perl RunTBA.pl -wkfolder /path/to/sing.maf_folder -gfolder /path/to/genomes_folder -kent /path/to/kentutil -tba /path/to/TBA -ref ref.fa
#
# Options: -wkfolder  [Path] Folder where the *.sing.maf files put                    [REQUIRED]
#          -gfolder   [Path] Folder where genomes put                                 [REQUIRED]
#          -kent      [Path] Kent' untilities                                         [OPTIONAL]
#          -tba       [Path] TBA program                                              [OPTIONAL]
#          -ref       [File] reference genome sequence                                [REQUIRED]
#          -help             print this information
#
####################################################################################
my ($wkfolder,$gfolder,$kent,$tba,$links,$ref,$Help);

GetOptions( 'wk=s' => \$wkfolder,
            'gfolder=s' => \$gfolder,
            'kent=s' => \$kent,
            'tba=s' => \$tba,
            'ref=s' => \$ref,
            'links=s' =>\$links,
            'help' => \$Help
            );

if ($Help){
print <<"END.";
Usage:  perl RunTBA.pl -wkfolder /path/to/sing.maf_folder -gfolder /path/to/genomes_folder -kent /path/to/kentutil -tba /path/to/TBA -ref ref.fa

Options: -wkfolder  [Path] Folder where the *.sing.maf files put                         [REQUIRED]
         -gfolder   [Path] Folder where genomes put                                      [REQUIRED]
         -kent      [Path] Kent' untilities                                              [OPTIONAL]
         -tba       [Path] TBA program                                                   [OPTIONAL] 
         -links     [File] Folder where the information for query chromosomes aligned to reference genome put [REQUIRED]
         -ref       [File] reference genome sequence                                     [REQUIRED]
         -help             print this information
END.
exit;
}

my $target_file = basename($ref);
my $target_path = dirname($ref);

unless (-e "$target_path/$target_file.sizes") {
        `${kent}faSize -detailed $ref > $target_path/$target_file.sizes`;
}

my @chr;
open SIZES, "$target_path/$target_file.sizes" or die "$!";
while (<SIZES>) {
chomp;
my @tmp = split (/\t/,$_);
push (@chr, $tmp[0]);
}
close SIZES;

foreach my $chromosome (@chr) {
unless (-e "$wkfolder/$chromosome") {
`mkdir $wkfolder/$chromosome`;
 }
&ExtractSeqByID ("$chromosome","$ref");
`cp $chromosome.tmp $wkfolder/$chromosome/$target_file`;
`rm $chromosome.tmp`;
}


`touch $wkfolder/splits.bed`;



opendir (DIR, $wkfolder ) || die "Error in opening dir $wkfolder\n";
my @pop;
while( (my $filename = readdir(DIR))) {
   if ($filename =~/(.*)\.(.*).sing.maf/) {
       print "Processing $1\t$2\t$filename\n";
 
my %links;
open LINKS, "$links/$1.$2.links" or die "$!";
while (<LINKS>) {
chomp;
my @tmp = split (/\t/,$_);
$links{$tmp[0]} = $tmp[1];
}
close LINKS;

      push (@pop, $2);
       `${kent}mafSplit splits.bed -byTarget -useFullSequenceName $1.$2.\. $wkfolder/$filename`;  
        foreach my $chrom (@chr) {  
          my @unit = split (/\./,$chrom);
          my @unit1= split (/\s/,$links{$unit[1]});
          foreach my $q_chr (@unit1) {
          &ExtractSeqByID ("$2.$q_chr","$gfolder/$2"); 
          `cat $2.$q_chr.tmp >> $2.tmp`;
          `rm $2.$q_chr.tmp`;      
          }

          system "cp $1.$2.$unit[1].maf $wkfolder/$chrom/$1.$2.sing.maf";
          system "rm $1.$2.$unit[1].maf";
          
              unless (-s "$wkfolder/$chrom/$2") {
                      system "cp $2.tmp $wkfolder/$chrom/$2";
                      system "rm $2.tmp";
              }
         }





        
    }  
}
closedir(DIR);

foreach my $ele (@chr) {
`echo '${tba}tba "($target_file @pop)" *.*.sing.maf $ele.tba.maf' > $ele.tba.sh`;
`cp $ele.tba.sh $ele/$ele.tba.sh`;
`rm $ele.tba.sh`;
}

########################################
############Sub function################
########################################

sub ExtractSeqByID {
my ($lookup,$genome) = @_; 
local $/ = "\n>";
open In, "$genome" or die "Couldn't open the genome file!";
open Out, ">$lookup.tmp" or die "$!";
while (my $seq = <In>) {
    chomp $seq;
    my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header 
    if ($id eq $lookup) {
         $seq=~ s/^>//;
    #    $seq =~ s/^>*.+\n//;  # remove FASTA header
    #   $seq =~ s/\n//g;  # remove endlines
         print Out ">$seq\n";    
         last;
       }
    }
}
