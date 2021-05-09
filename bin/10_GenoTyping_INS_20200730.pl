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
# Authors: Yi Liao (07/01/2020)
#
# Copyright (C) Not for commercial use
#
# Further computationally validated the merge raw SVs
#
# Prerequisite : TBA tools; MUSCLE
#
####################################################################################

my $time= strftime "%F", localtime;
my ($svraw,$svtype,$glst,$singmaf,$genomes,$stretcher,$tba,$refsize,$refname,$speciesname,$refseq,$OUTPUTdir,$extend,$head,$help);

GetOptions(
    "svraw=s"=>\$svraw,
    "svtype=s"=>\$svtype,
    "glst=s"=>\$glst,
    "singmaf=s"=>\$singmaf,
    "gseqs=s"=>\$genomes,
    "stretcher=s"=>\$stretcher,
    "tba=s"=>\$tba,
    "refsize=s"=>\$refsize,
    "refname=s"=>\$refname,
    "speciesname=s"=>\$speciesname,
    'refseq=s'=>\$refseq,
    "extend=s"=>\$extend,
    'outdir=s'=>\$OUTPUTdir,
    'head'=>\$head,
    'help'=>\$help,
);

############# Help information ##########
if ($help){
print <<"END.";
  Usage: perl $0 -svraw SV.txt -svtype INS -glst genome.lst -refsize B73.sizes -singmaf /path/to/sing/ -gseqs /path/tp/genomes -stretcher /path/to/stretcher -tba /path/to/tba -refname B73 -speciesname maize -outdir B73_tmp
  -svraw        raw merged SV file                                      [REQUIRED]
  -svtype       SV types [INS|DEL|CNV|INV]                              [REQUIRED]
  -gseqs        path to all genome sequences                            [REQUIRED]
  -singmaf      path to single coverage maf file                        [REQUIRED]
  -glst         a file with the list of genomes (one line one genome)   [REQUIRED]
  -refsize      reference sizes                                         [REQUIRED]
  -refname      reference name, e.g. hg38, IRGSP1                       [REQUIRED]
  -speciesname  species name, e.g. human, Drosophila, rice , maize      [REQUIRED]
  -tba          path for tba program                                    [REQUIRED]
  -outdir       define a output dir                                     [REQUIRED]
  -stretcher    path to stretcher                                       [REQUIRED]
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
############ Print VCF head ###########################
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
    `${tba}mafFind $singmaf/$reference.$g.sing.maf $beg $end > $OUTPUTdir/$svraw.$reference.$g.sing.maf`;
     my $liftover_coord = &liftOver("$OUTPUTdir/$svraw.$reference.$g.sing.maf");
     $liftover{$g} = $liftover_coord;
}

open IN, "$svraw" or die "$!";
while (<IN>) {
chomp;
my ($chr,$beg,$end,$len,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,8,9,10,11,12];
my $ref = (split(/\./,$chr))[0];
my $alt = (split(/\./,$qchr))[0];
my $flank;

if ($len < 1000) {
$flank=1000;
} elsif ($len >=1000 and $len < 5000) {
$flank = 5000;
} else {
$flank = 8000;
}

my $beg0;
my $end0;
my $lenextend;

if ($fiv >= $flank && $thr >= $flank) {
$beg0 = $beg - $flank;
$end0 = $end + $flank;
$lenextend = 2*$flank;
} elsif ($fiv >=$flank && $thr < $flank) {
$beg0 = $beg - $flank;
$end0 = $end + $thr - 5;
$lenextend = $flank + $thr - 5;
} elsif ($fiv < $flank && $thr >= $flank) {
$beg0 = $beg - $fiv + 5;
$end0 = $end + $flank;
$lenextend = $flank + $fiv -5;
} else {
$beg0 = $beg - $fiv + 5;
$end0 = $end + $thr - 5;
$lenextend = $fiv+$thr-10;
}

eval {
my $altseq;
my ($gap_s,$gap_e);

if (exists ${$liftover{$alt}}{$beg0} && exists ${$liftover{$alt}}{$end0}) {
  my @begflank= split (/_/,${$liftover{$alt}}{$beg0});
  my @endflank = split (/_/,${$liftover{$alt}}{$end0});
 if ($begflank[1]<$endflank[1] && ($begflank[0] eq $qchr) && ($endflank[0] eq $qchr)) {
    $altseq = &SubSeq("$genomes/$alt",$qchr,$begflank[1],$endflank[1],"plus");
    $gap_s = $qs-$begflank[1];
    $gap_e = $qe-$begflank[1];
  } elsif ($begflank[1]>$endflank[1] && ($begflank[0] eq $qchr) && ($endflank[0] eq $qchr)) {
    $altseq = &SubSeq("$genomes/$alt",$qchr,$endflank[1],$begflank[1],"minus");
    $gap_s = $qs-$endflank[1];
    $gap_e = $qe-$endflank[1];
  } else {
    print ERR "Reason 1: $_\n";
    next;
  }
} else {

    my $k;
    for ($k=$beg0;$k<=$beg;$k++) {
    if (exists ${$liftover{$alt}}{$k}) {
    $beg0 = $k;
    last;
    }
    }

    my $l;
    for ($l=$end0;$l>=$end;$l--) {
    if (exists ${$liftover{$alt}}{$l}) {
    $end0 = $l;
    last;
    }
    }

    if (exists ${$liftover{$alt}}{$beg0} && exists ${$liftover{$alt}}{$end0}) {
           my @begflank= split (/_/,${$liftover{$alt}}{$beg0});
           my @endflank = split (/_/,${$liftover{$alt}}{$end0});
        if ($begflank[1]<$endflank[1] && ($begflank[0] eq $qchr) && ($endflank[0] eq $qchr)) {
              $altseq = &SubSeq("$genomes/$alt",$qchr,$begflank[1],$endflank[1],"plus");
              $gap_s = $qs-$begflank[1];
              $gap_e = $qe-$begflank[1];
         } elsif ($begflank[1]>$endflank[1] && ($begflank[0] eq $qchr) && ($endflank[0] eq $qchr)) {
              $altseq = &SubSeq("$genomes/$alt",$qchr,$endflank[1],$begflank[1],"minus");
              $gap_s = $qs-$endflank[1];
              $gap_e = $qe-$endflank[1];
         } else {
              print ERR "Reason 2: $_\n";
              next;
         }
    } else {
              print ERR "Reason 3:$_\n";
              next;
    }
}


my $SVend = $beg + $len;
print VCF "$chr\t$beg\tSV$beg\tN\t\<$svtype\>\t22\tPASS\tSVTYPE=INS;END=$SVend;SVLEN=\-$len\tGT:GQ\t";

foreach my $ind (@Genome) {
     if ($ind eq $alt) {
           if ($ind eq $Genome[-1]) {
             print VCF "1|1:22";
           } else {
             print VCF "1|1:22\t";
           }
             next;
      }


  my $k;

  if (!exists ${$liftover{$ind}}{$beg0}) {
       for ($k=$beg0;$k<=$beg;$k++) {
        if (exists ${$liftover{$ind}}{$k}) {
           $beg0 =$k;
           last;
         }
       }
  }

  if (!exists ${$liftover{$ind}}{$end0}) {
       for ($k=$end0;$k>=$end;$k--) {
        if (exists ${$liftover{$ind}}{$k}) {
           $end0 =$k;
           last;
         }
       }
  }


  if (exists ${$liftover{$ind}}{$beg0} && exists ${$liftover{$ind}}{$end0}) {
    my @begflank1= split (/_/,${$liftover{$ind}}{$beg0});
    my @endflank1 = split (/_/,${$liftover{$ind}}{$end0});

    if ( abs ($begflank1[1]-$endflank1[1]) - $lenextend < 13 ) {
     if ($ind eq $Genome[-1]) {
      print VCF "0|0:22";
      } else {
      print VCF "0|0:22\t";
      }
      next;
    } elsif (abs ($begflank1[1]-$endflank1[1]) - $lenextend > 200000) {
     if ($ind eq $Genome[-1]) {
       print VCF "\.|\.:22";
       } else {
       print VCF "\.|\.:22\t";
       }
       next;
    }

    my $indseq;
    if ($begflank1[1]<$endflank1[1] && ($begflank1[0] eq $endflank1[0] )) {
    $indseq = &SubSeq("$genomes/$ind",$begflank1[0],$begflank1[1],$endflank1[1],"plus");
    } elsif ($begflank1[1]>$endflank1[1] && ($begflank1[0] eq $endflank1[0])) {
    $indseq = &SubSeq("$genomes/$ind",$begflank1[0],$endflank1[1],$begflank1[1],"minus");
    } else {
      
       if ($ind eq $Genome[-1]) {
         print VCF "\.|\.:22";
       } else {
         print VCF "\.|\.:22\t";
       }
     next;
    }

  my ($seq1,$seq2,$e_value) = &Realign($alt,$ind,$altseq,$indseq);

      if (defined ($seq1) && defined ($seq2) && defined ($e_value) ) {
  
                    if ($e_value>95) {
                      if ($ind eq $Genome[-1]) {
                             print VCF "1|1:22";
                      } else {
                             print VCF "1|1:22\t";
                      }
                     } else {
                               my $evaluation = &ParseMaf($seq1,$seq2,$gap_s,$gap_e);

                               if (defined ($evaluation) ) {
                                             if ($ind eq $Genome[-1]) {
                                                      print VCF "$evaluation|$evaluation:22";
                                             } else {
                                                      print VCF "$evaluation|$evaluation:22\t";
                                             }
                             } else {
                                            if ($ind eq $Genome[-1]) {
                                                      print VCF "\.|\.:22";
                                            } else {
                                                      print VCF "\.|\.:22\t";
                                            }
                                                 next;                           
                             }
                          
                   }
      } else {

        if ($ind eq $Genome[-1]) {
         print VCF "\.|\.:22";
         } else {
         print VCF "\.|\.:22\t";
         }
         next;
         
    }
  
} else {
 
  if ($ind eq $Genome[-1]) {
   print VCF "\.|\.:22";
   } else {
   print VCF "\.|\.:22\t";
   }
   next;
}
}
print VCF "\n";
};

if ($@) {
  print ERR "Failed\t$_\n";
}
}

`rm $OUTPUTdir/*.maf`;

#################
##### Module ####
#################
sub ParseMaf {
my ($sequenceA,$sequenceB,$svbeg,$svend) = @_;
my @seqA = split(//,$sequenceA);
my @seqB = split(//,$sequenceB);

my $i;
my $len_n=0;
my @gap_coord;
my $start;
my $end;
for ($i=0;$i<=$#seqA;$i++) {
 if ($seqA[$i]=~/[ATGCNatgcn]/) {
  $len_n++;
      if ($len_n==$svbeg) {
          $start = $i;
      }
      if ($len_n==$svend) {
          $end = $i;
      }
   }
}

my $j;
my $m=0;
my $n=0;
my $h=0;
for ($j=0;$j<$start;$j++) {
     if ($seqA[$j] =~/[ATGCNatgcn]/) {
         $n++;
        $seqA[$j] = uc ($seqA[$j]);
       if ($seqB[$j] =~/[ATGCNatgcn]/) {
           $h++;
           $seqB[$j] = uc ($seqB[$j]);
              if ($seqA[$j] eq $seqB[$j]) {
                 $m++;
              }
       }
    }
}

my $gap_identity5 = $m/($h+1);
my $gap_overlap5 = $h/($n+1);

my $q;
my $n1=0;
my $m1=0;
my $h1=0;
for ($q=$end+1;$q<=$#seqA;$q++) {
     if ($seqA[$q] =~/[ATGCNatgcn]/) {
         $n1++;
        $seqA[$q] = uc ($seqA[$q]);
       if ($seqB[$q] =~/[ATGCNatgcn]/) {
              $h1++;
              $seqB[$q] = uc ($seqB[$q]);
              if ($seqA[$q] eq $seqB[$q]) {
                 $m1++;
              }
           }
       }
   }
my $gap_identity3 = $m1/($h1+1);
my $gap_overlap3 = $h1/($n1+1);

my $n2=0;
my $m2=0;
my $y;
my $h2;
 for ($y=$start;$y<=$end;$y++) {
  if ($seqA[$y] =~/[ATGCNatgcn]/) {
    $n2++;
    $seqA[$y] = uc ($seqA[$y]);
    if ($seqB[$y] =~/[ATGCNatgcn]/) {
       $h2++;
       $seqB[$y] = uc ($seqB[$y]);
       if ($seqA[$y] eq $seqB[$y]) {
                    $m2++;
        }
     }
   }
 }
my $gap_region = $m2/($n2+1);

my $evalue;
   if ($gap_identity5>=0.75 && $gap_identity3 >=0.75 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region >= 0.4) {
       $evalue=1; ### similar to reference
    } elsif ($gap_identity5>=0.75 && $gap_identity3 >=0.75 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region < 0.1) {
       $evalue=0; ### support a deletion
    } elsif ($gap_region >= 0.9) {
      $evalue=1;
    } else {
      $evalue="\."; ### missing data
    }
return $evalue;
}

#############################################

sub stretcher2maf {
my ($aname,$aseq,$alen,$bname,$bseq,$blen) = @_;
open OUT, ">$OUTPUTdir/$aname.$alen.$bname.$blen.maf" or die "$!";
print OUT "a score=1000000.0\n";
print OUT "s $aname 0 $alen + $alen $aseq\n";
print OUT "s $bname 0 $blen + $blen $bseq\n";
}

##############################################

sub Realign {
my ($t,$q,$seq1,$seq2)=@_;
open SEQ1, ">$OUTPUTdir/$t.fa" or die "$!";
print SEQ1 ">$t\n$seq1\n";
open SEQ2, ">$OUTPUTdir/$q.fa" or die "$!";
print SEQ2 ">$q\n$seq2\n";
`${stretcher}stretcher -asequence $OUTPUTdir/$t.fa -bsequence $OUTPUTdir/$q.fa -outfile $OUTPUTdir/$t.$q.stretcher`;
my ($s1, $s2) = &ParseStretcherOutPut("$OUTPUTdir/$t.$q.stretcher",$t,$q);
my $eValue = &Parse_Stretcher("$OUTPUTdir/$t.$q.stretcher");
`rm $OUTPUTdir/$t.fa $OUTPUTdir/$q.fa $OUTPUTdir/$t.$q.stretcher`;
return ($s1,$s2,$eValue);
}
##############################################

sub ParseStretcherOutPut {
my ($stretcher,$t,$q) = @_;
local $/="\n";
open EMBOSS, "$stretcher" or die "$!";
my @ref;
my @query;
while (<EMBOSS>) {
chomp;
$_ =~ s/^\s+//;
if ($_=~/\#/) {
}elsif ($_=~/^$t/) {
my @tmp1 = split (/\s/,$_);
push (@ref,$tmp1[1]);
} elsif ($_=~/^$q/) {
my @tmp2 = split (/\s/,$_);
push (@query,$tmp2[1]);
 }
 }
 my $seq1 = join ("",@ref);
 my $seq2 = join ("",@query);
 return ($seq1,$seq2);
 }
#################################################
sub Parse_Stretcher {
my $stretcher = shift;
my $value;
local $/="\n";
open IN2, "$stretcher" or die "$!";
while (<IN2>) {
  if ($_=~/\# Identity:/) {
    $_=~/\# Identity:   (.*) \((.*)%\)/;
    $value = $2;
    last;
  }
}
return $value;
}
################################################## Extract sequence based on coordinate
sub SubSeq {
my ($seq,$chr,$beg,$end,$strand) = @_;
local $/="\n>";
open FA, "$seq" or die "$!";
my $seqrange;
   while (my $sequence = <FA>) {
         chomp $sequence;
         my ($id) = $sequence =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
   if ($id eq $chr) {
         $sequence =~ s/^>*.+\n//;  # remove FASTA header
         $sequence =~ s/\n//g;  # remove endlines
         my $len = $end - $beg + 1;
         $seqrange = substr $sequence, $beg -1, $len;
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
################################################# Reverse and complement the sequence
sub RevCom {
my $seq = shift;
my $revcomp = reverse $seq;
$revcomp =~ tr/ATGCNatgcn/TACGNtacgn/;
return $revcomp;
}
#################################################
sub Maf2fasta {
my($maf,$ref,$qs,$qe,$extend)=@_;
open MAF, "$maf" or die "Missing maf file!";
$extend ||=0;
my %coord_proxy;
my $beg_index=0;
my $end_index=0;
my $s=1;
my $n=0;
my $ref_end=0;
my $gap_ref = &GetRefCoord("$maf",$ref);
my @species;
$ref = (split(/\./,$ref))[0];
local $/="\n\n"; # record separator as a blank line when reading
my $alt_strand;
my $alt_size;
while (<MAF>) {
chomp;
next if ($_=~/eof/);
my  @temp = split(/\n/,$_);
$s = $s + $ref_end + ${$gap_ref}[$n];
foreach my $line (@temp) {
  if ($line =~/^s/) {
    my ($spec,$start,$extend1,$strand,$qsize,$seq) = (split(/\s+/,$line))[1,2,3,4,5,6];
    push (@species,$spec);
    $spec = (split(/\./,$spec))[0];
    if ($spec !~ /$ref/) {
    $alt_strand = $strand;
    $alt_size = $qsize;
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

  my @Arr = sort {$a<=>$b} (keys %{$coord_proxy{$ref}});

  if ($alt_strand =~/\+/) {
  my $begin = $qs-1 - $extend;
  my $end2 = $qe + $extend;
  my @all_matches = grep {${$coord_proxy{$ref}}{$_} eq $begin} keys %{$coord_proxy{$ref}};
  if ($#all_matches>-1) {
   $beg_index = $all_matches[0];
  } else {
   $beg_index = $Arr[0];
  }

  my @all_matches1 = grep {${$coord_proxy{$ref}}{$_} eq $end2} keys %{$coord_proxy{$ref}};
  if ($#all_matches1>-1){
    $end_index = $all_matches1[0];
  } else {
    $end_index = $Arr[-1];
  }
} elsif ($alt_strand =~/\-/) {
  my $begin = $qs-$extend;
  my $end2 = $qe+1+$extend;


my @all_matches = grep {${$coord_proxy{$ref}}{$_} eq $begin} keys %{$coord_proxy{$ref}};
 if ($#all_matches>-1) {
 $beg_index = $all_matches[0];
} else {
  $beg_index = $Arr[0];
}

my @all_matches1 = grep {${$coord_proxy{$ref}}{$_} eq $end2} keys %{$coord_proxy{$ref}};
if ($#all_matches1>-1) {
$end_index = $all_matches1[0];
} else {
  $end_index = $Arr[-1];
  }
}
  if ($alt_strand =~/\+/) {
  $alt_strand = "plus";
  } elsif ($alt_strand=~/\-/) {
  $alt_strand = "minus";
  }
  my $ref_hash = \%coord_proxy;
  my $ref_spec = \@species;
  return ($ref_hash,$beg_index,$end_index,$alt_strand,$ref_spec,$alt_size);
}

###############################################################
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

###############################################################
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
