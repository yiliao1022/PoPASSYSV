#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

####################################################################################
# GenometypingDeletions_Sing.pl
#
# Authors: Yi Liao (07/23/2020)
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

### read genoms/assemblies ID list from a input file (How many genomes are you processed)
open (GENOME,"$glst") or die "Genome/species names do not exist, $!";
my @Genome=<GENOME>;
close GENOME;
chomp @Genome;
my $reference=shift @Genome;
my $ref_glst = \@Genome;


open VCF, ">$OUTPUTdir/$svraw.vcf" or die "$!";
open IN, "$svraw" or die "$!";
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
###################### Main ###########################
#######################################################

### read genoms/assemblies ID list from a input file (How many genomes are you processed)
my $first=`head -1 $svraw`;
my $last=`tail -1 $svraw`;
my @beg0 = split (/\t/,$first);
my @end0 = split (/\t/,$last);
my $beg1 = $beg0[1]-$extend;
if ($beg1<0) {
$beg1=0;
}
my $end1 = $end0[2]+$extend;
my %liftover;

foreach my $g (@Genome) {
    `${tba}mafFind $singmaf/$reference.$g.sing.maf $beg1 $end1 > $OUTPUTdir/$svraw.$reference.$g.sing.maf`;
     my $liftover_coord = &liftOver("$OUTPUTdir/$svraw.$reference.$g.sing.maf");
     $liftover{$g} = $liftover_coord;
}

while (<IN>) {
chomp;
my ($chr,$beg,$end,$len,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,8,9,10,11,12];
my $ref = (split(/\./,$chr))[0];
my $alt = (split(/\./,$qchr))[0];
my $flank;

if ($len <= 100) {
  $flank = 100;
} elsif ($len > 100 and $len < 500 ) {
  $flank = 200;
} elsif ($len >= 500 && $len < 2000) {
  $flank = 500;
} elsif ($len >=2000 && $len < 5000) {
  $flank = 2000;
}  else {
  $flank = 5000;
}

my $beg0 = $beg - $fiv + 5;
my $end0 = $end + $thr - 5;

eval {
print VCF "$chr\t$beg\tSV$beg\tN\t\<$svtype\>\t22\tPASS\tSVTYPE=$svtype;END=$end;SVLEN=$len\tGT:GQ\t";
foreach my $ind (@Genome) {
   if (exists ${$liftover{$ind}}{$beg0} && exists ${$liftover{$ind}}{$end0}) {
     my @begflank1= split (/_/,${$liftover{$ind}}{$beg0});
     my @endflank1 = split (/_/,${$liftover{$ind}}{$end0});

     if ( abs ($begflank1[1]-$endflank1[1]) - ($fiv + $thr - 10) < 13 ) {
                  if ($ind eq $Genome[-1]) {
                        print VCF "1|1:22";
                  } else {
                        print VCF "1|1:22\t";
                  }
                  next;
     } elsif ( (abs ($begflank1[1]-$endflank1[1]) - ($fiv + $thr - 10)) > $len-13 && (abs ($begflank1[1]-$endflank1[1]) - ($fiv + $thr - 10)) < $len + 13) {
                 if ($ind eq $Genome[-1]) {
                      print VCF "0|0:22";
                 } else {
                      print VCF "0|0:22\t";
                 }
                    next;
     } elsif ((abs ($begflank1[1]-$endflank1[1]) - ($fiv + $thr - 10)) > 200000)  {
                if ($ind eq $Genome[-1]) {
                    print VCF "\.|\.:22";
                } else {
                    print VCF "\.|\.:22\t";
                }
                  next;

     } else {
       
       my $begin = $beg - $flank;
       my $endin = $end + $flank;

       `${tba}maf2fasta $refseq $OUTPUTdir/$svraw.$reference.$ind.sing.maf $begin $endin > $OUTPUTdir/TBA.$beg\_$end.fasta`;

             if (-z "$OUTPUTdir/TBA.$beg\_$end.fasta") {
                    system "rm $OUTPUTdir/TBA.$beg\_$end.fasta";
                    if ($ind eq $Genome[-1]) {
                      print VCF "\.|\.:22";
                    } else {
                      print VCF "\.|\.:22\t";
                    }
                    next;
             }

        my $evaluation = &ParseMaf2fasta("$OUTPUTdir/TBA.$beg\_$end.fasta",$flank,$ind,$len);

        if ($ind eq $Genome[-1]) {
        print VCF "$evaluation|$evaluation:22";
        } else {
        print VCF "$evaluation|$evaluation:22\t";
        }   
     }
   
   } else {
      my $begin = $beg - $flank;
      my $endin = $end + $flank;
     
      `${tba}maf2fasta $refseq $OUTPUTdir/$svraw.$reference.$ind.sing.maf $begin $endin > $OUTPUTdir/TBA.$beg\_$end.fasta`;

        if (-z "$OUTPUTdir/TBA.$beg\_$end.fasta") {
               system "rm $OUTPUTdir/TBA.$beg\_$end.fasta";
               if ($ind eq $Genome[-1]) {
                 print VCF "\.|\.:22";
               } else {
                 print VCF "\.|\.:22\t";
               }
               next;
        }

   my $evaluation = &ParseMaf2fasta("$OUTPUTdir/TBA.$beg\_$end.fasta",$flank,$ind,$len);

   if ($ind eq $Genome[-1]) {
   print VCF "$evaluation|$evaluation:22";
   } else {
   print VCF "$evaluation|$evaluation:22\t";
   }
} 
  `rm $OUTPUTdir/TBA.$beg\_$end.fasta`;
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
##################################################################################
sub ParseMaf2fasta10X {
my ($file,$extending,$query,$length,$seed) = @_;
my @data;
$/="\n";

open(DATA, "$file") or die "Couldn't open file sub maf2fasta file, $!";
while (<DATA>) {
chomp;
push (@data,$_);
}
close DATA;

my @seq_num = split (/\s/,$data[0]);
my @ref_name = split (/:/,$data[1]);
my $ref = $ref_name[0];
my @ref_seq = split(//,$data[3]);
my @query_seq = split(//,$data[4]);

my @ref_tens = unpack '(a$seed)*',$data[3];
my @query_tens = unpack '(a$seed)*',$data[4];

my @Gap_coord;
my $len_n=0;
my $t;
for ($t=0;$t<$seq_num[1];$t++) {
 if ($ref_seq[$t] =~/[ATGCNatgcn]/) {
  $len_n++;
  my $s;
  for ($s=$extending+1;$s<$extending+$length+1;$s++) {
  if ($len_n == $s) {
  push (@Gap_coord,$t);
   }
  }
 }
}
my $start = $Gap_coord[0];
my $end = $Gap_coord[-1];
my $start_int = int ($start/$seed);
my $end_int = int ($end/$seed);
my $total_int = int (($#ref_tens+1)/$seed);

my $j;
my $m=0;
my $n=0;
my $h=0;
for ($j=0;$j<$start_int;$j++) {
     if ($ref_tens[$j] =~/[ATGCNatgcn]/) {
         $n++;
        $ref_tens[$j] = uc ($ref_tens[$j]);
       if ($query_tens[$j] =~/[ATGCNatgcn]/) {
           $h++;
           $query_tens[$j] = uc ($query_tens[$j]);
              if ($ref_tens[$j] eq $query_tens[$j]) {
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
for ($q=$end_int;$q<=$total_int;$q++) {
        if ($ref_tens[$q] =~/[ATGCNatgcn]/) {
            $n1++;
           $ref_tens[$q] = uc ($ref_tens[$q]);
          if ($query_tens[$q] =~/[ATGCNatgcn]/) {
                 $h1++;
                 $query_tens[$q] = uc ($query_tens[$q]);
                 if ($ref_tens[$q] eq $query_tens[$q]) {
                    $m1++;
                 }
           }
        }
      }
my $gap_identity3 = $m1/($h1+1);
my $gap_overlap3 = $h1/($n1+1);

my $n2=0;
my $m2=0;
my $x;

for ($x=$start_int;$x<=$end_int;$x++) {
  $n2++;
  $ref_tens[$x] = uc ($ref_tens[$x]);
  $query_tens[$x] = uc ($query_tens[$x]);
  if ($ref_tens[$x]=~/[ATGCNatgcn]/ and ($ref_tens[$x] eq $query_tens[$x])) {
  $m2++;
}
}
my $gap_region = $m2/($n2+1);

my $evalue;
 #  print "IDENTITY: $gap_identity5\t$gap_identity3\t$gap_region\n";
   if ($gap_identity5>=0.9 && $gap_identity3 >=0.9 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region >= 0.5) {
       $evalue=0; ### similar to reference
    } elsif ($gap_identity5>=0.9 && $gap_identity3 >=0.9 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region < 0.1) {
       $evalue=1; ### support a deletion
    } else {
       $evalue="\."; ### missing data
    }
return $evalue;
}
##################################################################################
sub Stretcher {
my ($file,$target,$query)=@_;
my @lines;
local $/="\n";
open(IN1, "$file") or die "Couldn't open file sub maf2fasta file, $!";
while (<IN1>) {
chomp;
push (@lines,$_);
}
close IN1;
$lines[3]=~s/-//;
$lines[4]=~s/-//;

open SEQ1, ">$OUTPUTdir/$target.fa" or die "$!";
print SEQ1 ">$target\n$lines[3]\n";
open SEQ2, ">$OUTPUTdir/$query.fa" or die "$!";
print SEQ2 ">$query\n$lines[4]\n";
`${stretcher}stretcher -asequence $OUTPUTdir/$target.fa -bsequence $OUTPUTdir/$query.fa -outfile $OUTPUTdir/$target.$query.stretcher`;
my $e_evalue = &Parse_Stretcher("$OUTPUTdir/$target.$query.stretcher");
`rm $OUTPUTdir/$target.fa $OUTPUTdir/$query.fa $OUTPUTdir/$target.$query.stretcher`;
return $e_evalue;
}
##################################################################################
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
##################################################################################
sub ParseMaf2fasta {
my ($file,$extending,$query,$length) = @_;
my @data;
$/="\n";

open(DATA, "$file") or die "Couldn't open file sub maf2fasta file, $!";
while (<DATA>) {
chomp;
push (@data,$_);
}
close DATA;

my @seq_num = split (/\s/,$data[0]);
my @ref_name = split (/:/,$data[1]);
my $ref = $ref_name[0];
my @ref_seq = split(//,$data[3]);
my @query_seq = split(//,$data[4]);

my $len_n=0;
my $t;
my $start;
my $end;

for ($t=0;$t<$seq_num[1];$t++) {
 if ($ref_seq[$t] =~/[ATGCNatgcn]/) {
  $len_n++;
 if ($len_n==$extending+1) {
     $start = $t;
 }
 if ($len_n==$extending+$length) {
    $end = $t;
   }
  }
}

my $j;
my $m=0;
my $n=0;
my $h=0;
for ($j=0;$j<$start;$j++) {
     if ($ref_seq[$j] =~/[ATGCNatgcn]/) {
         $n++;
        $ref_seq[$j] = uc ($ref_seq[$j]);
       if ($query_seq[$j] =~/[ATGCNatgcn]/) {
           $h++;
           $query_seq[$j] = uc ($query_seq[$j]);
              if ($ref_seq[$j] eq $query_seq[$j]) {
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
for ($q=$end+1;$q<=$#ref_seq;$q++) {
        if ($ref_seq[$q] =~/[ATGCNatgcn]/) {
            $n1++;
           $ref_seq[$q] = uc ($ref_seq[$q]);
          if ($query_seq[$q] =~/[ATGCNatgcn]/) {
                 $h1++;
                 $query_seq[$q] = uc ($query_seq[$q]);
                 if ($ref_seq[$q] eq $query_seq[$q]) {
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
  if ($ref_seq[$y] =~/[ATGCNatgcn]/) {
    $n2++;
    $ref_seq[$y] = uc ($ref_seq[$y]);
    if ($query_seq[$y] =~/[ATGCNatgcn]/) {
       $h2++;
       $query_seq[$y] = uc ($query_seq[$y]);
       if ($ref_seq[$y] eq $query_seq[$y]) {
                    $m2++;
        }
     }
   }
 }

my $gap_region = $m2/($n2+1);

my $evalue;
 #  print "IDENTITY: $gap_identity5\t$gap_identity3\t$gap_region\n";
   if ($gap_identity5>=0.9 && $gap_identity3 >=0.9 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region >= 0.5) {
       $evalue=0; ### similar to reference
    } elsif ($gap_identity5>=0.9 && $gap_identity3 >=0.9 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region < 0.1) {
       $evalue=1; ### support a deletion
    } else {
       $evalue="\."; ### missing data
    }
return $evalue;
}
##################################################################################
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
##################################################################################
