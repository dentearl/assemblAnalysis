#!/bin/env perl
# A hack created by Hiram and 
# clumsily adapted by Dent.
# 11 Feb 2011
#
# For use with assemblathon analysis pipeline,
# this script takes a chrom.sizes (which is really contig
# sizes) and a .list file containing either lines of .2bit specs 
# or lines of .lst files and produces a .lift file.
#
##############################
use warnings;
use strict;
use File::Basename;

my $argc = scalar(@ARGV);

if ($argc != 2) {
    printf STDERR "usage: ./constructLiftFile.pl chrom.sizes target.list > target.lift\n";
    exit 255;
}

my %chromSize;
my $file = shift;
open (FH, "<$file") or die "can not read $file";
while (my $line = <FH>) {
    chomp $line;
    my ($chr, $size) = split('\s+', $line);
    $chromSize{$chr} = $size;
}
close (FH);

$file = shift;
open (FH, "<$file") or die "can not read $file";
# print a bogus line to prevent empty lift file freak-out by liftOver in the query stage.
print "-28659142	HSCHR6_MHC_SSTO	171115067	chr6_ssto_hap7	4928567\n";
while (my $line = <FH>) {
    chomp $line;
    next if($line =~ m/.lst/);
    my $twoBitSpec = basename($line);
    my ($twoBitFile, $chrom, $startEnd) = split(':', $twoBitSpec);
    my ($start, $end) = split('-', $startEnd);
    my $length = $end - $start;
    die "can not find chrom size for $chrom" if (! exists($chromSize{$chrom}));
    my $size = $chromSize{$chrom};
    # HEY!
    # This is a very specific regex, it corresponds to running partitionSequence.pl 10000000 10000
    # SO if we change the way we run partitionSequence.pl, we need to change the second regex here.
    if (($twoBitSpec =~ m/:0-/) && ($twoBitSpec !~ m/:0-.00.0000/)) {
        printf "%d\t%s\t%d\t%s\t%d\n", $start, $chrom, $length, $chrom, $size;
    } else {
        printf "%d\t%s:%d-%d\t%d\t%s\t%d\n", $start, $chrom, $start, $end, $length, $chrom, $size;
    }
}
close (FH);

__END__
# this is an example of the sort of line we expect:
# /hive/data/genomes/ce9/bed/testLastz/ce9.2bit:chrV:0-10010000
