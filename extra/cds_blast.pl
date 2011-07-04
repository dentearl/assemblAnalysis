#!/usr/bin/perl
# Script to calculate how many CDSs are present in an assembly
#
# Written by Keith Bradnam
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author: keith $
# Last updated on: $Date: 2011/05/03 18:40:29 $

use strict; use warnings;
use FAlite;
use Getopt::Long;

###############################################
# 
#  C o m m a n d   l i n e   o p t i o n s
#
###############################################

my $csv;   # produce CSV output file of results
my $verbose; 
GetOptions ("csv"     => \$csv,
			"verbose" => \$verbose);

# check we have a suitable input file
my $usage = "Usage: $0 [options] <cds file> <gzipped assembly file>
options:
        -csv         produce a CSV output file of all results
        -verbose     turn on extra output for each CDS match
";

die "$usage" unless (@ARGV == 2);

my ($cds, $ASSEMBLY) = @ARGV;
my ($assembly_ID) = $ASSEMBLY =~ m/(\w\d+)_scaffolds/;


# Read sequences from CDS file and store in hash
my %cds_sequences;

open(my $fh, $cds)  or die "Can't open $cds\n";
my $fasta = new FAlite($fh);
while (my $entry = $fasta->nextEntry) {
	my ($id) = $entry->def =~ /^>(\d+)\s+/;
	$cds_sequences{$id} = $entry->seq;
}

##################
# BLAST
##################

# make suitable output file name for blast & csv results
my $output = "$assembly_ID.blast.out";

# format BLAST databases if not already done
unless (-s "$ASSEMBLY.xni")  {system("xdformat -n -I $ASSEMBLY")  == 0 or die "Can't create blast database for $ASSEMBLY\n"}

#  Run BLAST, but only if no output file already exists
my $params = "-S 50 -M 1 -N -1 -Q 3 -R 3 -W 15 -mformat 2 -kap -errors";
unless (-s $output) {system("blastn $ASSEMBLY $cds $params -o $output") && die "Can't run blastn\n"}


# want matches at 95% identity
my $identity = 95;

# main loop over BLAST output
open(my $blast, "<", "$output") or die "Can't open $output\n";

while (<$blast>) {
	my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, $pct, $ppos,
		$qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se, $gr) = split;

	# skip HSP matches < 95% identity
	next unless $pct >= $identity;

	# need to get start/end coordinates and length of match, factoring in that start coordinate is $qe/$se 
	# if match is on negative strand
	my ($qstart, $qend) = ($qs, $qe);
       ($qstart, $qend) = ($qe, $qs) if ($qe < $qs);
	my $qlen = ($qend - $qstart) + 1;

	# now mask out regions of CDS sequence with '1' corresponding to matching region
	substr($cds_sequences{$qid}, $qstart-1,$qlen) = ("1" x $qlen);

}
close($blast);

# can now loop through all sequences in %cds_sequences and calculate final results
my ($transcriptome_length, $total_match_length, $cds_hits) = (0, 0, 0);

foreach my $cds (sort {$a cmp $b or $a <=> $b} keys %cds_sequences){
	my $bases = $cds_sequences{$cds} =~ tr/1/1/;
	$total_match_length += $bases;

	my $length = length($cds_sequences{$cds});
	$transcriptome_length += $length;
	
	my $percent = sprintf("%.2f", ($bases / $length) * 100);
	print "$cds $bases/$length ($percent%)\n" if ($verbose);
	
	# count a CDS if present if there is more than 95% of it present
	$cds_hits++ if $percent >= 95;
}

# basic output is just four values
my $percent = sprintf("%.2f", ($total_match_length / $transcriptome_length) * 100);
print "Assembly $assembly_ID) $cds_hits CDSs present out of 176, $total_match_length of $transcriptome_length transcriptome nt present ($percent%)\n";
	
# do we have CSV output to print
if($csv){
	my $csv_out;
	my $csv_file = "$assembly_ID.blast.csv";
	open($csv_out, ">", "$csv_file") or die "Can't open $csv_file\n";	
	print $csv_out "Assembly,Number of CDSs present in assembly,Number of transcriptome nt present in assembly, Percentage of transcriptome present in assembly\n";
	print $csv_out "$assembly_ID,$cds_hits,$total_match_length,$percent\n";

	close($csv_out);
}

exit;