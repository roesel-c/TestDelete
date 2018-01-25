#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;

# Load into a hash all 21-mers ending in GG from dmel-all-chromosome.r6.02.fasta
# where the key is the last 12 positions of the k-mer, and the value is the k-mer.
# Create a second hash to count how many times each 12-mer occurs in the genome.
# For each 12-mer that only occurs ONCE, the corresponding 21-mer is a potential CRISPR.
# Print the crisprs.fasta

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO;

#reading the fasta file using Bio::SeqIO
my $sequenceRef = Bio::SeqIO->new(
	-file   => "dmel-all-chromosome-r6.17.fasta",
	-format => 'fasta'
);

#hash to store occurrences of last 12 positions
my %last12Counts = ();

#hash to store kmers
my %kMerHash = ();

#loop through each sequence returned by Bio::SeqIO and get the 21mer substring
while ( my $seq_obj = $sequenceRef->next_seq ) {
	my $seq           = $seq_obj->seq;
	my $kmerSubstring = get21merSubstrings( \$seq );
}

sub get21merSubstrings {
	my ($sequence) = @_;

	#declare scalars to characterize sliding window
	#Set the size of the sliding window
	my $windowSize = 21;

	#Set the step size
	my $stepSize  = 1;
	my $seqLength = length($$sequence);

#for loop to increment the starting position of the sliding window
#starts at position zero; doesn't move past end of file; advance the window by step size
	for (
		my $windowStart = 0 ;
		$windowStart <= ( $seqLength - $windowSize ) ;
		$windowStart += $stepSize
	  )
	{

	   #Get a 21-mer substring from sequenceRef (two $ to deference reference to
	   #sequence string) starting at the window start for length $windowStart
		my $crisprSeq = substr( $$sequence, $windowStart, $windowSize );

#if the 21-mer ends in GG, create a hash with key=last 12 of k-mer and value is 21-mer
#Regex where $1 is the crispr, and $2 contains the last 12 of crispr.
		if ( $crisprSeq =~ /([ATGC]{9}([ATGC]{10}GG))$/ ) {

			#Put the crispr in the hash with last 12 as key, full 21 as value.
			$kMerHash{$2} = $1;
			$last12Counts{$2}++;

		}

	}
}

#Initialize the CRISPR count to zero
my $crisprCount = 0;

#writing crispr sequences to fasta file
my $crisprOutput =
  Bio::SeqIO->new( -file => '>crisprs1.fasta', -format => 'fasta' );

#Loop through the hash of last 12 counts
for my $last12Seq ( sort ( keys %last12Counts ) ) {

	#Check if count == 1 for this sequence
	if ( $last12Counts{$last12Seq} == 1 ) {

		#The last 12 seq of this CRISPR is unique in the genome.
		#Increment the CRISPR count.
		$crisprCount++;

		#Print the CRISPR in FASTA format.

		my $crispr_obj = Bio::Seq->new(
			-display_id => "crispr_$crisprCount",
			-desc       => 'CRISPR',
			-seq        => $kMerHash{$last12Seq}
		);

		$crisprOutput->write_seq($crispr_obj);

	}
}
