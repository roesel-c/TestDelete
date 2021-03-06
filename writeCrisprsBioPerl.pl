#!/usr/bin/perl
use warnings;
use strict;

use Bio::Seq;
use Bio::SeqIO;


#hash to store kmers
my %kMerHash = ();

#hash to store occurrences of last 12 positions
my %last12Counts = ();

my $seqio_obj = Bio::SeqIO->new(

	-file => 'dmel-all-chromosome-r6.17.fasta',

	#-file   => 'subset.fasta',
	-format => 'fasta'
);

while ( my $seq_obj = $seqio_obj->next_seq ) {
	my $seq = $seq_obj->seq;
	processSeq( \$seq );
}

sub processSeq {
	my ($seqRef) = @_;

	#declare scalars to characterize sliding window
	#Set the size of the sliding window
	my $windowSize = 21;

	#Set the step size
	my $stepSize  = 1;
	my $seqLength = length($$seqRef);

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
		my $crisprSeq = substr( $$seqRef, $windowStart, $windowSize );

#if the 21-mer ends in GG, create a hash with key=last 12 of k-mer and value is 21-mer
#Regex where $1 is the crispr, and $2 contains the last 12 of crispr.
		if ( $crisprSeq =~ /([ATGC]{9}([ATGC]{10}GG))$/ ) {

			#Put the crispr in the hash with last 12 as key, full 21 as value.
			$kMerHash{$2} = $1;
			$last12Counts{$2}++;

		}

	}
}
#open( FASTA_OUT, ">", 'crisprs1.fasta' ) or die $!;    #open or die
my $seqio_crisprs = Bio::SeqIO->new(
	-file   => '>crisprs1.fasta',
	-format => 'fasta'
);

#Initialize the CRISPR count to zero
my $crisprCount = 0;

#Loop through the hash of last 12 counts
for my $last12Seq ( sort ( keys %last12Counts ) ) {

	#Check if count == 1 for this sequence
	if ( $last12Counts{$last12Seq} == 1 ) {

		#The last 12 seq of this CRISPR is unique in the genome.
		#Increment the CRISPR count.
		$crisprCount++;
		my $seq_obj = Bio::Seq->new(
			-seq        => $kMerHash{$last12Seq},
			-desc       => "CRISPR",
			-display_id => "crispr_$crisprCount",
			-alphabet   => "dna"
		);
		$seqio_crisprs->write_seq($seq_obj);
	}
}
