# Counts total number of Sequences and Amino Acids in the file

use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => "file.fa", '-format' => 'Fasta');

my $totalAminoAcid = 0;
my $totalSequences = 0;

while(my $seq = $seqio->next_seq) {
  my $string = $seq->seq;
  print "Sequence = ", $seq->id, "->", length($string), "\n";
  
  $totalAminoAcid = $totalAminoAcid + length($string);
  $totalSequences = $totalSequences + 1;	
}
print "Total Sequences = ", $totalSequences,", Total Amino Acid Count = ", $totalAminoAcid , "\n";
