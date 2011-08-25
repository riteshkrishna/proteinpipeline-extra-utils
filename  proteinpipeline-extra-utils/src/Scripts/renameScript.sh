!/bin/sh
# This script is used for modifying the accession strings in the files produced
# by Augustus. This appends the Chromosome name ahead of the protein accession.
# Example - For Chr VIII, the protein ">g1.t" in protein file will become 
# ">TGME49_chrVIII.aa_g1.t". This will help in matching proteins to corresponding GFFs
 
FILE_DIR="*.aa"

for FILE_NAME in $FILE_DIR
do
	FASTA_INDICATOR=">"
	STR_TO_REPLACE="$FASTA_INDICATOR$FILE_NAME""_"
	
	echo "$STR_TO_REPLACE"
	sed  "s/$FASTA_INDICATOR/$STR_TO_REPLACE/g" $FILE_NAME > temp.tmp 
	mv temp.tmp $FILE_NAME
	echo "** DONE $FILE_NAME **"
done

