package uk.ac.liverpool.signalP;

import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;

public class FastaSequenceReader {

	String fastaFile;
	
	public FastaSequenceReader(String fastaFile){
		this.fastaFile = fastaFile;
	}
	
	/**
	 * Accept a list of Protein Accessions and retrieve the sequences.
	 * Return a Map of Accession and Sequences
	 * @param accessionList
	 * @return
	 */
	public HashMap<String, String> getSequence(ArrayList<String> accessionList) throws Exception{
		
		HashMap<String, String> seqMap = new HashMap<String, String>();
		
		FileInputStream inStream = new FileInputStream( fastaFile );
			
		
		FastaReader<ProteinSequence,AminoAcidCompound> fastaReader = 
			new FastaReader<ProteinSequence,AminoAcidCompound>(
					inStream, 
					new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), 
					new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		
		LinkedHashMap<String, ProteinSequence> b = fastaReader.process();

		for( String accession : accessionList){
			Iterator<String> seqIds = b.keySet().iterator();
			while(seqIds.hasNext()){
				String seq_accn = seqIds.next();
				if(seq_accn.contains(accession)){
					String proteinSeq = b.get(seq_accn).getSequenceAsString();
					seqMap.put(accession, proteinSeq);
					break;
				}
			}
		}
				
		return seqMap;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String inputFasta = "/Users/riteshk/Dropbox/PipelineResults/Toxo/Signal-P/GenesWithSignalPeptide_Source-ToxoDB.fasta";
		
		ArrayList<String> accessionList = new ArrayList<String>();
		accessionList.add("TGME49_036000");
		accessionList.add("TGME49_058980");
		accessionList.add("TGME49_044580");
		
		
		try{
			FastaSequenceReader fsr = new FastaSequenceReader(inputFasta);
			HashMap<String, String> coll = fsr.getSequence(accessionList);
			
			// Print
			Iterator<String> acc = coll.keySet().iterator();
			while(acc.hasNext()){
				String accn = acc.next();
				System.out.println(">" + accn);
				System.out.println(coll.get(accn));
			}
		}catch(Exception e){
			e.printStackTrace();
		}

	}

}
