package uk.ac.liverpool.gff;

/**
 * Create object from Peptide File for the purpose of mapping co-ordinates.
 * 
 * @author riteshk
 *
 */
public class ProteinPeptideObject {
	
	String accession;
	String peptideSequence;
	long start;
	long end;

	public ProteinPeptideObject(String accession,String peptideSequence,long start,long end) {
		this.accession = accession;
		this.peptideSequence = peptideSequence;
		this.start = start;
		this.end = end;
	}
	
	public String getAccession(){
		return this.accession;
	}
	
	public String getpeptideSequence(){
		return this.peptideSequence;
	}
	
	public long getstart(){
		return this.start;
	}
	
	public long getEnd(){
		return this.end;
	}
}
