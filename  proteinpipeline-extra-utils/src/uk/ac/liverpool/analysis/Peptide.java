package uk.ac.liverpool.analysis;

public class Peptide {
	
	public double fdr;
	public int start, end;
	public String specID, protAccn;
	
	public Peptide(double fdr, int start,int end,String specID, String protAccn){
		this.fdr = fdr;
		this.start = start;
		this.end = end;
		this.specID = specID;
		this.protAccn = protAccn;
	}
	
}
