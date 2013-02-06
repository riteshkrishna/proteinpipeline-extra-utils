package uk.ac.liverpool.analysis2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;

/**
 * This program maps 1% FDR filtered peptides on protein sequences obtained from Toxo ver 8.0 GFF.
 * The protein sequences are tryptically digested with 1 allowed cleavage and the peptides are searched
 * for their presence. If a peptide is found, the start and end locations are recored in a TSV file.
 * @author riteshk
 *
 */
public class MapPeptidesOnToxoVersion8 {
	String fastaFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/Version_8-mapping/ToxoDB-8.0_TgondiiME49.gff.fasta";
	String peptideFile  = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/Version_8-mapping/All-peptide-FDR-filtered-uniq-list.txt";
	String outputFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/Version_8-mapping/All-peptide-mappedOnProteins.txt";
	
	HashMap<String, ArrayList<String>> peptideProteinsMap = new HashMap<String, ArrayList<String>>();
	ArrayList<String> fdrFilteredPeptides = new ArrayList<String>();
	
	/*
	 * Read all the FDR filtered peptides from a file in an ArrayList
	 */
	public void readFdrFilteredPeptides(){
		try{
			BufferedReader in =  new BufferedReader(new FileReader(peptideFile));
			String line;
		
			while((line = in.readLine()) != null){
				fdrFilteredPeptides.add(line.trim());
			}
			
			in.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	
	public void searchPeptidesInFasta() throws Exception{
		
		readFdrFilteredPeptides();
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		
		// Digest Fasta (digestion criteria coded in ProcessFasta)
		ProcessFasta ps = new ProcessFasta();
    	ps.init(fastaFile);
    	
    	// Create a peptide -> [prot1, prot 2,..] map
    	Iterator<String> proteins = ps.accToPeptides.keySet().iterator();
    	while(proteins.hasNext()){
    		String protein = proteins.next();
    		String [] peptides = ps.accToPeptides.get(protein);
    		
    		for(String peptide : peptides){
    			
    			ArrayList<String> proteinsForThisPeptide;
    			if(peptideProteinsMap.containsKey(peptide))
    				proteinsForThisPeptide = peptideProteinsMap.get(peptide);
    			else
    				proteinsForThisPeptide = new ArrayList<String>();
    			
    			if(!proteinsForThisPeptide.contains(protein)){
					proteinsForThisPeptide.add(protein);
					peptideProteinsMap.put(peptide, proteinsForThisPeptide);
    			}
    		}
    	}
    	
    	// Search the location of each peptide in the proteins
    	int counter_notFoundPeptides = 0;
    	for(String peptide : fdrFilteredPeptides){
    		
    		if(!peptideProteinsMap.containsKey(peptide)){
    			counter_notFoundPeptides++;
    			System.out.println("Not found : " + counter_notFoundPeptides + "\t" + peptide);
    			continue;
    		}
    		
    		ArrayList<String> proteinsForThisPeptide = peptideProteinsMap.get(peptide);
    		
    		for(String protein : proteinsForThisPeptide){
    			String sequence = ps.accToSeq.get(protein);
    			
    			if(sequence.contains(peptide)){
    				int startIndex = sequence.indexOf(peptide) + 1; // To make it 1-based
    				int endIndex = startIndex + peptide.length() - 1; // To make it 1-based
    				
    				String toWrite = protein + "\t" + peptide + "\t" + startIndex +"\t" + endIndex + "\n";
    				out.write(toWrite);
    			}
    		}
    	}
    	
    	out.close();
	}
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String [] args) throws Exception{
		MapPeptidesOnToxoVersion8 map = new MapPeptidesOnToxoVersion8();
		map.searchPeptidesInFasta();
	}
}
