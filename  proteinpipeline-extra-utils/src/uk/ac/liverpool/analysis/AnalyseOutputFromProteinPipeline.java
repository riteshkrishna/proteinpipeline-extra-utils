package uk.ac.liverpool.analysis;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Vector;

/* This program is for parsing the output produced by Proteomics Pipeline for analysis of results. The output is produced
 * in the following manner with \t delimiter and 15 columns in total.
 * --
 Rnd1psu|NCLIV_010750	0.710180051	Spectrum5154 scans: 6944	FPALAPATAPPMGAPR	0.710180051	1580.79	GROUP_ID	1166	1181	0.092796095	## 15.9949_M:12	3	1580.825	t	/ifs/h/riteshk/ProteomicsResults/Neo/master/mainresult/dir37_8/Summary__129.txt
 Rnd1psu|NCLIV_010750	0.69254807	Spectrum2153 scans: 3840	CGDRSGVLAGQAPMSNVNR	0.69254807	1987.98	GROUP_ID	516	534	0.15909638		3	1987.93	o	/ifs/h/riteshk/ProteomicsResults/Neo/master/mainresult/dir40_9/Summary__450.txt
 Rnd1psu|NCLIV_010750	0.665853679	Spectrum759 scans: 2280	DKMSCAWCYALR	0.665853679	1575.64	GROUP_ID	929	940	0.1351814	## -999_M:3	2	1575.66	o	/ifs/h/riteshk/ProteomicsResults/Neo/master/mainresult/dir33_8/Summary__694.txt
 .......
 .......
 *--
 *
 * We want the following information -
 * 	- Map for Protein  = Accessions -> {peptide [], score}
 *  - Map for peptides = seq -> {FDR, start, end, spectrum, protein}
 *  - Map for Spectrum = specID -> {pep-seq}
 */

public class AnalyseOutputFromProteinPipeline {

	// Set according to the columns in the summary file
	final int protAccn_column 		= 0;
	final int protScore_column 		= 1;
	final int specId_column 		= 2;
	final int pepSeq_column 		= 3;
	final int fdr_column 			= 4;
	final int start_column 			= 7;
	final int end_column 			= 8;
	final int searchEngine_column 	= 13;
	
	String wholeSummaryFile;
	double fdrThreshold;
	String delimiter;
	
	HashMap <String, ArrayList<String>> proteinPeptideMap;
	HashMap <String, ArrayList<Double>> proteinScoreMap;
	HashMap <String, ArrayList<String>> spectrumPeptideSeqMap;
	HashMap <String, ArrayList<Peptide>> peptideMap;
	
	
	/**
	 * 
	 * @param pipelineSummaryFile - The tab delimited summary file from pipeline
	 * @param fdrThreshold - The FDR threshold for filtering of result
	 */
	AnalyseOutputFromProteinPipeline(String pipelineSummaryFile, double fdrThreshold, String delimiter) {
		wholeSummaryFile = new String(pipelineSummaryFile);
		this.fdrThreshold = fdrThreshold;
		this.delimiter = delimiter;
		
		proteinPeptideMap = new HashMap <String, ArrayList<String>>();
		proteinScoreMap = new HashMap <String, ArrayList<Double>>();
		spectrumPeptideSeqMap = new HashMap <String, ArrayList<String>>();
		peptideMap = new HashMap <String, ArrayList<Peptide>>();
	}
	
	/**
	 * 
	 */
	void createMaps(){
		
		Scanner scanner ;
		
		try{
			scanner = new Scanner(new FileReader(new File(wholeSummaryFile)));
			
			String prevProtein = new String();
			double prevProteinScore = 0.0;
			
			while(scanner.hasNextLine()){
				 String line = scanner.nextLine();
				 if(line.isEmpty())
					 continue;
				 String [] values = line.split(this.delimiter);
				 
				 try{
				 if(values.length < 14)
					 throw new Exception(" Less than 14 columns found in the line \n " + line);
				 }catch(Exception e){
					 System.out.println("Exception :: " + e.getMessage());
				 }
				 
				 String protAccn;
				 if(values[protAccn_column].trim().isEmpty())
					 protAccn = prevProtein;
				 else protAccn = values[protAccn_column].trim();
				 
				 double protScore;
				 if(values[protScore_column].trim().isEmpty())
					 protScore = prevProteinScore;
				 else protScore = Double.parseDouble(values[protScore_column].trim());
				 
				 String specID = values[specId_column].trim();
				 String pepSeq = values[pepSeq_column].trim();
				 int start = Integer.parseInt(values[start_column].trim());
				 int end = Integer.parseInt(values[end_column].trim());
				 double fdrScore = Double.parseDouble(values[fdr_column].trim());
				 
				 // Skip if FDR score is greater than threshold
				 if (fdrScore >= this.fdrThreshold)
					 continue;
				 
				 // Fill the maps
				 
				 // Protein - Peptide Map
				 ArrayList<String> peptides;
				 if(proteinPeptideMap.containsKey(protAccn))
					 peptides = proteinPeptideMap.get(protAccn);
				 else
					 peptides = new ArrayList<String>();
				 
				peptides.add(pepSeq);
				proteinPeptideMap.put(protAccn, peptides);
				  
				// Spectrum - Peptide Map
				ArrayList<String> pepSeqColl;
				if(spectrumPeptideSeqMap.containsKey(specID))
					pepSeqColl = spectrumPeptideSeqMap.get(specID);
				else
					pepSeqColl = new ArrayList<String>();
				
				pepSeqColl.add(pepSeq);
				spectrumPeptideSeqMap.put(specID, pepSeqColl);
				
				// Protein - score map
				ArrayList<Double> scores;
				if(proteinScoreMap.containsKey(protAccn))
					scores = proteinScoreMap.get(protAccn);
				else 
					scores = new ArrayList<Double>();
				
				scores.add(protScore);
				proteinScoreMap.put(protAccn, scores);
				
				// Peptide Map
				ArrayList<Peptide> peptideColl;
				if(peptideMap.containsKey(pepSeq))
					peptideColl = peptideMap.get(pepSeq);
				else 
					peptideColl = new ArrayList<Peptide>();
				
				Peptide p = new Peptide(fdrScore,start,end,specID,protAccn);
				
				peptideColl.add(p);
				peptideMap.put(pepSeq, peptideColl);
				
			}	
			
			scanner.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/************************ Query functions below ****************************/
	
	/**
	 * 
	 * @return Get the number of unique peptides
	 */
	public int getNumberOfUniquePeptides(){
		return peptideMap.size();
	}
	
	/**
	 * 
	 * @return Get the number of unique proteins
	 */
	public int getNumberOfUniqueProteins(){
		return proteinScoreMap.size();
	}
	
	/**
	 * Find total true protein accessions
	 * @param decoyString - String to identify decoy and true accessions
	 * @return - ArrayList of true protein accessions 
	 */
	public ArrayList<String> getTrueProteins(String decoyString){
		
		ArrayList<String> totalTrueProteins = new ArrayList<String>();
		
		Iterator<String> protAccns = proteinScoreMap.keySet().iterator();
		while(protAccns.hasNext()){
			String protAccn = protAccns.next();
			if(!protAccn.contains(decoyString))
				totalTrueProteins.add(protAccn);
		}
		
		return totalTrueProteins;
	}

	/**
	 * Find total decoy protein accessions
	 * @param decoyString - String to identify decoy and true accessions
	 * @return - ArrayList of true protein accessions 
	 */
	public ArrayList<String> getDecoyProteins(String decoyString){
		
		ArrayList<String> totalDecoyProteins = new ArrayList<String>();
		
		Iterator<String> protAccns = proteinScoreMap.keySet().iterator();
		while(protAccns.hasNext()){
			String protAccn = protAccns.next();
			if(protAccn.contains(decoyString))
				totalDecoyProteins.add(protAccn);
		}
		
		return totalDecoyProteins;
	}
	
	/**
	 * Get the Proteins with unique peptide sequences 
	 * @return - HashMap with Protein Accn and peptide sequence
	 */
	public HashMap<String,String> getProteinsWithUniquePeptideSequences(){
		HashMap<String,String> uniqueProtPeptide = new HashMap<String,String>();
		
		Iterator <String> protAccns = proteinPeptideMap.keySet().iterator();
		while(protAccns.hasNext()){
			String accn = protAccns.next();
			if(proteinPeptideMap.get(accn).size() == 1){
				uniqueProtPeptide.put(accn, proteinPeptideMap.get(accn).get(0));
			}
		}
		
		return uniqueProtPeptide;
	}
	
	/**
	 * 
	 * @return ArrayList of 3 elements -
	 * 		ArrayList[0] = ArrayList of peptide sequences
	 * 	    ArrayList[1] = ArrayList of protein accessions
	 *  	ArrayList[2] = ArrayList of start positions
	 *  
	 *  Example use - To access the first record, we can write -
	 *  
	 *   ArrayList<ArrayList<String>> nTerminals = getN_TerminalPeptides();
	 *   String pepSeqOfFirstRecord = nTerminals.get(0).get(0);
	 *   String protAccOfFirstRecord = nTerminals.get(1).get(0);
	 *   int startOfFirstRecord     = Integer.parseInt(nTerminals.get(2).get(0));
	 */
	public ArrayList<ArrayList<String>> getN_TerminalPeptides(){
		
		ArrayList<ArrayList<String>> nTerminals = new ArrayList<ArrayList<String>>();
		
		ArrayList<String> pepSeqCollection   = new ArrayList<String>();
		ArrayList<String> protAccnCollection = new ArrayList<String>();
		ArrayList<String> startCollection   = new ArrayList<String>();
		
		Iterator<String> pepKey = peptideMap.keySet().iterator();
		
		while(pepKey.hasNext()){	
			String pepSeq = pepKey.next();
			ArrayList<Peptide> pepColl = peptideMap.get(pepSeq);
			
			for (int i = 0; i < pepColl.size(); i++){
				if( (pepColl.get(i).start == 1 ) || (pepColl.get(i).start == 2 )){
					pepSeqCollection.add(pepSeq);
					protAccnCollection.add(pepColl.get(i).protAccn);
					startCollection.add(Integer.toString(pepColl.get(i).start));
				}
			}
		}
		
		nTerminals.add(pepSeqCollection);
		nTerminals.add(pepSeqCollection);
		nTerminals.add(startCollection);
		
		return nTerminals;
	}
	
	
	// Use peptideMap to find how many proteins share a peptide
	
	public void makeProteinGroups(){
		// TODO
	}
	
	
	/**
	 * The test function..
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception{
		
		String pipelineSummaryFile = args[0];
		double fdrThreshold = Double.parseDouble(args[1]);
		String delimiter = args[2];
		String decoyString = args[3];
		
		AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(pipelineSummaryFile, fdrThreshold, delimiter);
		ap.createMaps();
		
		/********* Queries ***************/
		
		System.out.println("Total unique peptides found = " + ap.getNumberOfUniquePeptides());
		
		// Get number of total proteins
		System.out.println( "Total Proteins Found = " + ap.getNumberOfUniqueProteins());
		
		// Get Number of True proteins
		System.out.println("Total TRUE Proteins Found = " + ap.getTrueProteins(decoyString).size());
		
		// Get Number of Decoy proteins
		System.out.println( "Total DECOY Proteins Found = " + ap.getDecoyProteins(decoyString).size());
		
		System.out.println("Total N-Terminals found = " + ap.getN_TerminalPeptides().get(0).size());
		
		System.out.println(" Total proteins with UNIQUE peptides found = " + ap.getProteinsWithUniquePeptideSequences().size());
		
		
		
	}
}
