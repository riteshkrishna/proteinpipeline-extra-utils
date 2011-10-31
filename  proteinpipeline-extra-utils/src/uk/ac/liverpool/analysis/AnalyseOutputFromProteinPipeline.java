package uk.ac.liverpool.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

import uk.ac.liverpool.inference.ProteinAmbiguityGrouping;

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
	String decoyString; 
	
	public HashMap <String, ArrayList<String>> proteinPeptideMap;
	public HashMap <String, ArrayList<Double>> proteinScoreMap;
	public HashMap <String, ArrayList<String>> spectrumPeptideSeqMap;
	public HashMap <String, ArrayList<Peptide>> peptideMap;
	
	/**
	 * 
	 * @param pipelineSummaryFile - Summary file for whole dataset
	 * @param fdrThreshold 		  - FDR to be used for analysis
	 * @param delimiter			  - \t or , or \w
	 * @param decoyString		  - Decoy identifier string	
	 */
	public AnalyseOutputFromProteinPipeline(String pipelineSummaryFile, double fdrThreshold, String delimiter,String decoyString) {
		wholeSummaryFile = new String(pipelineSummaryFile);
		this.fdrThreshold = fdrThreshold;
		this.delimiter = delimiter;
		this.decoyString = decoyString;
		
		proteinPeptideMap = new HashMap <String, ArrayList<String>>();
		proteinScoreMap = new HashMap <String, ArrayList<Double>>();
		spectrumPeptideSeqMap = new HashMap <String, ArrayList<String>>();
		peptideMap = new HashMap <String, ArrayList<Peptide>>();
	}
	
	/**
	 * 
	 */
	public void createMaps(){
		
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
				
				prevProtein = protAccn;
				prevProteinScore = protScore;
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
	 * @return - ArrayList of true protein accessions 
	 */
	public ArrayList<String> getTrueProteins(){
		
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
	 * @return - ArrayList of true protein accessions 
	 */
	public ArrayList<String> getDecoyProteins(){
		
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
	
	
	/**
	 * 
	 */
	
	public void makeProteinGroups(String pagListOutputFile){
		
		// Remove decoy proteins from PeptideMap
		HashMap<String,ArrayList<Peptide>> cleanPeptideMap = new HashMap<String,ArrayList<Peptide>>();
		
		Iterator <String> pepkey = this.peptideMap.keySet().iterator();
	    while(pepkey.hasNext()){
	    	String pepseq = pepkey.next();
	    	ArrayList<Peptide> pepColl = this.peptideMap.get(pepseq);
	    	
	    	ArrayList<Peptide> cleanPepColl = new ArrayList<Peptide>();
	    	HashMap<String,String> seenAccessions = new HashMap<String,String>();
	    	
	    	for(int i = 0 ; i < pepColl.size(); i++){
	    		String protAccn = pepColl.get(i).protAccn;
	    		if(!protAccn.contains(this.decoyString) && !seenAccessions.containsKey(protAccn)){
	    				seenAccessions.put(protAccn, "");
	    				cleanPepColl.add(pepColl.get(i));
	    		}	
	    	}
	    	if(!cleanPepColl.isEmpty())
	    		cleanPeptideMap.put(pepseq, cleanPepColl);
	    }
		////////////
		try{
			BufferedWriter out_pag = new BufferedWriter(new FileWriter(pagListOutputFile));
		
			ProteinAmbiguityGrouping pag = new ProteinAmbiguityGrouping(cleanPeptideMap,this.proteinPeptideMap);
			HashMap<Integer,ArrayList<String>> pagsFormed = pag.createAmbiguityGroup();
		
			out_pag.write("Total Pags found - " + pagsFormed.size());
			
			Iterator<Integer> pagKeys = pagsFormed.keySet().iterator();
			while(pagKeys.hasNext()){
				Integer key = pagKeys.next();
				ArrayList<String> values = pagsFormed.get(key);
			
				out_pag.write("\n Pag -> " + key.intValue() +"\n");
				for(int i=0;i<values.size(); i++){
					out_pag.write(values.get(i) + "\t");
				}
			}
			
			out_pag.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		// Dump the PeptideMap
		try {
		    BufferedWriter out = new BufferedWriter(new FileWriter(pagListOutputFile + "_Peptide_Dump.txt"));
		    
		    Iterator <String> pepList = cleanPeptideMap.keySet().iterator();
		    while(pepList.hasNext()){
		    	String pepseq = pepList.next();
		    	ArrayList<Peptide> pepColl = cleanPeptideMap.get(pepseq);
		    	if(pepColl == null){
		    		System.out.println("null ->" + pepseq);
		    		continue;
		    	}
		    	out.write("\n\nPepSeq = " + pepseq + "\n");
		    	for(int i=0 ; i < pepColl.size(); i++){
		    		out.write(pepColl.get(i).protAccn + "\t");
		    	}
		    }
		    out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	

	/**
	 * 
	 * @param trueProteins
	 * @param trueProteinListOutputFile
	 */
	public void writeProteinsToFile(ArrayList<String> trueProteins, String trueProteinListOutputFile){
		try{
			BufferedWriter out_prot = new BufferedWriter(new FileWriter(trueProteinListOutputFile));
			
			for(String protein : trueProteins){
				ArrayList<String> peptides = proteinPeptideMap.get(protein);
				
				out_prot.write("\n" + protein);
				for(int i = 0; i < peptides.size(); i++){
					String peptide = peptides.get(i);
					out_prot.write("\t" + peptide);
				}
			}
			
			out_prot.close();
		}catch(Exception e){
			
		}
	}
	
	/**
	 * Filter proteins based on the minimum number of peptides
	 * @param noOfPeptideForFilter
	 * @return
	 */
	public ArrayList<String> getTrueProteins_filteredOnNumberOfPeptides(int noOfPeptideForFilter){
		ArrayList<String> totalTrueProteinsFilteredOnPeptides = new ArrayList<String>();
		
		// First find out the true proteins
		ArrayList<String> totalTrueProteins = getTrueProteins();
		
		// Then filter them according to number of Peptides
		for(String protein : totalTrueProteins){
			if (proteinPeptideMap.get(protein).size() >= noOfPeptideForFilter)
				totalTrueProteinsFilteredOnPeptides.add(protein);
		}
		return totalTrueProteinsFilteredOnPeptides;
	}
	
	
	/**
	 * The test function..
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception{
		
		try{
			
		String pipelineSummaryFile = args[0];
		double fdrThreshold = Double.parseDouble(args[1]);
		String delimiter = args[2];
		String decoyString = args[3];
		
		int noOfPeptideThreshold = 2; // Used for filtering proteins based on number of peptides
		
		AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(pipelineSummaryFile, fdrThreshold, delimiter,decoyString);
		ap.createMaps();
		
		/******************************** Queries **************************************************/
		
		System.out.println("Total unique peptides found = " + ap.getNumberOfUniquePeptides());
		
		// Get number of total proteins
		System.out.println( "Total Proteins Found = " + ap.getNumberOfUniqueProteins());
		
		// Get Number of True proteins
		System.out.println("Total TRUE Proteins Found = " + ap.getTrueProteins().size());
		
		// Get Number of Proteins - with atleast 2 peptide criteria
		System.out.println(" Total TRUE Proteins with 2 PEPTIDE threshold = " + ap.getTrueProteins_filteredOnNumberOfPeptides(noOfPeptideThreshold).size());
		
		// Get Number of Decoy proteins
		System.out.println( "Total DECOY Proteins Found = " + ap.getDecoyProteins().size());
		
		System.out.println("Total N-Terminals found = " + ap.getN_TerminalPeptides().get(0).size());
		
		System.out.println(" Total proteins with UNIQUE peptides found = " + ap.getProteinsWithUniquePeptideSequences().size());
		
		
		/***************  Create PAGS **********************/
		String pagListOutputFile = "result/PAG_FOUND_" + Double.toString(fdrThreshold);
		ap.makeProteinGroups(pagListOutputFile);
		
		System.out.println("Done..");
		
		/************** Write details of true proteins ***************/
		String trueProteinListOutputFile = "result/TRUE_PROTEINS_" + Double.toString(fdrThreshold);
		ap.writeProteinsToFile(ap.getTrueProteins_filteredOnNumberOfPeptides(noOfPeptideThreshold),  trueProteinListOutputFile);
		
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
