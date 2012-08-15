package uk.ac.liverpool.analysis;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

/**
 * Take a set of peptide sequences and their protein accessions, along with the
 * list of proteoannotator result files. We can find out that how many spectrums
 * were identified corresponding to these peptides.
 * 
 *  Takes as input, the files producded from PeptideMappingOnChromosome.java
 *  with the following fields... 
 *  gb|TGME49_chrXI	820404	820452	SSVGIDVEDDIMDAFFK	gb|TGME49_chrXI.path1.gene136
    gb|TGME49_chrXI	1238985	1239045	RLCAELFALAIHVFGLLLLDK	gb|TGME49_chrXI.path1.gene205

 * @author riteshk
 *
 */
public class CountSpectrumsForGivenPeptides {
	
	// For Toxo dataset
	int totalFiles_toxo = 6;
	String [] summaryFiles_augustus_toxo = {
	                          "tmp/WholeSummary_1D_AG.txt",
	                          "tmp/WholeSummary_2D_AG.txt",
	                          "tmp/WholeSummary_mudpit_AG.txt",
	                          "tmp/WholeSummary_SFIF_AG.txt",
	                          "tmp/WholeSummaryDong_AG.txt",
	                          "tmp/WholeSummarySanya_AG.txt"
	                          };
	String [] summaryFiles_glimmer_toxo = {
							 "tmp/WholeSummary_sanya_GL.txt",
							 "tmp/WholeSummary_1D_GL.txt",
							 "tmp/WholeSummary_2D_GL.txt",
							 "tmp/WholeSummary_dong_GL.txt",
							 "tmp/WholeSummary_Mudpit_GL.txt",
							 "tmp/WholeSummary_SFIF_GL.txt"
            				 };

	
	// For N. Caninum dataset
	int totalFiles_neo = 2;
	String [] summaryFiles_augustus_neo = {
	                          "tmp/WholeSummaryDong_AG.txt",
	                          "tmp/WholeSummarySarah_AG.txt"
	                          };
	
	String [] summaryFiles_glimmer_neo = {
								"tmp/WholeSummary_Dong_GL.txt",
            					 "tmp/WholeSummary_sarah_GL.txt"
    							};


	String peptideList_augustus_toxo = "tmp/mapped_augustus_peptides.txt";
	String peptideList_glimmer_toxo  = "tmp/mapped_glimmer_peptides.txt";
	
	String peptideList_augustus_neo = "";
	String peptideList_glimmer_neo  = "";
	
	
	double fdrThreshold = 0.01;
	String delimiter = "\t";
	String decoyString = "Rev";
	
	HashMap [] PeptideMap_augutus_toxo  = new HashMap[totalFiles_toxo];
	HashMap [] PeptideMap_glimmer_toxo  = new HashMap[totalFiles_toxo];
	
	HashMap [] PeptideMap_augutus_neo  = new HashMap[totalFiles_neo];
	HashMap [] PeptideMap_glimmer_neo  = new HashMap[totalFiles_neo];
	
	/**
	 * 
	 */
	public void makeMaps_toxo(){
		for (int i= 0; i< totalFiles_toxo; i++){
			AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(summaryFiles_augustus_toxo[i], fdrThreshold, delimiter,decoyString);
			ap.createMaps();
			PeptideMap_augutus_toxo[i] = ap.peptideMap;
		}
		
		for (int i= 0; i< totalFiles_toxo; i++){
			AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(summaryFiles_glimmer_toxo[i], fdrThreshold, delimiter,decoyString);
			ap.createMaps();
			PeptideMap_glimmer_toxo[i] = ap.peptideMap;
		}
	}
	
	/**
	 * 
	 */
	public void makeMaps_neo(){
		for (int i= 0; i< totalFiles_neo; i++){
			AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(summaryFiles_augustus_neo[i], fdrThreshold, delimiter,decoyString);
			ap.createMaps();
			PeptideMap_augutus_neo[i] = ap.peptideMap;
		}
		
		for (int i= 0; i< totalFiles_neo; i++){
			AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(summaryFiles_glimmer_neo[i], fdrThreshold, delimiter,decoyString);
			ap.createMaps();
			PeptideMap_glimmer_neo[i] = ap.peptideMap;
		}
	}
	
	/**
	 * 
	 * @param peptideFile
	 * @return
	 */
	public ArrayList<String> readPeptideList(String peptideFile){
		ArrayList<String> pepSeqs = new ArrayList<String>();
		
		try{
			Scanner scanner = new Scanner(new FileReader(new File(peptideFile)));
			while(scanner.hasNextLine()){
				 String line = scanner.nextLine();
				 if(line.isEmpty())
					 continue;
				 String [] values = line.split("\t");
				 
				 if(values.length < 5){
					 System.out.println("Less than 5 columns.." + line);
				 }
				 String pepSeq = values[3].trim();
				 pepSeqs.add(pepSeq);
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}
		return pepSeqs;
	}
	
	/**
	 * 
	 * @param pepList
	 * @param peptideMap_model_species
	 * @return
	 */
	public HashMap<String,ArrayList<String>> findSpecCount(ArrayList<String> pepList,HashMap [] peptideMap_model_species){
		HashMap<String,ArrayList<String>> pepSpecMap = new HashMap<String,ArrayList<String>>();
		
		for(String pepseq : pepList){
			ArrayList<String> specColl = new ArrayList<String>();
			
			for(HashMap thisPepMap : peptideMap_model_species){
				if(thisPepMap.containsKey(pepseq)){
					ArrayList<Peptide> peptideColl = (ArrayList<Peptide>)thisPepMap.get(pepseq);
					for (Peptide p : peptideColl){
						String specID = p.specID.trim(); 
						if(!specColl.contains(specID)) // keep only unique IDs
							specColl.add(specID);
					}
				}
			}
			/*
			try{
			if(pepSpecMap.containsKey(pepseq))
			{
				System.out.println("Already present : " + pepseq);
				//System.in.read();
			}
			}catch(Exception e){
				
			}*/
			
			pepSpecMap.put(pepseq, specColl);
		}
		
		return pepSpecMap;
	}
	
	public void count_toxo(){
		makeMaps_toxo();
		
		ArrayList<String> pepList_augustus = readPeptideList(peptideList_augustus_toxo);
		HashMap<String,ArrayList<String>> pep_spec_augustus = findSpecCount(pepList_augustus,this.PeptideMap_augutus_toxo);
		
		ArrayList<String> pepList_glimmer = readPeptideList(peptideList_glimmer_toxo);
		HashMap<String,ArrayList<String>> pep_spec_glimmer = findSpecCount(pepList_glimmer,this.PeptideMap_glimmer_toxo);
		
		// Print..
		for(String pep : pep_spec_augustus.keySet()){
			System.out.println("Augustus" + "\t" + pep + "\t" + pep_spec_augustus.get(pep).size() + "\t" + pep_spec_augustus.get(pep).toString());
		}
		
		for(String pep : pep_spec_glimmer.keySet()){
			System.out.println("Glimmer" + "\t" + pep + "\t" + pep_spec_glimmer.get(pep).size() + "\t" + pep_spec_glimmer.get(pep).toString());
		}
	}
	
	/**
	 * Similar steps like toxo for neo to
	 */
	public void count_neo(){
		//TODO
		//makeMaps_neo();
	}
	
	/**
	 * 
	 */
	public  static void main(String []args) throws Exception{
		
		CountSpectrumsForGivenPeptides cs = new CountSpectrumsForGivenPeptides();
		cs.count_toxo();
	}
	
	
	
}
