package uk.ac.liverpool.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;


/**
 * The class should take all the summary files from all the experiments across all the models, and 
 * find out 
 * 	- total unique spectra-peptide match 
 * 
 * @author riteshk
 *
 */
public class AnalyseTotalSpectraAndPeptideAcrossModels {
	
	/*
	// For Toxo dataset
	int totalFiles = 18;
	String [] summaryFiles = {"tmp/WholeSumary_sanya_GL.txt",
	                          "tmp/WholeSummary_1D_GL.txt",
	                          "tmp/WholeSummary_1D_AG.txt",
	                          "tmp/WholeSummary_1D_OF.txt",
	                          "tmp/WholeSummary_2D_AG.txt",
	                          "tmp/WholeSummary_2D_GL.txt",
	                          "tmp/WholeSummary_2D_OF.txt",
	                          "tmp/WholeSummary_dong_GL.txt",
	                          "tmp/WholeSummary_Dong_OF.txt",
	                          "tmp/WholeSummary_mudpit_AG.txt",
	                          "tmp/WholeSummary_Mudpit_GL.txt",
	                          "tmp/WholeSummary_mudpit_OF.txt",
	                          "tmp/WholeSummary_SFIF_AG.txt",
	                          "tmp/WholeSummary_SFIF_GL.txt",
	                          "tmp/WholeSummary_SFIF_OF.txt",
	                          "tmp/WholeSummaryDong_AG.txt",
	                          "tmp/WholeSummarySanya_AG.txt",
	                          "tmp/WholeSummarySanya_OF.txt"};
	*/
	// For N. Caninum dataset
	int totalFiles = 6;
	String [] summaryFiles = {"tmp/WholeSummary_Dong_GL.txt",
	                          "tmp/WholeSummary_Dong_OF.txt",
	                          "tmp/WholeSummary_sarah_GL.txt",
	                          "tmp/WholeSummary_Sarah_OF.txt",
	                          "tmp/WholeSummaryDong_AG.txt",
	                          "tmp/WholeSummarySarah_AG.txt"};
	
	double fdrThreshold = 0.1;
	String delimiter = "\t";
	String decoyString = "Rev";
	HashMap [] spectrumPeptideSeqMap = new HashMap[totalFiles];
	public int totalPSMIncludeingRedundant = 0;
	
	public HashMap <String, String>globalPSMmap = new HashMap<String, String>();
	
	public AnalyseTotalSpectraAndPeptideAcrossModels(){
		for (int i= 0; i< totalFiles; i++){
			spectrumPeptideSeqMap[i] = new HashMap<String, ArrayList<String>>();
		}
	}
	
	public void countPSMsInAllFiles(){
		
			for (int i= 0; i< totalFiles; i++){
				AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(summaryFiles[i], fdrThreshold, delimiter,decoyString);
				ap.createMaps();
				spectrumPeptideSeqMap[i] = ap.spectrumPeptideSeqMap;
			}
			
			for (int i = 0; i < totalFiles; i++){
				
				Iterator<String> spectrums = spectrumPeptideSeqMap[i].keySet().iterator();
				
				while(spectrums.hasNext()){
					String specID = spectrums.next();
					ArrayList<String> peptides = (ArrayList<String>) spectrumPeptideSeqMap[i].get(specID);
					
					for (int k = 0; k < peptides.size(); k++){
						String key_made = specID + peptides.get(k).trim();
						
						totalPSMIncludeingRedundant++;
						
						if(!globalPSMmap.containsKey(key_made))
							globalPSMmap.put(key_made, specID);
					}
				}
				
			}
			
	}
	
	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception{
		AnalyseTotalSpectraAndPeptideAcrossModels ats = new AnalyseTotalSpectraAndPeptideAcrossModels();
		ats.countPSMsInAllFiles();
		System.out.println(" Total non-redundant size = " + ats.globalPSMmap.size());
		System.out.println(" Total PSM including redundant = " + ats.totalPSMIncludeingRedundant);
	}
	
	
}
