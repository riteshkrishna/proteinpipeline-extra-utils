package uk.ac.liverpool.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * This class is written to analyse peptides identified across different experiments for
 * a single gene model. Eg. - 1D, 2D, SFIF etc run against official gene models. We want
 * to know what are the total proteins and what are the non-redundant peptides for each 
 * proteins.
 * 
 * Example usage :
 * Official toxo data -
 * /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_SummaryFiles/WholeSummary_1D.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_SummaryFiles/WholeSummary_2D.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_SummaryFiles/WholeSummary_sf.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_SummaryFiles/WholeSummary_mudpitOfficial.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Official/Dong/WholeSummary_Dong_orbitrap_official.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Official/Sanya/WholeSummarySanyaOrbitrap_official.txt
 * 
 * Augustus - 
 * /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/WholeSummaryToxoAugustus_mudpit_solAndInsol.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/WholeSummary_SFIF_augustus.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/WholeSummary_2D_augustus.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/WholeSummary_1D_augustus.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Augustus/Dong/WholeSummaryToxoOrbitrapAugustus.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Augustus/Sanya/WholeSummarySanyaAugustus.txt
 * 
 * Glimmer -
 * /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Glimmer/Toxo_Glimmer/WholeSummary_Mudpit_glimmer.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Glimmer/Toxo_Glimmer/WholeSummary_SFIF_glimmer.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Glimmer/Toxo_Glimmer/WholeSummary_2D_glimmer.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Glimmer/Toxo_Glimmer/WholeSummary_1D_glimmer.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Glimmer/Sanya/WholeSumary_sanya_glimmer.txt /Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Glimmer/Dong/WholeSummary_dong_glimmer.txt  
 * @author riteshk
 *
 */
public class AnalysePeptidesAcrossExperiments {
	
	int noOfExperiments;
	String [] pipelineSummaryFiles;
	double fdrThreshold;
	String delimiter;
	String decoyString;
	
	HashMap [] proteinPeptideMaps;
	public HashMap<String, ArrayList<String>> allProteinsPeptides = new HashMap<String, ArrayList<String>>();
	
	public AnalysePeptidesAcrossExperiments(String [] summaryFiles, double fdrThreshold, String delimiter,String decoyString) {
		this.fdrThreshold = fdrThreshold;
		this.delimiter = delimiter;
		this.decoyString = decoyString;
		
		this.noOfExperiments = summaryFiles.length;
		
		pipelineSummaryFiles = new String[this.noOfExperiments];
		proteinPeptideMaps = new HashMap[this.noOfExperiments];
		
		for (int i=0; i< noOfExperiments; i++){
			pipelineSummaryFiles[i] = new String(summaryFiles[i]);
			proteinPeptideMaps[i] = new HashMap<String, ArrayList<String>>();
		}
	}
	
	/**
	 * Construct Protein-Peptide Map for all experiments by removing decoy entries
	 */
	void constructPeptideMaps(){
	
		for(int i = 0; i < noOfExperiments; i++){
			AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(pipelineSummaryFiles[i], fdrThreshold, delimiter,decoyString);
			ap.createMaps();
			
			HashMap<String, ArrayList<String>> protMap = new HashMap<String, ArrayList<String>>(ap.proteinPeptideMap);
			HashMap<String, ArrayList<String>> protMap_copy = new HashMap<String, ArrayList<String>>(ap.proteinPeptideMap);
			
			// Remove decoy accessions from the Map
			Iterator<String> accns = protMap.keySet().iterator();
			while(accns.hasNext()){
				String accn = accns.next();
				
				try{
				if(accn.contains(this.decoyString)){
					protMap_copy.remove(accn);
				}
				}catch(Exception e){
					e.printStackTrace();
				}
			}
			
			proteinPeptideMaps[i] = protMap_copy;
		}	
	}
	
	
	public void detectProteinsAndPeptides(){
		
		constructPeptideMaps();
		
		for(int i = 0 ; i < noOfExperiments; i++){
			HashMap<String, ArrayList<String>>  thisExp = proteinPeptideMaps[i];
			
			Iterator<String> proteins = thisExp.keySet().iterator();
			while(proteins.hasNext()){
				String accn = proteins.next();
				ArrayList<String> pepSeqs = thisExp.get(accn);
				
				if(allProteinsPeptides.containsKey(accn)){
					ArrayList<String> peptidesAlreadyFound = allProteinsPeptides.get(accn);
					
					// Since the peptides found by SE are already unique(seq + location wise),
					// so we don't need to compare just the sequences, and can add the whole seq to the list
					for(String pep : pepSeqs){
							peptidesAlreadyFound.add(pep);
					}
					allProteinsPeptides.put(accn, peptidesAlreadyFound);
				}else{
					allProteinsPeptides.put(accn, pepSeqs);
				}
			}
		}
		
	}
	
	/**
	 * 
	 * @param outFile
	 */
	public void writeProteinPeptideMapForAllExpInFile(String outFile){
		try{
			BufferedWriter out_prot = new BufferedWriter(new FileWriter(outFile));
			
			Iterator<String> proteins = allProteinsPeptides.keySet().iterator();
			while(proteins.hasNext()){
				String accn = proteins.next();
				ArrayList<String> peptides = allProteinsPeptides.get(accn);
				int totalPeptides = peptides.size();
				
				out_prot.write(accn + "\t" + totalPeptides + "\t" + peptides.toString() + "\n");
			}
			out_prot.close();
		}catch(Exception e){
			e.printStackTrace();
		}	
	}
	
	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception{
		
		// change parameters if needed....
		double fdrThreshold = 0.01;
		String delimiter = "\t";
		String decoyString = "Rnd";
		
		String [] summaryFiles = args;
		
		AnalysePeptidesAcrossExperiments aae = new AnalysePeptidesAcrossExperiments(summaryFiles,fdrThreshold, delimiter,decoyString);
		aae.detectProteinsAndPeptides();
		
		System.out.println(aae.allProteinsPeptides.toString());
		System.out.println(aae.allProteinsPeptides.size());
		
		aae.writeProteinPeptideMapForAllExpInFile("result/Protein-Peptide-Count.txt");
	}
}
