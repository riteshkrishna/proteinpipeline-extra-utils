package uk.ac.liverpool.gff;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import uk.ac.liverpool.analysis.AnalysePeptidesAcrossExperiments;

/**
 * Need this file for creating circos specific format. It will query a GFF for finding the co-ordinates of
 * proteins/genes and also query the protein-peptide output producded by AnalysePeptidesAcrossExperiments.java.
 * 
 * It will then produce a combined result as what was the count of peptides for each protein, and what was the 
 * location of those proteins on the chromosomes.
 * 
 * @author riteshk
 *
 */
public class GffAndPeptideFileCombo {

	public void process(String [] args) throws Exception{
		
		// change parameters if needed....
		double fdrThreshold = 0.01;
		String delimiter = "\t";
		String decoyString = "Rnd";
		
		// ** For Official model
		//String gffInput = "/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_GFF_produced/Toxo_1D_GFF_newoutput.gff";
		//String model = "official";
		
		// ** For Augustus model
		//String gffInput = "/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/Augustus_1D.gff";
		//String model = "augustus";
		
		// ** For Glimmer model
		String gffInput = "/Users/riteshk/Ritesh_Work/Toxo/Toxo_Predictions_Glimmer/gffwithfasta/GFF_Processed/All-Glimmer-ME49.gff";
		String model = "glimmer";
		
		String outputFile = "result/Gene-Peptide.txt";
		
		String [] summaryFiles = args;
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		
		AnalysePeptidesAcrossExperiments aae = new AnalysePeptidesAcrossExperiments(summaryFiles,fdrThreshold, delimiter,decoyString);
		aae.detectProteinsAndPeptides();
		  
		GFF3Routines gffHandle = new GFF3Routines();
		gffHandle.processGff(gffInput);
		
		HashMap<String, Gene_Information> geneRecords = new HashMap<String, Gene_Information>(gffHandle.geneRecords);
		HashMap<String, ArrayList<String>> allProteinsPeptides = new HashMap<String, ArrayList<String>>(aae.allProteinsPeptides);
		
		
		Iterator<String> proteins = allProteinsPeptides.keySet().iterator();
		while(proteins.hasNext()){
			String protAccn = proteins.next();
			
			int peptideCounts = allProteinsPeptides.get(protAccn).size();
			protAccn = parseProteinAccessionForGeneName(protAccn,model); // Needs to be parsed - apidb|cds_TGME49_113430-1
			
			if(protAccn != null){
				System.out.println(protAccn);
				Gene_Information geneObj = geneRecords.get(protAccn);
				String chromosome = geneObj.getSeqID();
				long start = geneObj.getStart();
				long end = geneObj.getEnd();
				
				out.write(chromosome + "\t" + start + "\t" + end + "\t" + peptideCounts + "\t" + protAccn + "\n");
			}
		}
		
		out.close();
	}
	
	/**
	 * 
	 * @param protein
	 * @return
	 */
	String parseProteinAccessionForGeneName(String protein,String model){
		String refinedAccn = null;
		
		 //Official models have accn like - apidb|cds_TGME49_113430-1, which need to be converted to apidb|TGME49_113430
		if(model.equalsIgnoreCase("official")){
			if(protein.contains("cds_")){
				protein = protein.replace("cds_","");
				protein = protein.replace("-1","");
				refinedAccn = protein;
			}
		}
		
		// Augustus models have TGME49_chrVIIa.aa_g226.t1, which need to converted to TGME49_chrVIIa.aa_g226
		if(model.equalsIgnoreCase("augustus")){
			if(protein.contains(".t1")){
				protein = protein.replace(".t1","");
				refinedAccn = protein;
			}
		}
		
		if(model.equalsIgnoreCase("glimmer")){
			refinedAccn = protein;
		}
		
		return refinedAccn;
	}
	
	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception{
		
		GffAndPeptideFileCombo gf = new GffAndPeptideFileCombo();
		gf.process(args);
		
	}
}
