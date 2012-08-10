package uk.ac.liverpool.analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

public class Classify_peptides_in_official_alternate {

	public HashSet<String> peptides_official = new HashSet<String>();
	public HashSet<String> peptides_augutus = new HashSet<String>();
	public HashSet<String> peptides_glimmer = new HashSet<String>();
	
	public HashSet<String> peptides_augustus_only = new HashSet<String>();
	public HashSet<String> peptides_glimmer_only = new HashSet<String>();
	public HashSet<String> peptides_augustus_glimmer_only = new HashSet<String>();
	
	public HashMap <String, ArrayList<Peptide>> official_peptideMap = new HashMap <String, ArrayList<Peptide>>();
	public HashMap <String, ArrayList<Peptide>> augustus_peptideMap = new HashMap <String, ArrayList<Peptide>>();
	public HashMap <String, ArrayList<Peptide>> glimmer_peptideMap = new HashMap <String, ArrayList<Peptide>>();
	
	/**
	 * 
	 * @param protein_Peptide_Count_File_official
	 * @param protein_Peptide_Count_File_augutus
	 * @param protein_Peptide_Count_File_glimmer
	 * @param trypticPeptides
	 */
	public void classifyPeptides(String protein_Peptide_Count_File_official, String protein_Peptide_Count_File_augutus,
			String protein_Peptide_Count_File_glimmer,HashSet<String> trypticPeptides){
		
		
		Process_output_from_AnalysePeptidesAcrossExperiments pa_o = new Process_output_from_AnalysePeptidesAcrossExperiments(protein_Peptide_Count_File_official);
		peptides_official = pa_o.getPeptideList();
		System.out.println("Total Official Peptides in file = " + peptides_official.size());	
		
		Process_output_from_AnalysePeptidesAcrossExperiments pa_a = new Process_output_from_AnalysePeptidesAcrossExperiments(protein_Peptide_Count_File_augutus);
		peptides_augutus = pa_a.getPeptideList();
		System.out.println("Total Augustus Peptides in file = " + peptides_augutus.size());
		
		Process_output_from_AnalysePeptidesAcrossExperiments pa_g = new Process_output_from_AnalysePeptidesAcrossExperiments(protein_Peptide_Count_File_glimmer);
		peptides_glimmer = pa_g.getPeptideList();
		System.out.println("Total Glimmer Peptides in file = " + peptides_glimmer.size());
		
	
		//--  All alternate peptides-------------------------------//
		// Add Augustus peptides
		Iterator<String> pep_itr = peptides_augutus.iterator();
		while(pep_itr.hasNext()){
			String pep = pep_itr.next();
			if(!peptides_official.contains(pep))
				peptides_augustus_glimmer_only.add(pep);
		}
		// Add Glimmer peptides
		pep_itr = peptides_glimmer.iterator();
		while(pep_itr.hasNext()){
			String pep = pep_itr.next();
			if(!peptides_official.contains(pep))
				peptides_augustus_glimmer_only.add(pep);
		}
		//---------------------------------------------------------//
		
		//----------------- Only Glimmer peptides-----------------//
		pep_itr = peptides_augustus_glimmer_only.iterator();
		while(pep_itr.hasNext()){
			String pep = pep_itr.next();
			if(!peptides_augutus.contains(pep))
				peptides_glimmer_only.add(pep);
		}

		//---------------------------------------------------------//
		
		//----------------- Only Augustus peptides-----------------//
		pep_itr = peptides_augustus_glimmer_only.iterator();
		while(pep_itr.hasNext()){
			String pep = pep_itr.next();
			if(!peptides_glimmer.contains(pep))
				peptides_augustus_only.add(pep);
		}
		//---------------------------------------------------------//
		
		//-------- Perform Tryptic cleaning here----------------------//
		//if any of the alternate peptide belongs to trypticPeptides, then ignore that
		HashSet<String> peptides_augustus_glimmer_only_copy = new HashSet<String>(peptides_augustus_glimmer_only);
		Iterator<String >pepIterator = peptides_augustus_glimmer_only_copy.iterator();
		while(pepIterator.hasNext()){
			String pep = pepIterator.next();
			if(trypticPeptides.contains(pep))
				peptides_augustus_glimmer_only.remove(pep);
		}
		
		HashSet<String> peptides_augustus_only_copy = new HashSet<String>(peptides_augustus_only);
		pepIterator = peptides_augustus_only_copy.iterator();
		while(pepIterator.hasNext()){
			String pep = pepIterator.next();
			if(trypticPeptides.contains(pep))
				peptides_augustus_only.remove(pep);
		}
		
		HashSet<String> peptides_glimmer_only_copy = new HashSet<String>(peptides_glimmer_only);
		pepIterator = peptides_glimmer_only_copy.iterator();
		while(pepIterator.hasNext()){
			String pep = pepIterator.next();
			if(trypticPeptides.contains(pep))
				peptides_glimmer_only.remove(pep);
		}
		//---------------------------------------------------------//
		
		System.out.println("Total Alternate Peptides = " + peptides_augustus_glimmer_only.size());
		System.out.println("Total Augustus only Peptides = " + peptides_augustus_only.size());
		System.out.println("Total Glimmer only  Peptides = " + peptides_glimmer_only.size());
	}
	
	/**
	 * 
	 * @param alternatePeptideFile
	 * @param officialPeptideFile
	 */
	public void writePeptidesToFile(String alternatePeptideFile, String officialPeptideFile){
		//---------------------------------------------------------//
		// Write to file
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(alternatePeptideFile));
			
			out.write("Unique Augustus peptides\n");
			Iterator<String>pep_itr = peptides_augustus_only.iterator();
			while(pep_itr.hasNext()){
				out.write(pep_itr.next() + "\n");
			}
			
			out.write("\n\n Unique Glimmer peptides\n");
			pep_itr = peptides_glimmer_only.iterator();
			while(pep_itr.hasNext()){
				out.write(pep_itr.next() + "\n");
			}
		
			out.close();
			
			BufferedWriter out_2 = new BufferedWriter(new FileWriter(officialPeptideFile));
			out_2.write("All Official peptides\n");
			pep_itr = peptides_official.iterator();
			while(pep_itr.hasNext()){
				out_2.write(pep_itr.next() + "\n");
			}
		
			out_2.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * The main function which does everything related to peptide classification
	 */
	public void perform_classification_Toxo()
	{
		String protein_Peptide_Count_File_official = new String("result/Protein-Peptide-Count-official-toxo.txt");
		String protein_Peptide_Count_File_augutus = new String("result/Protein-Peptide-Count-Augustus-toxo.txt");
		String protein_Peptide_Count_File_glimmer = new String("result/Protein-Peptide-Count-Glimmer-toxo.txt");
		
		String trypticFile = "result/Digested-Peptides-Toxo.txt";
		TrypticPeptideReader tr = new TrypticPeptideReader(trypticFile);
		HashSet<String> trypticPeptides = tr.getTrypticPeptides();
		
		classifyPeptides(protein_Peptide_Count_File_official, protein_Peptide_Count_File_augutus, protein_Peptide_Count_File_glimmer,trypticPeptides);
		
		String alternatePeptideFile = "result/Toxo-Peptide-Augustus-Glimmer.txt";
		String officialPeptideFile = "result/Toxo-Peptide-Official.txt";
		
		writePeptidesToFile(alternatePeptideFile, officialPeptideFile);
		
	}
	
	public void perform_classification_Neo()
	{
		String protein_Peptide_Count_File_official = new String("result/Protein-Peptide-Count-official-neo.txt");
		String protein_Peptide_Count_File_augutus = new String("result/Protein-Peptide-Count-augustus-neo.txt");
		String protein_Peptide_Count_File_glimmer = new String("result/Protein-Peptide-Count-glimmer-neo.txt");
		
		String trypticFile = "result/Digested-Peptides-Neo.txt";
		TrypticPeptideReader tr = new TrypticPeptideReader(trypticFile);
		HashSet<String> trypticPeptides = tr.getTrypticPeptides();
		
		classifyPeptides(protein_Peptide_Count_File_official, protein_Peptide_Count_File_augutus, protein_Peptide_Count_File_glimmer,trypticPeptides);
		
		String alternatePeptideFile = "result/Neo-Peptide-Augustus-Glimmer.txt";
		String officialPeptideFile = "result/Neo-Peptide-Official.txt";
		
		writePeptidesToFile(alternatePeptideFile, officialPeptideFile);
		
	}
	
	/**
	 * This is for preparing peptide maps
	 * 
	 * @param summaryFiles
	 */
	public HashMap <String, ArrayList<Peptide>> obtainInformationAboutPeptides(String [] summaryFiles,double fdrThreshold,
			String delimiter,String decoyString){
		
		AnalysePeptidesAcrossExperiments aae = new AnalysePeptidesAcrossExperiments(summaryFiles,fdrThreshold, delimiter,decoyString);
		//aae.detectProteinsAndPeptides();
		
		// create peptide maps
		aae.createCompletePeptideMapForAllExperiments();
		
		return aae.allPeptideMap; 
	}
	
	/**
	 * Perform SETTINGS and create peptide maps
	 */
	public void obtainPeptideForAllModels(){
		double fdrThreshold = 0.01;
		String delimiter = "\t";
		String decoyString = "Rnd";
		
		String [] summaryFiles_official = {"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_SummaryFiles/WholeSummary_1D.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_SummaryFiles/WholeSummary_2D.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_SummaryFiles/WholeSummary_sf.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_official_Results/ToxoOfficial_SummaryFiles/WholeSummary_mudpitOfficial.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Official/Dong/WholeSummary_Dong_orbitrap_official.txt" ,
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Official/Sanya/WholeSummarySanyaOrbitrap_official.txt"};
		
		this.official_peptideMap = obtainInformationAboutPeptides(summaryFiles_official,fdrThreshold, delimiter,decoyString);
		
		String [] summaryFiles_augustus = {"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/WholeSummaryToxoAugustus_mudpit_solAndInsol.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/WholeSummary_SFIF_augustus.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/WholeSummary_2D_augustus.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Augustus/Toxo_Augustus/WholeSummary_1D_augustus.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Augustus/Dong/WholeSummaryToxoOrbitrapAugustus.txt",
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Augustus/Sanya/WholeSummarySanyaAugustus.txt"};
		this.augustus_peptideMap = obtainInformationAboutPeptides(summaryFiles_augustus,fdrThreshold, delimiter,decoyString);
		
		String [] summaryFiles_glimmer = {"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Glimmer/Toxo_Glimmer/WholeSummary_Mudpit_glimmer.txt" ,
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Glimmer/Toxo_Glimmer/WholeSummary_SFIF_glimmer.txt" ,
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Glimmer/Toxo_Glimmer/WholeSummary_2D_glimmer.txt" ,
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Glimmer/Toxo_Glimmer/WholeSummary_1D_glimmer.txt" ,
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Glimmer/Sanya/WholeSumary_sanya_glimmer.txt" ,
				"/Users/riteshk/Ritesh_Work/Dropbox-PipelineResultsPvtCopy/Toxo/Toxo_Orbitrap/Glimmer/Dong/WholeSummary_dong_glimmer.txt"};
		this.glimmer_peptideMap = obtainInformationAboutPeptides(summaryFiles_glimmer,fdrThreshold, delimiter,decoyString);
		
	}
	
	/**
	 * 
	 * @param peptideSet
	 * @param peptideMap
	 * @param outputFile
	 */
	public void obtainLocationsForPeptides(HashSet<String> peptideSet, HashMap <String, ArrayList<Peptide>> peptideMap, 
			String outputFile){
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
			
			Iterator<String> pepSeqs  = peptideSet.iterator();
			while(pepSeqs.hasNext()){
				String pepseq = pepSeqs.next();
				
				if(!peptideMap.containsKey(pepseq))
					continue;
				
				ArrayList<Peptide> pepColl = peptideMap.get(pepseq);
				HashMap<String,Integer> uniqueProtAccn = new HashMap<String,Integer>();
				for(int i = 0 ; i < pepColl.size(); i++){
					uniqueProtAccn.put(pepColl.get(i).protAccn,i);
				}
				
				Iterator<String> protAccnColl = uniqueProtAccn.keySet().iterator();
				while(protAccnColl.hasNext()){
					String protAccn = protAccnColl.next();
					int index = uniqueProtAccn.get(protAccn);
					int pep_start = pepColl.get(index).start;
					int pep_end = pepColl.get(index).end;
					
					out.write(protAccn + "\t" + pepseq + "\t" + pep_start + "\t" + pep_end + "\n");
				}
			}
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * Do all the setting in the perform_classification() 
	 */
	public static void  main(String [] args){
		
		//Classify_peptides_in_official_alternate cp= new Classify_peptides_in_official_alternate();
		
		//System.out.println("Toxo -" + "\n");
		//cp.perform_classification_Toxo();
		
		/* ------------------------------------------------------------------ */
		// -- Only for Toxo where we want to create circos plots ----//
		// create peptide maps
		//cp.obtainPeptideForAllModels();	
		
		// For Augustus
		//String outputFile_augustus = "result/Toxo-Augustus-PeptideMap.txt";
		//cp.obtainLocationsForPeptides(cp.peptides_augustus_only, cp.augustus_peptideMap, outputFile_augustus);
		
		// For Glimmer
		//String outputFile_glimmer = "result/Toxo-Glimmer-PeptideMap.txt";
		//cp.obtainLocationsForPeptides(cp.peptides_glimmer_only, cp.glimmer_peptideMap, outputFile_glimmer);
		/* ------------------------------------------------------------------ */
		
		// For Neo Data - Comment out the above lines and run
		Classify_peptides_in_official_alternate cp_neo= new Classify_peptides_in_official_alternate();
		System.out.println("Neo -" + "\n");
		cp_neo.perform_classification_Neo();
		
	}
	
	
}
