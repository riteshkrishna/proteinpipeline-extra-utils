package uk.ac.liverpool.signalP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.Set;

/**
 * We can download a TSV file from ToxoDB.org containing the list of genes with Signal-P predictions.
 * We will need to parse the file and create appropriate FASTA files which can be used for querying
 * for proteomics evidences.
 *  
 *  Files needed for processing -
 *  	- downloaded TSV file with following columns -
 *  		[Gene ID]	[Gene Group]	[Organism]	[Genomic Location]	[Product Description]	[SignalP Scores]	[SignalP Peptide]
 *  	- FASTA sequence file of the genes listed in the Excel file
 * @author riteshk
 *
 */
public class ParseSignalPResultFromToxoDB {

	String signalPInputFile;
	String fastaInputFile;
	
	String delimiter = "\t"; // Assuming this is a TSV file
	
	public ParseSignalPResultFromToxoDB(String signalPFileName, String fastaFile){
		try{
			this.signalPInputFile = new String(signalPFileName);
			this.fastaInputFile = new String(fastaFile);
		}
		catch(NullPointerException e){
			System.out.println("Problem with file paths");
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	
	public void processSignalP_File(String outputFasta_NN, String outputFasta_HMM, double score_threshold){
		Scanner scanner ;
		
		try{
			
			BufferedWriter out_nn = new BufferedWriter(new FileWriter(outputFasta_NN));
			BufferedWriter out_hmm = new BufferedWriter(new FileWriter(outputFasta_HMM));
			
			HashMap<String,String> geneMap_NN = new HashMap<String,String>();
			HashMap<String,String> geneMap_HMM = new HashMap<String,String>();
			
			scanner = new Scanner(new FileReader(new File(signalPInputFile)));
			
			int lineCounter = 0;
			
			String [] header = new String[7];
			
			int index_GeneID = -1;
			int index_SignalPScore = -1;
			int index_SignalPPeptide = -1;
			
			while(scanner.hasNextLine()){
				
				String line = scanner.nextLine();
				 if(line.isEmpty())
					 continue;
				 
				 lineCounter++;
				 String [] values = line.split(this.delimiter);
				 
				 try{
				 if(values.length < 7)
					 throw new Exception(" Less than 7 columns found in the line \n " + line);
				 }catch(Exception e){
					 System.out.println("Exception :: " + e.getMessage());
				 }
				 
				 // Fill the header
				 if(lineCounter == 1){
					 
					 System.arraycopy(values, 0, header, 0, values.length);
					 
					 for(int i = 0 ;i < header.length; i++){
						 if(header[i].contains("[Gene ID]"))
							 index_GeneID = i;
						 if(header[i].contains("[SignalP Scores]"))
							 index_SignalPScore = i;
						 if(header[i].contains("[SignalP Peptide]"))
							 index_SignalPPeptide = i;
					 }
					 
					 try{
					 if(!(index_GeneID != -1) && (index_SignalPPeptide != -1) && (index_SignalPScore != -1))
						 throw new Exception("Header not initialized properly..");
					 }catch(Exception e){
						 System.out.println("Exception :: " + e.getMessage());
						 System.exit(0);
					 }
					 
					 continue;
				 }
				 
				 String gene_id = values[index_GeneID];
				 String signalp_score = values[index_SignalPScore];
				 String signal_peptide = values[index_SignalPPeptide];
				 
				 double nn_score = get_nn_score(signalp_score);
				 double hmm_score = get_hmm_score(signalp_score);
				 String nn_peptide = get_nn_peptide(signal_peptide);
				 String hmm_peptide = get_hmm_peptide(signal_peptide);
				
				 if(nn_score >= score_threshold){
					geneMap_NN.put(gene_id,nn_peptide); 
					System.out.println("Gene ID = " + gene_id + " , NN = " + nn_score + " , NN Pep = " + nn_peptide);
				 }
				 if(hmm_score >= score_threshold){
					 geneMap_HMM.put(gene_id,hmm_peptide);
					 System.out.println("Gene ID = " + gene_id + " , HMM = " + hmm_score + " , HMM Pep = " + hmm_peptide);
				 }
			}
			
			
			// Prepare ArrayList of NN genes
			Set<String> nn_keys = geneMap_NN.keySet();
			LinkedList<String> nn_list = new LinkedList<String>(nn_keys);
			ArrayList<String> genes_nn = new ArrayList<String>(nn_list);
			// Prepare ArrayList of HMM genes
			Set<String> hmm_keys = geneMap_HMM.keySet();
			LinkedList<String> hmm_list = new LinkedList<String>(hmm_keys);
			ArrayList<String> genes_hmm = new ArrayList<String>(hmm_list);
			
			/*
			ArrayList<String> genes_nn = new ArrayList<String>();
			ArrayList<String> genes_hmm = new ArrayList<String>();
			
			Iterator<String> nn_keys = geneMap_NN.keySet().iterator();
			while(nn_keys.hasNext())
				genes_nn.add(nn_keys.next());
			
			Iterator<String> hmm_keys = geneMap_HMM.keySet().iterator();
			while(hmm_keys.hasNext())
				genes_hmm.add(hmm_keys.next());
			*/
			
			// Get Protein Sequences from Fasta File
			FastaSequenceReader fsr = new FastaSequenceReader(this.fastaInputFile);
			HashMap<String, String> collection_NN = fsr.getSequence(genes_nn);
			HashMap<String, String> collection_HMM = fsr.getSequence(genes_hmm);
			
			System.out.println(" Total NN Protein queried = " + genes_nn.size() + " Total HMM Protein queried = " + genes_hmm.size());
			System.out.println(" Total NN Protein found = " + collection_NN.size() + " Total HMM Protein found = " + collection_HMM.size());
			System.in.read();
			
			// Remove signal-peptides from Protein sequences - For NN
			for(String gene : genes_nn){
				String peptide = geneMap_NN.get(gene).trim();
				String protein = collection_NN.get(gene);
				if(protein != null){
					String truncatedProtein = protein.replaceAll(peptide, "");
					out_nn.write(">" + gene);
					out_nn.write("\n" + truncatedProtein + "\n");
				}//else
					//System.out.println("No Protein sequence found for -" + gene);
			}
			
			// Remove signal-peptides from Protein sequences - For HMM
			for(String gene : genes_hmm){
				String peptide = geneMap_HMM.get(gene).trim();
				String protein = collection_HMM.get(gene);
				if(protein != null){
					String truncatedProtein = protein.replaceAll(peptide, "");
					out_hmm.write(">" + gene);
					out_hmm.write("\n" + truncatedProtein + "\n");
				}//else
					//System.out.println("No Protein sequence found for -" + gene);
			}
			
			
			out_nn.close();
			out_hmm.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	/**
	 * Parse the string like - "NN Sum: 1, NN D: .39, HMM Prob: .55"
	 * and extract NN D score, here 0.39
	 * @param signalp_score
	 * @return
	 */
	double get_nn_score(String signalp_score){
		double nn_score = 0;
		
		signalp_score = signalp_score.replaceAll("\"", "");
		String [] scores = signalp_score.split(",");
		if(scores[1] != null){
			String [] NN_D_score = scores[1].split(":");
			if(NN_D_score[1] != null)
				nn_score = Double.parseDouble(NN_D_score[1]);
			else
				System.out.println("No NN-D Score found for - " + signalp_score);
		}
		return nn_score;
	}
	

	/**
	 * Parse the string like - "NN Sum: 1, NN D: .39, HMM Prob: .55"
	 * and extract HMM Prob score, here 0.55
	 * @param signalp_score
	 * @return
	 */
	double get_hmm_score(String signalp_score){
		double hmm_score = 0;
		
		signalp_score = signalp_score.replaceAll("\"", "");
		String [] scores = signalp_score.split(",");
		
		if(scores[2] != null){
			String [] HMM_score = scores[2].split(":");
			if(HMM_score[1] != null)
				hmm_score = Double.parseDouble(HMM_score[1]);
			else
				System.out.println("No HMM Score found for - " + signalp_score);
		}
		return hmm_score;
	}
	
	/**
	 * Extract NN peptide from an entry like this -
	 * 	"HMM: MFLLLDPVSPRHSLVSPSSFRNSLFLLSPLFPLLSCLSPLSSV, NN: MFLLLDPVSPRHSLVSPSSFRNSLFLLSPLFPLLSCL"
	 * @param signal_peptide
	 * @return
	 */
	String get_nn_peptide(String signal_peptide){
		String nn_peptide = new String();
		
		signal_peptide = signal_peptide.replaceAll("\"", "");
		
		String []peptides = signal_peptide.split(",");
		
		if(peptides[1] != null){
			String [] nn_peps = peptides[1].split(":");
			if(nn_peps[1] != null)
				nn_peptide = nn_peps[1];
			else
				System.out.println("No NN Peptide found for -" + signal_peptide);
		}
		
		return nn_peptide;
	}
	
	/**
	 * Extract HMM peptide from an entry like this -
	 * 	"HMM: MFLLLDPVSPRHSLVSPSSFRNSLFLLSPLFPLLSCLSPLSSV, NN: MFLLLDPVSPRHSLVSPSSFRNSLFLLSPLFPLLSCL"
	 * @param signal_peptide
	 * @return
	 */
	
	String get_hmm_peptide(String signal_peptide){
		String hmm_peptide = new String();
		
		signal_peptide = signal_peptide.replaceAll("\"", "");
		
		String []peptides = signal_peptide.split(",");
		
		if(peptides[0] != null){
			String [] hmm_peps = peptides[0].split(":");
			if(hmm_peps[1] != null)
				hmm_peptide = hmm_peps[1];
			else
				System.out.println("No NN Peptide found for -" + signal_peptide);
		}
		return hmm_peptide;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String signalPFileName = args[0];
		String fastaFile = args[1];
		String outputFasta_NN = args[2]; 
		String outputFasta_HMM = args[3];
		double score_threshold = Double.parseDouble(args[4]);
		
		ParseSignalPResultFromToxoDB ptb = new ParseSignalPResultFromToxoDB(signalPFileName,fastaFile);
		ptb.processSignalP_File(outputFasta_NN, outputFasta_HMM,score_threshold);
	}

}
