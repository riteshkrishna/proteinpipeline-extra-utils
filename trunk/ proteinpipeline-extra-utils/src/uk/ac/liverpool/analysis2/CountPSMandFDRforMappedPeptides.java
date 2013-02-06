package uk.ac.liverpool.analysis2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

/**
 * This has come as a final requirement for uploading the data on ToxoDB browser. We have the
 * peptides mapped on Genome using the MapPeptidesOnToxoVersion8 and PeptideMappingOnChromosome,
 * but we want to include some additional information like FDR score and PSM count.
 * 
 * So, I have again requeried Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/All-WholeSummary.txt
 * and have extracted columns having specID, peptide-seq and FDR scores using the cut command.
 * The file was further filtered for the peptides with <= 0.01 FDR score, and the file produced is -
 * /Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/Version_8-mapping/All-WholeSummary-Column3-4-5Only-filteredOn_1Percent.txt
 * 
 * In the working directory -/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/Version_8-mapping/
 * 
 * We can scan through All-peptide-mappedOnGenome.txt and for each peptide, we can count the PSMs and
 * the best FDR score from  All-WholeSummary-Column3-4-5Only-filteredOn_1Percent.txt. We can just append 
 * these numbers next to the values in  All-peptide-mappedOnGenome.txt and that should serve the purpose.
 * 
 * @author riteshk
 *
 */
public class CountPSMandFDRforMappedPeptides {
	
	
	public static void main(String [] args) throws Exception{
		
		String genomeFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/Version_8-mapping/All-peptide-mappedOnGenome.txt";
		String summaryFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/Version_8-mapping/All-WholeSummary-Column3-4-5Only-filteredOn_1Percent.txt";
		String outputFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/Version_8-mapping/All-peptide-mappedOnGenome-PSM-FDR.txt";
		
		
		//String summaryFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/All-WholeSummary.txt";
		
		double fdr_threshold = 0.01;
		
		HashMap<String, ArrayList<String>> map_peptide_specList = new HashMap<String, ArrayList<String>>();
		HashMap<String, Double> map_peptide_bestFDR = new HashMap<String, Double>();
		HashMap<String, Integer> map_peptide_specCount = new HashMap<String, Integer>();
		
		try{
			// Create the Peptide, FDR, Spec_count map
			Scanner scanner_summary = new Scanner(new File(summaryFile));
			while(scanner_summary.hasNextLine()){
				String line = scanner_summary.nextLine();
				
				String [] values = line.split("\t");
				
				//System.out.println("Size : " + values.length +" Line : " + line);
				
				String specID = values[0];
				String peptide = values[1];
				double fdr = Double.parseDouble(values[2]);
		
				
				if(fdr <= fdr_threshold){
					
					if(map_peptide_specList.containsKey(peptide)){
						ArrayList<String> specs = map_peptide_specList.get(peptide);
						if(!specs.contains(specID)){
							specs.add(specID);
							map_peptide_specList.put(peptide,specs);
						}
					}else{
						ArrayList<String> specs = new ArrayList<String>();
						specs.add(specID);
						map_peptide_specList.put(peptide,specs);
					}
					
					if(map_peptide_bestFDR.containsKey(peptide)){
						double best_fdr = map_peptide_bestFDR.get(peptide);
						if(best_fdr > fdr)
							map_peptide_bestFDR.put(peptide, fdr);
					}else{
						map_peptide_bestFDR.put(peptide, fdr);
					}	
				}
				
			}
			// Count the specID for eac peptide
			Iterator<String> peptides = map_peptide_specList.keySet().iterator();
			while(peptides.hasNext()){
				String peptide = peptides.next();
				map_peptide_specCount.put(peptide, map_peptide_specList.get(peptide).size());
			}
			
			scanner_summary.close();
			
			BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
			
			// Process the genome file
			Scanner scanner = new Scanner(new File(genomeFile));
			while(scanner.hasNextLine()){
				String line = scanner.nextLine();
				
				String [] values = line.split("\t", -1);
				String peptide = values[3];
				
				try{
					double best_fdr = map_peptide_bestFDR.get(peptide);
					int total_spec = map_peptide_specCount.get(peptide);
					
					out.write(line + "\t" + best_fdr + "\t" + total_spec + "\n");
					
				}catch(Exception e){
					System.out.println("Not found for : " + peptide);
					e.printStackTrace();
				}
			}
			
			out.close();
			scanner.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
