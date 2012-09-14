package uk.ac.liverpool.analysis2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Scanner;

import uk.ac.liverpool.gff.GFF3Routines;
import uk.ac.liverpool.gff.Gene_Information;

/**
 * This file takes a list of unique protein accessions identified during post-processing
 * in a single column text file. It also takes the corresponding GFF file. It will create
 * a TSV file with 4 columns with genomic locations listed for the proteins provided -
 * File format produced - 
 * [Chr name] [start] [end] [Protein Accn]
 * 
 * This file can be directly ported to circos for 'text Track' plotting.
 * @author riteshk
 *
 */
public class CreateProteinTrackFileForCircos {
	
	/**
	 * 
	 * @param queryProteins
	 * @param gffFile
	 * @param circosFile
	 */
	public void createCircosFile(String queryProteins, String gffFile, String circosFile){
		try{
			
			BufferedWriter out = new BufferedWriter(new FileWriter(circosFile));
			
			// Get the info about all the genes
			GFF3Routines gffHandle = new GFF3Routines();
			gffHandle.processGff(gffFile);
			HashMap<String, Gene_Information> geneRecords = new HashMap<String, Gene_Information>(gffHandle.geneRecords);
			
			
			Scanner scanner = new Scanner(new FileReader(new File(queryProteins)));
			while(scanner.hasNextLine()){
				String line = scanner.nextLine();
				if(line.startsWith("#"))
					continue;
				
				// Adjustment in the name due to different scheme used in GFF
				String proteinToQueryToSearch = line.replace("gb", "apidb");
				
				if(geneRecords.containsKey(proteinToQueryToSearch)){
					Gene_Information gene = geneRecords.get(proteinToQueryToSearch);
					String chr = gene.getSeqID();
					long start = gene.getStart();
					long end = gene.getEnd();
					
					// Naming convention for circos
					chr = chr.substring(chr.indexOf("_")+1, chr.length());
					out.write(chr + "\t" + start + "\t" + end +  "\t" + line +"\n");
				}
			}
			
			scanner.close();
			out.close();
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
		String queryProteins = "result2/identified-official-proteins.txt";
		String gffFile = "result2/Toxo_1D_OFF-modified.gff";
		String circosFile = "result2/circos-officialprotein-track.txt";
		
		CreateProteinTrackFileForCircos cs = new CreateProteinTrackFileForCircos();
		cs.createCircosFile(queryProteins, gffFile, circosFile);
	}

}
