package uk.ac.liverpool.gff;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Scanner;

/**
 * Program to map peptides on chromosome co-ordinates. We can use this file for circos drawings for
 * peptides identified from alternate models, and N-terminal discoveries.
 * 
 * Take a TSV file with fields [prot-accession peptide-seq start end], produced by 
 * Classify_peptides_in_official_alternate.java 
 * 
 * Take the corresponding GFF3 file.
 * Map and produce a file with the fields -
 * [chromosome start(on chromosome) end(on chromosome) peptide-seq prot-accession]
 * 
 * @author riteshk
 *
 */
public class PeptideMappingOnChromosome {

	String gff;
	String peptideFile;
	String delimiter;
	
	public PeptideMappingOnChromosome(String gffFile, String peptideFile, String delimiter) {
		this.gff = new String(gffFile);
		this.peptideFile = new String(peptideFile);
		this.delimiter = delimiter;
	}

	 
	 public void performMapping(String outputFile) throws Exception{
		 
		 GFFRoutinesForPeptideMapping gffHandle  = new GFFRoutinesForPeptideMapping(this.gff);
		 gffHandle.processGffFile();
		 
		 BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		   
		 Scanner scanner ;
		 try{
				scanner = new Scanner(new FileReader(new File(peptideFile)));
				
				while(scanner.hasNextLine()){
					 String line = scanner.nextLine();
					 if(line.isEmpty())
						 continue;
					 String [] values = line.split(this.delimiter);
					 
					 try{
						 if(values.length < 4)
							 throw new Exception(" Less than 4 columns found in the line \n " + line);
					  }catch(Exception e){
							 System.out.println("Exception :: " + e.getMessage());
					  }
					  
					  String accession = values[0];
					  String peptideSequence = values[1];
					  long start = Long.parseLong(values[2]);
					  long end  = Long.parseLong(values[3]);
					  
					  ProteinPeptideObject p = new ProteinPeptideObject(accession, peptideSequence,start, end);
					  String[]  co_ords = new String[2];
					  long start_map,end_map;
					  
					  try{
						  co_ords = gffHandle.mapToGff(p);
					  
						  // Take care of the negative strand reporting in GFF
						  start_map = Long.parseLong(co_ords[0]);
						  end_map  =  Long.parseLong(co_ords[1]);
						  if(end_map < start_map){
							  long tmp = start_map;
							  start_map = end_map;
							  end_map = tmp;
					  		}
					  }catch(Exception e){
						  System.out.println("Exception - " + accession +" \t"+ peptideSequence +" \t"+ start +" \t"+ end);
						  e.printStackTrace();
						  //System.exit(0);
						  continue;
					  }
					  out.write(co_ords[2] + "\t" + start_map + "\t" + end_map + "\t" + peptideSequence + "\t" + accession + "\n");
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
		 
		 /*
	     String gffFile = "tmp/new-Glimmer-ME49.gff"; 
	     String peptideFile = "tmp/all-PSM-glimmer-filtered.txt"; 
	     String delimiter = "\t";
	     String outFile = "tmp/all-PSM-glimmer-filtered-mapped.txt";
	     */
		 /*
	     String gffFile = "tmp/Augustus_1D.gff"; 
	     String peptideFile = "tmp/all-PSM-augustus-filtered.txt"; 
	     String delimiter = "\t";
	     String outFile = "tmp/all-PSM-augustus-filtered-mapped.txt";
	     */
		 
	     String gffFile = "tmp/Toxo_1D_OFF-modified.gff"; 
	     String peptideFile = "result2/all-PSM-official-filtered.txt"; 
	     String delimiter = "\t";
	     String outFile = "result2/all-PSM-official-filtered-mapped.txt";
	     
		 
		 PeptideMappingOnChromosome pc = new PeptideMappingOnChromosome(gffFile, peptideFile, delimiter);
		 pc.performMapping(outFile);
		 
	 }
	 
}
