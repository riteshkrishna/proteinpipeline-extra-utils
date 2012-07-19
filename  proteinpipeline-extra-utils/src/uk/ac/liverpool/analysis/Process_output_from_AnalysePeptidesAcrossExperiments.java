package uk.ac.liverpool.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Scanner;
import java.util.StringTokenizer;

/**
 * This program analyses output file produced by AnalysePeptidesAcrossExperiments.java. It creates
 * a HashSet of total peptides identified, so we can count that how many are official peptides and
 * how many are alternates.
 * 
 * @author riteshk
 *
 */
public class Process_output_from_AnalysePeptidesAcrossExperiments {

	String inputFile;
	HashSet<String> peptideList;
	
	String delimiter = "\t"; // Protein-Peptide-Count files are delimited by "\t"
	
	/**
	 * 
	 * @param protein_Peptide_Count_File
	 */
	public Process_output_from_AnalysePeptidesAcrossExperiments(String protein_Peptide_Count_File) {
		
		inputFile = protein_Peptide_Count_File;
		peptideList= new HashSet<String>();
		
	}
	
	
	void createPeptideList(){
		Scanner scanner ;
		
		try{
			scanner = new Scanner(new FileReader(new File(inputFile)));
			
			while(scanner.hasNextLine()){
				String line = scanner.nextLine();
				
				if(line.isEmpty())
					 continue;
				
				 String [] values = line.split(this.delimiter);
				 if(values.length != 3)
					 continue;
				 
				 String peptide_for_thisprotein = values[2];
				 peptide_for_thisprotein = peptide_for_thisprotein.replace("[", "");
				 peptide_for_thisprotein = peptide_for_thisprotein.replace("]", "");
				 
				 StringTokenizer stok = new StringTokenizer(peptide_for_thisprotein,",");
				 while(stok.hasMoreTokens()){
					 this.peptideList.add(stok.nextToken().trim());
				 }
			}
			scanner.close();
		}catch(Exception e){
			
		}
	}
	
	/**
	 * 
	 * @return
	 */
	public HashSet<String> getPeptideList(){
		createPeptideList();
		return this.peptideList;
	}
	
	
}
