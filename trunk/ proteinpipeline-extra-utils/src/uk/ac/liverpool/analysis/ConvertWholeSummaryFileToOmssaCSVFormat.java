package uk.ac.liverpool.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

/**
 * This program is written for converting WholeSummary_ pipeline output files to a OMSSA CVS output format.
 * The output file will go through other routine for protein inference.
 *
 *   The OMSSA headers are -
 *   0. Spectrum number	 
 *   1. Filename/id	 
 *   2. Peptide	 
 *   3. E-value	 
 *   4. Mass	 
 *   5. gi	 
 *   6. Accession	 
 *   7. Start	 
 *   8. Stop	 
 *   9. Defline	 
 *   10. Mods	 
 *   11. Charge	 
 *   12. Theo Mass	 
 *   13. P-value	 
 *   14. NIST score
 * @author riteshk
 *
 */

public class ConvertWholeSummaryFileToOmssaCSVFormat {
	
	// Set according to the columns in the summary file
	final int protAccn_column 		= 0;
	final int protScore_column 		= 1;
	final int specId_column 		= 2;
	final int pepSeq_column 		= 3;
	final int fdr_column 			= 4;
	final int exp_mass_column		= 5;
	final int group_id_column		= 6;
	final int start_column 			= 7;
	final int end_column 			= 8;
	final int simple_fdr_column	    = 9; 
	final int mod_column			= 10;
	final int charge_column 		= 11;	
	final int theo_mass_column      = 12;
	final int searchEngine_column 	= 13;
	
	String wholeSummaryFile;
	double fdrThreshold;
	String delimiter;
	String decoyString; 
	String omssaLikeOutputFile;
	
	/**
	 * 
	 * @param pipelineSummaryFile - Summary file for whole dataset
	 * @param fdrThreshold 		  - FDR to be used for analysis
	 * @param delimiter			  - \t or , or \w
	 * @param decoyString		  - Decoy identifier string	
	 * @param omssaLikeOutputFile - Omssa like output file 
	 */
	public ConvertWholeSummaryFileToOmssaCSVFormat(String pipelineSummaryFile, double fdrThreshold, String delimiter,String decoyString, String omssaLikeOutputFile) {
		wholeSummaryFile = new String(pipelineSummaryFile);
		this.fdrThreshold = fdrThreshold;
		this.delimiter = delimiter;
		this.decoyString = decoyString;
		this.omssaLikeOutputFile = new String(omssaLikeOutputFile);
	}
	
	/**
	 * 
	 */
	public void createOmssaLikeFile(){
		
		Scanner scanner ;
		
		try{
			BufferedWriter out_prot = new BufferedWriter(new FileWriter(this.omssaLikeOutputFile));
			
			scanner = new Scanner(new FileReader(new File(wholeSummaryFile)));
			
			String prevProtein = new String();
			double prevProteinScore = 0.0;
			
			int counter = 1;
			while(scanner.hasNextLine()){
				 String line = scanner.nextLine();
				 if(line.isEmpty())
					 continue;
				 String [] values = line.split(this.delimiter);
				 
				 try{
				 if(values.length < 14)
					 throw new Exception(" Less than 14 columns found in the line \n " + line);
				 }catch(Exception e){
					 System.out.println("Exception :: " + e.getMessage());
				 }
				 
				 String protAccn;
				 if(values[protAccn_column].trim().isEmpty())
					 protAccn = prevProtein;
				 else protAccn = values[protAccn_column].trim();
				 
				 double protScore;
				 if(values[protScore_column].trim().isEmpty())
					 protScore = prevProteinScore;
				 else protScore = Double.parseDouble(values[protScore_column].trim());
				 
				 String specID = values[specId_column].trim();
				 String pepSeq = values[pepSeq_column].trim();
				 int start = Integer.parseInt(values[start_column].trim());
				 int end = Integer.parseInt(values[end_column].trim());
				 double fdrScore = Double.parseDouble(values[fdr_column].trim());
				 
				 double expMass = Double.parseDouble(values[exp_mass_column].trim());
				 String groupID = values[group_id_column].trim();
				 double theoMass = Double.parseDouble(values[theo_mass_column].trim());
				 int charge = Integer.parseInt(values[charge_column].trim());
				 String mods = values[mod_column].trim();
				 double simple_fdr = Double.parseDouble(values[simple_fdr_column].trim());
				 
				 
				 
				 // Skip if FDR score is greater than threshold
				 //if (fdrScore >= this.fdrThreshold)
				//	 continue;
				 // Skip if protein is a decoy
				 //if(protAccn.contains(this.decoyString))
				//	 continue;
				 
				 
				 // else write in omssa format...
				 int gi = 0;
				 int nistScore = 0;
				 String omssaString = counter + "," + specID + "," + pepSeq + "," + fdrScore + "," 
				 					  + expMass + "," + gi + "," + protAccn + "," + start + "," 
				 					  + end + "," + protAccn + "," + mods + "," + charge + "," 
				 					  + theoMass + "," + simple_fdr + "," + nistScore + "\n";
				 
				 out_prot.write(omssaString);
				 
				 counter++  ;
				 prevProtein = protAccn;
				 prevProteinScore = protScore;
			}	
			
			scanner.close();
			
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
		
		String pipelineSummaryFile = "tmp/WholeSummary_pH3-10.txt"; 
		double fdrThreshold = 0.05; 
		String delimiter = "\t";
		String decoyString = "Rnd";
		String omssaLikeOutputFile = "tmp/omssaLike-pH3-10.txt";
		
		ConvertWholeSummaryFileToOmssaCSVFormat cv = new ConvertWholeSummaryFileToOmssaCSVFormat(pipelineSummaryFile, fdrThreshold, delimiter, decoyString, omssaLikeOutputFile);
		cv.createOmssaLikeFile();
		
	}

}
