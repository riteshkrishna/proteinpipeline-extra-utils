package uk.ac.liverpool.analysis2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Scanner;

/**
 * The Toxo data was searched against a concatenated database of different models. We need to
 * first combine all the WholeSummary files and then split the records in the merged file into
 * different files according to the gene models they belong to. This way, each produced file can
 * be mapped to the corresponding GFF using the ProteoAnnotator. 
 * 
 * The concatenated file is available at - /Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/All-WholeSummary.txt
 *
 * @author riteshk
 *
 */
public class SplitWholeSummaryFileForModels {
	
	// Set according to the columns in the summary file
	final int protAccn_column 		= 0;

	final String A_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/all-Summary-official.txt";
	final String B_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/all-Summary-augustus.txt";
	final String C_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/all-Summary-geneMark.txt";
	final String D_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/all-Summary-glimmer.txt";
	final String E_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/all-Summary-ORF.txt";
	final String Rev_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/all-Summary-Rev.txt";
	
	final String delimiter = "\t";
	String summaryFile;
	
	public SplitWholeSummaryFileForModels(String summaryFile){
		this.summaryFile = new String(summaryFile);
	}
	
	public void performParsing(){
		try{
			BufferedWriter out_A = new BufferedWriter(new FileWriter(A_file));
			BufferedWriter out_B = new BufferedWriter(new FileWriter(B_file));
			BufferedWriter out_C = new BufferedWriter(new FileWriter(C_file));
			BufferedWriter out_D = new BufferedWriter(new FileWriter(D_file));
			BufferedWriter out_E = new BufferedWriter(new FileWriter(E_file));
			BufferedWriter out_Rev = new BufferedWriter(new FileWriter(Rev_file));
			
			String prevProtein = new String();
			Scanner scanner = new Scanner(new FileReader(new File(this.summaryFile)));
			
			while(scanner.hasNextLine()){
				String line = scanner.nextLine();
				if(line.isEmpty())
					 continue;
				 String [] values = line.split(this.delimiter);
				 
				 String protAccn;
				 if(values[protAccn_column].trim().isEmpty())
					 protAccn = prevProtein;
				 else protAccn = values[protAccn_column].trim();
				
				 String underscore = "_";
				 String model = protAccn.substring(0,protAccn.indexOf(underscore));
				 
				 line = line.replace(model+"_", "");
				 
				 if(model.equals("A")){
						out_A.write(line + "\n");
					}else if(model.equals("B")){
						out_B.write(line + "\n");
					}else if (model.equals("C")){
						out_C.write(line + "\n");
					}else if(model.equals("D")){
						out_D.write(line + "\n");
					}else if(model.equals("E")){
						out_E.write(line + "\n");
					}else if(model.equals("Rev")){
						out_Rev.write(line + "\n");
					}
				 
				 prevProtein = protAccn;
			}
			
			scanner.close();
			out_A.close();
			out_B.close();
			out_C.close();
			out_D.close();
			out_E.close();
			out_Rev.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 * 
	 * @param args
	 */
	public static void main(String [] args){
		String inputSummaryFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/All-WholeSummary.txt";
		SplitWholeSummaryFileForModels sp = new SplitWholeSummaryFileForModels(inputSummaryFile);
		sp.performParsing();
	}
}
