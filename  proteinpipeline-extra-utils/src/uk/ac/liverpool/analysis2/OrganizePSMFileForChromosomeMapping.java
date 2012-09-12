package uk.ac.liverpool.analysis2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Scanner;

/**
 * This class is written to parse the PSM file produced after running the Post-processing jar
 * on the search results.
 * 
 * We can create a single CSV file by combining identifications from all the experiments, by keeping the 
 * *_processed_PSM.csv files in a single folder and running the following script -
 * ******
 	!/bin/sh
	FILE_DIR="*.csv"
	for FILE_NAME in $FILE_DIR
	do
		cut -d, -f11,15 $FILE_NAME >> all-PSM.txt
	done
 * *******
 * The script will create a single file will all the information we need for creating TSV files for
 * each gene models. Before processing the file, remove the repeated header lines composed of strings
 * "Sequence, proteinacc_start_stop_pre_post_;"
 * 
 * After running the script, a number of files will be created. Each file then needs to be filtered to keep
 * only the unique records - this filtering can  be done in Excel. The filtered files should be ready to be
 * sent to the mapping program.
 * 
 * @author riteshk
 */

public class OrganizePSMFileForChromosomeMapping {
		final String file_delimiter = "\t";
		final String field_delimiter = ";";
		
		final String A_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/postprocess/All-PSM-files/All-PSM_models/all-PSM-official.txt";
		final String B_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/postprocess/All-PSM-files/All-PSM_models/all-PSM-augustus.txt";
		final String C_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/postprocess/All-PSM-files/All-PSM_models/all-PSM-geneMark.txt";
		final String D_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/postprocess/All-PSM-files/All-PSM_models/all-PSM-glimmer.txt";
		//final String E_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/postprocess/All-PSM-files/All-PSM_models/all-PSM-E.txt";
		final String Rev_file = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/postprocess/All-PSM-files/All-PSM_models/all-PSM-Rev.txt";
		
		public void performParsing(String psmFile){
			try{
				BufferedWriter out_A = new BufferedWriter(new FileWriter(A_file));
				BufferedWriter out_B = new BufferedWriter(new FileWriter(B_file));
				BufferedWriter out_C = new BufferedWriter(new FileWriter(C_file));
				BufferedWriter out_D = new BufferedWriter(new FileWriter(D_file));
				//BufferedWriter out_E = new BufferedWriter(new FileWriter(E_file));
				BufferedWriter out_Rev = new BufferedWriter(new FileWriter(Rev_file));
				
				Scanner scanner = new Scanner(new FileReader(new File(psmFile)));
				
				while(scanner.hasNextLine()){
					String line = scanner.nextLine();
					
					String [] pepAndAccn = line.split(this.file_delimiter);
					if(pepAndAccn.length < 2){
						System.out.println(pepAndAccn[0]);
						continue;
					}
					String pepSeq 		 = pepAndAccn[0];
					String accessions 	 = pepAndAccn[1];
					
					String underscore = "_";
					String last_underscore = "_null_null";
					
					if(accessions != null){
						String [] protAccns = accessions.split(this.field_delimiter);
						for(String prot : protAccns){
							String model = prot.substring(0,prot.indexOf(underscore));
							
							int index_last_underscore = prot.indexOf(last_underscore);
							String prot_sub_withend = prot.substring(0,index_last_underscore);
							int end_index = prot_sub_withend.lastIndexOf(underscore);
							
							String prot_sub_withstart = prot_sub_withend.substring(0,end_index);
							int start_index = prot_sub_withstart.lastIndexOf(underscore);
							
							String start = prot.substring(start_index + 1,end_index);
							String end  = prot.substring(end_index + 1,index_last_underscore);
							
							String protein = prot.substring(prot.indexOf(underscore) + 1, start_index);
							
							String protInfo = protein + "\t" + pepSeq + "\t" + start + "\t" + end;
							
							if(model.equals("A")){
								out_A.write(protInfo + "\n");
							}else if(model.equals("B")){
								out_B.write(protInfo + "\n");
							}else if (model.equals("C")){
								out_C.write(protInfo + "\n");
							}else if(model.equals("D")){
								out_D.write(protInfo + "\n");
							}else if(model.equals("E")){
							//	out_E.write(protInfo + "\n");
							}else if(model.equals("Rev")){
								out_Rev.write(protInfo + "\n");
							}
								
							System.out.println("\t" + model + "\t" + protInfo);
						}
					}
					
					//System.in.read();
				}
				scanner.close();
				out_A.close();
				out_B.close();
				out_C.close();
				out_D.close();
				//out_E.close();
				out_Rev.close();
				
			}catch(Exception e){
				e.printStackTrace();
			}
			
		}
		
		/**
		 * 
		 * @param args
		 * @throws Exception
		 */
		public static void main(String [] args) throws Exception {
			
			String inputFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/postprocess/All-PSM-files/all-PSM.txt";
			OrganizePSMFileForChromosomeMapping oc = new OrganizePSMFileForChromosomeMapping();
			oc.performParsing(inputFile);
			
		}
	
}
