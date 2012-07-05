package uk.ac.liverpool.signalP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Vector;

/**

I have summarised results from all the experiments in two files with filter set to 1% threshold -
	- All-semi-argc-results-threshold-01.t txt (for Semi-Arg-C)
	- All-full-argc-results-threshold-01.txt     (for full Arg-C)
	
 * The files need to be processed so we can identify the official n-termials, signal-peptides and alternate-n-terminals.
 * The following code does the above and create an output file with the above classes. 
	


 * @author riteshk
 *
 */
public class Querying_Gian_argc_searches {
	
	String decoy_identifier = "Rev";
	
	String gian_search_file;
	
	
	public class Groups{
		public String specId;
		public String pepseq;
		public String accn;
		public int start;
	};
	
	HashMap<String, ArrayList<Groups>> records;
	
	public HashMap<String, ArrayList<Groups>> official_nTermnal;
	public HashMap<String, ArrayList<Groups>> alternate_nTermnal;
	public HashMap<String, ArrayList<Groups>> signalPeptides;
	
	public Querying_Gian_argc_searches(String gian_search_file){
		this.gian_search_file = new String(gian_search_file);
		this.records = new HashMap<String,ArrayList<Groups>>();
		
		this.official_nTermnal = new HashMap<String,ArrayList<Groups>>();
		this.alternate_nTermnal = new HashMap<String,ArrayList<Groups>>();
		this.signalPeptides = new HashMap<String,ArrayList<Groups>>();
	}
	
	
	public void processSearchFile() throws Exception{
		
		Scanner scanner = new Scanner(new FileReader(new File(this.gian_search_file)));
		
		while(scanner.hasNextLine()){
			
			String line = scanner.nextLine();
			 if(line.isEmpty())
				 continue;
			 
			 String [] values = line.split("\t");
			 
			 try{
			 if(values.length < 13)
				 throw new Exception(" Less than 13 columns found in the line \n " + line);
			 }catch(Exception e){
				 System.out.println("Exception :: " + e.getMessage());
			 }
			 
			 String specId = values[1];
			 String pepSeq = values[2];
			 String accn = values[6];
			 int start = Integer.parseInt(values[7]);
			 
			 // Ignore the Decoy hits
			 if(accn.contains(decoy_identifier))
				 continue;
			 
			 String key = specId.trim() + pepSeq.trim();
			 
			 Groups g = new Groups();
			 g.pepseq = pepSeq;
			 g.specId = specId;
			 g.accn   = accn;
			 g.start  = start;
			 
			 if(!records.containsKey(key)){
				 ArrayList<Groups> g_g = new ArrayList<Groups>();
				 g_g.add(g);
				 records.put(key, g_g);
			 }else{
				 ArrayList<Groups> g_g = records.get(key);
				 g_g.add(g);
				 records.put(key, g_g);
			 }
				 
		}
		
		scanner.close();
		
		// After we have read the relevant information from the file, we do the processing...
		classify();
	}
		
	public void classify(){
		// The protein accessions in the file have the following strings in the beginning, except in SP which is at the end.
		String official_identifier = "gb|";
		String augustus_identifier = "augustus";
		String glimmer_identifier = "Glimmer";
		String signalp_identifier = ".SP";
		
		Iterator<String> keys = records.keySet().iterator();
		while(keys.hasNext()){
			String key = keys.next();
			
			Vector<Integer> start_locs = new Vector<Integer>();
			Vector<String> accessions = new Vector<String>();
			
			ArrayList<Groups> groupList = records.get(key);
			for(Groups g : groupList){
				start_locs.add(g.start);
				accessions.add(g.accn);
			}
			

			// Process only if any of the start locations is 1 or 2
			if((start_locs.contains(1) || start_locs.contains(2)) == false){
				continue;
			}else{
				
				// find the accessions with 1 or 2
				Vector<String> interestingAccessions = new Vector();
				Vector<Integer> interseting_indices = new Vector();
				
				for(int i = 0; i < start_locs.size(); i++){
					if((start_locs.get(i) == 1) || (start_locs.get(i) == 2)){
						String accn = accessions.get(i);
						interestingAccessions.add(accn);
						interseting_indices.add(i);
					}
				}
				
				boolean foundOfficial = false;
				for(String thisaccn : interestingAccessions){
					if(thisaccn.startsWith(official_identifier)){
						
						boolean sigP = false;
						for(String searchSigp : interestingAccessions){
							if(searchSigp.contains(signalp_identifier)){
								sigP = true;
								break;
							}
						}
						if(sigP == true)
							// We are keeping only one copy from multiple protein-peptide combo, hash collision not handled
							signalPeptides.put(thisaccn, groupList);  
						else
							official_nTermnal.put(thisaccn, groupList);
						
						foundOfficial = true;
					}	
				}
				
				// if no official accn was found in the list, then it's alternate or sigP only
				if(foundOfficial == false){
					for(String alternateaccn : interestingAccessions){
						// If alternate
						if(alternateaccn.startsWith(augustus_identifier) || alternateaccn.startsWith(glimmer_identifier)){
							alternate_nTermnal.put(alternateaccn, groupList);
							break;
						}
						// if signal-P 
						if(alternateaccn.contains(signalp_identifier)){
							signalPeptides.put(alternateaccn, groupList);
							break;
						}
					}
				}
				
			}
		} // end of while
		
	}
	
	
	
	public void print_hash(HashMap<String, ArrayList<Groups>> results, String name, BufferedWriter writer){
		try{
		writer.write(" Printing " + name + "\n" );
		Iterator<String> accns = results.keySet().iterator();
		while(accns.hasNext()){
			String accn = accns.next();
			ArrayList<Groups> gList = results.get(accn);
			for(Groups g : gList){
				writer.write(accn + "\t"+ g.accn +"\t" + g.pepseq + "\t" + g.start + "\n");
			}
			writer.write("\n");
		}
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	/*
	 * 
	 */
	public static void main(String [] args) throws Exception{
		String gian_search_file = "tmp/All-full-argc-results-threshold-01.txt";
		Querying_Gian_argc_searches qc = new Querying_Gian_argc_searches(gian_search_file);
		qc.processSearchFile();
		
		 BufferedWriter writer = new BufferedWriter(new FileWriter("tmp/gian_analysis.txt"));
		 
		qc.print_hash(qc.official_nTermnal,"Official N-Terminals", writer);
		qc.print_hash(qc.alternate_nTermnal, "Alternate N-Terminals",writer);
		qc.print_hash(qc.signalPeptides, "Signal-Peptides",writer);
		
		writer.close();
		 
	}

}
