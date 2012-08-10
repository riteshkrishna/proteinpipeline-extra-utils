package uk.ac.liverpool.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Scanner;

/**
 * Written to analyse the N-termnial outputs producded by the pipeline. Input file for analysis are
 * in TSV format with -
 *           [Accn]              [PepSeq]              [start-location]
 *  TGME49_chrXI.aa_g64.t1	    AMKYVAAYLMVVLSGTDAPTKK	     2
	TGME49_chrXII.aa_g290.t1	AASGHPIPELGEFIIANK	         2
	TGME49_chrIV.aa_g129.t1	    MLANKLGIQDVGAQLTGK	         1
 * The example files are located in - 
 * /Users/riteshk/Google Drive/Research/Pipeline_paper/Supp-For-Paper/Toxo/N-Terminal
 * 
 * We need two files, one for official n-terminal, and the other for alternate n-terminals
 * @author riteshk
 *
 */
public class AnalyseNTerminalFiles {
	
	String officialFile;
	String alternateFile;
	
	String delimiter = "\t";
	
	public HashMap<String,String> official_nterminals = new HashMap<String,String>();
	public HashMap<String,String> alternate_nterminals = new HashMap<String,String>();
	
	public AnalyseNTerminalFiles(String officialFile,String alternateFile){
		this.officialFile = officialFile;
		this.alternateFile = alternateFile;
	}
	
	public void makeMaps(){
		
		Scanner scanner_off ;
		Scanner scanner_alt ;
		try{
			// Official map
			scanner_off = new Scanner(new FileReader(new File(this.officialFile)));
			while(scanner_off.hasNextLine()){
				String line = scanner_off.nextLine();
				
				if(line.isEmpty())
					 continue;
				
				 String [] values = line.split(this.delimiter);
				 if(values.length != 3)
					 continue;
				 
				 String accn = values[0].trim();
				 String peptideSeq = values[1].trim();
				 
				 official_nterminals.put(peptideSeq, accn);
			}
			scanner_off.close();
			
			// Alternate map
			scanner_alt = new Scanner(new FileReader(new File(this.alternateFile)));
			while(scanner_alt.hasNextLine()){
				String line = scanner_alt.nextLine();
				
				if(line.isEmpty())
					 continue;
				
				 String [] values = line.split(this.delimiter);
				 if(values.length != 3)
					 continue;
				 
				 String accn = values[0].trim();
				 String peptideSeq = values[1].trim();
				 
				 alternate_nterminals.put(peptideSeq, accn);
			}
			scanner_alt.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	/**
	 * Compare alternate and official maps and remove the sequences from alternate map if they are present in official map
	 */
	public void compareMaps(){
		
		HashMap<String,String> alternate_nterminals_copy = new HashMap<String,String>(alternate_nterminals);
		
		Iterator<String> altSeqs  = alternate_nterminals_copy.keySet().iterator();
		
		while(altSeqs.hasNext()){
			String seq  = altSeqs.next();
			if(this.official_nterminals.containsKey(seq))
				this.alternate_nterminals.remove(seq);
		}
	}
	
	/**
	 * To filter alternate N-terminals from official theoretical tryptic-peptides
	 * @param trypticFile
	 */
	public void filter_alternate_nterminals_with_trypticpeptides(String trypticFile){
		
		TrypticPeptideReader tr = new TrypticPeptideReader(trypticFile);
		HashSet<String> trypticPeptides = tr.getTrypticPeptides();

		HashMap<String,String> alternate_nterminals_copy = new HashMap<String,String>(this.alternate_nterminals);
		Iterator<String> altSeqs  = alternate_nterminals_copy.keySet().iterator();
		
		while(altSeqs.hasNext()){
			String seq  = altSeqs.next();
			if(trypticPeptides.contains(seq))
				this.alternate_nterminals.remove(seq);
		}
	}
	/**
	 * 
	 * @param args
	 */
	public static void main(String [] args){
		
		String officialFile = "tmp/All-official-N-terminal-neo.txt";
		String alternateFile = "tmp/All-Alternate-N-terminal-neo.txt";
		String nTerminal_outputFile = "tmp/N-terminal-analysis-neo.txt";
		String trypticFile = "result/Digested-Peptides-Neo.txt";
		
		/*
		String officialFile = "tmp/All-official-N-terminal-toxo.txt";
		String alternateFile = "tmp/All-Alternate-N-terminal-toxo.txt";
		String nTerminal_outputFile = "tmp/N-terminal-analysis-toxo.txt";
		String trypticFile = "result/Digested-Peptides-Toxo.txt";
		*/
		
		AnalyseNTerminalFiles as = new AnalyseNTerminalFiles(officialFile,alternateFile);
		as.makeMaps();
		as.compareMaps();
		as.filter_alternate_nterminals_with_trypticpeptides(trypticFile);
		
		HashMap<String,String> official_nterminals_final = new HashMap<String,String>(as.official_nterminals);
		HashMap<String,String> alternate_nterminals_final = new HashMap<String,String>(as.alternate_nterminals);
		
		System.out.println("Total official N-terminals = " + official_nterminals_final.size());
		System.out.println("Total alternate N-terminals = " + alternate_nterminals_final.size());
		
		
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(nTerminal_outputFile));
			out.write("\n Official N-terminal candidates \n");
			Iterator<String> offSeqs  = official_nterminals_final.keySet().iterator();
			while(offSeqs.hasNext()){
				String pep = offSeqs.next();
				out.write("\n" + pep + "\t" + official_nterminals_final.get(pep));
			}
			
			out.write("\n Alternate N-terminal candidates \n");
			Iterator<String> altSeqs  = alternate_nterminals_final.keySet().iterator();
			while(altSeqs.hasNext()){
				String pep = altSeqs.next();
				out.write("\n" + pep + "\t" + alternate_nterminals_final.get(pep));
			}
			out.close();
		}catch(Exception e){
			
		}
	}

}
