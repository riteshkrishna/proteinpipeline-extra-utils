package uk.ac.liverpool.analysis;
/**
 * Take a file with each line as a peptide string and returns the contents in a HashSet
 */
import java.io.File;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Scanner;

public class TrypticPeptideReader {
	
	String trypticFile;
	public HashSet<String> trypticPeptides;
	
	public TrypticPeptideReader(String peptideFile){
		this.trypticFile = peptideFile;
		trypticPeptides = new HashSet<String>();
	}
	
	void createHash(){
		Scanner scanner;
		try{
			scanner = new Scanner(new FileReader(new File(this.trypticFile)));
			
			while(scanner.hasNextLine()){
				String line = scanner.nextLine();	
				if(line.isEmpty())
					 continue;
				trypticPeptides.add(line.trim());
			}
			scanner.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 * @return
	 */
	public HashSet<String> getTrypticPeptides(){
		createHash();
		return trypticPeptides;
	}
	
	
	/**
	 * Test function
	 * @param args
	 */
	public static void main(String [] args){
		String trypticFile = "tmp/Digested-Peptides-Toxo.txt";
		TrypticPeptideReader tr = new TrypticPeptideReader(trypticFile);
		System.out.println("Total peptides in this file  = " + tr.getTrypticPeptides().size());
	}
}
