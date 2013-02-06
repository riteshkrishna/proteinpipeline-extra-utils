package uk.ac.liverpool.analysis2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Scanner;

/* We have cut two columns (peptide and FDR score) from All-WholeSummary.txt file. We want to 
 * filter peptides which passed the FDR threshold. This doesn't remove the duplicate peptides.
 * Use uniq on command line to find the unique list.
 * 
 * Command used to obtain the two column file -
 * cut -f4,5 All-WholeSummary.txt  >> All-peptide-FDR.txt
*/
public class FilterPeptideOnFDRScore {
	
	public static void main (String [] args) throws Exception{
		
		String peptideFDRColumnFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/All-peptide-FDR.txt";
		double fdrThreshold = 0.01;
		String outFile = "/Users/riteshk/Ritesh_Work/ProteoAnnotator-ToxoResults/SearchResults/All-peptide-FDR-filtered.txt";
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
		
		Scanner scanner = new Scanner(new File(peptideFDRColumnFile));
		
		while(scanner.hasNextLine()){
			String line = scanner.nextLine();
			String [] parts = line.split("\t");
			
			String peptide = parts[0];
			double fdr = Double.parseDouble(parts[1]);
			
			if(fdr <= fdrThreshold){
				out.write(peptide + "\t" + fdr + "\n");
			}
		}
		
		scanner.close();
		out.close();
	}
}
