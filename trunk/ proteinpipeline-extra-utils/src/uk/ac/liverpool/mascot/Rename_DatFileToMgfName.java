package uk.ac.liverpool.mascot;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

/**
 * Mascot has it's own way of naming the result files produced after the search.
 * The output files are named as .dat (e.g.- F010723.dat), which may be produced
 * from a .mgf (e.g. - Cryptosporidium_whole_cell_lysate_1D_gel_slices_Cp1D63_merge.mgf).
 * 
 *  To make this correspondence apparent, we need to prase the .dat file and retrieve the 
 *  line with FILE attribute , e.g. -
 *  "FILE=C:\Cryptosporidium_whole_cell_lysate_1D_gel_slices_Cp1D63_merge.mgf"
 *  and rename the .dat accordingly, so we know that which output file corresponds to which
 *  input file. 
 *  
 *  @author Ritesh Krishna
 */
public class Rename_DatFileToMgfName {
	
	HashMap<String, String> indexer;
	
	/**
	 * 
	 * @param mainDirectory
	 * @return
	 * @throws Exception
	 */
	File [] getNamesOfDatFiles(String mainDirectory) throws Exception{
        
        class DatFilter implements FilenameFilter {
            public boolean accept(File dir, String name) {
                return (name.endsWith(".dat"));
            }
        }
        
        FilenameFilter fileFilter = new DatFilter();
        File dir = new File(mainDirectory);
        File[] files = dir.listFiles(fileFilter);
        
        return files;
	} 
	
	/**
	 * Read the dat files one by one and find the corresponding mgf names.
	 * Set this information in the indexer HashMap.
	 * 
	 * @param datfiles
	 */
	void findMgfNamesForDatFiles(File[] datfiles) throws Exception{
		
		for(int i = 0; i < datfiles.length; i++){
			
			 Scanner scanner = new Scanner(new FileReader(datfiles[i]));
             
			 try{
				 while(scanner.hasNextLine()){
					 String line = scanner.nextLine();
					 if(line.contains("FILE=")){
						 String[] keyVal = line.split("=");
						 if(keyVal[0].contains("FILE"))
							 indexer.put(datfiles[i].getName(), keyVal[1]);
						 break;
					 }
				 }
             }finally{
                     scanner.close();
             }
		}
	}
	
	/**
	 * Rename the files. 
	 */
	void renameDatFilesToMgfNames(){
		
		Iterator<String> datfiles = indexer.keySet().iterator();
		while(datfiles.hasNext()){
			String file = datfiles.next();
			File f = new File(file);
			f.renameTo(new File(indexer.get(file)));
		}
	}
	
	/**
	 * 
	 * @param mainDirectory
	 */
	void flowOfCommands(String mainDirectory) throws Exception{
		File []  datfiles = getNamesOfDatFiles(mainDirectory);
		findMgfNamesForDatFiles(datfiles);
		renameDatFilesToMgfNames();
	}
	
	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args)throws Exception{
		
		String mainDirectory = args[0];
		
		Rename_DatFileToMgfName rd = new Rename_DatFileToMgfName();
		rd.flowOfCommands(mainDirectory);
		
	}
	
}
