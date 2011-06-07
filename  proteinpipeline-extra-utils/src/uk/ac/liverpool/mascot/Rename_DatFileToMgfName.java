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
	
	HashMap<File, String> indexer = new HashMap<File, String>();
	
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
			while(scanner.hasNextLine()){
				 String line = scanner.nextLine();
				 if(line.contains("FILE=")){
					 String[] keyVal = line.split("=");
					 if(keyVal[0].contains("FILE")){
						 
						 // Parse the path in keyVal[1] and extract the file name
						 String pathToParse = keyVal[1];
						 String fileName1 = "";
                         int lastSlashIndex = pathToParse.lastIndexOf("\\");
                         if(lastSlashIndex > -1){
                                 fileName1 = pathToParse.substring(lastSlashIndex + 1);
                         }else {
                                 fileName1 = pathToParse;
                         }
                         // Deal with the extension
                         int extDotIndex = fileName1.indexOf(".");
                         if(extDotIndex > -1){
                             fileName1 = fileName1.substring(0,extDotIndex);
                         }
                         else{
                             fileName1 = fileName1;
                         }
                         // Add to the HashMap with correct path
                         indexer.put(datfiles[i], datfiles[i].getParent().concat("\\" + fileName1 + ".dat"));
                         break;
					 }
				}
			}        
			scanner.close();
		}
		
	}
	
	/**
	 * Rename the files. 
	 */
	void renameDatFilesToMgfNames(){
		
		Iterator<File> datfiles = indexer.keySet().iterator();
		while(datfiles.hasNext()){
			File file = datfiles.next();
			boolean result = file.renameTo(new File(indexer.get(file)));
			System.out.println("Success - " + result);
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
