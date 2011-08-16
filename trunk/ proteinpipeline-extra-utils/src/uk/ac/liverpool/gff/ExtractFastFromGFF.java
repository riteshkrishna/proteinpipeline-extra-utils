/**
 * Reads a GFF file and creates Protein files for each Chromosome
 */
package uk.ac.liverpool.gff;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
public class ExtractFastFromGFF {
	
	HashMap<String, ArrayList<CDS_Information>> cdsRecords;
	File temp;
	
	HashMap<String, BufferedWriter> chromosomeFastaMap; 
	
	public ExtractFastFromGFF(){
		cdsRecords = new HashMap<String, ArrayList<CDS_Information>>();
		chromosomeFastaMap = new HashMap<String, BufferedWriter>();
	}
	
	/**
	 * 
	 * @param gffFile
	 * @throws Exception
	 */
	public void processGffFileAndCreateFasta(String gffFile) throws Exception{
		
		try{
			BufferedReader in =  new BufferedReader(new FileReader(gffFile));
			String line;
			boolean fastaRegionFlag = false;
			
			// flag for creating chromosome files..
			boolean fileCreatorFlag = false;
			
		    // Variables needed during processing of ##FASTA file
		    String accn = new String();
			String seq = new String();
			int fastaCount = 0;
			
			while((line = in.readLine()) != null){
				
				//System.out.println(line);
				
				if(line.startsWith("##")){
					if(line.equals("##FASTA")){
						fastaRegionFlag = true;
					}
				}
				
				if(fastaRegionFlag == false){
					
					String [] records = line.split("\t",-1);
					if(records.length != 9){
						System.out.println("Skipping... - " + line);
						continue;
					}
					
					String seqId     	= records[0];
					String source 	 	= records[1]; 
					String type 	 	= records[2];
					String start 	 	= records[3];
					String end 		 	= records[4];
					String score 		= records[5];
					String strand 		= records[6];
					String phase 		= records[7];
					String attribute 	= records[8];
					
					if(type.trim().matches("CDS")){
						long startPos;
						long endPos;
						try{
							startPos = Long.parseLong(start);
							endPos = Long.parseLong(end);
							//A simple check...
							if(endPos < startPos){
								System.out.println("End position = " + endPos + " greater than start position " + startPos+ " in Gff file");
								throw new Exception();
							}
						}catch(NumberFormatException e){
							System.out.println("Problem in processing the start and end fields of Line - \n" + line + "\n...exiting");
							e.printStackTrace();
							throw e;
						}
						
						CDS_Information cdsObj = new CDS_Information(seqId, source, startPos, endPos, strand, phase, attribute);
						String cdsId = parseAttributeFieldToExtractId(attribute); 
						
						ArrayList<CDS_Information> cdsColl = cdsRecords.get(cdsId);
						
						if(cdsColl == null){
							cdsColl = new ArrayList<CDS_Information>();
						}

						cdsColl.add(cdsObj);
						this.cdsRecords.put(cdsId, cdsColl);
					}	
					
				}else{
					
					// Now we are in the FASTA section, so we should have all the cds information
					System.out.println("Entering the ##FASTA section...total CDS found - " + cdsRecords.size() );
					
					// Create separate fasta file for each chromosome, execute only once
					if(fileCreatorFlag == false){
						HashMap<String, String> fastaFileHash = createFastaFileForEachChromosome();
						
						Iterator<String>chromosomes = fastaFileHash.keySet().iterator();
						while(chromosomes.hasNext()){
							String chrm = chromosomes.next();
							String fileName = fastaFileHash.get(chrm);
							
							try{
							FileWriter fstream_chrm = new FileWriter(fileName);
							BufferedWriter out_chrm = new BufferedWriter(fstream_chrm);
							
							this.chromosomeFastaMap.put(chrm, out_chrm);
							
							}catch(Exception e){
								System.out.println("Problem creating Chromosome file : " + chrm);
								e.printStackTrace();
							}
						}
						
						fileCreatorFlag = true;
					}
					
					if(cdsRecords.size() == 0){
						System.out.println("No CDS field found.....exiting");
						throw new Exception();
					}
					
					// If all well, then process fasta - skip the line having ##FASTA 
					if(line.equals("##FASTA")){
						continue;
					}
				
					// Write all other lines
					if(line.contains(">")){
						boolean cdsAccnFound = verifyNeededAccnOrNot(accn);
						
						// Write to file only if it is a cds Accession
						if(cdsAccnFound){
							String content = "\n" + accn + "\n" + seq;
							
							String refindedCDS = new String(accn);
							if(refindedCDS.contains(">"))
								refindedCDS = refindedCDS.replace(">", "").trim(); 
							String chrmForThisCds = this.cdsRecords.get(refindedCDS).get(0).seqid;
							
							chromosomeFastaMap.get(chrmForThisCds).write(content);
							
							fastaCount++;
						}
						accn = line;
						seq = "";
					}else{
						seq = seq.concat(line);
					}
					
				}
				
			} // end of while
			
			System.out.println("Out of while loop");
			
			// Signal error if no ##FASTA is encountered...
			if(fastaRegionFlag == false){
				System.out.println("No ##FASTA directive encountered.....Quitting");
				throw new Exception();
			}
			
			// Write the last fasta sequence in the file
			if(line == null){
				boolean cdsAccnFound = verifyNeededAccnOrNot(accn);
				if(cdsAccnFound){
					String content = "\n" + accn + "\n" + seq;
					
					String refindedCDS = new String(accn);
					if(refindedCDS.contains(">"))
						refindedCDS = refindedCDS.replace(">", "").trim(); 
					String chrmForThisCds = this.cdsRecords.get(refindedCDS).get(0).seqid;
					
					chromosomeFastaMap.get(chrmForThisCds).write(content);
					fastaCount++;
					}
			}
			
			in.close();
			
			// close the Chromosome files
			Iterator<String>chromosomes = this.chromosomeFastaMap.keySet().iterator();
			while(chromosomes.hasNext()){
				String chrm = chromosomes.next();
				chromosomeFastaMap.get(chrm).close();
			}
			
			
		}catch(Exception e){
			throw e;
		}
		
	}
	
	/**
	 * find total chromosomes and create fasta file names for each of them
	 * @return
	 */
	public HashMap<String, String> createFastaFileForEachChromosome(){
		HashMap<String, String> fastaFileHash = new HashMap<String, String>();
		
		// Find unique chromosomes
		Iterator<String> cdsKeys = this.cdsRecords.keySet().iterator();
		while(cdsKeys.hasNext()){
			String cdsId = cdsKeys.next();
			String seqID = this.cdsRecords.get(cdsId).get(0).seqid; // All CDS are on the same seq
			
			if(!fastaFileHash.containsKey(seqID)){
				fastaFileHash.put(seqID, "result/Chromosome_" + seqID);
			}
		}
		
		return fastaFileHash;
	}
	
	/**
	 * 
	 * @param attribute
	 * @return
	 */
	String parseAttributeFieldToExtractId (String attribute) throws Exception{
		String id = new String();
		
		String [] splitOnColon = attribute.split(";",-1);
		if(splitOnColon.length == 0){
			String [] keyVal = attribute.split("=",-1);
			if(keyVal.length == 0){
				System.out.println("No ID found for the CDS entry");
			}else if (keyVal[0].compareToIgnoreCase("ID")  < 0){
				System.out.println("No ID found for the CDS entry");
			}else{
				id  = keyVal[1];
			}
		}else{
			boolean found = false;
			for(int i = 0 ; i< splitOnColon.length; i++){
				
				String subAttr = splitOnColon[i];
				String [] keyVal = subAttr.split("=",-1);
				if(keyVal.length == 0){
					System.out.println("No ID found for the CDS entry");
				}else if (keyVal[0].compareToIgnoreCase("ID")  < 0){
					System.out.println("No ID found for the CDS entry");
				}else{
					id = keyVal[1];
					found = true;
				}
				if(found)
					break;
			}
			
			if(!found){
				System.out.println("No ID found for the CDS entry");
				throw new Exception();
			}
			
		}
		
		return id;
	}
	
    /**
	 * 
	 * @param accnToCheck
	 * @return
	 */
    boolean verifyNeededAccnOrNot(String accnToCheck){
    	boolean found = false;
    	
    	// Remove the ">" from accnToCheck
    	if(accnToCheck.contains(">"))
    		accnToCheck = accnToCheck.replace(">", "").trim(); 
    	
    	found = cdsRecords.containsKey(accnToCheck);
    	
    	return found;
    }
    
    /**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception{
		
		ExtractFastFromGFF vgf = new ExtractFastFromGFF();
		
		String gffFile = args[0];
		
		// GFF Reading part here...
		if(args.length == 1){	
			System.out.println("Creating Fasta file from the GFF.....");
			 vgf.processGffFileAndCreateFasta(gffFile);
			 return;
		}	
	}

	
	
}