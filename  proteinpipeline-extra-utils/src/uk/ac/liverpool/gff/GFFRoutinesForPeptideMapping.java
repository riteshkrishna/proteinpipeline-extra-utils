package uk.ac.liverpool.gff;

/**
 * Program to map peptide sequences on genome co-ordinates based on the start and end locations 
 * on protein sequences. 
 * 
 * Used in PeptideMappingOnChromosome.java.
 *   
 */
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class GFFRoutinesForPeptideMapping {
		
		HashMap<String, ArrayList<CDS_Information>> cdsRecords;
		String gffFile;
		
		public GFFRoutinesForPeptideMapping(String gffFile){
			cdsRecords = new HashMap<String, ArrayList<CDS_Information>>();
			this.gffFile = gffFile;
		}
		
		/*
		 * Process GFF and create CDS records
		 */
		public void processGffFile() throws Exception{
			
			try{
				BufferedReader in =  new BufferedReader(new FileReader(gffFile));
				String line;
				boolean fastaRegionFlag = false;
				
				while((line = in.readLine()) != null){
						
					if(line.startsWith("##")){
						if(line.equals("##FASTA")){
							fastaRegionFlag = true;
						}
					}
					
					if(fastaRegionFlag == false){
						
						String [] records = line.split("\t",-1);
						if(records.length != 9){
							//System.out.println("Skipping... - " + line);
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
						
					}		
					
				} // end of while
				
				
				in.close();
	
			}catch(Exception e){
				throw e;
			}
			
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
	     * Map the co-ordinates
	     * Returns a String array - idx=0=seq-id, idx=1=start and idx=2=end
	     */
		public String[] mapToGff(ProteinPeptideObject pr){
	    	String[] gffMapping = new String[3];
	    	
	    	// Get accession from the protein object
	    	String accession = pr.getAccession();
	    	
	    	// return if no key found
	    	if(!cdsRecords.containsKey(accession))
	    		return null;
	    	
	    	//...find the cds this range falls into...
	    	ArrayList<CDS_Information> cdsColl = cdsRecords.get(accession);
	    			
	    	long[] mapped_cords = determineTheLocationOfSeqOnCds(pr, cdsColl);
	    	
	    	gffMapping[0] = Long.toString(mapped_cords[0]);
	    	gffMapping[1] = Long.toString(mapped_cords[1]);
	    	gffMapping[2] = cdsRecords.get(accession).get(0).seqid;
	    	
	    	return gffMapping;
	    }
	    
	    /**
	     * 
	     * @param pr
	     * @param cdsColl
	     * @return
	     */
	    long[] determineTheLocationOfSeqOnCds(ProteinPeptideObject pr, ArrayList<CDS_Information> cdsColl){
	    	long[] gffEntry = new long[2];
	    	
	    	// Get locations from the protein object
	    	long start = pr.getstart();
	    	long end = pr.getEnd();
	  
	    	// sort the CDS collection according to the strand
	    	ArrayList<CDS_Information> sortedCDS = sortCDSAccordingToStartPosition(cdsColl);
	
	    	long mapped_start = getMappedCordinates(start, sortedCDS);
	    	long mapped_end   = getMappedCordinates(end, sortedCDS);
	    	
	    	gffEntry[0] = mapped_start;
	    	gffEntry[1] = mapped_end;
	    	
	    	return gffEntry;
	    }
	    

	    /**
	     * 
	     * @param number
	     * @param cdsColl
	     * @return
	     */
	    long getMappedCordinates(long number, ArrayList<CDS_Information> cdsColl){
	    	long mappedCord = 0;
	    	
	    	// compute the cumm array
	    	long [] cummArray = cumulativeStartPositions(cdsColl);
	    	//long number_toMap = (number - 1) * 3;
	    	long number_toMap = number * 3;
	    	
	    	int idx = determineTheIndexInCummulativeArray(number_toMap, cummArray);
	    	
	    	long shift = cummArray[idx] - number_toMap;
	    	
	    	if(cdsColl.get(0).getStrand().contains("+")){
	    		long end_cds = cdsColl.get(idx).getEnd();
	    		mappedCord = end_cds - shift;
	    	}else{
	    		long start_cds = cdsColl.get(idx).getStart();
	    		mappedCord = start_cds + shift;
	    	}
	    	
	    	return mappedCord;
	    }
	    
	    /**
	     * Sort the CDS record, with the ascending or descending order of the start position. If the CDS is on +ve strand
	     * then sort in ascending order, otherwise, sort in descending order.
	     * 
	     * @param cdsColl
	     * @return
	     */
	    ArrayList<CDS_Information> sortCDSAccordingToStartPosition(ArrayList<CDS_Information> cdsCollection){
	    	
	    	ArrayList<CDS_Information> cdsColl = new ArrayList<CDS_Information>(cdsCollection); 
	    	
	    	String strand = cdsCollection.get(0).getStrand();
	    	boolean sortAscending = false;
	    	
	    	if(strand.contains("+"))
	    		sortAscending = true;
	    			
	    	// Ascending sort...
	    	for(int i =  0; i < cdsColl.size() - 1; i++){
	    		for(int j = i+1; j < cdsColl.size(); j++){
	    			CDS_Information cds_i = cdsColl.get(i);
	    			CDS_Information cds_j = cdsColl.get(j);
	    			if(cds_i.getStart() > cds_j.getStart()){
	    				CDS_Information temp = cds_i;
	    				cdsColl.set(i, cds_j);
	    				cdsColl.set(j, temp);
	    			}
	    		}
	    	}
	    	
	    	// If we need it in descending order, then reverse the sorted array
	    	if(!sortAscending){
	    		for(int i = cdsColl.size() - 1; i >= 0; i--){
	    			int j = cdsColl.size() - 1 - i;
	    			if(j < i){
	    				CDS_Information temp = cdsColl.get(i);
	    				cdsColl.set(i, cdsColl.get(j));
	    				cdsColl.set(j, temp);
	    			}
	    		}
	    	}
	    	
	    	return cdsColl;
	    	
	    }
	    
	    /**
	     * 
	     * @param cdsCollection
	     * @return
	     */
	    long [] cumulativeStartPositions(ArrayList<CDS_Information> cdsCollection){
	    	long [] cummStartPosition = new long[cdsCollection.size()];
	    	
	    	// compute Diff = end - start
	    	for(int i = 0 ; i < cdsCollection.size() ; i++){
	    		//cummStartPosition[i] = cdsCollection.get(i).getEnd() - cdsCollection.get(i).getStart();
	    		cummStartPosition[i] = cdsCollection.get(i).getEnd() - cdsCollection.get(i).getStart() + 1;
	    	}
	    	
	    	// Form cumulative array 
	    	for(int i = 1; i < cummStartPosition.length ; i++){
	    		cummStartPosition[i] = cummStartPosition[i] + cummStartPosition[i-1];
	    	}
	    	
	    	return cummStartPosition;
	    }
	    
	    /**
	     * 
	     * @param number
	     * @param cummArray
	     * @return
	     */
	    int determineTheIndexInCummulativeArray(long number, long[] cummArray){
	    	int idx = -1;
	    	
	    	try{
	    		
	    		for(int i = 0; i < cummArray.length; i++){
	    			if(number <= cummArray[i]){
	    				idx = i;
	    				break;
	    			}
	    		}
	    	
	    	}catch(Exception e){
	    		System.out.println(e.getMessage());
	    		e.printStackTrace();
	    	}
	    	
	    	return idx;
	    }
	    
			
}