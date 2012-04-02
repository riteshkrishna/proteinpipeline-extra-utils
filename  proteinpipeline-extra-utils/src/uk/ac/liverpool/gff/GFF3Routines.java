package uk.ac.liverpool.gff;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Set of some general purpose functions to extract information from a GFF3 file.
 * We had bit and pieces of code here and there, but as we are dealing increasingly with
 * GFF3 format, I decided to put some general purpose routines in this class.
 * 
 * There is a Fasta relate routine already in uk.ac.liverpool.gff.ExtractFastFromGFF class
 * 
 * @author riteshk
 *
 */
public class GFF3Routines {

	HashMap<String, ArrayList<CDS_Information>> cdsRecords;
	HashMap<String, Gene_Information> geneRecords;
	HashMap<String, ArrayList<Peptide_Information>> peptideRecords;
	
	public GFF3Routines(){
		cdsRecords = new HashMap<String, ArrayList<CDS_Information>>();
		geneRecords = new HashMap<String, Gene_Information> ();
		peptideRecords = new HashMap<String, ArrayList<Peptide_Information>> ();
	}
	
	/**
	 * 
	 * @param gffFile
	 * @throws Exception
	 */
	public void processGff(String gffFile) throws Exception{
		
		try{
			BufferedReader in =  new BufferedReader(new FileReader(gffFile));
			String line;
			boolean fastaRegionFlag = false;
			
			
			while((line = in.readLine()) != null){
				
				//System.out.println(line);
				
				if(line.startsWith("##")){
					if(line.equals("##FASTA")){
						fastaRegionFlag = true;
						break;
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
					
					// Create CDS map
					if(type.trim().matches("CDS"))
						process_cds(line, seqId, source, startPos, endPos, strand, phase, attribute);
					
					// Create gene map
					if(type.trim().matches("gene"))
						process_genes(line, seqId, source, startPos, endPos, strand, phase, attribute);
					
					// Create peptide map
					if(type.trim().matches("peptide"))
						process_peptides(line, seqId, source, startPos, endPos, strand, phase, attribute);
					
					
				}else{
					// Now we are in the FASTA section, so we should have all the cds information
					System.out.println("Entering the ##FASTA section...total CDS found - " + cdsRecords.size() );
					break;
				}
				
			} // end of while
			
			
			in.close();
					
		}catch(Exception e){
			throw e;
		}
	}
	
	/*
	 * Process peptide fields in GFF
	 */
	void process_peptides(String line, String seqId, String source, long start,long end,String strand,String phase,String attribute)
	throws Exception{
		
		Peptide_Information peptideObj = new Peptide_Information(seqId, source, start, end, strand, phase, attribute);
		String pepId = parseAttributeFieldToExtractId(attribute);
		
		ArrayList<Peptide_Information> pepColl = peptideRecords.get(pepId);
		if(pepColl == null){
			pepColl = new ArrayList<Peptide_Information>();
		}

		pepColl.add(peptideObj);
		this.peptideRecords.put(pepId, pepColl);
	}
	
	/*
	 * Process gene fields in GFF
	 */
	void process_genes(String line, String seqId, String source, long start,long end,String strand,String phase,String attribute)
	throws Exception{
		
		Gene_Information geneObj = new Gene_Information(seqId, source, start, end, strand, phase, attribute);
		String geneId = parseAttributeFieldToExtractId(attribute); 
		this.geneRecords.put(geneId, geneObj);
	}
	
	
	/*
	 * Process CDS fields in the GFF
	 */
	void process_cds(String line, String seqId, String source, long start,long end,String strand,String phase,String attribute)
			throws Exception{
			
			CDS_Information cdsObj = new CDS_Information(seqId, source, start, end, strand, phase, attribute);
			String cdsId = parseAttributeFieldToExtractId(attribute); 
		
			ArrayList<CDS_Information> cdsColl = cdsRecords.get(cdsId);
		
			if(cdsColl == null){
				cdsColl = new ArrayList<CDS_Information>();
			}

			cdsColl.add(cdsObj);
			this.cdsRecords.put(cdsId, cdsColl);
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
	 * @param outputFileName
	 * @throws Exception
	 */
   public void writeGeneInfoInFile(String outputFile) throws Exception{
	   BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
	   
	   Iterator<String> geneList = geneRecords.keySet().iterator();
	   while(geneList.hasNext()){
		   String geneId = geneList.next();
		   Gene_Information geneObj = geneRecords.get(geneId);
		   String chromosome = geneObj.getSeqID();
		   long start = geneObj.getStart();
		   long end = geneObj.getEnd();
		   
		   out.write("\n" + chromosome + " " + start + " " + end + " " + geneId);
	   }
	   out.close();
   }
   
   /**
    * 
    * @param outputFile
    * @throws Exception
    */
   public void writePeptideInfoInFile(String outputFile) throws Exception{
	   BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
	   
	   Iterator<String> peptideList = peptideRecords.keySet().iterator();
	   while(peptideList.hasNext()){
		   String peptideId = peptideList.next();
		   ArrayList<Peptide_Information> peptideColl = peptideRecords.get(peptideId);
		   for(Peptide_Information pepObj : peptideColl){
			   String chromosome = pepObj.getSeqID();
			   long start = pepObj.getStart();
			   long end = pepObj.getEnd();
			   
			   out.write("\n" + chromosome + " " + start + " " + end + " " + peptideId);
		   }
		   
	   }
	   out.close();
   }
   
   /**
    * 
    * @param args
    */
	public static void main(String[] args) throws Exception{
		String gffInput = args[0];  
		
		GFF3Routines gffHandle = new GFF3Routines();
		gffHandle.processGff(gffInput);
		
		//String geneOutput = "result/GeneInfo.txt";
		//gffHandle.writeGeneInfoInFile(geneOutput);
		
		String allPeptideOutput = "result/PeptideInfo.txt";
		gffHandle.writePeptideInfoInFile(allPeptideOutput);
		
	}

}
