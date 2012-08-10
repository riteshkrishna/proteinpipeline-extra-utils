package uk.ac.liverpool.gff;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;

/**
 * We want to find the nearest official genes to the alternate peptides. This can help us browse the
 * annotation using gene browsers, as we will know that which official genes are closer to the peptides
 * which were not present in the official annotation.
 * 
 * @author riteshk
 *
 */
public class FindNearestGenesForPeptides {
	
	public HashMap<String, Gene_Information> geneRecords;
	public HashMap<String, HashMap<String, Gene_Information>> chrBuckets;
	
	public FindNearestGenesForPeptides(){
		geneRecords = new HashMap<String, Gene_Information> ();
		chrBuckets = new HashMap<String, HashMap<String, Gene_Information>>();
	}
	
	/**
	 * Read all the genes related information from the GFF
	 * @param gffInput
	 */
	public void getGeneRecordsFromGff(String gffInput){
		try{
			GFF3Routines gffHandle = new GFF3Routines();
			gffHandle.processGff(gffInput);
		
			geneRecords = gffHandle.geneRecords;
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 */
	public void divideGeneRecordsInChromosomes(){
		
		Iterator<String> geneAccn = geneRecords.keySet().iterator();
		while(geneAccn.hasNext()){
			String thisaccn = geneAccn.next();
			String thischr = geneRecords.get(thisaccn).getSeqID().trim();
			
			// Ignoring the peculiar chromosomes from Toxo gff files
			if(thischr.contains("apidb|DS"))
				continue;
			
			if(chrBuckets.containsKey(thischr)){
				HashMap<String, Gene_Information> contained = chrBuckets.get(thischr);
				contained.put(thisaccn, geneRecords.get(thisaccn));
				chrBuckets.put(thischr, contained);
			}else{
				HashMap<String, Gene_Information> contained = new HashMap<String, Gene_Information>();
				contained.put(thisaccn, geneRecords.get(thisaccn));
				chrBuckets.put(thischr, contained);
			}
		}
	}
	
	/**
	 * Sort gene list in each chr bucket according to their start position on the chr
	 */
	public void sortGeneListInEachChrBucket(){
		
		Iterator<String> chrs = chrBuckets.keySet().iterator();
		while(chrs.hasNext()){
			String chr = chrs.next();
			HashMap<String, Gene_Information> geneColln = chrBuckets.get(chr);
			geneColln = sortHashMap(geneColln);
			chrBuckets.put(chr, geneColln);
		}
	}
	
	private HashMap<String, Gene_Information> sortHashMap(HashMap<String, Gene_Information> input){
	    Map<String, Gene_Information> tempMap = new HashMap<String, Gene_Information>();
	    
	    for (String geneAccn : input.keySet()){
	        tempMap.put(geneAccn,input.get(geneAccn));
	    }

	    List<String> mapKeys = new ArrayList<String>();
	    List<Long> mapValues = new ArrayList<Long>();
	    
	    Iterator<String> keys = tempMap.keySet().iterator();
	    while(keys.hasNext()){
	    	String key = keys.next();
	    	long start = geneRecords.get(key).getStart(); // Pick it from the other ds, contains the same info.
	    	mapKeys.add(key);
	    	mapValues.add(start);
	    }
	    
	    HashMap<String, Gene_Information> sortedMap = new LinkedHashMap<String, Gene_Information>();
	    TreeSet<Long> sortedSet = new TreeSet<Long>(mapValues);
	    Object[] sortedArray = sortedSet.toArray();
	    int size = sortedArray.length;
	    for (int i=0; i<size; i++){
	        sortedMap.put(mapKeys.get(mapValues.indexOf(sortedArray[i])), 
	                      geneRecords.get(mapKeys.get(mapValues.indexOf(sortedArray[i]))));
	    }
	    return sortedMap;
	}
	
	/**
	 * 
	 * @param peptideInput
	 * @param outputFile
	 */
	public void mapFromPeptideFile(String peptideInput, String outputFile){
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
			  
			BufferedReader in =  new BufferedReader(new FileReader(peptideInput));
			
			String line;
			
			while((line = in.readLine()) != null){
				
				String [] records = line.split("\t",-1);
				
				if(records.length != 5){
					System.out.println("Skipping... - " + line);
					continue;
				}
				
				String pep_chr = records[0];
				long pep_start = Long.parseLong(records[1]);
				long pep_end = Long.parseLong(records[2]);
				String pep_seq = records[3];
				String prot_accn = records[4];
				
				// now find the chr, and the proper index....
				String chrToSearchFor = pep_chr.substring(pep_chr.indexOf("chr"), pep_chr.length());
				chrToSearchFor = "apidb|TGME49_" + chrToSearchFor;
				
				System.out.println("Chr to search : " + chrToSearchFor + " Pep Seq : " + pep_seq);
				
				HashMap<String, Gene_Information> geneColl = chrBuckets.get(chrToSearchFor);
				
				
				
				String prevAccn = new String();
				long prevstart = 0;
				long prevend = 0;
			    for(String accn : geneColl.keySet()){
			    	long genestart = geneColl.get(accn).start;
			    	long geneend   = geneColl.get(accn).end;
			    	
			    	if(pep_end <= geneend){
			    		if(pep_start >= genestart){
			    			// report - in the middle of cds
			    			out.write(chrToSearchFor + "\t" + pep_start + "\t" + pep_end + "\t" + pep_seq + "\t"
			    					+ prot_accn + "\t" + "Within" + "\t" + accn + "("+ genestart + "," + geneend +")" + "\n");
			    			break;
			    		}else{
			    			// report - in the middle of cds
			    			out.write(chrToSearchFor + "\t" + pep_start + "\t" + pep_end + "\t" + pep_seq + "\t"
			    					+ prot_accn + "\t" + "Middle" + "\t" + accn + "("+ genestart + "," + geneend +")" + "\t" + prevAccn + "("+ prevstart + "," + prevend +")"+ "\n");
			    			break;
			    		}
			    	}
			    	
			    	
			    	prevAccn = accn;
			    	prevstart = genestart;
			    	prevend = geneend;
			    }
			}

		    in.close();
		    out.close();

		}catch(Exception e){
				e.printStackTrace();
		}
	}
	
	/**
	 * Step by step instructions...
	 * @param gffInput
	 */
	public void performStepsForGff(String gffInput, String peptideInput,String outputFile){
		getGeneRecordsFromGff(gffInput);
		divideGeneRecordsInChromosomes();
		sortGeneListInEachChrBucket();
		
		mapFromPeptideFile(peptideInput,outputFile);
	}
	
	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception{
		
		String gffInput = "tmp/Toxo_1D_OFF.gff";
		String peptideInput = "tmp/mapped_glimmer_peptides.txt";
		String outputFile = "tmp/whereabout-in-gff-glimmer.txt";
		
		FindNearestGenesForPeptides fg = new FindNearestGenesForPeptides();
		fg.performStepsForGff(gffInput,peptideInput,outputFile);
		
		//System.out.println(fg.geneRecords.size());
		//System.out.println(fg.chrBuckets.size() + "\t" + fg.chrBuckets.keySet().toString());
		
		//HashMap<String, Gene_Information> trialColl = fg.chrBuckets.get(fg.chrBuckets.keySet().iterator().next());
		//for(String accn : trialColl.keySet())
		//	System.out.println(accn + "\t" + trialColl.get(accn).start + "\t" + trialColl.get(accn).end);
		
	}
}
