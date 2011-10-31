package uk.ac.liverpool.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class DiscoverPeptidesAcrossDiffernetModels {
	
	public HashMap <String, ArrayList<Peptide>> peptideMap_official;
	public HashMap <String, ArrayList<Peptide>> peptideMap_augustus;
	public HashMap <String, ArrayList<Peptide>> peptideMap_glimmer;
	
	// Common peptides across models
	public HashMap<String, ArrayList<ArrayList<Peptide>>> bucket_official_augustus
						= new HashMap<String, ArrayList<ArrayList<Peptide>>>() ;
	public HashMap<String, ArrayList<ArrayList<Peptide>>> bucket_official_glimmer
						= new HashMap<String, ArrayList<ArrayList<Peptide>>>() ;
	public HashMap<String, ArrayList<ArrayList<Peptide>>> bucket_agustus_glimmer
						= new HashMap<String, ArrayList<ArrayList<Peptide>>>() ;
	public HashMap<String, ArrayList<ArrayList<Peptide>>> bucket_official_augustus_glimmer
						= new HashMap<String, ArrayList<ArrayList<Peptide>>>() ;
	
	// Unique peptides in each model
	public HashMap<String, ArrayList<Peptide>> bucket_unique_official
			= new HashMap<String, ArrayList<Peptide>>();
	public HashMap<String, ArrayList<Peptide>> bucket_unique_augustus
			= new HashMap<String, ArrayList<Peptide>>();
	public HashMap<String, ArrayList<Peptide>> bucket_unique_glimmer
			= new HashMap<String, ArrayList<Peptide>>();
	
	/**
	 * All the parameters are needed for calls to AnalyseOutputFromProteinPipeline 
	 * @param pipelineSummaryFile
	 * @param fdrThreshold
	 * @param delimiter
	 * @param decoyString
	 * @return PeptideMap datastructure
	 */
	HashMap <String, ArrayList<Peptide>> createMaps(String pipelineSummaryFile, double fdrThreshold, String delimiter,String decoyString){
		
		AnalyseOutputFromProteinPipeline ap = new AnalyseOutputFromProteinPipeline(pipelineSummaryFile, fdrThreshold, delimiter,decoyString);
		ap.createMaps();
		
		return ap.peptideMap;
	}
	
	/**
	 * 
	 * @param summaryFile_official
	 * @param summaryFile_augustus
	 * @param summaryFile_glimmer
	 * @param fdrThreshold
	 * @param delimiter
	 * @param decoyString
	 */
	void findCommonPeptidesAcrossModels(String summaryFile_official, String summaryFile_augustus,String summaryFile_glimmer,
			double fdrThreshold, String delimiter,String decoyString){
		
		// Find Peptide related information from all three models
		this.peptideMap_official = createMaps(summaryFile_official, fdrThreshold, delimiter, decoyString);
		this.peptideMap_augustus = createMaps(summaryFile_augustus, fdrThreshold, delimiter, decoyString);
		this.peptideMap_glimmer = createMaps(summaryFile_glimmer, fdrThreshold, delimiter, decoyString);
		
		
		Iterator<String> peptide_official = peptideMap_official.keySet().iterator();
		while(peptide_official.hasNext()){
			String peptideSeq_official = peptide_official.next();
			
			// For Official + Augustus 
			if(peptideMap_augustus.containsKey(peptideSeq_official)){
				ArrayList<ArrayList<Peptide>> pepColl = new ArrayList<ArrayList<Peptide>>();
				
				pepColl.add(peptideMap_official.get(peptideSeq_official)); // 0 - Official
				pepColl.add(peptideMap_augustus.get(peptideSeq_official)); // 1 -Augustus
				
				bucket_official_augustus.put(peptideSeq_official, pepColl);
			}
			
			// For Official + Glimmer
			if(peptideMap_glimmer.containsKey(peptideSeq_official)){
				ArrayList<ArrayList<Peptide>> pepColl = new ArrayList<ArrayList<Peptide>>();
				
				pepColl.add(peptideMap_official.get(peptideSeq_official)); // 0 - Official
				pepColl.add(peptideMap_glimmer.get(peptideSeq_official));  // 1 -Glimmer
				
				bucket_official_glimmer.put(peptideSeq_official, pepColl);
			}
			
			// For Official + Augustus + Glimmer
			if(peptideMap_augustus.containsKey(peptideSeq_official) && peptideMap_glimmer.containsKey(peptideSeq_official)){
				ArrayList<ArrayList<Peptide>> pepColl = new ArrayList<ArrayList<Peptide>>();
				
				pepColl.add(peptideMap_official.get(peptideSeq_official)); // 0 - Official
				pepColl.add(peptideMap_augustus.get(peptideSeq_official)); // 1 -Augustus
				pepColl.add(peptideMap_glimmer.get(peptideSeq_official));  // 2 -Glimmer
				
				bucket_official_augustus.put(peptideSeq_official, pepColl);
			}
			
		}
		
		// For Augustus + Glimmer
		Iterator<String> peptide_augustus = peptideMap_augustus.keySet().iterator();
		while(peptide_augustus.hasNext()){
			String peptideSeq_augustus = peptide_augustus.next();
			if(peptideMap_glimmer.containsKey(peptideSeq_augustus)){
				ArrayList<ArrayList<Peptide>> pepColl = new ArrayList<ArrayList<Peptide>>();
				
				pepColl.add(peptideMap_augustus.get(peptideSeq_augustus)); // 0 - Augustus
				pepColl.add(peptideMap_glimmer.get(peptideSeq_augustus));  // 1 -Glimmer
				
				bucket_official_glimmer.put(peptideSeq_augustus, pepColl);
			}
		}
		
		/**************************************************************************************/
		// Unique peptide finding...
		
		Iterator<String> peptide_official_unique = peptideMap_official.keySet().iterator();
		while(peptide_official_unique.hasNext()){
			String peptideSeq = peptide_official_unique.next();
			if(!(peptideMap_augustus.containsKey(peptideSeq) || peptideMap_glimmer.containsKey(peptideSeq))){
				bucket_unique_official.put(peptideSeq, peptideMap_official.get(peptideSeq));
			}
		}
		
		Iterator<String> peptide_augustus_unique = peptideMap_augustus.keySet().iterator();
		while(peptide_augustus_unique.hasNext()){
			String peptideSeq = peptide_augustus_unique.next();
			if(!(peptideMap_official.containsKey(peptideSeq) || peptideMap_glimmer.containsKey(peptideSeq))){
				bucket_unique_augustus.put(peptideSeq, peptideMap_augustus.get(peptideSeq));
			}
		}
		
		Iterator<String> peptide_glimmer_unique = peptideMap_glimmer.keySet().iterator();
		while(peptide_glimmer_unique.hasNext()){
			String peptideSeq = peptide_glimmer_unique.next();
			if(!(peptideMap_official.containsKey(peptideSeq) || peptideMap_augustus.containsKey(peptideSeq))){
				bucket_unique_glimmer.put(peptideSeq, peptideMap_glimmer.get(peptideSeq));
			}
		}
		
	}
	
	/**
	 * 
	 * @param bucket
	 */
	void analyzeBuckets(HashMap<String, ArrayList<ArrayList<Peptide>>> bucket){
		int totalEntries = bucket.size();
		System.out.println(" Total common peptides in bucket = " + totalEntries);
	}
	
	void analyzeBuckets_unique(HashMap<String, ArrayList<Peptide>> bucket){
		int totalEntries = bucket.size();
		System.out.println(" Total common peptides in bucket = " + totalEntries);
	}
	
	
	public static void main(String [] args){
		
		String summaryFile_official = args[0];
		String summaryFile_augustus = args[1];
		String summaryFile_glimmer = args[2];
		double fdrThreshold = Double.parseDouble(args[3]);
		String delimiter = args[4];
		String decoyString = args[5];
		
		DiscoverPeptidesAcrossDiffernetModels dp = new DiscoverPeptidesAcrossDiffernetModels();
		dp.findCommonPeptidesAcrossModels(summaryFile_official, summaryFile_augustus, summaryFile_glimmer, fdrThreshold, delimiter, decoyString);
		
		System.out.println("\n Bucket - Official + Augustus");
		dp.analyzeBuckets(dp.bucket_official_augustus);
		
		System.out.println("\n Bucket - Official + Glimmer");
		dp.analyzeBuckets(dp.bucket_official_glimmer);
		
		System.out.println("\n Bucket - Augustus + Glimmer");
		dp.analyzeBuckets(dp.bucket_agustus_glimmer);
		
		System.out.println("\n Bucket - Official + Augustus + Glimmer");
		dp.analyzeBuckets(dp.bucket_official_augustus_glimmer);
		
		System.out.println("\n Bucket - Unique Official ");
		dp.analyzeBuckets_unique(dp.bucket_unique_official);
		
		System.out.println("\n Bucket - Unique Augustus ");
		dp.analyzeBuckets_unique(dp.bucket_unique_augustus);
		
		System.out.println("\n Bucket - Unique Glimmer ");
		dp.analyzeBuckets_unique(dp.bucket_unique_glimmer);
		
	}
}