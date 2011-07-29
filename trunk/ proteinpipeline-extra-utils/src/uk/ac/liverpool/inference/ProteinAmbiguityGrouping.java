package uk.ac.liverpool.inference;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import uk.ac.liverpool.analysis.Peptide;

public class ProteinAmbiguityGrouping {
	
	HashMap<String, Integer> protein_pagnumber;
	HashMap<Integer,ArrayList<String>> pags;
	
	HashMap <String, ArrayList<Peptide>> peptideMap;
	HashMap <String, ArrayList<String>> proteinPeptideMap;
	
	int pagCounter = 0; 
	
	/**
	 * 
	 * @param peptideMap
	 * @param proteinPeptideMap
	 */
	public ProteinAmbiguityGrouping(HashMap <String, ArrayList<Peptide>> peptideMap, 
			HashMap <String, ArrayList<String>> proteinPeptideMap){
		
		this.peptideMap = peptideMap;
		this.proteinPeptideMap = proteinPeptideMap;
		this.protein_pagnumber = new HashMap<String, Integer>();
		this.pags = new HashMap<Integer,ArrayList<String>>();
	}
	
	/**
	 * 
	 * @return
	 */
	public HashMap<Integer,ArrayList<String>> createAmbiguityGroup(){
		
		Iterator<String> peptides = peptideMap.keySet().iterator();
		
		HashMap<String,String> proteinsWithSinglePeptides = new HashMap<String,String>();
		
		// Iterate for each peptide 
		while(peptides.hasNext()){
			
			String peptideSeq = peptides.next();
			
			// Get the name of proteins containing this peptide
			ArrayList<Peptide> associatedPeptides = peptideMap.get(peptideSeq);
			ArrayList<String> associatedProteins  = new ArrayList<String>();
			for(int i = 0; i < associatedPeptides.size(); i++)
				associatedProteins.add(associatedPeptides.get(i).protAccn);
			
			// Get the name of peptides, the associated proteins are containing
			HashMap<String, ArrayList<String>> proteinsAndRelatedPeptides = new HashMap<String, ArrayList<String>>();
			for(int i = 0; i < associatedProteins.size(); i++){
				String protAccn = associatedProteins.get(i);
				ArrayList<String> relatedPeptide = proteinPeptideMap.get(protAccn);
				
				// Proteins with single peptide should not be considered
				if(relatedPeptide.size() == 1){
					proteinsWithSinglePeptides.put(protAccn, peptideSeq);
					continue;
				}
				
				proteinsAndRelatedPeptides.put(protAccn, relatedPeptide);
			}	
			
			if(proteinsAndRelatedPeptides.size() <= 1)
				continue;
			
			ArrayList<ArrayList<String>> groups = formGroups(proteinsAndRelatedPeptides);
			
			if(groups != null){
				
				for(int i = 0; i < groups.size(); i++){
					
					String prot1 = groups.get(i).get(0);
					String prot2 = groups.get(i).get(1);
					
					int pag_idx = -1;
					if(protein_pagnumber.containsKey(prot1))
						pag_idx = protein_pagnumber.get(prot1);
					else if(protein_pagnumber.containsKey(prot2))
						pag_idx = protein_pagnumber.get(prot2);
					
					ArrayList<String>existingPag;
					if(pag_idx == -1){
						existingPag = new ArrayList<String>();
						pag_idx = ++pagCounter;
					}else{
						existingPag = pags.get(new Integer(pag_idx));
					}
					if(!existingPag.contains(prot1))
						existingPag.add(prot1);
					if(!existingPag.contains(prot2))
						existingPag.add(prot2);
					pags.put(new Integer(pag_idx), existingPag);
					protein_pagnumber.put(prot1, pag_idx);
					protein_pagnumber.put(prot2, pag_idx);
						
				}
			}
			
		} // End of while loop
		
		return pags;
	}
	
	/**
	 * 
	 * @param proteinsAndRelatedPeptides
	 * @return
	 */
	ArrayList<ArrayList<String>> formGroups(HashMap<String, ArrayList<String>>  proteinsAndRelatedPeptides){
		
		ArrayList<ArrayList<String>> groups = new ArrayList<ArrayList<String>>();
		
		Iterator<String> proteins = proteinsAndRelatedPeptides.keySet().iterator();
		
		ArrayList<ArrayList<String>> peptidesForProteins = new ArrayList<ArrayList<String>>();
		ArrayList<String> protAccns = new ArrayList<String>();
		
		while(proteins.hasNext()){
			String protein = proteins.next();
			protAccns.add(protein);
			peptidesForProteins.add(proteinsAndRelatedPeptides.get(protein));
		}
		
		for(int i = 0 ; i < peptidesForProteins.size() - 1; i++){
			ArrayList<String> set1 = new ArrayList<String>(peptidesForProteins.get(i));
			ArrayList<String> set2 = new ArrayList<String>(peptidesForProteins.get(i+1));
			
			ArrayList<String> copy_set1 = new ArrayList<String>(set1);
			
			copy_set1.retainAll(set2); // Take set intersection
			
			if((copy_set1.size() == set1.size()) || (copy_set1.size() == set2.size())){
				// Group exists
				ArrayList<String> g = new ArrayList<String>();
				g.add(protAccns.get(i));
				g.add(protAccns.get(i+1));
				groups.add(g); // Make one group
			}	
		}
		
		if(groups.isEmpty())
			return null;
		
		return groups;
	}
	
	
	/**
	 * The test function..
	 * @param args
	 */
	public static void main(String [] args){
		
	}

}
