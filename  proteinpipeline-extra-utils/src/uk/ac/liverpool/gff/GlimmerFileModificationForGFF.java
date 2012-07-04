package uk.ac.liverpool.gff;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

/**
 * There are some incompatibilities with the GFF produced by Glimmer. The entries are as following -
 * 
gb|TGME49_chrIa	GlimmerHMM	gene	35888	49132	.	-	.	ID=gb|TGME49_chrIa.path1.gene4;Name=gb|TGME49_chrIa.path1.gene4
gb|TGME49_chrIa	GlimmerHMM	CDS	35888	41860	.	-	0	ID=gb|TGME49_chrIa.cds4.1;Parent=gb|TGME49_chrIa.path1.gene4;Name=gb|TGME49_chrIa.path1.gene4;Note=final-exon
gb|TGME49_chrIa	GlimmerHMM	CDS	46551	46703	.	-	0	ID=gb|TGME49_chrIa.cds4.2;Parent=gb|TGME49_chrIa.path1.gene4;Name=gb|TGME49_chrIa.path1.gene4;Note=internal-exon

We want to change the naming convention in the following way -
gb|TGME49_chrIa	GlimmerHMM	gene	35888	49132	.	-	.	ID=gb|TGME49_chrIa.path1.gene4.t1;Name=gb|TGME49_chrIa.path1.gene4.t1
gb|TGME49_chrIa	GlimmerHMM	CDS	35888	41860	.	-	0	ID=gb|TGME49_chrIa.path1.gene4;Parent=gb|TGME49_chrIa.path1.gene4.t1;Name=gb|TGME49_chrIa.path1.gene4.t1;Note=final-exon
gb|TGME49_chrIa	GlimmerHMM	CDS	46551	46703	.	-	0	ID=gb|TGME49_chrIa.path1.gene4;Parent=gb|TGME49_chrIa.path1.gene4.t1;Name=gb|TGME49_chrIa.path1.gene4.t1;Note=internal-exon

 This means that we replace the CDS-ID with gene-ID, and modify the Gene-ID by adding a '.ti' to it. '.ti' is nothing but a way to make the gene-ID 
 distinct from the children CDS-IDs.
 * 
 * @author riteshk
 *
 */
public class GlimmerFileModificationForGFF {
	
	String glimmerFile;
	String modifiedFile;
	
	public GlimmerFileModificationForGFF(String glimmerFile, String modifiedFile) {
		this.glimmerFile = new String(glimmerFile);
		this.modifiedFile = new String(modifiedFile);
	}

	public void process(){
		try{
			GFF3Routines gffHandle = new GFF3Routines();
			gffHandle.processGff(glimmerFile);
			
			HashMap<String, ArrayList<CDS_Information>> cdsRecords = gffHandle.cdsRecords;
			HashMap<String, Gene_Information> geneRecords = gffHandle.geneRecords;
			
			 // based on ID key
			 Iterator<String> geneIds = geneRecords.keySet().iterator();
			 Iterator <String> cdsIds = cdsRecords.keySet().iterator();
			 
			 // Read the whole file text 
			 String text = new Scanner( new File(glimmerFile) ).useDelimiter("\\A").next();
			 
			 while(geneIds.hasNext()){
				 String geneId = geneIds.next();
				 String new_geneId = geneId + ".t1";				 
				 text = text.replace(geneId, new_geneId);
			 }
			 
			 while(cdsIds.hasNext()){
				 String cdsId = cdsIds.next();
				 
				// something here-gb|TGME49_chrIa.cds4.1 ->gb|TGME49_chrIa.path1.gene4
				 String pattern = cdsId.substring(cdsId.lastIndexOf("cds")+3, cdsId.length()-1); //4.1
				 String number = pattern.substring(0,pattern.indexOf(".")); //4
				  
				 String preCds_tag  = cdsId.substring(0, cdsId.lastIndexOf("cds"));
				 String new_cds_id = preCds_tag + "path1.gene" + number;
				 
				 text = text.replace(cdsId, new_cds_id);
				 
			 }
			 
			 BufferedWriter writer = new BufferedWriter(new FileWriter(modifiedFile));
			 writer.write(text);
			 writer.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception{
		String glimmer_gff = "tmp/All-Glimmer-ME49.gff";
		String new_glimmer = "tmp/new-Glimmer-ME49.gff";
		
		GlimmerFileModificationForGFF glm = new GlimmerFileModificationForGFF(glimmer_gff, new_glimmer);
		glm.process();
		
	}
}
