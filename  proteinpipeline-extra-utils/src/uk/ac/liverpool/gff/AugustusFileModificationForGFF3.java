package uk.ac.liverpool.gff;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

/**
 * Code for modifying the GFF files produced by Augustus.
 * 
 * This code should be run for each Chromosome prediction produced by Augustus. Each run
 * should refer to the chromosome's corresponding GFF. All the produced GFF should be merged
 * in the end. Finally, the combined FASTA file can be concatenated to the merged GFF to produce
 * a final file.  
 * 
 * 
 * This program takes Augustus' native  GFF and FASTA file, and produces a GFF3 version of it.
 * A sample of a Augustus GFF is -
 * 
	# start gene g57
	TGME49_chrXII	AUGUSTUS	gene	587471	588066	0.01	+	.	g57
	TGME49_chrXII	AUGUSTUS	transcript	587471	588066	0.01	+	.	g57.t1
	TGME49_chrXII	AUGUSTUS	tss	587471	587471	.	+	.	transcript_id "g57.t1"; gene_id "g57";
	TGME49_chrXII	AUGUSTUS	5'-UTR	587471	587555	0.18	+	.	transcript_id "g57.t1"; gene_id "g57";
	TGME49_chrXII	AUGUSTUS	start_codon	587556	587558	.	+	0	transcript_id "g57.t1"; gene_id "g57";
	TGME49_chrXII	AUGUSTUS	CDS	587556	587852	0.97	+	0	transcript_id "g57.t1"; gene_id "g57";
	TGME49_chrXII	AUGUSTUS	stop_codon	587850	587852	.	+	0	transcript_id "g57.t1"; gene_id "g57";
	TGME49_chrXII	AUGUSTUS	3'-UTR	587853	588066	0.03	+	.	transcript_id "g57.t1"; gene_id "g57";
	TGME49_chrXII	AUGUSTUS	tts	588066	588066	.	+	.	transcript_id "g57.t1"; gene_id "g57";
	# protein sequence = [MAPKKVTKKGTEGKKKRAKKDPNAPKKPLSSYMFFAKDKRAEILKKQPTLKSDIGKVGKMIGEEWAKLSSSQKMTYQK
	# KAEQEKIRYQREMSLYNKKK]
	# end gene g57
 	
 	And the corresponding FASTA sequence is -
 	>TGME49_chrXII.aa_g57.t1
	MAPKKVTKKGTEGKKKRAKKDPNAPKKPLSSYMFFAKDKRAEILKKQPTLKSDIGKVGKMIGEEWAKLSSSQKMTYQKKAEQEKIRYQREMSLYNKKK
 * 
 * @author riteshk
 *
 */
public class AugustusFileModificationForGFF3 {

	String augustusGFF;
	String outputGFF;
	
	public AugustusFileModificationForGFF3(String inputAugustusGff, String gffToProduce){
		this.augustusGFF   = inputAugustusGff;
		this.outputGFF     = gffToProduce;
	}
	
	void performConversion(){
		
		try{
			BufferedReader in =  new BufferedReader(new FileReader(this.augustusGFF));
			String line;
			
			BufferedWriter out = new BufferedWriter(new FileWriter(this.outputGFF));
			
			while((line = in.readLine()) != null){
				
				if(line.startsWith("##")){
					continue;
				}
				
				String [] records = line.split("\t",-1);
				
				if(records.length != 9){
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
				
				if(type.trim().matches("GENE") || type.trim().matches("gene")){
					// The genes have to be renamed as well as shown in the example -
					// gene 57 in chr XII is named as - TGME49_chrXII.aa_g57
					String gene_attr = "ID=" + seqId + ".aa_" +attribute;
					
					String revised_line = seqId + "\t" + source + "\t" + type + "\t" +
					start + "\t" + end + "\t" + score + "\t" +
					strand + "\t" + phase + "\t" + gene_attr;

					out.write(revised_line + "\n");
				}
				if(type.trim().matches("CDS") || type.trim().matches("cds")){
					// replace transcript_id "g1.t1"; gene_id "g1";
					//line.replaceAll("transcript_id \"", "ID=");
					//line.replaceAll("gene_id \"", "parent=");
					//line.replaceAll("\";",";" );
					String [] attrs = attribute.split(";",-1);
					for(String thisAttr : attrs){
						if(thisAttr.trim().isEmpty())
							continue;
						String [] keyval = thisAttr.trim().split(" "); 
						// Remove "" from the values
						if(keyval[1] != null && keyval[1].length() > 0 && keyval[1].contains("\"")){
							keyval[1] = keyval[1].substring(keyval[1].indexOf("\""), keyval[1].lastIndexOf("\""));
						}
						
						if(thisAttr.contains("gene_id")) {
							attribute = attribute.replaceAll(keyval[1], seqId.concat(".aa_" + keyval[1]));
							attribute = attribute.replaceAll("gene_id", "parent=");	
						}else if(thisAttr.contains("transcript_id")){
							attribute = attribute.replaceAll("transcript_id", "ID=");
							//attribute = attribute.replaceAll(keyval[1], seqId.concat(".aa_" + keyval[1]));
						} 
					}
					// Remove all "
					attribute = attribute.replaceAll("\"","");
					
					String revised_line = seqId + "\t" + source + "\t" + type + "\t" +
											start + "\t" + end + "\t" + score + "\t" +
											strand + "\t" + phase + "\t" + attribute;
					
					System.out.println(revised_line + "\n");
					
					out.write(revised_line + "\n");
				}
				
			} // end of while
			
			in.close();
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String inputAugustusGff = args[0]; 		// Input : Augustus Native GFF
		String gffToProduce = args[1];			// Output : GFF3 file with FASTA sequences
		
		AugustusFileModificationForGFF3 ag = new AugustusFileModificationForGFF3(inputAugustusGff, gffToProduce);
		ag.performConversion();

	}

}
