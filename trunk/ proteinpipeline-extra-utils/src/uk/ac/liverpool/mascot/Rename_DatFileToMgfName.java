package uk.ac.liverpool.mascot;

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

}
