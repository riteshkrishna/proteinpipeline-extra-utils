
package uk.ac.liverpool.analysis2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 *
 * @author jonesar
 */
public class ProcessFasta {
    
    
    private String enzymeRegex = "(?<=[KR])(?!P)";                  //defaults to trypsin
    
    public HashMap<String, String> accToSeq = new HashMap();
    public HashMap<String, String> accToDefline = new HashMap();
    public HashMap<String, String[]> accToPeptides = new HashMap();
    public ArrayList<String> allProteinIDs = new ArrayList();
    
    private HashMap<String,Double> aaMap = new HashMap();
    
    
    private String accessionRegex = "";     
    private int missedCleavages = 1;                                //Not currently implemented?
    
    private HashMap<String, Double> aminoAcidMap = new HashMap();
 
    //private Double crossLinkerMass= 156.078644;       //These need to be parameterized and sent by the calling code
    //private String crossLinkerVarResidue = "K";       //These need to be parameterized and sent by the calling code 
    //private String crossLinkerName = "dss";
     
    /*
     * E.g. for disulphide analysis, the residue is C and not yet sure what the mass shift is
     * 
     * 
     */
    private Double crossLinkerMass= -2.0;       //These need to be parameterized and sent by the calling code
    private String crossLinkedResidue = "C";       //These need to be parameterized and sent by the calling code 
    private String crossLinkerName = "_ss_";
    
 
    
    private static Double HMASS =1.007825;
    private static Double OMASS = 15.994915;
    
    private HashMap <String,Double> varMods;    //hash map of var mods to be applied  ;key = "one letter code amino acids plus [ or ] for termini" value = double mass shift
    private HashMap <String,Double> fixedMods;  //hash map of fixed mods to be applied ;key = "one letter code amino acids plus [ or ] for termini" value = double mass shift
    
    
    private HashMap <String,Integer> varModStringToUniqueID = new HashMap();    //hashmap storing a unique ID for each var mod
    private HashMap <Integer,String> varModUniqueIDToString = new HashMap();    //reverse mapping
    
    /* 
     * Empty constructor
     */
    public ProcessFasta(){
    	
    }
    
    /*
     * Constructor for testing only
     */
    public ProcessFasta(String fastaFile, Double crosslinkMass, String crossLinkResidue, String crossLinkName, int missCleave, HashMap varModifications, HashMap fixedModifications){
        

        crossLinkerMass  = crosslinkMass;
        crossLinkedResidue = crossLinkResidue;
        crossLinkerName = crossLinkName;
        missedCleavages = missCleave;
        varMods = varModifications;
        
        //Populate the hashmap giving unique IDs for all variable mods
        Integer id = 1;
        for(String varMod : varMods.keySet()){
            varModStringToUniqueID.put(varMod, id);
            varModUniqueIDToString.put(id, varMod);
            id++;            
        }
        
        fixedMods = fixedModifications;    
        init(fastaFile);
    }
    
    
    public void init(String fastaFile){
        fillAAMap();
        /*
        //Handle fixed mods
        for(String fixedMod : fixedMods.keySet()){
            addFixedModsToAAMap(fixedMod,fixedMods.get(fixedMod));            
        }
        */
        readFasta(fastaFile);
        
    }
    
    /**
     * RK adds this 
     * The main program
     * @param args
     */
    public static void main(String [] args){
    	
    	String fastaFile = "/Users/riteshk/Ritesh_Work/Toxo/Toxo_Predictions_Augustus/TGME49_chrXII.aa";
    	
    	ProcessFasta ps = new ProcessFasta();
    	ps.init(fastaFile);
    	
    	System.out.println("Total proteins : " + ps.accToPeptides.size());
    	
    	String testKey = "TGME49_chrXII.aa_g1.t1";
    	String [] peptides  = ps.accToPeptides.get(testKey);
    	for (String peptide : peptides){
    		System.out.print( "\t" + peptide);
    	}
    	
    } 
    
    
    private void readFasta(String inputFasta){
        
        try{
            
            BufferedReader in = new BufferedReader(new FileReader(inputFasta));
            String line;
            
            String currSequence = "";
            String currProtAcc = null;
            String currDefline = null; 
            
            int recordCounter = 0;
            while ((line = in.readLine()) != null) {         
                line = line.replaceAll("\n","");
                line = line.replaceAll("\r","");
                      
                if(line.contains(">")){                
                    //Insert previous into hash and reset
                    if(recordCounter != 0){              
                        currSequence = currSequence.replaceAll(" ","");
                        allProteinIDs.add(currProtAcc);                     
                       
                        String[] peptides = currSequence.split(enzymeRegex);
                        addPeptidesToMap(currProtAcc,peptides);
                        
                        accToSeq.put(currProtAcc, currSequence);
                        //System.out.println("Inserting:" + currProtAcc + "_" + currSequence);
                        accToDefline.put(currProtAcc, currDefline);
                        //System.out.println("Inserting2:" + currProtAcc + "_" + currDefline);
                        
                        currSequence = "";
                    }   
                    
                    
                    line = line.replaceAll(">", "");
                    
                    /* RK
                    int splitPos = line.indexOf(accessionRegex);
                    if(splitPos!=-1){
                        currProtAcc = line.substring(0,splitPos);
                        currDefline = line.substring(splitPos+1);
                    }
                    else{
                        System.out.println("Regular expression not found for split: " + line);
                        System.exit(0);
                    }
                    */
                    currProtAcc = line; // RK
                    currDefline = line; //RK
                    
                    
                    recordCounter++;
                }
                else{
                    currSequence += line;
                }
            }
            //handle last
            currSequence = currSequence.replaceAll(" ","");
            accToSeq.put(currProtAcc, currSequence);
            accToDefline.put(currProtAcc, currDefline);
            
            String[] peptides = currSequence.split(enzymeRegex);
            addPeptidesToMap(currProtAcc,peptides);
            //System.out.println("Put last:" + currProtAcc + " seq: " + currSequence);
            allProteinIDs.add(currProtAcc);
            
            //Close the input stream
            in.close();
        }
        catch (IOException e){
            System.err.println("Error: " + e.getMessage());
        }        
    }
    
    
    private void addPeptidesToMap(String protAcc, String[] peptides){
        
        ArrayList<String> tempPeptides = new ArrayList();
        
        if(missedCleavages ==0){
            accToPeptides.put(protAcc,peptides);
        }
        else{
            
            for(int i=0;i<peptides.length;i++){
                String outerPep = peptides[i];
                tempPeptides.add(outerPep);
                
                //TODO - test missed cleaves code!
                String tempPeptide="";
                for(int j =i+1; j<=i+missedCleavages;j++){           
                    
                    if(j<peptides.length){
                        tempPeptide+= peptides[j];
                        String mc = outerPep+tempPeptide;
                        tempPeptides.add(mc);
                        //System.out.println("Adding mc:" + mc );
                    }  
                }
            }
            
            String[] newPeptides = new String[tempPeptides.size()];
            for(int i =0; i<tempPeptides.size();i++){
               newPeptides[i] = tempPeptides.get(i);
            }
            
            accToPeptides.put(protAcc,newPeptides);
            
        }
        
    }
    
    
    private void fillAAMap(){
        
        aaMap.put("A",71.037114);
        aaMap.put("R",156.101111);
        aaMap.put("N",114.042927);
        aaMap.put("D",115.026943);
        aaMap.put("C",103.009185);
        aaMap.put("E",129.042593);
        aaMap.put("Q",128.058578);
        aaMap.put("G",57.021464);
        aaMap.put("H",137.058912);
        aaMap.put("I",113.084064);
        aaMap.put("L",113.084064);
        aaMap.put("K",128.094963);
        aaMap.put("M",131.040485);
        aaMap.put("F",147.068414);
        aaMap.put("P",97.052764);
        aaMap.put("S",87.032028);
        aaMap.put("T",101.047679);
        aaMap.put("U",150.95363);
        aaMap.put("W",186.079313);
        aaMap.put("Y",163.06332);
        aaMap.put("V",99.068414);
        aaMap.put("[",0.0);
        aaMap.put("]",0.0);
        //aaMap.put("_",0.0);             //special character inserted into crosslinked peptides
        /*        
        A=71.037114
        B=114.534940
        C=160.030649
        D=115.026943
        E=129.042593
        F=147.068414
        G=57.021464
        H=137.058912
        I=113.084064
        J=0.000000
        K=128.094963
        L=113.084064
        M=131.040485
        N=114.042927
        O=0.000000
        P=97.052764
        Q=128.058578
        R=156.101111
        S=87.032028
        T=101.047679
        U=150.953630
        V=99.068414
        W=186.079313
        X=111.000000
        Y=163.063329
        Z=128.550590
        Hydrogen=1.007825
        Carbon=12.000000
        Nitrogen=14.003074
        Oxygen=15.994915
        Electron=0.000549
        C_term=17.002740
        N_term=1.007825
        */
        
    }
    
    private void addFixedModsToAAMap(String aa, double modMass){
        
        double mass = aaMap.get(aa);
        double newMass = mass+modMass;
        aaMap.put(aa,newMass);
        System.out.println("Changed mass of " + aa +" to " + newMass);
    }
    
    private double getPepMr(String peptidePlusModString){
        double mr = 0.0;
        //double hmass = 1.007825;
        
        peptidePlusModString = peptidePlusModString.toUpperCase();
        peptidePlusModString = peptidePlusModString.trim();
        System.out.print("Peptide:" + peptidePlusModString);
        
        String[] temp = peptidePlusModString.split("_");
        String peptide = "[" + temp[0] + "]";
        String modString = temp[1];
        
        String[] aas = peptide.split("");
        //double pepMass = 0.0;
        for(String aa : aas){
            if(aa.length()==1){  
                if(aaMap.get(aa) != null){
                    mr+= aaMap.get(aa);
                }
                else{
                    //System.out.print("i:" + aa);
                    //We need to ignore modifications
                }
            }            
        }     
        
        String[] mods = modString.split("");
        System.out.println("modString:" + modString);
        //Now loop through the var mod string
       for(String mod : mods){           
            if(!mod.equals("0") && !mod.equals("")){  
                String varMod = varModUniqueIDToString.get(Integer.parseInt(mod));
                System.out.print("mod" + mod + " varMmod:" + varMod);
                Double varModMass = varMods.get(varMod);
                System.out.println(" modmass:" + varModMass);
                mr+=varModMass;
            }            
        }  
        
        
        mr +=  2* HMASS + OMASS;
        System.out.println( " mr: " + mr); 
        
        return mr;
        
    }
    
   /*
    * Returns a hashmap of all peptide sequences (including potentially crosslinked), along with their mass
    * HashMap String = peptide sequence (plus special characters for crosslinker)
    * HashMap Double = corresponding mass
    * 
    */
    public HashMap<String,Double> getAllPeptideMasses(){
       
       HashMap<String,Double> pepMasses = new HashMap();
       
       /*All proteinIDs is an arraylist of all proteins seen (tested with only one protein so far
       * Use this to retrieve each array of peptides in turn
       * 
       * Code as written does not search for crosslinked peptides in different proteins, only for peptides within one protein
       * 
       */
       for(String protAcc : allProteinIDs){
           
           //System.out.println("Processing: " + protAcc);           
           String[] peptides = accToPeptides.get(protAcc);
           
           
           for(int i=0;i<peptides.length;i++){
               String tempPep = peptides[i]; //Add characters for termini       
               tempPep = makeUnmodifiedModString(tempPep);  //Add the empty mod string  
               
               Double tempPepMass = getPepMr(tempPep);

               pepMasses.put(tempPep, tempPepMass);     //Add unmodified peptide to the hash
               
                //Handle var mods               
               for(String varMod : varMods.keySet()){
                   
                   Double varModMass = varMods.get(varMod);
                   
                   
                   Set<Integer> modPositions = new HashSet();
                   
                   int pos = tempPep.indexOf(varMod);
                   while(pos !=-1){                    
                       modPositions.add(pos);
                       pos = tempPep.indexOf(varMod,pos+1);  
                   }                
               
              
                   //Generate x * y combinations of var mods - TODO code needs checking it is correct
                   //If a mod is inserted in the peptide string, it is prefixed with "+" so it can be easily recognized later on
                   //TODO - This code is not correct; if there are > 2 residues that could be modified, it does not capture all combintations...
                   String[] temp = tempPep.split("_");                  
                   

                   ArrayList<String> modPeptides = new ArrayList();
                    for (Set<Integer> modSet : powerSet(modPositions)) {
                        if(modSet.size()>0){ //not interested in empty set
                            //System.out.println("tempPep" + modSet);
                            
                            
                            //tempPep = temp[0] + "_" +insertMod(temp[1],varMod,modSet);
                            modPeptides= insertMod(temp[0],temp[1],varMod,modSet);
                        }
                    }
                    
                    for(String modPep : modPeptides){
                       pepMasses.put(modPep, getPepMr(modPep)); 
                    }
                   
                    /* Old code now not used, probably can be deleted after testing
                    ArrayList<String> allModStrings = getModString(temp[0]);                      
                   for(int x = 0; x <modPositions.size();x++){
                       int modPos = modPositions.get(x);
                       
                       String[] temp = tempPep.split("_");
                       tempPep = temp[0] + "_" +insertMod(temp[1],varMod,modPos);
                       String modTempPep = tempPep;
                       Double modTempPepMass = getPepMr(modTempPep); 
                       pepMasses.put(modTempPep, modTempPepMass);        
                       
                       for(int y = x+1; y <modPositions.size();y++){
                           modPos = modPositions.get(y);
                           modTempPep = tempPep.subSequence(0, modPos+1) + "+"+varModMass + tempPep.subSequence(modPos+1,tempPep.length());
                           modTempPepMass += varModMass; 
                           pepMasses.put(modTempPep, modTempPepMass);
                       }
                   }                   
                   */ 
               }
           }
           
       }
       
        /* TODO
* This code adds the crosslinker as a variable modification along the peptide sequence
* This code is not correct for several reasons:
* - May be other variable modifications, we wish to apply e.g. Oxidation on M
* - This only linearly adds the var mod - in fact, it should (ultimately) be added in all combintatorial positions - not strictly needed until we look at MS2 data
* 
* Fundamental error:
* - In the dss example, if K is at position (say), peptide A 4, 5 and 8  and peptide B at 6 and 9, following logic should apply:
*   - dss can also be a dead-end variable modification at positions not-crosslinked
*   - current implementation of code, can apply the dss reagent as an additional mass on both peptides at the same site i.e. add the mass twice. This is not correct.
*   - the mass should only be added once.
* 
* 
* 
*/
       
       //Now this code only adds potentially crosslinked peptides to the hashmap if they both contain the crosslinked residue (unmodifed)
       Object[] peptides = pepMasses.keySet().toArray();
       for(int i=0;i<peptides.length;i++){
           
            String outerPepPlusModString = (String)peptides[i];
            String[] temp = outerPepPlusModString.split("_");
            String outerPep = temp[0];
            String outerPepModString = temp[1];
            Double outerPepMass = pepMasses.get(outerPep);
        
            System.out.println("outerpepmass:" + outerPepMass);
            
            ArrayList<Integer> xlOuterPositions = new ArrayList();
                   
            int pos = outerPep.indexOf(crossLinkedResidue);
            System.out.println("Processing outerPep:" + outerPep + " modString:" + outerPepModString);
 
            
            while(pos !=-1){
                
                //System.out.println("Pos:" + pos);
                
                if(pos != outerPep.length()-1){     
                        //This code serves two purposes:
                            //1) precludes crosslinks on C-terminal amino acid - assumption that cleavage and crosslinking cannot occur on same amino acid - TODO check this!
                            //2) then allows the check for the position not being modified
                     xlOuterPositions.add(pos);
                     //TODO insert code to check whether this position is also modified

                }
                pos = outerPep.indexOf(crossLinkedResidue,pos+1);
            }                

            //Generate x * y combinations of potential crosslinked sites - TODO code needs checking it is correct
            //If a mod is inserted in the peptide string, it is prefixed with "-" so it can be easily recognized later on
            for(int x = 0; x <xlOuterPositions.size();x++){
               int xlPos = xlOuterPositions.get(x);
               String newOuterPep = outerPep.subSequence(0, xlPos+1) + crossLinkerName + outerPep.subSequence(xlPos+1,outerPep.length());                
               
               for(int j=i+1;j<peptides.length;j++){
                   String innerPep = (String)peptides[j];
                   
                   ArrayList<Integer> xlInnerPositions = new ArrayList();
                   
                    pos = innerPep.indexOf(crossLinkedResidue);
                    while(pos !=-1){
                        if(pos != innerPep.length()-1){ 
                            xlInnerPositions.add(pos);
                            //TODO insert code to check whether this position is also modified

                        }
                        pos = innerPep.indexOf(crossLinkedResidue,pos+1);
                    }                
                   
                    for(int y = 0; y <xlInnerPositions.size();y++){
                       int xlInnerPos = xlInnerPositions.get(y);
                       String newInnerPep = innerPep.subSequence(0, xlInnerPos+1) + crossLinkerName + innerPep.subSequence(xlInnerPos+1,innerPep.length());
                       Double innerPepMass = getPepMr(innerPep);
                       System.out.println("innerpepmass:" + innerPepMass);
                       String crosslinkedPep = newOuterPep + "----" + newInnerPep;
                       Double crosslinkedMass = outerPepMass + innerPepMass + crossLinkerMass;  //TODO - this should only be valid IF the two peptides contain the crosslinked residue
                       pepMasses.put(crosslinkedPep, crosslinkedMass);
                   }
               }
               
           }
       }

       
       return pepMasses;
       
   }
    
    //Powerset code from http://stackoverflow.com/questions/1670862/obtaining-powerset-of-a-set-in-java
    public static <T> Set<Set<T>> powerSet(Set<T> originalSet) {
        Set<Set<T>> sets = new HashSet<Set<T>>();
        if (originalSet.isEmpty()) {
            sets.add(new HashSet<T>());
            return sets;
        }
        List<T> list = new ArrayList<T>(originalSet);
        T head = list.get(0);
        Set<T> rest = new HashSet<T>(list.subList(1, list.size())); 
        for (Set<T> set : powerSet(rest)) {
            Set<T> newSet = new HashSet<T>();
            newSet.add(head);
            newSet.addAll(set);
            sets.add(newSet);
            sets.add(set);
        }		
        return sets;
    }
  
    /*
    * Helper method to turn PEPTIDER into PEPTIDER_0000000000 where each 0 corresponds to unmodified n-terminus (pos 0), unmodified each amino acid and unmodified N-terminuis
    */
   private String makeUnmodifiedModString(String peptide){
       
       String peptideWithModString = peptide+"_";
       
       for(int i=0; i<peptide.length()+2;i++){
           peptideWithModString+= 0;
       }
       
       return peptideWithModString;
   }
   
   /*
    * Helper method to insert a mod into a modstring i.e. PMER_000000 to PMER_001000 for Oxidation on M, where oxidation on M has ID = 1
    */
   private ArrayList<String> insertMod(String peptide, String modString, String varMod, Set<Integer> modSet){
       
       ArrayList<String> peptides = new ArrayList();
       
       System.out.print("Mod string: " +  modString + " for peptide: " + peptide + " modSet:" + modSet);
       String newPeptidePlusModString = "";
       
       for(int pos : modSet){
           pos++;       //need to increment the position, since the modString has N-terminal mod at the beginning
           int modID = varModStringToUniqueID.get(varMod);        
           modString = modString.subSequence(0, pos-1) + ""+modID + modString.subSequence(pos+1,modString.length());
                        
       }
       newPeptidePlusModString += peptide + "_" + modString;
       System.out.println(" to: " +  newPeptidePlusModString);
       peptides.add(newPeptidePlusModString);
       
       return peptides;
   }
    
}
