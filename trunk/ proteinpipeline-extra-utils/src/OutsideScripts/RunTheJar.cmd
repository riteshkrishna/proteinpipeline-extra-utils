REM This calls the multipleSearch routine for Mascot produced mzIdentMLs.
REM The original multiple search code is modified to produce peptide seq as well and runs only on a SINGLE search engine
REM The decoy is hard coded to be "Rnd" 
REM copy the jar and the cmd file in the input dir and run from there

for %%f in (*.mzid) do java -jar multipleSearchForMascot.jar C:\Ritesh_Work\From_Beta_UPenn\20110405_MascotResults_converted\%%f mascot C:\Ritesh_Work\From_Beta_UPenn\output\%%f.txt 3

pause