!/bin/bash

## Script to remove comma from MGF files
## Change the path of the folder and run the script

FILES=/Users/riteshk/Ritesh_Work/TestSpace/tmp/*.mgf

for file in $FILES
    do
      tr , ' ' < $file 	> /tmp/tempfile.tmp
      mv /tmp/tempfile.tmp $file
      echo "Modified: " $file
    done

echo " *** Done! *** "