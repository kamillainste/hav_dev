nextclade dataset get --name  community/masphl-bioinformatics/hav/vp1-2b-junction --output-dir data/vp1-2b-junction
nextclade dataset get --name  community/masphl-bioinformatics/hav/whole-genome --output-dir data/whole-genome 

nextclade run --include-reference --input-dataset data/vp1-2b-junction/ --output-all output_vp1-2b-junction/ data/Sekvens* 
# sliter med IIIa

nextclade run --include-reference --input-dataset data/whole-genome/ --output-all output_whole-genome/ data/Sekvens* 
# for mange blir Ib

nextclade run data/Sekvens* --input-dataset data/3a_AJ299464.3/ --output-all tmp
# fungerer for IIIa, men ikke for alle genotyper



# Make the algorithm MORE tolerant of mismatches
nextclade run data/Sekvens* \
  --input-dataset data/vp1-2b-junction/ \
  --penalty-mismatch -1 \
  --score-match 2 \
  --output-all tmp

#Try these adjustments (from most to least permissive):
    #Most permissive: --penalty-mismatch -1 --score-match 2
    #Moderate: --penalty-mismatch -2 --score-match 3
    #Stricter: --penalty-mismatch -5 --score-match 5

# ikke testet noe av dette enda


nextclade run data/Sekvens* \
  --input-dataset data/vp1-2b-junction/ \
  --alignment-preset high-diversity \
  --output-all output_vp1-2b-junction_high-diversity/

#Fungerer!
#Jobber videre med nextclade_vp1_2b_alginment.sh


nextclade_vp1_2b_alignment.sh ./data/
# videre arbeid: få result.tsv til å ha heading "strain" istedenfor "seqName", 
  # og fjern index kolonnen slik at strain blir første kolonne
  # dette for å kunne legge den over treet i auspice og bruke til å fargelegge
  
