nextclade dataset get --name  community/masphl-bioinformatics/hav/vp1-2b-junction --output-dir data/vp1-2b-junction
nextclade dataset get --name  community/masphl-bioinformatics/hav/whole-genome --output-dir data/whole-genome 

nextclade run --include-reference --input-dataset data/vp1-2b-junction/ --output-all output_vp1-2b-junction/ data/Sekvens* 
# sliter med IIIa

nextclade run --include-reference --input-dataset data/whole-genome/ --output-all output_whole-genome/ data/Sekvens* 
# for mange blir Ib

nextclade run data/Sekvens* --input-dataset data/3a_AJ299464.3/ --output-all tmp
# fungerer for IIIa, men ikke for alle genotyper


nextclade run data/Sekvens* \
  --input-dataset data/vp1-2b-junction/ \
  --alignment-preset high-diversity \
  --output-all output_vp1-2b-junction_high-diversity/

#Fungerer!
#Jobber videre med nextclade_vp1_2b_alginment.sh


./nextclade_vp1_2b_alignment.sh ./data/
# videre arbeid: få result.tsv til å ha heading "strain" istedenfor "seqName", 
  # og fjern index kolonnen slik at strain blir første kolonne
  # dette for å kunne legge den over treet i auspice og bruke til å fargelegge
  # kan vente med dette til de genotype-spesifikke analysene
  
  # ved å velge node-type i auspice får man alle sekvensene man la til i analysen markert
  # dersom vi ikke har andre sekvenser som også legges til i analysen er jo det dekkende for 
    ## å se sine prøver, men er jo fint å få med litt mer metadata også

# et problem med å bruke dette data-settet er at referansen er kortere enn mange av sekvensene våre, så det blir mye "missing data" i analysen. 
  # men hvis vi bare skal bruke det til å bestemme genotype, så gjør ikke det noe




nextclade run data/Sekvens* \
  --input-dataset data/IIIA_AJ299464.3/ \
  --alignment-preset high-diversity \
  --output-all output_IIIA_AJ299464.3_high-diversity/
  