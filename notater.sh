nextclade dataset get --name  community/masphl-bioinformatics/hav/vp1-2b-junction --output-dir data/vp1-2b-junction
nextclade dataset get --name  community/masphl-bioinformatics/hav/whole-genome --output-dir data/whole-genome 

nextclade run --include-reference --input-dataset data/vp1-2b-junction/ --output-all output_vp1-2b-junction/ data/Sekvens* 
# sliter med IIIa

nextclade run --include-reference --input-dataset data/whole-genome/ --output-all output_whole-genome/ data/Sekvens* 
# for mange blir Ib

nextclade run data/Sekvens* --input-dataset data/3a_AJ299464.3/ --output-all tmp
# fungerer for IIIa, men ikke for alle genotyper
