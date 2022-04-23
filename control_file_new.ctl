### This is an example of control file. GriffinDetector v1.1

SPECIES
# this part include two columns which is split by ":", the first column is species name, the second column is the annotated protein sequences in fasta format. the path of the protein sequences should be the absolute path
jap:/home/User_ID_Path/Data_Path/O.sativa.ssp.japonica/O.sativa.ssp.japonica.maker.proteins.fasta
indica:/home/User_ID_Path/Data_Path/O.sativa.ssp.indica/O.sativa.ssp.indica.maker.proteins.fasta
rufi:/home/User_ID_Path/Data_Path/O.rufipogon/O.rufipogon.maker.protein.fasta 
glab:/home/User_ID_Path/Data_Path/O.glaberrima/O.glaberrima.maker.proteins.fasta
barthii:/home/User_ID_Path/Data_Path/MAKER_annotations_O.barthii/O.barthii.maker.proteins.fasta
longi:/home/User_ID_Path/Data_Path/O.longistaminta/O.longistaminta.maker.proteins.fasta 
punctata:/home/User_ID_Path/Data_Path/O.punctata/O.punctata.maker.proteins.fasta 
brachy:/home/User_ID_Path/Data_Path/O.brachyantha/O.brachyantha.maker.proteins.fasta 
leersia:/home/User_ID_Path/Data_Path/Leersia/leersia.maker.proteins.fasta
tair:/home/User_ID_Path/data-TAIR/TAIR10_pep 
SPECIES

SPEGFFS
jap:/home/User_ID_Path/Data_Path/O.sativa.ssp.japonica/O.sativa.ssp.japonica.gff
indica:/home/User_ID_Path/Data_Path/O.sativa.ssp.indica/O.sativa.ssp.indica.gff
rufi:/home/User_ID_Path/Data_Path/O.rufipogon/O.rufipogon.maker.gff
glab:/home/User_ID_Path/Data_Path/O.glaberrima/O.glaberrima.maker.gff
barthii:/home/User_ID_Path/Data_Path/MAKER_annotations_O.barthii/O.barthii.maker.gff
longi:/home/User_ID_Path/Data_Path/O.longistaminta/O.longistaminta.maker.gff
punctata:/home/User_ID_Path/Data_Path/O.punctata/O.punctata.maker.gff
brachy:/home/User_ID_Path/Data_Path/O.brachyantha/O.brachyantha.maker.gff
leersia:/home/User_ID_Path/Data_Path/Leersia/leersia.maker.gff
tair:/home/User_ID_Path/data-TAIR/TAIR10_GFF3_genes.gff 
SPEGFFS


TREE:((((((((jap,indica),rufi),(glab,barthii)),longi),punctata),brachy),leersia),tair);
#this is the phylogeny tree among these species


INGROUPS:jap,glab,barthii,indica
#these are the species that you are focus on, it can be more than one in total, and cannot also be the outgroup, and finally, the chimeric genes will be given out based on these species

OUTGROUPS:tair,leersia,brachy
#these species are the outgroup,they are fixed.

Makeblastdb:/home/User_ID_Path/software/ncbi-blast-2.2.28+/bin/makeblastdb
Blastp:/home/User_ID_Path/software/ncbi-blast-2.2.28+/bin/blastp 

SKIP:Y
#skip the data prepare steps, if it was set as Y, then, the OUTFILE option must be given.

OUTFILE:blast_plus_db
