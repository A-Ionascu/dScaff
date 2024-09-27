#!/bin/bash

# Input files: 
# -a (assembly) = draft assembly in FASTA format
# -g (genes) = gene.fna file with gene sequences from reference genome
# -d (dataset) = ncbi_dataset.tsv containing all genes in reference genome


# Define the help method
get_help() {
  echo " "
  echo "Usage: $(basename $0) [-h] [-a FILE] [-g FILE] [-d FILE]"
  echo " "
  echo "Options:"
  echo "  -h, --help        Display this help message"
  echo "  -a, --assembly    Draft assembly in FASTA format"
  echo "  -g, --genes       gene.fna file with gene sequences from reference genome"
  echo "  -d, --dataset     ncbi_dataset.tsv containing all genes in reference genome"
  echo " "
  echo "Use files from working directory or provide absolute path to input files."
  echo " "
  exit 0
}

# Check if no arguments provided
if [ $# -eq 0 ]; then
  echo " "
  echo "Error: No arguments provided."
  echo " "
  get_help
fi

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  get_help
fi

while [ ! -z "$1" ]; do
  case "$1" in
     --assembly|-a)
         shift
         echo " "
         echo "Input draft assembly is: $1"
         #echo " "
         assembly=$1
         if [ -z "$1" ]; then
            echo "Error: Missing argument for draft assembly file."
            echo " "
            exit 1
         fi         
         ;;
     --genes|-g)
         shift
         echo " "
         echo "Input reference gene sequences is: $1"
         #echo " "
         genes=$1
         if [ -z "$1" ]; then
            echo "Error: Missing argument for reference genes sequences."
            echo " "
            exit 1
         fi  
         ;;
     --dataset|-d)
        shift
        echo " "
        echo "Input genes dataset is: $1"
        #echo " "
        dataset=$1 
        if [ -z "$1" ]; then
            echo "Error: Missing argument for genes dataset."
            echo " "
            exit 1
         fi         
        ;;
     *)
        get_help
        ;;
  esac
shift
done

if [ "$assembly" == "" ]; then
    echo "Error: Missing input draft assembly file"
    echo " "
    exit 1
fi
if [ "$genes" == "" ]; then
    echo "Error: Missing input reference genes file"
    echo " "
    exit 1
fi
if [ "$dataset" == "" ]; then
    echo "Error: Missing input genes dataset file"
    echo " "
    exit 1
fi


#while getopts a:g:d: flag
#do
#    case "${flag}" in
#        a) assembly=${OPTARG};;
#        g) genes=${OPTARG};;
#        d) dataset=${OPTARG};;
#    esac
#done


START=$(date +%s)

###
# First part --> filter genes by chromosomes
###
CWD="$(pwd)"

echo " "
echo "Identifing chromosomes ..."
echo " "
awk -F "\t" 'NR>1{print $4}' $dataset | sort | uniq > chromosomes.txt

sed -i '/^$/d' chromosomes.txt

subdir="$(basename -- ${assembly%%.*})"
mkdir $subdir

echo "Selecting genes for each chromosome ..."
echo " "
cat chromosomes.txt | while read line 
do 
awk -v line=$line -F "\t" '{if ($4 == line) print $0}' $dataset > ./$subdir/$line.tsv
done

mv chromosomes.txt $subdir


echo "Filtering genes datasets ..."
echo " "
cp gene_filtering.R $subdir
cd $subdir
Rscript gene_filtering.R

 


###
# Second part --> blastn between each gene of interest and draft assembly
###

cd $CWD

echo "Preparing BLAST database ..."
mkdir assembly_database
cp $assembly assembly_database
cd assembly_database
makeblastdb -in *.fasta -dbtype nucl
cd $CWD/$subdir/

echo "Performing BLAST for selected genes ..."
echo " "
ls -d ./*/ | while read line
do 

#[[ ! -d "line" ]] && continue
 
cd $line

echo "Entered " $line

awk -F "," 'NR > 1 {print $2":"$3"-"$4}' *_filtered.csv > gene_headers.txt
sed -i 's/"//g' gene_headers.txt

cd $CWD/$subdir
grep ">" $genes | sed 's/>//g' > headers_list_db.txt

cd $line
cat gene_headers.txt | while read gene; do grep $gene ../headers_list_db.txt >> gene_headers_db.txt; done


seqtk subseq $genes gene_headers_db.txt > genes_of_interest.fasta

mkdir genes

awk -F ">| " '/^>/ {s=$2".gene.fasta"}; {print > s}' genes_of_interest.fasta

mv *.gene.fasta ./genes

threads=$(nproc --all)


for i in genes/*.gene.fasta; do blastn -db $CWD/assembly_database/*.fasta -query $i -out ${i%.gene.fasta}".lucru.csv" -outfmt 6 -num_threads $threads; done

for i in genes/*.gene.fasta; do lung=$( echo $i | sed -e '1d' $i | wc -c ); length=$( expr $lung - 1 ); awk -v Length=$length -F "\t" '{ FS = OFS = "\t" } {print $0,Length}' ${i%.gene.fasta}".lucru.csv" > ${i%.gene.fasta}".csv"; done

rm -r genes/*.gene.fasta genes/*.lucru.csv


cd $CWD/$subdir
rm headers_list_db.txt

done



###
# Third part --> map draft assembly contigs using the genes of interest
###
echo " "
echo "Indexing and mapping contigs ..."
ls -d ./*/ | while read line
do 

#[[ ! -d "line" ]] && continue
 
cd $CWD/$subdir/$line

#echo "Entered " $line

mv *_distances_filtered.csv genes_filtered.csv

cd genes
find . -name '*' -size 0 -print0 | xargs -0 rm 2> /dev/null


done

cd $CWD
cp contigs_mapping.R $subdir
cd $subdir

Rscript contigs_mapping.R 2> /dev/null

############################
cd $CWD
grep ">" $assembly > headers_assembly.txt
############################
cd $subdir

ls -d ./*/ | while read line
do
cd $line

mkdir tmp
mv *.csv tmp
mv *.txt tmp
mv *.fasta tmp

#######################################
if [ -d "chromosome" ]; then
cd chromosome/
##################################################################################################################
sed -i 's/"//g' *selected_contigs.csv 
sed -i 's/,/ /g' *selected_contigs.csv
#Pentru fiecare fisier care se termina cu selected_contigs.csv, il parcurg cu while, apoi extrag headerele de interes in baza fisierului cu all headers.
for i in *selected_contigs.csv; do cat $i | while read contig index; do grep $contig $CWD/headers_assembly.txt >> ${i%_selected_contigs.csv}".selected.headers.txt"; done; done

#Adaug o coloana care contine cromozomul si scaffoldul.
for i in *selected_contigs.csv; do awk -v loc="${i%_selected_contigs.csv}" '{print $0, loc}' $i >> contigs_chr_loc.txt; done

#Elimin ">".
sed -i 's/>//g' *selected.headers.txt

#Extrag contigurile, dar nu sunt in ordine.
for i in *.selected.headers.txt; do seqtk subseq $CWD/$assembly $i > ${i%.selected.headers.txt}".new1.assembly.fasta"; done

#Sparg fisierul fasta ce contine contiguri.
awk -F ">| " '/^>/ {s=$2".fna"}; {print > s}' *.new1.assembly.fasta

#Am ordonat contigurile
cat contigs_chr_loc.txt | while read contig loc chr; do for i in *.fna; do if [[ ${i%.fna} == $contig ]]; then cat $i >> $chr".dScaff.assembly1.fasta"; fi; done; done
rm -r *.fna
rm -r *.new1.assembly.fasta

for chr in *.dScaff.assembly1.fasta
do
awk '/^>/{gsub(/^>/,">"i++" ");}1' i=1 $chr > ${chr%.dScaff.assembly1.fasta}"_dScaff_assembly.fasta"
sed -i "s/>/>${chr%.dScaff.assembly1.fasta} /g" ${chr%.dScaff.assembly1.fasta}"_dScaff_assembly.fasta"
done
rm -r contigs_chr_loc.txt
rm -r *.dScaff.assembly1.fasta
##################################################################################################################

fi
#######################################
#######################################
if [ -d "scaffolds" ]; then
cd scaffolds
ls -d ./*/ | while read scaff
do
cd $scaff

##################################################################################################################
sed -i 's/"//g' *selected_contigs.csv 
sed -i 's/,/ /g' *selected_contigs.csv
#Pentru fiecare fisier care se termina cu selected_contigs.csv, il parcurg cu while, apoi extrag headerele de interes in baza fisierului cu all headers.
for i in *selected_contigs.csv; do cat $i | while read contig index; do grep $contig $CWD/headers_assembly.txt >> ${i%_selected_contigs.csv}".selected.headers.txt"; done; done

#Adaug o coloana care contine cromozomul si scaffoldul.
for i in *selected_contigs.csv; do awk -v loc="${i%_selected_contigs.csv}" '{print $0, loc}' $i >> contigs_chr_loc.txt; done

#Elimin ">".
sed -i 's/>//g' *selected.headers.txt

#Extrag contigurile, dar nu sunt in ordine.
for i in *.selected.headers.txt; do seqtk subseq $CWD/$assembly $i > ${i%.selected.headers.txt}".new1.assembly.fasta"; done

#Sparg fisierul fasta ce contine contiguri.
awk -F ">| " '/^>/ {s=$2".fna"}; {print > s}' *.new1.assembly.fasta

#Am ordonat contigurile
cat contigs_chr_loc.txt | while read contig loc chr; do for i in *.fna; do if [[ ${i%.fna} == $contig ]]; then cat $i >> $chr".dScaff.assembly1.fasta"; fi; done; done
rm -r *.fna
rm -r *.new1.assembly.fasta

for chr in *.dScaff.assembly1.fasta
do
awk '/^>/{gsub(/^>/,">"i++" ");}1' i=1 $chr > ${chr%.dScaff.assembly1.fasta}"_dScaff_assembly.fasta"
sed -i "s/>/>${chr%.dScaff.assembly1.fasta} /g" ${chr%.dScaff.assembly1.fasta}"_dScaff_assembly.fasta"
done
rm -r contigs_chr_loc.txt
rm -r *.dScaff.assembly1.fasta
##################################################################################################################
cd -
done
fi
#######################################
cd $CWD/$subdir


done

rm gene_filtering.R
rm contigs_mapping.R
rm chromosomes.txt

cd $CWD
mv assembly_database $subdir

echo " "
echo "Finished !"
### End

END=$(date +%s)
DIFF=$(( $END - $START ))
time=$(expr $DIFF / 60)
echo "dScaff ran for $DIFF seconds"
echo "That's around $time minutes"

# CWD este dScaff folder
# subdir este folder cu numele asamblarii
# line este folder cu numele chromozomului

