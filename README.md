# dScaff
Scaffolding strategy for draft assemblies (contigs) based on reference assembly using either genes queries or ranked queries.

<div align="justify">

#
### Strategy:

Digital Scaffolding (**dScaff**) aims improve new draft assemblies of organisms that have an NCBI reference assembly. This strategy is most valuable when performing *de novo* assemblies of local species. Two strategies are available, one using the gene annotations for the reference assembly and one using consecutive reference fragments obtained with the SubSequencesExtractor complementary script.

#
#
### Dependencies:

**dScaff** was implemented for use on Debian-based systems. Before running the main script:

+ Download the files in this repository.
+ Run the requirements.sh script in order to download R, seqtk and ncbi-blast+ packages. 
+ Open a bash terminal, enter R environment with the *R* command. Inside the R environment run the *install.packges("dplyr")* command and type *yes* for each threshold. This step was required during testing on a Linux Mint 21.3 system, for locally installing R packages. Other operating systems might not require this step.
+ Make sure the main script is executable.
+ Keep the *dScaff.sh* , *query_filtering.R* and *contigs_mapping.R* files in the same directory.

#
#
### Running dScaff:

Information about **dScaff** can be obtained by running the script with the *-h* or *--help* command. The script should be run from its main directory with either *bash* , *sh* or *./* prefixes. 

**dScaff** takes as arguments:  
-a, --assembly            Draft assembly in FASTA format\
-q, --query               gene.fna file with gene sequences from reference genome or ranked queries.fasta output of ranked queries SubSequencesExtractor.sh script\
-d, --dataset             ncbi_dataset.tsv containing all genes in reference genome\
-gq, --gene_queries       perform gene queries strategy\
-rq, --ranked_queries     perform ranked queries strategy


It is recommended to use the absolute paths for each file. The script will only work if all four arguments are provided. The file names should not contain dots in the name.

#
#
### Contact and feedback:

For offering feedback or any type of inquires about the application, please contact us at **a.ionascu20@s.bio.unibuc.com**.  


</div>


