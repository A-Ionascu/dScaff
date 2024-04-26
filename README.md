# dScaff
Scaffolding strategy for draft assembled contigs based on reference genes annotations

<div align="justify">

#
#
### Strategy:

Digital Scaffolding (**dScaff**) aims improve new draft assemblies of organisms that have a better reference assebly in with annotated genes in various databases. This strategy is most valuable when performing *de novo* assemblies of local species. The inplemented strategy aims to use the draft assembly, sequences of all genes from the reference genome and a table of their annotation in order to enhance scaffolding posibilities.

**dScaff** recognizez the chromosomes associated to annotated genes in the reference genome and creates separate working directories. The script makes a subselection of genes that are more than 15k nt apart in the reference genome and performs BLAST between the sequences of these genes and the draft assembly. The BLAST results for each chromosome or for each scaffold within a chromosome in order to produce an indexed table and a map of draft assembly contigs in reference genome order. 



#
#
### Dependencies:

**dScaff** was implemented for use on Debian-based systems. Before running the main script:

+ Run the requirements.sh script in order to download R, seqtk and ncbi-blast+ packages. 
+ Open a bash terminal, enter R environment with the *R* command. Inside the R environment run the *install.packges("dplyr")* command and type *yes* for each threshold. This step was required during testing on a Linux Mint 21.3 system, for locally installing R packages. Other operating systems might not require this step.
+ Make sure the main script is executable.

#
#
### Running dScaff:

Information about **dScaff** can be obtained by running the script with the *-h* or *--help* command. The script should be run from its main directory with either *bash* , *sh* or *./* prefixes. 

**dScaff** takes as arguments:  
-a, --assembly    Draft assembly in FASTA format
-g, --genes       gene.fna file with gene sequences from reference genome
-d, --dataset     ncbi_dataset.tsv containing all genes in reference genome

It is recommended to use the absolute paths for each file. The script will only work if all three arguments are provided.

#
#
### Contact and feedback:

For offering feedback or any type of inquires about the application, please contact us at **a.ionascu20@s.bio.unibuc.com**.  


</div>
