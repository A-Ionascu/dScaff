

###############################################################################

# This script is intended to work within a bash pipeline 
#   with the purpose of identifying genes from the reference genome
#   that are more than 15000 nucleotides apart within a unique scaffold

# The script should be placed in a new directory
#   that contains only the .csv files to be analysed.

# The script creates a directory for results
#   in which 2 .csv output files are created: unfiltered and filtered

# Call the script in bash with the commands

#   chmod +x distante_genice_10apr2024.R
#   and
#   Rscript distante_genice_10apr2024.R



# List packages to check installation
packages <- c("dplyr","readr","utils")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages]) }
# Packages loading
suppressWarnings(suppressMessages(invisible(lapply(packages, library, character.only = TRUE))))

setwd(paste(getwd(),sep=""))



# read the files in input directory
all.files <- list.files(getwd())
tsv.files <- grep(".tsv",all.files,value=T)


for(f in tsv.files){


  # import chromosome table
  
  #chr_full <- read.csv(f, header = FALSE)
  
  
  chr_full <- readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE)
  
  
  chr_full <- as.data.frame(chr_full)
  
  
  # make selection
  chr <- chr_full[ , 1:5]

  # write columns name
  colnames(chr) <- c("scaff", "start", "stop", "chr", "strand")

  # add empty columns
  chr[,6:11] <- rep(NA,nrow(chr))
  colnames(chr)[6:11] <- c("g1", "g2", "result","locus","gene","full_id")

  chr$result[1] <- paste("good")
  chr$locus[1] <- paste( chr_full$X7[1] )
  chr$gene[1] <- paste( chr_full$X6[1] )
  chr$full_id[1] <- paste(chr_full$X1[1],":",chr_full$X2[1],"-",chr_full$X3[1],
                          sep="")

  i <- 1 # pentru stop
  j <- 2 # pentru start
  
  while(i < nrow(chr)){ 
    
    if(chr$scaff[i] != chr$scaff[j]){
      
      chr$result[j] <- paste("good")
      chr$locus[j] <- paste( chr_full$X7[j] )
      chr$gene[j] <- paste( chr_full$X6[j] )
      chr$full_id[j] <- paste(chr_full$X1[j],":",chr_full$X2[j],"-",chr_full$X3[j],
                              sep="")
      i <- j
      j <- j+1
      next }
    
    else{
      
      distance <- abs(chr$stop[i] - chr$start[j]) # calculate distance
    
    
      if( distance <= 15000 ){ # distance condition bad
        j <- j + 1 # jump to next gene to compare
        
        if(j > nrow(chr)){ # do not go over table end
          break }
        else{next} } # restart while loop
  
      else { # distance condition good
        chr$result[j] <- paste("good") # write info
        chr$g1[j] <- paste(i) # write row of gene 1 (stop)
        chr$g2[j] <- paste(j) # write row of gene 2 (start)
        chr$locus[j] <- paste(chr_full$X7[j]) # write locus name
        chr$gene[j] <- paste( chr_full$X6[j] ) # write gene name
        chr$full_id[j] <- paste(chr_full$X1[j],":",chr_full$X2[j],"-",chr_full$X3[j],
                                sep="") # write full gene id from ref
    
        i <- j # jump to the gene that satisfied condition
        j <- i+1 # compare to next gene
    
        if(j > nrow(chr)){ # do not go over table end
          break }
        else{next} # restart while loop
    
      }
    }
  }

  # filter only the good results 
  chr_result <- chr %>%
    filter(result == "good")

  # create a directory for results
  mainDir <- getwd()
  
  chrDir <- gsub(".tsv", "_chr", f)
  #subDir <- "results_gene_interes"

  ifelse(!dir.exists(file.path(mainDir, chrDir)), 
          dir.create(file.path(mainDir, chrDir)), FALSE)
  
  # output to results directory
  setwd(file.path(mainDir, chrDir))

  # output unfiltered table (chr)
  unfiltered_output <- gsub( ".tsv", "_hits_distances_unfiltered.csv", f)
  write.csv(chr, paste(unfiltered_output))


  # output filtered table (chr_results)
  filtered_output <- gsub( ".tsv", "_hits_distances_filtered.csv", f)
  write.csv(chr_result, paste(filtered_output))

  # go back to input file
  setwd(mainDir)

}






