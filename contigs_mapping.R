

# List packages to check installation
packages <- c("dplyr","readr","utils","tibble")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
suppressWarnings(suppressMessages(invisible(lapply(packages, library, character.only = TRUE))))



mainDir <- (paste(getwd(),sep=""))

setwd(mainDir)
directories <- list.dirs(recursive=FALSE)

for(d in directories){
  dir <- gsub("./","",d)
  path <- file.path(mainDir,dir,"genes")
  setwd(path)
  

  all.files <- list.files(getwd())
  csv.files <- grep(".csv",all.files,value=T)
  
  
  contigs_of_interest <- data.frame(matrix(NA,nrow=0,ncol=13))
  colnames(contigs_of_interest) <- c("query_name", "subject_name", "identity", "length",
                                     "mismatch", "gap", "query_start", "query_end", 
                                     "subject_start", "subject_end", "E_value", 
                                     "bit_score","gene_length")
  
  for(f in csv.files){
    
    data <- readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE)
    
    data <- as.data.frame(data)
    
    colnames(data) <- c("query_name", "subject_name", "identity", "length",
                        "mismatch", "gap", "query_start", "query_end", 
                        "subject_start", "subject_end", "E_value", 
                        "bit_score","gene_length")



  if(data$gene_length[1] <= 2500){
    next
  }
  
  else{
    
    if(any(data$length >= data$gene_length[1] * 0.5)){
      hgl <- data$gene_length[1] * 0.5
      aliniament <- data[which(data$length >= hgl),]
      contigs_of_interest <- rbind(contigs_of_interest, aliniament)  
      contigs_of_interest <- contigs_of_interest[!duplicated(contigs_of_interest),]
      next }
    
    else{
      if(any(data$length >= data$gene_length[1] * 0.4)){
        hgl <- data$gene_length[1] * 0.4
        aliniament <- data[which(data$length >= hgl),]
        contigs_of_interest <- rbind(contigs_of_interest, aliniament)
        contigs_of_interest <- contigs_of_interest[!duplicated(contigs_of_interest),]
        next }
      
      else{
        if(any(data$length >= data$gene_length[1] * 0.3)){
          hgl <- data$gene_length[1] * 0.3
          aliniament <- data[which(data$length >= hgl),]
          contigs_of_interest <- rbind(contigs_of_interest, aliniament)  
          contigs_of_interest <- contigs_of_interest[!duplicated(contigs_of_interest),]
          next }
        
        #  else{
        #    if(any(tmp_ctg$length >= data$gene_length[1] * 0.25)){
        #      hgl <- data$gene_length[1] * 0.25
        #      aliniament <- data[which(data$length >= hgl),]
        #      contigs_of_interest <- rbind(contigs_of_interest, aliniament)  }
        #  }
      }
    }
    
  }
    
  }
  

  setwd(file.path(mainDir,dir))

  genes_of_interest <- read.csv("genes_filtered.csv")

  contigs_of_interest <- add_column(contigs_of_interest, genomic_start=NA, .before="E_value")
  contigs_of_interest <- add_column(contigs_of_interest, genomic_end=NA, .before="E_value")
  contigs_of_interest <- add_column(contigs_of_interest, ref_scaff=NA, .after="gene_length")
  contigs_of_interest <- add_column(contigs_of_interest, scaff_start=NA, .after="ref_scaff")
  contigs_of_interest <- add_column(contigs_of_interest, scaff_stop=NA, .after="scaff_start")
  
  for(i in seq(1,nrow(contigs_of_interest),1)){
  
  contigs_of_interest$genomic_start[i] <- genes_of_interest$start[which(genes_of_interest$full_id == contigs_of_interest$query_name[i])]
  contigs_of_interest$genomic_end[i] <- genes_of_interest$stop[which(genes_of_interest$full_id == contigs_of_interest$query_name[i])]
  contigs_of_interest$ref_scaff[i] <- genes_of_interest$scaff[which(genes_of_interest$full_id == contigs_of_interest$query_name[i])]
  contigs_of_interest$scaff_start[i] <- genes_of_interest$start[which(genes_of_interest$full_id == contigs_of_interest$query_name[i])]
  contigs_of_interest$scaff_stop[i] <- genes_of_interest$stop[which(genes_of_interest$full_id == contigs_of_interest$query_name[i])]
  
  }

  contigs_of_interest <- arrange(contigs_of_interest, genomic_start)


  write.csv(contigs_of_interest, "contigs_of_interest_all_scaffolds.csv")

################


  
  if(length(unique(contigs_of_interest$ref_scaff)) > 1){
    
    setwd(file.path(mainDir,dir))

    for(s in unique(contigs_of_interest$ref_scaff)){
     
      ifelse(!dir.exists(file.path(mainDir, dir, "scaffolds")), 
            dir.create(file.path(mainDir, dir, "scaffolds")), FALSE)
      setwd(file.path(mainDir,dir,"scaffolds"))

      ifelse(!dir.exists(file.path(mainDir, dir,"scaffolds", s)), 
             dir.create(file.path(mainDir, dir,"scaffolds", s)), FALSE)
      setwd(file.path(mainDir,dir,"scaffolds",s))
      
      scaff_int <- contigs_of_interest %>%
        filter(ref_scaff == s)
      
      map <- data.frame(matrix(NA,
                               nrow = length(unique(scaff_int$subject_name)),
                               ncol = length(unique(scaff_int$query_name))))
      
      colnames(map) <- unique(scaff_int$query_name)
      rownames(map) <- unique(scaff_int$subject_name)
      
      
      for(i in seq(1,ncol(map),1)){
        
        gene <- scaff_int %>%
          filter(scaff_int$query_name == colnames(map)[i])
        
        contiguri <- unique(gene$subject_name)
        
        for(j in contiguri){
          
          rand <- which(rownames(map) == j)
          map[rand,i] <- c("______")
        }
        
      }
      
      write.csv(map, paste(s,"mapped_contigs.csv"), na="")
      
      
      #############
      
      contig_intex <- data.frame(matrix(NA, nrow=nrow(map), ncol=7))  
      colnames(contig_intex) <- c("contigs","genes_hit",
                                  "contig_start","contig_end",
                                  "genomic_start","genomic_end","scaffold")
      contig_intex$contigs <- rownames(map)
      
      for(i in seq(1,nrow(map),1)){
        
        hits <- which(!is.na(map[i,]))
        contig_intex$genes_hit[i] <- length(hits)
        contig_intex$scaffold[i] <- colnames(map)[i]
        
        tmp_contig <- contigs_of_interest %>%
          filter(subject_name == rownames(map)[i])
        
        if(min(tmp_contig$subject_start) < min(tmp_contig$subject_end)){
          contig_intex$contig_start[i] <- min(tmp_contig$subject_start)
          contig_intex$contig_end[i] <- max(tmp_contig$subject_end) }
        else{
          contig_intex$contig_start[i] <- max(tmp_contig$subject_start)
          contig_intex$contig_end[i] <- min(tmp_contig$subject_end) }
        
        if(min(tmp_contig$genomic_start) < min(tmp_contig$genomic_end)){
          contig_intex$genomic_start[i] <- min(tmp_contig$genomic_start)
          contig_intex$genomic_end[i] <- max(tmp_contig$genomic_end) }
        else{
          contig_intex$genomic_start[i] <- max(tmp_contig$genomic_start)
          contig_intex$genomic_end[i] <- min(tmp_contig$genomic_end) }
      }
      
      write.csv(contig_intex, paste(s,"indexed_contigs.csv"))
      
    }

  }

  if(length(unique(contigs_of_interest$ref_scaff)) == 1){
    
    setwd(file.path(mainDir,dir))
    ifelse(!dir.exists(file.path(mainDir, dir, "chromosome")), 
            dir.create(file.path(mainDir, dir, "chromosome")), FALSE)
    setwd(file.path(mainDir,dir,"chromosome"))
    
    map <- data.frame(matrix(NA,
                             nrow = length(unique(contigs_of_interest$subject_name)),
                             ncol = length(unique(contigs_of_interest$query_name))))
    
    colnames(map) <- unique(contigs_of_interest$query_name)
    rownames(map) <- unique(contigs_of_interest$subject_name)
    
    
    for(i in seq(1,ncol(map),1)){
      
      gene <- contigs_of_interest %>%
        filter(contigs_of_interest$query_name == colnames(map)[i])
      
      contiguri <- unique(gene$subject_name)
      
      for(j in contiguri){
        
        rand <- which(rownames(map) == j)
        map[rand,i] <- c("______")
      }
      
    }
    
    write.csv(map, "mapped_contigs.csv", na="")
    
    
    #############
    
    contig_intex <- data.frame(matrix(NA, nrow=nrow(map), ncol=6))  
    colnames(contig_intex) <- c("contigs","genes_hit",
                                "contig_start","contig_end",
                                "genomic_start","genomic_end")
    contig_intex$contigs <- rownames(map)
    
    for(i in seq(1,nrow(map),1)){
      
      hits <- which(!is.na(map[i,]))
      contig_intex$genes_hit[i] <- length(hits)
      
      tmp_contig <- contigs_of_interest %>%
        filter(subject_name == rownames(map)[i])
      
      if(min(tmp_contig$subject_start) < min(tmp_contig$subject_end)){
        contig_intex$contig_start[i] <- min(tmp_contig$subject_start)
        contig_intex$contig_end[i] <- max(tmp_contig$subject_end) }
      else{
        contig_intex$contig_start[i] <- max(tmp_contig$subject_start)
        contig_intex$contig_end[i] <- min(tmp_contig$subject_end) }
      
      if(min(tmp_contig$genomic_start) < min(tmp_contig$genomic_end)){
        contig_intex$genomic_start[i] <- min(tmp_contig$genomic_start)
        contig_intex$genomic_end[i] <- max(tmp_contig$genomic_end) }
      else{
        contig_intex$genomic_start[i] <- max(tmp_contig$genomic_start)
        contig_intex$genomic_end[i] <- min(tmp_contig$genomic_end) }
    }
    
    write.csv(contig_intex, "indexed_contigs.csv")
    
    
  }
    
    
  }
  
  
  
  
  
  
  
  
############

