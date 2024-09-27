

# List packages to check installation
packages <- c("dplyr","readr","utils","tibble","ggplot2","ggrepel","scales")
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

  complete_blast_results <- data.frame(matrix(NA, nrow=1, ncol=1))
  
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
      
      scaff_int <- scaff
      
      #########################################################################
      #########################################################################
      #########################################################################
      #########################################################################
      
      sample_cns <- structure(list(
        gene = scaff$subject_name, 
        chromosome = scaff$query_name, 
        start = c(scaff$genomic_start), end = c(scaff$genomic_end), 
        cn = c(rep(1L,nrow(scaff))), CNA = c(rep("gain",nrow(scaff))),
        size = c(rep(max(scaff$scaff_stop), nrow(scaff) ))
      ),
      
      .Names = c("gene", "chromosome", "start", "end", "cn", "CNA","size"),
      row.names = c(NA,nrow(scaff)), class = "data.frame" )
      
      
      ##############################################################################
      
      
      sample_cns <- sample_cns %>% add_column(query_code = rep(NA,nrow(sample_cns)), .after = "chromosome")
      sample_cns$query_code <- sample_cns$chromosome
      sample_cns <- sample_cns %>% add_column(query_coordinates = rep(NA,nrow(sample_cns)), .after = "query_code")
      sample_cns$query_coordinates <- sample_cns$chromosome
      
      sample_cns$query_code <- gsub("\\:.*","",sample_cns$query_code)
      sample_cns$query_coordinates <- gsub(".*:","",sample_cns$query_coordinates)
      
      sample_cns <- sample_cns %>% add_column(fragments = rep(NA,nrow(sample_cns)), .after = "size")
      
      n <- c(1)
      edited_sample_cns <- sample_cns[-(1:nrow(sample_cns)),]
      
      
      
      for(i in unique(sample_cns$gene)){
        
        contig <- sample_cns %>%
          filter(gene == i)
        
        if(nrow(contig) > 1){
          
          rows <- which(sample_cns$gene == contig$gene[1])
          
          outliers <- boxplot.stats(rows)$out
          while(length(outliers) != 0){
            for(o in seq(1,length(outliers),1)){
              remove_row <- which(rows == outliers[o])
              contig <- contig[-remove_row,]
              rows <- rows[-remove_row]
            }
            outliers <- boxplot.stats(rows)$out
          }
          
          
          if(nrow(contig) > 1){
            if(any(diff(rows) > 10)){
              
              breaks <- c(0, which(diff(rows) > 10), length(rows))  
              breaks_list <- sapply(seq(length(breaks) - 1),function(i) rows[(breaks[i] + 1):breaks[i+1]])
              
              if(typeof(breaks_list) == "list"){
                
                if(any(lengths(breaks_list) >= 1)){ #### aici era 3
                  
                  for(sequence in seq(1,length(breaks_list),1)){
                    
                    if(length(breaks_list[[sequence]]) >= 1){ #### aici era 3
                      partial_contig <- sample_cns[breaks_list[[sequence]],]
                      edited_sample_cns[n,] <- partial_contig[1,]
                      edited_sample_cns$start[n] <- min(partial_contig$start)
                      edited_sample_cns$end[n] <- max(partial_contig$end)
                      edited_sample_cns$fragments[n] <- nrow(partial_contig)
                      n <- n+1
                    }
                  }
                  contig <- contig[0,]
                } else { contig <- contig[0,] }
              } else { contig <- contig[0,] }
              
            } }
          
          #contig_filtered <- sample_cns[breaks_list[[which(lengths(breaks_list) == max(lengths(breaks_list)))]],]
          
          #Breaks <- c(0, which(diff(vec2) >= 10), length(vec2))
          #sapply(seq(length(Breaks) - 1),function(i) vec2[(Breaks[i] + 1):Breaks[i+1]])
          #a[[which(lengths(a) == max(lengths(a)))]]
        }
        
        if(nrow(contig) != 0){
          edited_sample_cns[n,] <- contig[1,]
          edited_sample_cns$start[n] <- min(contig$start)
          edited_sample_cns$end[n] <- max(contig$end)
          edited_sample_cns$fragments[n] <- nrow(contig)
          n <- n+1 }
      }
      
      edited_sample_cns <- na.omit(edited_sample_cns)
      
      edited_sample_cns <- edited_sample_cns %>% arrange(start)
      
      sample_cns <- edited_sample_cns
      
      sample_cns <- sample_cns[!duplicated(sample_cns), ]
      
      ###############################################################################
      ###############################################################################
      
      
      ####
      # new chrom sizes
      ##
      chrom_sizes <- structure(list(
        chromosome = unique(sample_cns$gene), 
        size = rep(c(max(sample_cns$size)),length(unique(sample_cns$gene)))),
        
        .Names = c("chromosome", "size"),
        class = "data.frame", row.names= c(NA,length(unique(sample_cns$gene))))
      ######
      
      
      
      chrom_order <- c(unique(chrom_sizes$chromosome))
      chrom_key <- setNames(object = as.character(c(seq(1,length(unique(sample_cns$gene)),1))),
                            nm = chrom_order)
      
      
      # chrom_order <- factor(x = chrom_order, levels = rev(chrom_order)) # reverse
      chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
      
      chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                            levels = chrom_order)
      sample_cns[["gene"]] <- factor(x = sample_cns[["gene"]], 
                                     levels = chrom_order)
      
      group.colors <- c(gain = "red", loss = "blue")
      
      
      
      ggplot(data = chrom_sizes) + 
        # base rectangles for the chroms, with numeric value for each chrom on the x-axis
        geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                      xmax = as.numeric(chromosome) + 0.2, 
                      ymax = size, ymin = 0), 
                  # "colour" este pentru border
                  colour="white", fill = "white") + 
        # rotate the plot 90 degrees
        coord_flip() +
        # black & white color theme 
        theme(#axis.text.x = element_text(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        # give the appearance of a discrete axis with chrom labels
        scale_x_discrete(name = "Contigs", limits = names(chrom_key)) +
        # add bands for CNA value
        geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                         xmax = as.numeric(gene) + 0.2, 
                                         ymax = end, ymin = start, fill = CNA)) + 
        scale_fill_manual(values = group.colors) +
        # add 'gain' gene markers
        #geom_text_repel(data = subset(sample_cns, sample_cns$CNA == "gain"), 
        #                aes(x = gene, y = start, label = gene), 
        #                color = "red", show.legend = FALSE) +
        # add 'loss' gene markers
        #geom_text_repel(data = subset(sample_cns, sample_cns$CNA == "loss"), 
        #                aes(x = gene, y = start, label = gene ), 
        #                color = "blue", show.legend = FALSE) +
        ggtitle("Contigs mapping on scaffold") +
        # supress scientific notation on the y-axis
        # scale_y_continuous(labels = comma) +
        scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = 250000),1)) +
        ylab("Scaffold (bp)")  +
        theme(legend.position = "none")
      
      #scale_x_continuous(name=paste("Scaffold",s), limits = c(0, max(chrom_sizes$size)*1.05))
      
      
      ggsave("scaffold_all_contigs.jpeg", plot=last_plot(), limitsize = FALSE)
      
      
      #####################################
      ### Adding blue lines
      
      for(i in seq(1,nrow(sample_cns),1)){
        included_rows <- intersect(which(sample_cns$start < sample_cns$start[i]),
                                   which(sample_cns$end > sample_cns$end[i]))
        if(length(included_rows) != 0){
          sample_cns$CNA[i] <- "loss"  }
        if(sample_cns$fragments[i] < 3){
          sample_cns$CNA[i] <- "loss" }
        
      }
      
      sample_cns_unfiltered <- sample_cns
      
      chrom_sizes <- structure(list(
        chromosome = unique(sample_cns$gene), 
        size = rep(c(max(sample_cns$size)),length(unique(sample_cns$gene)))),
        .Names = c("chromosome", "size"),
        class = "data.frame", row.names= c(NA,length(unique(sample_cns$gene))))
      chrom_order <- c(unique(chrom_sizes$chromosome))
      chrom_key <- setNames(object = as.character(c(seq(1,length(unique(sample_cns$gene)),1))),
                            nm = chrom_order)
      chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
      chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                            levels = chrom_order)
      sample_cns[["gene"]] <- factor(x = sample_cns[["gene"]], 
                                     levels = chrom_order)
      group.colors <- c(gain = "red", loss = "blue")
      ggplot(data = chrom_sizes) + 
        geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                      xmax = as.numeric(chromosome) + 0.2, 
                      ymax = size, ymin = 0), 
                  colour="white", fill = "white") + 
        coord_flip() +
        theme(#axis.text.x = element_text(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        scale_x_discrete(name = "Contigs", limits = names(chrom_key)) +
        geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                         xmax = as.numeric(gene) + 0.2, 
                                         ymax = end, ymin = start, fill = CNA)) + 
        scale_fill_manual(values = group.colors) +
        ggtitle("Contigs mapping on scaffold") +
        scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = 250000),1)) +
        ylab("Scaffold (bp)")  +
        theme(legend.position = "none")
      
      
      ggsave("scaffold_all_contigs_red_blue.jpeg", plot=last_plot(), limitsize = FALSE)
      
      ###############################################################################
      ###############################################################################
      # Filtering out blue contigs
      
      sample_cns <- sample_cns %>% filter(CNA == "gain")
      
      
      chrom_sizes <- structure(list(
        chromosome = unique(sample_cns$gene), 
        size = rep(c(max(sample_cns$size)),length(unique(sample_cns$gene)))),
        .Names = c("chromosome", "size"),
        class = "data.frame", row.names= c(NA,length(unique(sample_cns$gene))))
      chrom_order <- c(unique(chrom_sizes$chromosome))
      chrom_key <- setNames(object = as.character(c(seq(1,length(unique(sample_cns$gene)),1))),
                            nm = chrom_order)
      chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
      chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                            levels = chrom_order)
      sample_cns[["gene"]] <- factor(x = sample_cns[["gene"]], 
                                     levels = chrom_order)
      group.colors <- c(gain = "red", loss = "blue")
      ggplot(data = chrom_sizes) + 
        geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                      xmax = as.numeric(chromosome) + 0.2, 
                      ymax = size, ymin = 0), 
                  colour="white", fill = "white") + 
        coord_flip() +
        theme(#axis.text.x = element_text(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        scale_x_discrete(name = "Contigs", limits = names(chrom_key)) +
        geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                         xmax = as.numeric(gene) + 0.2, 
                                         ymax = end, ymin = start, fill = CNA)) + 
        scale_fill_manual(values = group.colors) +
        ggtitle("Contigs mapping on scaffold") +
        scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = 250000),1)) +
        ylab("Scaffold (bp)")  +
        theme(legend.position = "none")
      
      
      ggsave("minimal_scaffold.jpeg", plot=last_plot(), limitsize = FALSE)
      
      
      
      ##################################################
      
      #setwd(file.path(mainDir,dir))
      # output list of red contigs for creating fasta
      selected_contigs <- data.frame(matrix(NA,nrow=length(unique(sample_cns$gene)),ncol=2))
      selected_contigs[,1] <- unique(sample_cns$gene)
      selected_contigs[,2] <- seq(1,length(unique(sample_cns$gene)),1)
      write.table(selected_contigs, file=paste(dir,"scaffold",s,"selected_contigs.csv",sep="_"),
                  sep=",",row.names = F, col.names = FALSE)
      
      
      # output table of all contigs of interest
      write.table(sample_cns_unfiltered, 
                  file=paste(dir,"scaffold",s,"all_contigs.csv",sep="_"),
                  sep=",",row.names = F, col.names = T)
      write.table(sample_cns, 
                  file=paste(dir,"scaffold",s,"minimal_contigs.csv",sep="_"),
                  sep=",",row.names = F, col.names = T)
      
      
      # table included
      sample_cns_included_contigs <- sample_cns_unfiltered %>% filter(CNA == "loss")
      included_contigs_vector <- sample_cns_included_contigs$gene
    
      maximum_included <- c()
      for(i in seq(1,nrow(sample_cns),1)){
        included_rows_vector <- intersect(which(sample_cns_included_contigs$start >= sample_cns$start[i]), 
                                          which(sample_cns_included_contigs$end <= sample_cns$end[i]))
        maximum_included <- append(maximum_included,length(included_rows_vector)) }
      
      included_contigs <- data.frame(matrix(NA,nrow=max(maximum_included), ncol=nrow(sample_cns)))
      colnames(included_contigs) <- as.character(sample_cns$gene)
      for(i in seq(1,nrow(sample_cns),1)){
        included_rows <- intersect(which(sample_cns_included_contigs$start >= sample_cns$start[i]), 
                                   which(sample_cns_included_contigs$end <= sample_cns$end[i]))
        if(length(included_rows) != 0){
          included_contigs[1:length(included_rows),i] <- as.character(sample_cns_included_contigs$gene[included_rows])
        } }
      write.table(included_contigs, 
                  file=paste(dir,"scaffold",s,"included_contigs.csv",sep="_"), 
                  sep=",", row.names = F, na="")
      
      #
      #setwd(file.path(mainDir,dir,"scaffolds",s))
      
      #########################################################################
      #########################################################################
      #########################################################################
      #########################################################################
      #########################################################################
      
      
      ###
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
      ###
      
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
    
    ###########################################################################
    ###########################################################################
    
    scaff <- contigs_of_interest
    
    ##############################################################################
    
    sample_cns <- structure(list(
      gene = scaff$subject_name, 
      chromosome = scaff$query_name, 
      start = c(scaff$genomic_start), end = c(scaff$genomic_end), 
      cn = c(rep(1L,nrow(scaff))), CNA = c(rep("gain",nrow(scaff))),
      size = c(rep(max(scaff$scaff_stop), nrow(scaff) ))
    ),
    
    .Names = c("gene", "chromosome", "start", "end", "cn", "CNA","size"),
    row.names = c(NA,nrow(scaff)), class = "data.frame" )
    
    
    ##############################################################################
    
    
    sample_cns <- sample_cns %>% add_column(query_code = rep(NA,nrow(sample_cns)), .after = "chromosome")
    sample_cns$query_code <- sample_cns$chromosome
    sample_cns <- sample_cns %>% add_column(query_coordinates = rep(NA,nrow(sample_cns)), .after = "query_code")
    sample_cns$query_coordinates <- sample_cns$chromosome
    
    sample_cns$query_code <- gsub("\\:.*","",sample_cns$query_code)
    sample_cns$query_coordinates <- gsub(".*:","",sample_cns$query_coordinates)
    
    sample_cns <- sample_cns %>% add_column(fragments = rep(NA,nrow(sample_cns)), .after = "size")
    
    n <- c(1)
    edited_sample_cns <- sample_cns[-(1:nrow(sample_cns)),]
    
    
    
    for(i in unique(sample_cns$gene)){
      
      contig <- sample_cns %>%
        filter(gene == i)
      
      if(nrow(contig) > 1){
        
        rows <- which(sample_cns$gene == contig$gene[1])
        
        outliers <- boxplot.stats(rows)$out
        while(length(outliers) != 0){
          for(o in seq(1,length(outliers),1)){
            remove_row <- which(rows == outliers[o])
            contig <- contig[-remove_row,]
            rows <- rows[-remove_row]
          }
          outliers <- boxplot.stats(rows)$out
        }
        
        
        if(nrow(contig) > 1){
          if(any(diff(rows) > 10)){
            
            breaks <- c(0, which(diff(rows) > 10), length(rows))  
            breaks_list <- sapply(seq(length(breaks) - 1),function(i) rows[(breaks[i] + 1):breaks[i+1]])
            
            if(typeof(breaks_list) == "list"){
              
              if(any(lengths(breaks_list) >= 1)){ #### aici era 3
                
                for(sequence in seq(1,length(breaks_list),1)){
                  
                  if(length(breaks_list[[sequence]]) >= 1){ #### aici era 3
                    partial_contig <- sample_cns[breaks_list[[sequence]],]
                    edited_sample_cns[n,] <- partial_contig[1,]
                    edited_sample_cns$start[n] <- min(partial_contig$start)
                    edited_sample_cns$end[n] <- max(partial_contig$end)
                    edited_sample_cns$fragments[n] <- nrow(partial_contig)
                    n <- n+1
                  }
                }
                contig <- contig[0,]
              } else { contig <- contig[0,] }
            } else { contig <- contig[0,] }
            
          } }
        
        #contig_filtered <- sample_cns[breaks_list[[which(lengths(breaks_list) == max(lengths(breaks_list)))]],]
        
        #Breaks <- c(0, which(diff(vec2) >= 10), length(vec2))
        #sapply(seq(length(Breaks) - 1),function(i) vec2[(Breaks[i] + 1):Breaks[i+1]])
        #a[[which(lengths(a) == max(lengths(a)))]]
      }
      
      if(nrow(contig) != 0){
        edited_sample_cns[n,] <- contig[1,]
        edited_sample_cns$start[n] <- min(contig$start)
        edited_sample_cns$end[n] <- max(contig$end)
        edited_sample_cns$fragments[n] <- nrow(contig)
        n <- n+1 }
    }
    
    edited_sample_cns <- na.omit(edited_sample_cns)
    
    edited_sample_cns <- edited_sample_cns %>% arrange(start)
    
    sample_cns <- edited_sample_cns
    
    sample_cns <- sample_cns[!duplicated(sample_cns), ]
    
    ###############################################################################
    ###############################################################################
    
    
    ####
    # new chrom sizes
    ##
    chrom_sizes <- structure(list(
      chromosome = unique(sample_cns$gene), 
      size = rep(c(max(sample_cns$size)),length(unique(sample_cns$gene)))),
      
      .Names = c("chromosome", "size"),
      class = "data.frame", row.names= c(NA,length(unique(sample_cns$gene))))
    ######
    
    
    
    chrom_order <- c(unique(chrom_sizes$chromosome))
    chrom_key <- setNames(object = as.character(c(seq(1,length(unique(sample_cns$gene)),1))),
                          nm = chrom_order)
    
    
    # chrom_order <- factor(x = chrom_order, levels = rev(chrom_order)) # reverse
    chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
    
    chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                          levels = chrom_order)
    sample_cns[["gene"]] <- factor(x = sample_cns[["gene"]], 
                                   levels = chrom_order)
    
    group.colors <- c(gain = "red", loss = "blue")
    
    
    
    ggplot(data = chrom_sizes) + 
      # base rectangles for the chroms, with numeric value for each chrom on the x-axis
      geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                    xmax = as.numeric(chromosome) + 0.2, 
                    ymax = size, ymin = 0), 
                # "colour" este pentru border
                colour="white", fill = "white") + 
      # rotate the plot 90 degrees
      coord_flip() +
      # black & white color theme 
      theme(#axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      # give the appearance of a discrete axis with chrom labels
      scale_x_discrete(name = "Contigs", limits = names(chrom_key)) +
      # add bands for CNA value
      geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                       xmax = as.numeric(gene) + 0.2, 
                                       ymax = end, ymin = start, fill = CNA)) + 
      scale_fill_manual(values = group.colors) +
      # add 'gain' gene markers
      #geom_text_repel(data = subset(sample_cns, sample_cns$CNA == "gain"), 
      #                aes(x = gene, y = start, label = gene), 
      #                color = "red", show.legend = FALSE) +
      # add 'loss' gene markers
      #geom_text_repel(data = subset(sample_cns, sample_cns$CNA == "loss"), 
      #                aes(x = gene, y = start, label = gene ), 
      #                color = "blue", show.legend = FALSE) +
      ggtitle("Contigs mapping on scaffold") +
      # supress scientific notation on the y-axis
      # scale_y_continuous(labels = comma) +
      scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = 250000),1)) +
      ylab("Scaffold (bp)")  +
      theme(legend.position = "none")
    
    #scale_x_continuous(name=paste("Scaffold",s), limits = c(0, max(chrom_sizes$size)*1.05))
    
    
    ggsave("scaffold_all_contigs.jpeg", plot=last_plot(), limitsize = FALSE)
    
    
    #####################################
    ### Adding blue lines
    
    for(i in seq(1,nrow(sample_cns),1)){
      included_rows <- intersect(which(sample_cns$start < sample_cns$start[i]),
                                 which(sample_cns$end > sample_cns$end[i]))
      if(length(included_rows) != 0){
        sample_cns$CNA[i] <- "loss"  }
      if(sample_cns$fragments[i] < 3){
        sample_cns$CNA[i] <- "loss" }
      
    }
    
    sample_cns_unfiltered <- sample_cns
    
    chrom_sizes <- structure(list(
      chromosome = unique(sample_cns$gene), 
      size = rep(c(max(sample_cns$size)),length(unique(sample_cns$gene)))),
      .Names = c("chromosome", "size"),
      class = "data.frame", row.names= c(NA,length(unique(sample_cns$gene))))
    chrom_order <- c(unique(chrom_sizes$chromosome))
    chrom_key <- setNames(object = as.character(c(seq(1,length(unique(sample_cns$gene)),1))),
                          nm = chrom_order)
    chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
    chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                          levels = chrom_order)
    sample_cns[["gene"]] <- factor(x = sample_cns[["gene"]], 
                                   levels = chrom_order)
    group.colors <- c(gain = "red", loss = "blue")
    ggplot(data = chrom_sizes) + 
      geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                    xmax = as.numeric(chromosome) + 0.2, 
                    ymax = size, ymin = 0), 
                colour="white", fill = "white") + 
      coord_flip() +
      theme(#axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_x_discrete(name = "Contigs", limits = names(chrom_key)) +
      geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                       xmax = as.numeric(gene) + 0.2, 
                                       ymax = end, ymin = start, fill = CNA)) + 
      scale_fill_manual(values = group.colors) +
      ggtitle("Contigs mapping on scaffold") +
      scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = 250000),1)) +
      ylab("Scaffold (bp)")  +
      theme(legend.position = "none")
    
    
    ggsave("scaffold_all_contigs_red_blue.jpeg", plot=last_plot(), limitsize = FALSE)
    
    ###############################################################################
    ###############################################################################
    # Filtering out blue contigs
    
    sample_cns <- sample_cns %>% filter(CNA == "gain")
    
    
    chrom_sizes <- structure(list(
      chromosome = unique(sample_cns$gene), 
      size = rep(c(max(sample_cns$size)),length(unique(sample_cns$gene)))),
      .Names = c("chromosome", "size"),
      class = "data.frame", row.names= c(NA,length(unique(sample_cns$gene))))
    chrom_order <- c(unique(chrom_sizes$chromosome))
    chrom_key <- setNames(object = as.character(c(seq(1,length(unique(sample_cns$gene)),1))),
                          nm = chrom_order)
    chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
    chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                          levels = chrom_order)
    sample_cns[["gene"]] <- factor(x = sample_cns[["gene"]], 
                                   levels = chrom_order)
    group.colors <- c(gain = "red", loss = "blue")
    ggplot(data = chrom_sizes) + 
      geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                    xmax = as.numeric(chromosome) + 0.2, 
                    ymax = size, ymin = 0), 
                colour="white", fill = "white") + 
      coord_flip() +
      theme(#axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_x_discrete(name = "Contigs", limits = names(chrom_key)) +
      geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                       xmax = as.numeric(gene) + 0.2, 
                                       ymax = end, ymin = start, fill = CNA)) + 
      scale_fill_manual(values = group.colors) +
      ggtitle("Contigs mapping on scaffold") +
      scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = 250000),1)) +
      ylab("Scaffold (bp)")  +
      theme(legend.position = "none")
    
    
    ggsave("minimal_scaffold.jpeg", plot=last_plot(), limitsize = FALSE)
    
    ###########################################################################
    #setwd(file.path(mainDir,dir))
    #
    # output list of red contigs for creating fasta
    selected_contigs <- data.frame(matrix(NA,nrow=length(unique(sample_cns$gene)),ncol=2))
    selected_contigs[,1] <- unique(sample_cns$gene)
    selected_contigs[,2] <- seq(1,length(unique(sample_cns$gene)),1)
    write.table(selected_contigs, file=paste(dir,"chromosome","selected_contigs.csv",sep="_"),
                sep=",",row.names = F, col.names = FALSE)
    
    
    # output table of all contigs of interest
    write.table(sample_cns_unfiltered, 
                file=paste(dir,"chromosome","all_contigs.csv",sep="_"),
                sep=",",row.names = F, col.names = T)
    write.table(sample_cns, 
                file=paste(dir,"chromosome","minimal_contigs.csv",sep="_"),
                sep=",",row.names = F, col.names = T)
    
    
    # table included
    sample_cns_included_contigs <- sample_cns_unfiltered %>% filter(CNA == "loss")
    included_contigs_vector <- sample_cns_included_contigs$gene
    
    maximum_included <- c()
    for(i in seq(1,nrow(sample_cns),1)){
      included_rows_vector <- intersect(which(sample_cns_included_contigs$start >= sample_cns$start[i]), 
                                        which(sample_cns_included_contigs$end <= sample_cns$end[i]))
      maximum_included <- append(maximum_included,length(included_rows_vector)) }
    
    included_contigs <- data.frame(matrix(NA,nrow=max(maximum_included), ncol=nrow(sample_cns)))
    colnames(included_contigs) <- as.character(sample_cns$gene)
    for(i in seq(1,nrow(sample_cns),1)){
      included_rows <- intersect(which(sample_cns_included_contigs$start >= sample_cns$start[i]), 
                                 which(sample_cns_included_contigs$end <= sample_cns$end[i]))
      if(length(included_rows) != 0){
        included_contigs[1:length(included_rows),i] <- as.character(sample_cns_included_contigs$gene[included_rows])
      } }
    write.table(included_contigs, 
                file=paste(dir,"chromosome","included_contigs.csv",sep="_"), 
                sep=",", row.names = F, na="")
    
    
    #
    ###########################################################################
    ###########################################################################
    
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

