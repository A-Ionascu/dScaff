# List packages to check installation
packages <- c("ggplot2", "dplyr","tidyr","shiny","shinythemes","DT",
              "gridExtra","ggthemes","nortest",
              "shinycssloaders","shinybusy","shinycustomloader",
              "fontawesome","RColorBrewer","bslib",
              "htmltools",
              
              "readr","utils","tibble","ggrepel","scales"
              )
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


if(interactive()){
  
  
  # just make ui nice
  
  ui <- shinyUI( navbarPage( 
    #theme = shinytheme("cerulean"),
    theme = bs_theme(version = 5, bootswatch = "cosmo",
                     base_font = "Times New Roman", heading_font = "Times New Roman"),
    title="dScaff",
    
    main_page <- tabPanel(title=HTML(paste0("App")), value="ica",
                          titlePanel(HTML(paste0("dScaff"))),
                          
                          sidebarLayout(
                            sidebarPanel(#style = "position:fixed;width:20%;",
                              title="File input",
                              fileInput("file", "Choose .csv file",
                                        accept = c('.csv','.xls', ".xlsx"), buttonLabel=list(icon("upload"), paste("Browse"))),
                              width = 3,
                              br(),
                              selectInput("scaffolds","Choose scaffold  (if using scaffolds based reference)",choices = c()),
                              br(),
                              numericInput("missing_hits",label=HTML(paste0("Minimum missing hits within contig fragments:")),value = c(10)),
                              actionButton("update","Update results", icon=icon("redo")),
                              br(),br(),br(),
                              numericInput("contig_fragments","Minimum hits for contig fragments:",value = c(3)),
                              actionButton("update2","Update results",  icon=icon("redo")),
                              br(),br(),br(),
                              numericInput("red_blue_threshold","Minimum hits for red contig fragments:",value = c(3)),
                              actionButton("update3","Update results",  icon=icon("redo")),
                              br(),br(),br(),
                              selectInput(inputId = "plot_orientation", label="Choose orientation of contigs progression on Y-axis:", 
                                          choices = list("From down to up" = "down_to_up", "From up to down" = "up_to_down"),
                                          selected = "down_to_up"),
                              br()
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel(
                                  title="Plots", 
                                  value="plots",
                                  br(),
                                  h4("Minimal scaffold"),
                                  br(),
                                  downloadButton("download_minimal_contigs_plot", "Export minimal scaffold plot"),
                                  br(),br(),
                                  textOutput("error_message_plot"),
                                  withLoader(plotOutput("minimal_contigs_plot", height = 700), type="html", loader="dnaspin"),
                                  
                                  
                                  
                                  br(),br(),br(),
                                  h4("All contigs (red/blue)"),
                                  br(),
                                  downloadButton("download_red_blue_contigs_plot", "Export processed contigs plot (red/blue)"),
                                  withLoader(plotOutput("red_blue_contigs_plot", height = 700), type="html", loader="dnaspin"),
                                  
                                  
                                  br(),br(),br(),
                                  h4("Processed contigs"),
                                  br(),
                                  downloadButton("download_all_contigs_plot", "Export processed contigs plot"),
                                  withLoader(plotOutput("all_contigs_plot", height = 700), type="html", loader="dnaspin"),
                                  
                                  
                                  br(),br(),br(),
                                  h4("Unprocessed contigs"),
                                  br(),
                                  downloadButton("download_unprocessed_contigs_plot", "Export unprocessed contigs plot"),
                                  withLoader(plotOutput("unprocessed_contigs_plot", height = 700), type="html", loader="dnaspin"),
                                  
                                  
                                  br(),br(),br(),br()),
                                
                                tabPanel(title="Minimal contigs", value="minimal_contigs",
                                         br(),
                                         h4("Minimal contigs table"),
                                         br(),
                                         downloadButton("download_minimal_contigs_table", "Export minimal contigs table"),
                                         withLoader(DT::dataTableOutput('minimal_contigs_table'), type="html", loader="dnaspin"),
                                         textOutput("error_message_table"),
                                        
                                         br(),br(),br(),br()),
                                tabPanel(title="All contigs", value="all_contigs",
                                         br(),
                                         h4("All contigs table"),
                                         br(),
                                         downloadButton("download_all_contigs_table", "Export all contigs table"),
                                         withLoader(DT::dataTableOutput('all_contigs_table'), type="html", loader="dnaspin"),
                                         br(),br(),br(),br()),
                                tabPanel(title="Included contigs table", value="included_contigs",
                                         br(),
                                         h4("Included contigs"),
                                         br(),
                                         downloadButton("download_included_contigs_table", "Export included contigs table"),
                                         withLoader(DT::dataTableOutput('included_contigs_table'), type="html", loader="dnaspin"),
                                         br(),br(),br(),br()),
                                tabPanel(title="Unprocessed contigs", value="unprocessed_contigs",
                                         br(),
                                         h4("Unprocessed contigs table"),
                                         br(),
                                         downloadButton("download_unprocessed_contigs_table", "Export unprocessed contigs table"),
                                         withLoader(DT::dataTableOutput('unprocessed_contigs_table'), type="html", loader="dnaspin"),
                                         br(),br(),br(),br()),
                                tabPanel(title="Original table", value="original",
                                         br(),
                                         h4("Original table"),
                                         br(),
                                         withLoader(DT::dataTableOutput('original_table'), type="html", loader="dnaspin"),
                                         br(),br(),br(),br()),
                                tabPanel(title="Orientation tables", value="orientation",
                                         br(),
                                         h4("Orientation tables"),
                                         br(),
                                         downloadButton("download_orientation_table", "Export all contigs orientation table"),
                                         br(),br(),
                                         downloadButton("download_merged_orientation_table", "Export contigs of interest orientation table"),
                                         br(),br(),
                                         downloadButton("download_undecided_orientation_table", "Export contigs with undecided orientation table"),
                                         #withLoader(DT::dataTableOutput('orientation_table'), type="html", loader="dnaspin"),
                                         br(),br(),br(),br()),
                                
                              ))
                          
                          )                     
    ),
    about_page <- tabPanel(title="About",
                           titlePanel="About",
                           h3(tags$a(href="https://github.com/DL-UB/dScaff","dScaff"), " is a scaffolding strategy for draft assemblies (contigs) based on reference assembly using either 
                              genes queries or ranked queries."),
                           br(),
                           h4("Digital Scaffolding (dScaff) aims improve new draft assemblies of organisms that have a better 
                              reference assembly in the NCBI databases. 
                              This strategy is most valuable when performing de novo assemblies of local species. 
                              The inplemented strategy aims to use the draft assembly, sequences of all genes from the reference genome 
                              and a table of their annotation in order to enhance scaffolding posibilities. Alternatively, the ranked queries 
                              strategy is suitable for reference assemblies without annotated genes and can be implemented using the ",
                              tags$a(href="https://github.com/DL-UB/SubSequencesExtractor","SubSequenceExtractor")," complementary script."),
                           br(),
                           h4("dScaff ICA allows for image creation and easy modifications to images without running the main dScaff script.
                              The app uses dScaff output (contigs_of_interest.csv from tmp folder) and
                              creates the images similar to the dScaff script. The user can adjust parameters and modify the images on demand.
                              Additionally, complementary tables are recreated outside the dScaff procedure."),
                           br(),
                           h4("All images and tables can be exported from the app.")
                           )
    
    ))
    
  server <- function(input, output, session){ 
    
    myData <- reactive({
      inFile <- input$file
      if(is.null(inFile)) return(NULL)
      data_user <- read.csv(inFile$datapath, header = TRUE)
      
      data_user <- data_user[,-19]
      data_user <- data_user[,-18]
      data_user <- data_user[,-1]
      colnames(data_user) <- c("Query","Contig","Identity","Alignment length","Mismatches","Gaps",
                               "Query start","Query end","Subject start","Subject end",
                               "Genomic start","Genomic end","E value","Bit score",
                               "Query length","Reference","Orientation","Contig length")
      
      data_user
    })
    
    output$original_table <- DT::renderDataTable( 
      DT::datatable(myData(),
                    filter="top", 
                    options = list(pageLength = -1,
                              lengthMenu = list(c(10,50,100,-1),c('10','50','100','All')),
                              columnDefs = list(list(className = 'dt-center', targets = 0:7)),
                              dom = "ft")))
    
    observeEvent(input$file, {
      req(myData())
      updateSelectInput(
        session,
        "scaffolds",
        choices = unique(myData()$Reference),
        selected = unique(myData()$Reference)[1])
    })
    
    # Ask for minimum missing hits within contig fragments
    missing_hits_reactive <- reactive({
      if(input$update == 0){
        missing_hits <- c(10)
        return(missing_hits) }
      
      else{
        a <- eventReactive(input$update, {
          missing_hits_updated <- input$missing_hits
          missing_hits_updated }) }
      
      if(is.null(a())){
        missing_hits <- c(10) }
      else{missing_hits <- a()}
      
      missing_hits
    })
    
    
    # Ask for minimum hits for contig fragments
    contig_fragments_reactive <- reactive({
      if(input$update2 == 0){
        contig_fragments <- c(3)
        return(contig_fragments) }
      
      else{
        a <- eventReactive(input$update2, {
          contig_fragments_updated <- input$contig_fragments
          contig_fragments_updated }) }
      
      if(is.null(a())){
        contig_fragments <- c(3) }
      else{contig_fragments <- a()}
      
      contig_fragments
    })
    
    
    # Ask for minimum hits for red contig fragments
    red_contig_fragments_reactive <- reactive({
      if(input$update3 == 0){
        red_contig_fragments <- c(3)
        return(red_contig_fragments) }
      
      else{
        a <- eventReactive(input$update3, {
          red_contig_fragments_updated <- input$red_blue_threshold
          red_contig_fragments_updated }) }
      
      if(is.null(a())){
        red_contig_fragments <- c(3) }
      else{red_contig_fragments <- a()}
      
      red_contig_fragments
    })
    
    
############################################################################### 
###############################################################################
###############################################################################
###############################################################################
###############################################################################
    
    
    unprocessed_hits_reactive <- reactive({
      req(input$file)
      req(input$scaffolds)
      inFile <- input$file
      data_raw <- read.csv(inFile$datapath, header = TRUE)
      scaff <- as.data.frame(data_raw)
      scaff <- scaff %>% filter(ref_scaff == input$scaffolds)
      #colnames(data) <- c("group", "case", "type", "gene", "RB", "RT", "Ct")
      
      sample_cns <- structure(list(
        gene = scaff$subject_name, 
        chromosome = scaff$query_name, 
        start = c(scaff$genomic_start), end = c(scaff$genomic_end), 
        cn = c(rep(1L,nrow(scaff))), CNA = c(rep("gain",nrow(scaff))),
        size = c(rep(max(scaff$scaff_stop), nrow(scaff) ))
      ),
      
      .Names = c("gene", "chromosome", "start", "end", "cn", "CNA","size"),
      row.names = c(NA,nrow(scaff)), class = "data.frame" )
      
      sample_cns <- sample_cns %>% add_column(query_code = rep(NA,nrow(sample_cns)), .after = "chromosome")
      sample_cns$query_code <- sample_cns$chromosome
      sample_cns <- sample_cns %>% add_column(query_coordinates = rep(NA,nrow(sample_cns)), .after = "query_code")
      sample_cns$query_coordinates <- sample_cns$chromosome
      sample_cns$query_code <- gsub("\\:.*","",sample_cns$query_code)
      sample_cns$query_coordinates <- gsub(".*:","",sample_cns$query_coordinates)
      sample_cns <- sample_cns %>% add_column(fragments = rep(NA,nrow(sample_cns)), .after = "size")
      sample_cns <- sample_cns %>% add_column(ctg_size = rep(NA,nrow(sample_cns)), .after = "end")
      for(csize in seq(1,nrow(sample_cns),1)){
        sample_cns$ctg_size[csize] <- (sample_cns$end[csize] - sample_cns$start[csize]) + 1
      }
      
      sample_cns
    })
    
###############################################################################
###############################################################################    
    
    unprocessed_hits_ggplot <- reactive({
      req(input$plot_orientation)
      
      sample_cns <- unprocessed_hits_reactive()
      sample_cns <- data.frame(sample_cns)
      
      if(input$plot_orientation == "up_to_down"){
        sample_cns <- sample_cns %>% arrange(desc(start)) 
      }
      
      chrom_sizes <- structure(list(
        chromosome = unique(sample_cns$gene), 
        size = rep(c(max(sample_cns$size)),length(unique(sample_cns$gene)))),
        .Names = c("chromosome", "size"),
        class = "data.frame", row.names= c(NA,length(unique(sample_cns$gene))))
      chrom_order <- c(unique(chrom_sizes$chromosome))
      chrom_key <- setNames(object = as.character(c(seq(1,length(unique(sample_cns$gene)),1))),nm = chrom_order)
      
      # chrom_order <- factor(x = chrom_order, levels = rev(chrom_order)) # reverse
      chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
      chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], levels = chrom_order)
      sample_cns[["gene"]] <- factor(x = sample_cns[["gene"]], levels = chrom_order)
      group.colors <- c(gain = "red", loss = "blue")
      
      
      unprocessed_contigs_ggplot_sub <- ggplot(data = chrom_sizes) + 
        geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                      xmax = as.numeric(chromosome) + 0.2, 
                      ymax = size, ymin = 0), 
                  colour="white", fill = "white") + 
        coord_flip() + 
        theme( 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_x_discrete(name = "Contigs", #limits = names(chrom_key)
                         ) +
        geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                         xmax = as.numeric(gene) + 0.2, 
                                         ymax = end, ymin = start, fill = CNA)) + 
        scale_fill_manual(values = group.colors) +
        #ggtitle("Contigs mapping on scaffold") +
        scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = max(sample_cns$size)/50),1)) +
        ylab("Scaffold (bp)")  +
        theme(legend.position = "none")
      
      
      unprocessed_contigs_ggplot_sub
    })

    
    output$unprocessed_contigs_plot <- renderPlot({ unprocessed_hits_ggplot() })
    
###############################################################################
############################################################################### 
###############################################################################
###############################################################################
###############################################################################
############################################################################### 
    
    unprocessed_hits_reactive_edited <- reactive({
      sample_cns <- unprocessed_hits_reactive()
      sample_cns <- data.frame(sample_cns)
      sample_cns <- sample_cns[,-11]
      sample_cns <- sample_cns[,-9]
      sample_cns <- sample_cns[,-8]
      colnames(sample_cns) <- c("Contig","Query","Query code","Query coordinates","Query start","Query end","Query size","Subject size")
      
      sample_cns
    })
    
    output$unprocessed_contigs_table <- renderDataTable(
      unprocessed_hits_reactive_edited(), filter="top", escape = FALSE, 
      options = list(lengthChange = FALSE, searching=TRUE, paging=FALSE,
                     columnDefs = list(list(className = 'dt-center', targets = 2:3)),dom="t"))
    
    
    output$download_unprocessed_contigs_plot <- downloadHandler(
      filename = function(){paste("unprocessed_contigs_plot_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
      content = function(file){ ggsave(file, plot=unprocessed_hits_ggplot(), dpi = 600, width = 16, height=10) } )
    
    output$download_unprocessed_contigs_table <- downloadHandler(
      filename = function(){paste("uprocessed_contigs_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(unprocessed_hits_reactive_edited(), fname) })
    
    
###############################################################################
############################################################################### 
###############################################################################
###############################################################################
###############################################################################
###############################################################################
    
    table_orientation_reactive <- reactive({
      req(input$file)
      req(input$missing_hits)
      req(input$contig_fragments)
      req(input$scaffolds)
      inFile <- input$file
      data_raw <- read.csv(inFile$datapath, header = TRUE)
      scaff <- as.data.frame(data_raw)
      scaff <- scaff %>% filter(ref_scaff == input$scaffolds)
      
      #########
      sample_cns <- structure(list(
        gene = scaff$subject_name, 
        chromosome = scaff$query_name, 
        start = c(scaff$genomic_start), end = c(scaff$genomic_end), 
        cn = c(rep(1L,nrow(scaff))), CNA = c(rep("gain",nrow(scaff))),
        size = c(rep(max(scaff$scaff_stop), nrow(scaff))),
        orientation = scaff$orientation,
        contig_start = scaff$subject_start,
        contig_end = scaff$subject_end,
        contig_length = scaff$contig_length
      ),
      
      .Names = c("gene", "chromosome", "start", "end", "cn", "CNA","size","orientation","contig_start","contig_end","contig_length"),
      row.names = c(NA,nrow(scaff)), class = "data.frame" )
      
      sample_cns <- sample_cns %>% add_column(query_code = rep(NA,nrow(sample_cns)), .after = "chromosome")
      sample_cns$query_code <- sample_cns$chromosome
      sample_cns <- sample_cns %>% add_column(query_coordinates = rep(NA,nrow(sample_cns)), .after = "query_code")
      sample_cns$query_coordinates <- sample_cns$chromosome
      sample_cns$query_code <- gsub("\\:.*","",sample_cns$query_code)
      sample_cns$query_coordinates <- gsub(".*:","",sample_cns$query_coordinates)
      sample_cns <- sample_cns %>% add_column(fragments = rep(NA,nrow(sample_cns)), .after = "size")
      
      #########
      table_orientation <- data.frame(matrix(NA,nrow=0,ncol=4))
      colnames(table_orientation) <- c("Contig","Plus_fragments","Minus_fragments","Contig_lengths")
      or <- c(1)
      for(i in unique(sample_cns$gene)){
        contig <- sample_cns %>% filter(gene == i)
        plus_count <- sum(contig$orientation == "plus")
        minus_count <- sum(contig$orientation == "minus")
        table_orientation[or,1] <- i
        table_orientation[or,2] <- plus_count
        table_orientation[or,3] <- minus_count
        table_orientation[or,4] <- sample_cns$contig_length[which(sample_cns$gene == table_orientation[or,1])[1]]
        or <- or+1
      }
      
      table_orientation 
      
    })
    
    table_orientation_reactive_edited <- reactive({
      table_orientation_reactive_edited <- table_orientation_reactive()
      table_orientation_reactive_edited <- data.frame(table_orientation_reactive_edited)
      colnames(table_orientation_reactive_edited) <- c("Contig","Fragments on plus","Fragments on minus","Contig length")
      table_orientation_reactive_edited
      })
    
    output$download_orientation_table <- downloadHandler(
      filename = function(){paste("orientation_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(table_orientation_reactive_edited(), fname) }) 
    
    
#    output$orientation_table <- renderDataTable(
#      table_orientation_reactive_edited(), filter="top", escape = FALSE, 
#      options = list(lengthChange = FALSE, searching=TRUE, paging=FALSE,
#                     columnDefs = list(list(className = 'dt-center', targets = 1:4)),dom="t"))
    
    
################################################################################    
############################################################################### 
###############################################################################
###############################################################################
###############################################################################
###############################################################################  
      
    all_contigs_reactive <- reactive({
      req(input$file)
      req(input$missing_hits)
      req(input$contig_fragments)
      req(input$scaffolds)
      inFile <- input$file
      data_raw <- read.csv(inFile$datapath, header = TRUE)
      scaff <- as.data.frame(data_raw)
      scaff <- scaff %>% filter(ref_scaff == input$scaffolds)
    
      sample_cns <- structure(list(
        gene = scaff$subject_name, 
        chromosome = scaff$query_name, 
        start = c(scaff$genomic_start), end = c(scaff$genomic_end), 
        cn = c(rep(1L,nrow(scaff))), CNA = c(rep("gain",nrow(scaff))),
        size = c(rep(max(scaff$scaff_stop), nrow(scaff))),
        orientation = scaff$orientation,
        contig_start = scaff$subject_start,
        contig_end = scaff$subject_end,
        contig_length = scaff$contig_length
      ),
      
      .Names = c("gene", "chromosome", "start", "end", "cn", "CNA","size","orientation","contig_start","contig_end","contig_length"),
      row.names = c(NA,nrow(scaff)), class = "data.frame" )
      
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
          
          #rows <- which(sample_cns$gene == contig$gene[1])
          
          rows <- c()
          for(r in seq(1,nrow(contig),1)){
            rows <- append(rows,which(unique(sample_cns$chromosome) == contig$chromosome[r]))
          }
          
          outliers <- boxplot.stats(rows)$out
          while(length(outliers) != 0){
            #for(o in seq(1,length(outliers),1)){
            for(o in unique(outliers)){
              remove_row <- which(rows == o) #which(rows == outliers[o])
              contig <- contig[-remove_row,]
              rows <- rows[-remove_row]
            }
            outliers <- boxplot.stats(rows)$out
          }
          
          
          if(nrow(contig) > 1){
            if(any(diff(rows) >= missing_hits_reactive())){
              
              breaks <- c(0, which(diff(rows) >= missing_hits_reactive()), length(rows))  
              breaks_list <- sapply(seq(length(breaks) - 1),function(i) rows[(breaks[i] + 1):breaks[i+1]])
              
              if(typeof(breaks_list) == "list"){
                
                if(any(lengths(breaks_list) >= contig_fragments_reactive())){ 
                  
                  for(sequence in seq(1,length(breaks_list),1)){
                    
                    if(length(breaks_list[[sequence]]) >= contig_fragments_reactive()){
                      #partial_contig <- sample_cns[breaks_list[[sequence]],]
                      partial_contig <- contig[which(contig$chromosome %in% unique(sample_cns$chromosome)[breaks_list[[sequence]]]),]
                      edited_sample_cns[n,] <- partial_contig[1,]
                      edited_sample_cns$start[n] <- min(partial_contig$start)
                      edited_sample_cns$end[n] <- max(partial_contig$end)
                      #edited_sample_cns$fragments[n] <- nrow(partial_contig)
                      edited_sample_cns$fragments[n] <- length(unique(partial_contig$chromosome))
                      ## orientation
                      orientation <- c(partial_contig$orientation)
                      plus_count <- sum(orientation == "plus")
                      minus_count <- sum(orientation == "minus")
                      if(plus_count > minus_count){
                        edited_sample_cns$orientation[n] <- "plus"
                        edited_sample_cns$contig_start[n] <- min(partial_contig$contig_start)
                        edited_sample_cns$contig_end[n] <- max(partial_contig$contig_end) }
                      if(plus_count < minus_count){
                        edited_sample_cns$orientation[n] <- "minus"
                        edited_sample_cns$contig_start[n] <- max(partial_contig$contig_start)
                        edited_sample_cns$contig_end[n] <- min(partial_contig$contig_end)  }
                      if(plus_count == minus_count){
                        edited_sample_cns$orientation[n] <- "undecided"
                        edited_sample_cns$contig_start[n] <- min(partial_contig$contig_start)
                        edited_sample_cns$contig_end[n] <- max(partial_contig$contig_end)}
                      ##
                      n <- n+1
                    }
                  }
                  contig <- contig[0,]
                } else { contig <- contig[0,] }
              } else { contig <- contig[0,] }
              
            } }

        }
        
        if(nrow(contig) != 0){
          edited_sample_cns[n,] <- contig[1,]
          edited_sample_cns$start[n] <- min(contig$start)
          edited_sample_cns$end[n] <- max(contig$end)
          #edited_sample_cns$fragments[n] <- nrow(contig)
          edited_sample_cns$fragments[n] <- length(unique(contig$chromosome))
          orientation <- c(contig$orientation)
          plus_count <- sum(orientation == "plus")
          minus_count <- sum(orientation == "minus")
          if(plus_count > minus_count){
            edited_sample_cns$orientation[n] <- "plus"
            edited_sample_cns$contig_start[n] <- min(contig$contig_start)
            edited_sample_cns$contig_end[n] <- max(contig$contig_end)  }
          if(plus_count < minus_count){
            edited_sample_cns$orientation[n] <- "minus"
            edited_sample_cns$contig_start[n] <- max(contig$contig_start)
            edited_sample_cns$contig_end[n] <- min(contig$contig_end)   }
          if(plus_count == minus_count){
            edited_sample_cns$orientation[n] <- "undecided"
            edited_sample_cns$contig_start[n] <- min(contig$contig_start)
            edited_sample_cns$contig_end[n] <- max(contig$contig_end)   }
          n <- n+1 }
      }
      
      edited_sample_cns <- na.omit(edited_sample_cns)
      
      edited_sample_cns <- edited_sample_cns %>% arrange(start)
      
      sample_cns <- edited_sample_cns
      
      sample_cns <- sample_cns[!duplicated(sample_cns), ]
      
    })
    
    
    output$unprocessed_contigs_table <- renderDataTable(
      unprocessed_hits_reactive_edited(), filter="top", escape = FALSE, 
      options = list(lengthChange = FALSE, searching=TRUE, paging=FALSE,
                     columnDefs = list(list(className = 'dt-center', targets = 2:3)),dom="t"))
    
    
    output$download_unprocessed_contigs_plot <- downloadHandler(
      filename = function(){paste("unprocessed_contigs_plot_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
      content = function(file){ ggsave(file, plot=unprocessed_hits_ggplot(), dpi = 600, width = 16, height=10) } )
    
    output$download_unprocessed_contigs_table <- downloadHandler(
      filename = function(){paste("uprocessed_contigs_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(unprocessed_hits_reactive_edited(), fname) })
    
    
############################################################################### 
###############################################################################
###############################################################################
###############################################################################
###############################################################################     

        
    all_contigs_ggplot <- reactive({
      req(input$plot_orientation)
      
      sample_cns <- all_contigs_reactive()
      sample_cns <- data.frame(sample_cns)
      
      if(input$plot_orientation == "up_to_down"){
        sample_cns <- sample_cns %>% arrange(desc(start)) 
      }
      
      chrom_sizes <- structure(list(
        chromosome = unique(sample_cns$gene), 
        size = rep(c(max(sample_cns$size)),length(unique(sample_cns$gene)))),
        .Names = c("chromosome", "size"),
        class = "data.frame", row.names= c(NA,length(unique(sample_cns$gene))))
      chrom_order <- c(unique(chrom_sizes$chromosome))
      chrom_key <- setNames(object = as.character(c(seq(1,length(unique(sample_cns$gene)),1))),nm = chrom_order)
      
      # chrom_order <- factor(x = chrom_order, levels = rev(chrom_order)) # reverse
      chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
      chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], levels = chrom_order)
      sample_cns[["gene"]] <- factor(x = sample_cns[["gene"]], levels = chrom_order)
      group.colors <- c(gain = "red", loss = "blue")
      
      
      all_contigs_ggplot_sub <- ggplot(data = chrom_sizes) + 
        geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                      xmax = as.numeric(chromosome) + 0.2, 
                      ymax = size, ymin = 0), 
                  colour="white", fill = "white") + 
        coord_flip() + 
        theme( 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_x_discrete(name = "Contigs", #limits = names(chrom_key)
        ) +
        geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                         xmax = as.numeric(gene) + 0.2, 
                                         ymax = end, ymin = start, fill = CNA)) + 
        scale_fill_manual(values = group.colors) +
        #ggtitle("Contigs mapping on scaffold") +
        scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = max(sample_cns$size)/50),1)) +
        ylab("Scaffold (bp)")  +
        theme(legend.position = "none")
      
      
      all_contigs_ggplot_sub
      
      
    })
    
    output$all_contigs_plot <- renderPlot({ all_contigs_ggplot() })
    
    
    output$download_all_contigs_plot <- downloadHandler(
      filename = function(){paste("processed_contigs_plot_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
      content = function(file){ ggsave(file, plot=all_contigs_ggplot(), dpi = 600, width = 16, height=10) } )
    
    
    
    
###########################################################################
############################################################################### 
###############################################################################
###############################################################################
###############################################################################
###############################################################################
    
    red_blue_contigs_reactive <- reactive({
      req(input$file)
      req(input$red_blue_threshold)
      sample_cns <- all_contigs_reactive()
      sample_cns <- data.frame(sample_cns)
      
      gene_frequency <- sample_cns %>% count(gene)
      
      gene_frequency$gene <- factor(gene_frequency$gene)
      sample_cns$gene <- factor(sample_cns$gene)
      
      
      for(i in seq(1,nrow(sample_cns),1)){
        included_rows <- intersect(which(sample_cns$start < sample_cns$start[i]),
                                   which(sample_cns$end > sample_cns$end[i]))
        if(length(included_rows) != 0){
          sample_cns$CNA[i] <- "loss"  }
        
        semi_included_rows_equal_start <- intersect(which(sample_cns$start == sample_cns$start[i]),
                                                    which(sample_cns$end > sample_cns$end[i]))
        if(length(semi_included_rows_equal_start) != 0){
          sample_cns$CNA[i] <- "loss"  }
        semi_included_rows_equal_end <- intersect(which(sample_cns$start < sample_cns$start[i]),
                                                  which(sample_cns$end == sample_cns$end[i]))
        if(length(semi_included_rows_equal_end) != 0){
          sample_cns$CNA[i] <- "loss"  }
        
        
        if(sample_cns$fragments[i] <= red_contig_fragments_reactive() ){
          sample_cns$CNA[i] <- "loss" }
        
        repeated_rows <- intersect(which(sample_cns$start == sample_cns$start[i]), 
                                   which(sample_cns$end == sample_cns$end[i]))
        if(length(repeated_rows) != 1){
          repeated_rows_frequency <- gene_frequency %>% 
            filter(gene_frequency$gene %in% levels(droplevels(sample_cns$gene[repeated_rows])))
          row_of_interest <- repeated_rows_frequency[which(repeated_rows_frequency$n == min(repeated_rows_frequency$n), arr.ind = TRUE)[1],]
          rows_to_eliminate <- repeated_rows[-which(repeated_rows_frequency$gene == row_of_interest[,1])]
          if(length(rows_to_eliminate)  != 0){
            for(rr in rows_to_eliminate){
              sample_cns$CNA[rr] <- "loss"
            }}}
      }
      
      sample_cns
      
    })
    
    
###############################################################################
###############################################################################    
    
    red_blue_contigs_ggplot <- reactive({
      req(input$plot_orientation)
      
      sample_cns <- red_blue_contigs_reactive()
      sample_cns <- data.frame(sample_cns)
      
      if(input$plot_orientation == "up_to_down"){
        sample_cns <- sample_cns %>% arrange(desc(start)) 
      }
      
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
      
      
      
      red_blue_contigs_ggplot_sub <- ggplot(data = chrom_sizes) + 
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
        scale_x_discrete(name = "Contigs", #limits = names(chrom_key)
        ) +
        geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                         xmax = as.numeric(gene) + 0.2, 
                                         ymax = end, ymin = start, fill = CNA)) + 
        scale_fill_manual(values = group.colors) +
        #ggtitle("Contigs mapping on scaffold") +
        scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = max(sample_cns$size)/50),1)) +
        ylab("Scaffold (bp)")  +
        theme(legend.position = "none")
      
      red_blue_contigs_ggplot_sub
      
    })
    
    output$red_blue_contigs_plot <- renderPlot({ red_blue_contigs_ggplot() })
    
    
    red_blue_contigs_reactive_edited <- reactive({
      sample_cns <- red_blue_contigs_reactive()
      sample_cns <- data.frame(sample_cns)
      sample_cns <- sample_cns %>% add_column(query_size = rep(NA,nrow(sample_cns)), .after = "end")
      sample_cns <- sample_cns %>% add_column(ctg_alg_len = rep(NA,nrow(sample_cns)), .after = "contig_end")
      for(csize in seq(1,nrow(sample_cns),1)){
        sample_cns$query_size[csize] <- (sample_cns$end[csize] - sample_cns$start[csize]) + 1 
        sample_cns$ctg_alg_len[csize] <- abs(sample_cns$contig_end[csize] - sample_cns$contig_start[1] +1) }
      #sample_cns <- sample_cns[,-10]
      #sample_cns <- sample_cns[,-8]
      sample_cns <- sample_cns[,-8]
      colnames(sample_cns) <- c("Contig","Query","Query code","Query coordinates",
                                "Query start","Query end","Alignment query size",
                                "Filter","Reference size","Hit fragments","Orientation",
                                "Contig start","Contig end","Alignment contig size","Contig size")
      
      sample_cns <- sample_cns %>% 
        mutate(Filter = ifelse(as.character(Filter) == "gain", "red", as.character(Filter)))
      sample_cns <- sample_cns %>% 
        mutate(Filter = ifelse(as.character(Filter) == "loss", "blue", as.character(Filter)))
      
      sample_cns
    })
    
    output$all_contigs_table <- renderDataTable(
      red_blue_contigs_reactive_edited(), filter="top", escape = FALSE, 
      options = list(lengthChange = FALSE, searching=TRUE, paging=FALSE,
                     columnDefs = list(list(className = 'dt-center', targets = 2:3)),dom="t"))
    
    output$download_red_blue_contigs_plot <- downloadHandler(
      filename = function(){paste("all_red_blue_contigs_plot_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
      content = function(file){ ggsave(file, plot=red_blue_contigs_ggplot(), dpi = 600, width = 16, height=10) } )
    
    output$download_all_contigs_table <- downloadHandler(
      filename = function(){paste("processed_contigs_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(red_blue_contigs_reactive_edited(), fname) }) 
    
    
###########################################################################
############################################################################### 
###############################################################################
###############################################################################
###############################################################################
###############################################################################
    
    
    minimal_contigs_reactive <- reactive({
      req(input$file)
      
      sample_cns <- red_blue_contigs_reactive()
      sample_cns <- data.frame(sample_cns)
      sample_cns <- sample_cns %>% filter(CNA == "gain")
      
      for(t in unique(sample_cns$gene)){
        contig <- sample_cns %>% filter(gene == t)
        if(nrow(contig) > 1){
          
          plus_subtable <-  contig %>% filter(orientation == "plus")
          minus_subtable <-  contig %>% filter(orientation == "minus")
          
          if(sum(plus_subtable$fragments) > sum(minus_subtable$fragments)){
            sample_cns$orientation[which(sample_cns$gene == t)] <- c("plus") }
          if(sum(plus_subtable$fragments) < sum(minus_subtable$fragments)){
            sample_cns$orientation[which(sample_cns$gene == t)] <- c("minus") }
          if(sum(plus_subtable$fragments) == sum(minus_subtable$fragments)){
            sample_cns$orientation[which(sample_cns$gene == t)] <- c("undecided") }
        }
      }
      
      sample_cns
      
    })
    
    minimal_contigs_ggplot <- reactive({
      req(input$plot_orientation)
      
      sample_cns <- minimal_contigs_reactive()
      sample_cns <- data.frame(sample_cns)
      
      coverage_lengths <- c()
      i <- 1
      while(i < nrow(sample_cns)){
        if(sample_cns$start[i+1] <= sample_cns$end[i]){
          for(j in seq(i+1, nrow(sample_cns),1)){
            if(j == nrow(sample_cns)){
              local_length <- sample_cns$end[j] - sample_cns$start[i]
              coverage_lengths <- append(coverage_lengths, local_length)
              i <- j+1
              break }
            else{
              if(sample_cns$start[j+1] <= sample_cns$end[j]){ next }
              else{
                local_length <- sample_cns$end[j] - sample_cns$start[i]
                coverage_lengths <- append(coverage_lengths, local_length)
                i <- j+1
                break }
            } } }
        else{
          local_length <- sample_cns$end[i] - sample_cns$start[i]
          coverage_lengths <- append(coverage_lengths, local_length)
          i <- i + 1
        } }
      coverage <- sum(coverage_lengths) / sample_cns$size[1]
      
      if(input$plot_orientation == "up_to_down"){
        sample_cns <- sample_cns %>% arrange(desc(start)) 
      }
      
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
      
      minimal_contigs_ggplot_sub <- ggplot(data = chrom_sizes) + 
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
        scale_x_discrete(name = "Contigs", limits = names(chrom_key)
        ) +
        geom_rect(data = sample_cns, aes(xmin = as.numeric(gene) - 0.2, 
                                         xmax = as.numeric(gene) + 0.2, 
                                         ymax = end, ymin = start, fill = CNA)) + 
        scale_fill_manual(values = group.colors) +
        ggtitle(paste("Contigs mapping with coverage =",round(coverage,5))) +
        scale_y_continuous(breaks =  round(seq(0, max(sample_cns$size), by = max(sample_cns$size)/50),1)) +
        ylab("Scaffold (bp)")  +
        theme(legend.position = "none")
      
      minimal_contigs_ggplot_sub
      
    })
    
    observeEvent(input$scaffolds, {
      #minimal_contigs_reactive <- minimal_contigs_reactive()
      if(nrow(minimal_contigs_reactive()) > 0){
      output$minimal_contigs_plot <- renderPlot({ minimal_contigs_ggplot() })
      output$error_message_plot <- renderText({ "" })
      } else {
        output$minimal_contigs_plot <- renderPlot({ NULL # Clear the table output 
          }) 
        output$error_message_plot <- renderText({ "Plot could not be created because no contigs passed the conditional arguments." })
      }
    })
    
    minimal_contigs_reactive_edited <- reactive({
      sample_cns <- minimal_contigs_reactive()
      if(nrow(sample_cns) > 0){
      sample_cns <- data.frame(sample_cns)
      sample_cns <- sample_cns %>% add_column(ctg_size = rep(NA,nrow(sample_cns)), .after = "end")
      for(csize in seq(1,nrow(sample_cns),1)){
        sample_cns$ctg_size[csize] <- (sample_cns$end[csize] - sample_cns$start[csize]) + 1 }
      #sample_cns <- sample_cns[,-10]
      sample_cns <- sample_cns[,-9]
      sample_cns <- sample_cns[,-8]
      sample_cns <- sample_cns[,-4]
      colnames(sample_cns) <- c("Contig","Query","Query code",
                                "Query start","Query end","Alignments lengths",
                                "Subject size","Hit fragments","Orientation",
                                "Contig start","Contig end","Contig length")
      }
      sample_cns
    })
    
    observeEvent(input$scaffolds, {
      #minimal_contigs_reactive_edited <- minimal_contigs_reactive_edited()
    
      if(nrow(minimal_contigs_reactive_edited()) > 0 && !is.null(minimal_contigs_reactive_edited())){
          output$minimal_contigs_table <- renderDataTable(
          minimal_contigs_reactive_edited(), filter="top", escape = FALSE, 
          options = list(lengthChange = FALSE, searching=TRUE, paging=FALSE,
                     columnDefs = list(list(className = 'dt-center', targets = 2:3)),dom="t"))
          output$error_message_table <- renderText({ "" })
      } else {
        output$minimal_contigs_table <- renderDataTable( NULL )
        output$error_message_table <- renderText({ "Table could not be created because no contigs passed the conditional arguments." })
      }
    })

    
#    output$minimal_contigs_table <- renderDataTable(
#      minimal_contigs_reactive_edited(), filter="top", escape = FALSE, 
#                options = list(lengthChange = FALSE, searching=TRUE, paging=FALSE,
#                               columnDefs = list(list(className = 'dt-center', targets = 2:3)),dom="t"))
    
    output$download_minimal_contigs_plot <- downloadHandler(
      filename = function(){paste("minimal_contigs_plot_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
      content = function(file){ ggsave(file, plot=minimal_contigs_ggplot(), dpi = 600, width = 16, height=10) } )
    
    output$download_minimal_contigs_table <- downloadHandler(
      filename = function(){paste("minimal_contigs_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(minimal_contigs_reactive_edited(), fname) }) 
    
    
    orientation_merged_reactive <- reactive({
      sample_cns <- minimal_contigs_reactive()
      table_orientation <- table_orientation_reactive()
      
      merged_orientation <- subset(table_orientation, Contig %in% sample_cns$gene)
      merged_orientation
    })
    
    orientation_merged_reactive_edited <- reactive({
      merged_orientation_edited <- orientation_merged_reactive()
      colnames(merged_orientation_edited) <- c("Contig","Plus fragments","Minus fragments","Contig length")
      merged_orientation_edited
    })
    
    output$download_merged_orientation_table <- downloadHandler(
      filename = function(){paste("merged_orientation_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(orientation_merged_reactive_edited(), fname) }) 
    
    
    orientation_undecided_reactive <- reactive({
      merged_orientation <- orientation_merged_reactive()
      undecided_orientation <- subset(merged_orientation, Plus_fragments != 0 & Minus_fragments != 0)
      colnames(undecided_orientation) <- c("Contig","Plus fragments","Minus fragments","Contig length")
      undecided_orientation
    })
    
    output$download_undecided_orientation_table <- downloadHandler(
      filename = function(){paste("undecided_orientation_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(orientation_undecided_reactive(), fname) }) 
    
    
    
###########################################################################
############################################################################## 
###############################################################################
###############################################################################
###############################################################################
###############################################################################
    
    
    included_contigs_reactive <- reactive({
      
      req(input$file)
      sample_cns <- red_blue_contigs_reactive()
      sample_cns <- data.frame(sample_cns)
      
      sample_cns_included_contigs <- sample_cns %>% filter(CNA == "loss")
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
      
      included_contigs
      
      
    })
    
    
    output$included_contigs_table <- renderDataTable(
      included_contigs_reactive(), escape = FALSE, 
      options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                     columnDefs = list(list(className = 'dt-center', targets = 2:3)),dom="t"))
    
    output$download_included_contigs_table <- downloadHandler(
      filename = function(){paste("included_contigs_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(included_contigs_reactive(), fname) }) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
      
    }
  
  
  
  
  
  
  
  
}




# shiny app
shinyApp(ui=ui, server=server)
