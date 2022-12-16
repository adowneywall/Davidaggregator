#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)
library(shiny)
library(shinyjs)

###################################################
################ Functions ########################
###################################################

## Function for spliting term strings 
# pos = 1 - The term ID
# pos = 2 - The term name
term_split <- function(x,pos=1){
  if(pos==1){
    y <- ifelse(isTRUE(grepl('GO:',x)),unlist(strsplit(x,split = '~'))[1],unlist(strsplit(x,split = '~|:'))[1])
  }
  if(pos==2){
    y <- ifelse(isTRUE(grepl('GO:',x)),unlist(strsplit(x,split = '~'))[2],unlist(strsplit(x,split = '~|:'))[2])
  }
  return(y)
}

## Removes white spaces from strings within a vector
wsr <- function(x){
  gsub("\\s", "",x)
}

## Wrapper function for processing an individual file
file_process <- function(input){
  name <- input$name
  file <- input$datapath
  delim='[_.]'
  dat <- read_lines(file)
  if(length(dat) > 1){
    GeneSetID <- name#str_split(basename(file),delim,simplify = T)[pos]
    
    ## Select cluster and enrichment information ##
    dat_annotation <- grep('Annotation',dat,value = T)
    cluster <- data.frame(str_split(dat_annotation,'\t| ',simplify = T))
    cluster_label <- data.frame(Annotation_Cluster=as.numeric(cluster$X3),Cluster_Score=as.numeric(cluster$X6))
    
    ## Processing category / term data
    c_sum <- clusters_process(dat,file)
    
    y <- cbind(GeneSetID,cluster_label,c_sum)
    
    y <- y %>%
      mutate(All_Unique_Terms_Name_Copy=All_Unique_Terms_Name) %>%
      relocate(All_Unique_Terms_Name_Copy,.after = All_Unique_Terms_ID) %>%
      relocate(Top_Genes,.after = Cluster_Score) %>%
      relocate(Top_Count,.after = Cluster_Score) %>%
      relocate(Top_FDR,.after = Cluster_Score) %>%
      relocate(All_Unique_Terms_Name,.after = Cluster_Score) %>%
      relocate(Top_Term,.after = Cluster_Score)
  }else{y <- NULL}
  
  return(y)
}

## Process each DAVID enrichment cluster
clusters_process <- function(dat,file){
  
  block_start <- grep('Annotation',dat)
  block_end <- which(dat == "")-2
  
  header <- c("Category","Term","Count","Percent","PValue","Genes","List_Total",
              "Pop_Hits","Pop_Total","Fold_Enrichment","Bonferroni","Benjamini","FDR")
  top_header <- paste0("Top_",header)
  
  for(i in 1:length(block_start)){
    out <- read_tsv(file = file,skip = block_start[i],n_max = block_end[i] - block_start[i])
    # Top process for each enrichment cluster   
    top_process <- out[1,] 
    # Relabel
    colnames(top_process) <- top_header
    top_process <- top_process %>%
      mutate(Top_Term_ID = term_split(Top_Term,pos=1),
             Top_Term_Name = str_to_sentence(term_split(Top_Term,pos=2))) %>%
      relocate(Top_Term_ID,.after=Top_Term) %>%
      relocate(Top_Term_Name,.after=Top_Term_ID)
    
    ## Add summary of terms within cluster ##
    top_process$All_Terms <- paste(out$Term,collapse=',')
    top_process$All_Terms_ID <- paste(sapply(out$Term,term_split,pos=1),collapse=',')
    top_process$All_Terms_Name <- paste(str_to_sentence(sapply(out$Term,term_split,pos=2)),collapse=',')
    top_process$All_Terms_Count <- length(out$Term)
    top_process$All_Unique_Terms_ID <- paste(unique(toupper(sapply(out$Term,term_split,pos=1))),collapse=',')
    top_process$All_Unique_Terms_Name <- paste(unique(str_to_sentence(sapply(out$Term,term_split,pos=2))),collapse=',')
    top_process$All_Unique_Terms_Name_Count <- length(unique(str_to_sentence(sapply(out$Term,term_split,pos=2))))
    
    ## Add summary of genes within cluster ##
    # wsr = white space remove (function above)
    unique_genes <- unique(wsr(unlist(strsplit(out$Genes,','))))
    
    top_process$All_Unique_Genes_Count <- length(unique_genes)
    top_process$All_Unique_Genes <- paste(unique_genes,collapse = ',')
    
    if(i == 1){
      y <- top_process 
    }else{
      y <- rbind(y,top_process)
    }
  }
  return(y)
}


# Define UI for application that draws a histogram
ui <- fluidPage(
  useShinyjs(),
  # Application title
  titlePanel("DAVID file aggregator"),
  
  ## change color for error messages
  tags$head(
    tags$style(HTML("
      #error {
        color: #8B0000;
        white-space:pre-wrap;
      }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "david_file",buttonLabel = "Upload",label="Upload DAVID .txt file(s)",multiple = TRUE),
      helpText(HTML("<i>Note 1: Files should be in the default format prodivded by the DAVID enrichment analysis tool <i/>")),
      helpText(HTML("<i>Note 2: Select multiple files by using either control or shift click <i/>")),
      tags$a(href="https://david.ncifcrf.gov/tools.jsp", 
             "Link to DAVID webpage"),
      helpText(""),
      actionButton(inputId = "go", label = "Convert"),
      
      conditionalPanel(
        "false", # always hide the download button
        downloadButton(outputId = "downloadData", "Download tab-delimited (.tsv) file")
      )
    ),
    
    mainPanel(
      verbatimTextOutput("error")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  onSessionEnded(function() {
    list.dirs(recursive=FALSE) %>% 
      keep(function(x) grepl("^./tmp", x)) %>% 
      unlink(recursive = TRUE)
  })
  
  resulted_zip <- reactiveVal()
  
  observeEvent(input$go,{
    # show error if files in not presented
    output$error <- renderText({
      validate(
        need(input$david_file != '', 'Please provide file(s) exported from DAVID')
      )
    })
    
    main_dir <- getwd()
    
    ## field should not be empty
    shiny::req(input$david_file)
    
    ## get list of filenames
    input_files <- input$david_file %>% 
      select(name, datapath)
    
    # Diagnostic code 
    # path_files <- '~/Desktop/WaxmanLab/2022_DAVIDEnrichmentOfHubLncRNAs/DAVIDEnrichmentFolders/DAVID_Chow'
    # input_files <- data.frame(name=list.files(path_files,full.names = F),datapath = list.files(path_files,full.names = T))

    ## create temporary directory
    zipdir <- tempfile("tmp", tmpdir = "./")
    if (!dir.exists(zipdir)) {
      dir.create(zipdir)  
    }
    setwd(zipdir)
    
    # conversion
    withProgress(message = "Processing files",
                 value = 0,
                 {
                   n <- nrow(input_files)
                   for(i in 1:n){
                     if(i == 1){
                       file_process(input_files[i,]) %>%
                         write_tsv(file="aggregated_file.tsv",col_names = T)
                     }else{
                       file_process(input_files[i,]) %>%
                         write_tsv(file="aggregated_file.tsv",col_names = F,append = T)
                     }
                     incProgress(1/n)
                   }

                   setProgress(detail = "Zipping files...")
                   zip("data.zip","aggregated_file.tsv")
                 }
    )
    
    resulted_zip(paste0(zipdir,"/data.zip"))
    setwd(main_dir)
    
    runjs("$('#downloadData')[0].click();")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("David_Aggregation_", format(Sys.time(), "%a_%b_%d_%Y_%I-%M%p"), ".zip", sep="")
    },
    content = function(file) {
      file.copy(resulted_zip(), file)
    },
    contentType = "application/zip")
}

# Run the application 
shinyApp(ui = ui, server = server)