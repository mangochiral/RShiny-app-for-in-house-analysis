library(shiny)
library(markdown)
library(readr)
library(stringr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(qvalue)
library(DESeq2)
library(genefilter)
library(tidyverse)
library(gridExtra)
library(vegan)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPpeakAnno)
library(GGally)
library(ChIPseeker)
library(clusterProfiler)
library(org.Hs.eg.db)
library(vsn)
library(sva)
library(rGREAT)

# Increase max upload size to 300MB
options(shiny.maxRequestSize = 300 * 1024^2)

# Define UI for data upload app
ui <- fluidPage(
  
  # App title
  titlePanel("Analyzing bulk ATAC data using DESeq2"),
  
  # Navbar with one tabPanel
  navbarPage("Navbar",
             tabPanel("Data",
                      sidebarLayout(
                        sidebarPanel(
                          fileInput("ataqrawcounts", 
                                    "Choose ATAC Counts TSV File",
                                    multiple = FALSE,
                                    accept = c(".csv", ".tsv")),
                          fileInput("design", 
                                    "Choose Experiment Design TSV File",
                                    multiple = FALSE,
                                    accept = c(".csv", ".tsv")),
                          tags$hr(),
                          radioButtons("disp", "Display",
                                       choices = c(Head = "head",
                                                   All = "all"),
                                       selected = "head")
                        ),
                        mainPanel(
                          tableOutput("counts"),
                          tableOutput("design")
                        )
                      )
             ),
             tabPanel("Exploratory Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("plottype", "Choose Transformation",
                                       c("Variance Stabilizing (VST)" = "vsd",
                                         "Regularized Log (rlog)" = "rld",
                                         "Z-score Normalization" = "norm"),
                                       selected = "norm"),
                          # Button
                          downloadButton("download", "Download")
                          
                        ),
                        mainPanel(
                          plotOutput("plot")
                        )
                      )
             )
  )
)

# Define server logic
server <- function(input, output) {
  
  # --- Reactive function to read ATAC counts file --- #
  countsData <- reactive({
    req(input$ataqrawcounts)
    # read_tsv does not use `header/sep/quote`, so remove or use read_delim if you need custom separators.
    df <- read_tsv(input$ataqrawcounts$datapath) 
    df
  })
  
  # --- Reactive function to read design file --- #
  designData <- reactive({
    req(input$design)
    df <- read_tsv(input$design$datapath) 
    df
  })
  
  # --- Show small preview of the uploaded tables --- #
  output$counts <- renderTable({
    df <- countsData()
    if (input$disp == "head") {
      return(head(df))
    } else {
      return(df)
    }
  })
  
  output$design <- renderTable({
    df <- designData()
    if (input$disp == "head") {
      return(head(df))
    } else {
      return(df)
    }
  })

  # --- Reactive function to make DESeq object --- #
  deATAC <- reactive({
    req(countsData(), designData())
    # Helper to turn the raw table into a proper count matrix
    makeDesq2counts <- function(count_input){
      count_input <- as.data.frame(count_input)
      # Use chr, start, end as rownames
      rownames(count_input) <- paste(count_input$chr, 
                                     count_input$start, 
                                     count_input$end, 
                                     sep = "-")
      # Remove the first 3 columns
      count_input <- count_input[,-c(1,2,3)]
      # Convert to matrix
      count_input <- as.matrix(count_input)
      # Filter out rows whose max count is <= 10
      count_input <- subset(count_input, apply(count_input, 1, max) > 10)
      count_input
    }
    # Read data and prepare
    dfCountData <- makeDesq2counts(countsData())
    dfConditions <- designData()
    
    # Create DESeq2 object
    dseq <- DESeqDataSetFromMatrix(countData = dfCountData,
                                  colData = dfConditions,
                                  design = ~ Group)    
    DESeq(dseq)
  })
  # Precompute transformations
  vsd <- reactive({
    req(deATAC())
    vst(deATAC(), nsub = round(nrow(deATAC()) * 0.8), blind = FALSE)
  })
  
  rld <- reactive({
    req(deATAC())
    rlog(deATAC(), blind = FALSE)
  })  
  
  znorm <- reactive({
    req(deATAC())
    dfNormalizedCount <- as.data.frame(counts(deATAC(), normalized = T))
    
  })
  
  znorm_matrix <- reactive({
    req(znorm())
    # Filter low-variance rows using genefilter::varFilter
    varFilter(as.matrix(znorm()), var.cutoff = 0.2, filterByQuantile = TRUE)
  })
  
  pcaObj <- reactive({
    req(znorm_matrix())
    # PCA on transposed data
    prcomp(t(znorm_matrix()), center = TRUE, scale. = TRUE)
  })
  
  # --- Main PCA plot --- #
  output$plot <- renderPlot({

    
    # Switch based on user choice
    choice <- input$plottype
    
    if (choice == "vsd") {
      # Quick PCA from DESeq2
      DESeq2::plotPCA(vsd(), intgroup = "Group")
      
    } else if (choice == "rld") {
      # Quick PCA from DESeq2
      DESeq2::plotPCA(rld(), intgroup = "Group")
      
    } else {
      pcaRes <- pcaObj()
      percentVar <- (pcaRes$sdev^2 / sum(pcaRes$sdev^2)) * 100
      dfConditions <- designData()
      # Basic ggplot2-based PCA
      qplot(pcaRes$x[,1], pcaRes$x[,2],
            colour = dfConditions$Group,
            xlab = paste0("PC1: ", round(percentVar[1], 1), "% variance"),
            ylab = paste0("PC2: ", round(percentVar[2], 1), "% variance"),
            size = I(3)) +
        scale_color_discrete(name = "") +
        theme_classic()
    }
  })
  #Download plots
  output$downloadData <- downloadHandler(
    filename = function() {
     paste0(input$plottype, "_Pca.pdf")
    },
    content = function(file) {
      # You can open a graphic device or use ggsave.
      # For the built-in DESeq2::plotPCA, we often do:
      
      choice <- input$plottype
      if (choice == "vsd") {
        # Generate the PCA plot object from VSD
        p <- DESeq2::plotPCA(vsd(), intgroup = "Group")
        # Save with ggsave
        ggsave(file, plot = p, width = 8, height = 6)
        
      } else if (choice == "rld") {
        p <- DESeq2::plotPCA(rld(), intgroup = "Group")
        ggsave(file, plot = p, width = 8, height = 6)
        
      } else {
        pcaRes <- pcaObj()
        percentVar <- (pcaRes$sdev^2 / sum(pcaRes$sdev^2)) * 100
        dfConditions <- designData()
        
        p <- qplot(pcaRes$x[,1], pcaRes$x[,2],
                   colour = dfConditions$Group,
                   xlab = paste0("PC1: ", round(percentVar[1], 1), "% variance"),
                   ylab = paste0("PC2: ", round(percentVar[2], 1), "% variance"),
                   size = I(3)) +
          scale_color_discrete(name = "") +
          theme_classic()
        
        ggsave(file, plot = p, width = 8, height = 6)
      }
    }
  )
  
}

# Create Shiny app
shinyApp(ui, server)
