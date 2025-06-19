# app/app.R

library(shiny)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)

# -----------------------------------------------
# ✅ Load precomputed data
# -----------------------------------------------
counts <- read.csv("../data/raw/count_matrix_with_symbols_clean.csv", row.names = 1, check.names = FALSE)
degs <- read.csv("../data/raw/DEG_results_mapped.csv")
gene_list <- unique(degs$SYMBOL)

# 📦 Sample metadata
metadata <- data.frame(
  sample = colnames(counts),
  group = ifelse(grepl("01A", colnames(counts)), "Tumor", "Normal")
)
metadata$group <- factor(metadata$group, levels = c("Normal", "Tumor"))

# -----------------------------------------------
# ⚙️ GO Enrichment Function
# -----------------------------------------------
get_go_plot <- function(gene_symbols) {
  entrez_ids <- bitr(gene_symbols,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  dotplot(ego, showCategory = 15) + ggtitle("GO Biological Process Enrichment")
}

# -----------------------------------------------
# ⚙️ KEGG Enrichment Function
# -----------------------------------------------
get_kegg_plot <- function(gene_symbols) {
  entrez_ids <- bitr(gene_symbols,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                      organism = 'hsa',
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)
  dotplot(ekegg, showCategory = 15) + ggtitle("KEGG Pathway Enrichment")
}

# -----------------------------------------------
# 🖼️ UI
# -----------------------------------------------
ui <- fluidPage(
  titlePanel("TCGA BRCA Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_gene", "Select a Gene:", choices = gene_list, selected = "TP53")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Boxplot", 
                 plotOutput("boxplot"), 
                 downloadButton("download_boxplot", "Download Boxplot")
        ),
        tabPanel("Volcano Plot", 
                 plotOutput("volcano"),
                 downloadButton("download_volcano", "Download Volcano Plot")
        ),
        tabPanel("GO Enrichment", 
                 plotOutput("go_plot"),
                 downloadButton("download_go", "Download GO Plot")
        ),
        tabPanel("KEGG Enrichment", 
                 plotOutput("kegg_plot"),
                 downloadButton("download_kegg", "Download KEGG Plot")
        ),
        tabPanel("About the Data", 
                 uiOutput("about_tab"),
                 br(),
                 downloadButton("download_metadata", "Download Metadata CSV")
        )
      )
    )
  ),
  
  tags$hr(),
  div(style = "text-align: center; font-size: 12px; color: gray;",
      "Developed by Prince Mishra © 2025 | For academic use only"
  )
)

# -----------------------------------------------
# 🧠 Server
# -----------------------------------------------
server <- function(input, output) {
  
  # 🎯 Boxplot Reactive
  boxplot_reactive <- reactive({
    gene <- input$selected_gene
    if (!(gene %in% rownames(counts))) return(NULL)
    
    gene_data <- data.frame(
      expression = as.numeric(counts[gene, ]),
      sample = colnames(counts),
      group = metadata$group
    )
    
    ggplot(gene_data, aes(x = group, y = expression, fill = group)) +
      geom_boxplot() +
      labs(title = paste("Expression of", gene),
           y = "Normalized Expression", x = "") +
      theme_minimal()
  })
  
  output$boxplot <- renderPlot({
    p <- boxplot_reactive()
    if (is.null(p)) {
      plot.new()
      title("Selected gene not found in count matrix")
    } else {
      p
    }
  })
  
  output$download_boxplot <- downloadHandler(
    filename = function() { paste0(input$selected_gene, "_boxplot.png") },
    content = function(file) {
      png(file)
      print(boxplot_reactive())
      dev.off()
    }
  )
  
  # 🎯 Volcano Plot Reactive
  volcano_plot <- reactive({
    EnhancedVolcano(degs,
                    lab = degs$SYMBOL,
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title = 'Volcano Plot: Tumor vs Normal',
                    pCutoff = 0.05,
                    FCcutoff = 1)
  })
  
  output$volcano <- renderPlot({
    volcano_plot()
  })
  
  output$download_volcano <- downloadHandler(
    filename = "volcano_plot.png",
    content = function(file) {
      png(file, width = 1000, height = 800)
      print(volcano_plot())
      dev.off()
    }
  )
  
  # 🧬 GO Plot
  output$go_plot <- renderPlot({
    sig_genes <- degs$SYMBOL[degs$pvalue < 0.05 & abs(degs$log2FoldChange) > 1]
    get_go_plot(sig_genes)
  })
  
  output$download_go <- downloadHandler(
    filename = "go_enrichment.png",
    content = function(file) {
      png(file, width = 1000, height = 800)
      sig_genes <- degs$SYMBOL[degs$pvalue < 0.05 & abs(degs$log2FoldChange) > 1]
      p <- get_go_plot(sig_genes)
      print(p)
      dev.off()
    }
  )
  
  # 🧬 KEGG Plot
  output$kegg_plot <- renderPlot({
    sig_genes <- degs$SYMBOL[degs$pvalue < 0.05 & abs(degs$log2FoldChange) > 1]
    get_kegg_plot(sig_genes)
  })
  
  output$download_kegg <- downloadHandler(
    filename = "kegg_enrichment.png",
    content = function(file) {
      png(file, width = 1000, height = 800)
      sig_genes <- degs$SYMBOL[degs$pvalue < 0.05 & abs(degs$log2FoldChange) > 1]
      p <- get_kegg_plot(sig_genes)
      print(p)
      dev.off()
    }
  )
  
  # 📤 Download Metadata
  output$download_metadata <- downloadHandler(
    filename = "metadata.csv",
    content = function(file) {
      write.csv(metadata, file, row.names = FALSE)
    }
  )
  
  # 📘 About Tab
  output$about_tab <- renderUI({
    tagList(
      h3("🧾 About the Dataset"),
      p("This app analyzes RNA-Seq data from the TCGA Breast Cancer (BRCA) dataset."),
      
      h4("📂 Sample Information"),
      tags$ul(
        tags$li("Samples were downloaded via ", code("TCGAbiolinks"), " and include primary tumor and normal breast tissue."),
        tags$li("For demonstration purposes, a subset of samples was used.")
      ),
      
      h4("🧬 Gene Expression & DEG Analysis"),
      p("DESeq2 was used to identify differentially expressed genes. Genes with adjusted p-value < 0.05 and |log2FC| > 1 were considered significant."),
      
      h4("🔍 Enrichment Analysis"),
      p("GO and KEGG enrichment was performed using ", code("clusterProfiler"), "."),
      
      h4("⚠️ Disclaimer"),
      p("This app is for educational and exploratory use only. Results may vary with different datasets or sample sizes.")
    )
  })
}

# -----------------------------------------------
# 🚀 Run the App
# -----------------------------------------------
shinyApp(ui = ui, server = server)
