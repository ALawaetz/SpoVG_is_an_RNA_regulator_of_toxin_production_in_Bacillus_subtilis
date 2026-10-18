#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")

### Make custom GO terms from SubtiWiki to input ClusterProfiler
# Load necessary libraries
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# Load your custom gene-to-category data
gene2category <- read.csv("/Users/andreaslawaetz/Streptomyces/Bacillus/SpoVG_paper/GSEA/GeneID_Category_All.csv")
colnames(gene2category) <- c("Category", "GeneID")

# Create a TERM2GENE data frame
TERM2GENE <- gene2category %>%
  rename(term = Category, gene = GeneID)


### GSEA plots
################################################################################

run_enrichment <- function(gene_list_file, plot_title = "Enrichment analysis") {
  
  # Read gene list
  genes_of_interest <- read.table(
    gene_list_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )$V1
  
  # Remove missing, empty, and duplicate gene IDs
  genes_of_interest <- unique(genes_of_interest)
  genes_of_interest <- genes_of_interest[
    !is.na(genes_of_interest) & genes_of_interest != ""
  ]
  
  # Perform enrichment analysis
  go_enrichment <- enricher(
    gene = genes_of_interest,
    TERM2GENE = TERM2GENE
  )
  
  # Create dot plot
  enrichment_plot <- dotplot(
    go_enrichment,
    showCategory = 20
  ) +
    ggplot2::labs(title = plot_title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  # Display plot
  print(enrichment_plot)
  
  # Return results and plot
  return(list(
    enrichment = go_enrichment,
    plot = enrichment_plot
  ))
}

results_CRAC_M9 <- run_enrichment(
  gene_list_file = "/Users/andreaslawaetz/Streptomyces/Bacillus/GitHub_Thesis-main/SpoVG_chapter/Supplementary_tables/S1_SpoVG_M9_CRAC_thesis_table_BSU_numbers_targets.txt",
  plot_title = "CRAC (M9)"
)

results_CRAC_LB <- run_enrichment(
  gene_list_file = "/Users/andreaslawaetz/Streptomyces/Bacillus/GitHub_Thesis-main/SpoVG_chapter/Supplementary_tables/S2_SpoVG_LB_CRAC_thesis_table.txt_BSU_numbers_targets.txt",
  plot_title = "CRAC (LB)"
)

results_RNAseq_M9 <- run_enrichment(
  gene_list_file = "/Users/andreaslawaetz/Streptomyces/Bacillus/GitHub_Thesis-main/SpoVG_chapter/Supplementary_tables/S3_spoVG_M9_RNAseq_thesis_table.txt_BSU_numbers_DEG.txt",
  plot_title = "RNAseq (M9)"
)

results_RNAseq_LB <- run_enrichment(
  gene_list_file = "/Users/andreaslawaetz/Streptomyces/Bacillus/GitHub_Thesis-main/SpoVG_chapter/Supplementary_tables/S4_spoVG_LB_RNAseq_thesis_table.txt_BSU_numbers_DEG.txt",
  plot_title = "RNAseq (LB)"
)




################################################################################

