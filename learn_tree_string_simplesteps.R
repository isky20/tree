library(ape)
library(dendextend)
library(tidyverse)
library(ggtree)
library(gridExtra)
library(ggtreeExtra)
#  term_id	Term	source	intersection_size	Genes	FDR	term_size	query_size
# File should have :
## Fold.Enrichment (should be calculated)
## Genes
## Term

tree_string_simple <- function(Table_S6.csv, size, offset) {
  #OPEN FILE
  title <- gsub("_list _enriched.csv", "",Table_S6.csv )
  ref <- read.csv(Table_S6.csv)
  #select top2 fold enrichment
  df <- ref %>% 
    arrange(desc(Fold.Enrichment)) %>% 
    slice(1:20)ì

  
  # build matrix between genes and terms
  all_genes <- unique(unlist(strsplit(df$Genes, ",")))

  binary_matrix <- matrix(0, nrow = nrow(df), ncol = length(all_genes), dimnames = list(df$Term, all_genes)) 
  for (i in seq_along(df$Genes)) {
    gene_list <- unlist(strsplit(df$Genes[i], "(,| )"))
    binary_matrix[df$Term[i], intersect(all_genes, gene_list)] <- 1 #df$term
  }
  binary_position_matrix <- (binary_matrix == 1)
  # Compute Jaccard similarity matrix
  matrix <- dist(binary_position_matrix,  method = "binary")
  hclust_object <- hclust(as.dist(matrix), method = "ward.D2" )
  # Create a dendrogram object
  dendrogram_object <- as.dendrogram(hclust_object)
  newick <- as.phylo(dendrogram_object)
  # make dendrogram
  p1 <- ggtree(newick, branch.length="none", layout="circular") + 
    geom_tiplab(align=TRUE, linesize = 0, size = size, offset = offset) + ggtitle(title)+
    theme(plot.title = element_text(hjust = 0))+ xlim(0, 20)
  # make bar
  
  
  p2 <- p1 + 
    geom_fruit(
      data = df, geom = geom_col, mapping = aes(y=Term, x=  Fold.Enrichment ))+
    theme(plot.margin = margin(2, 2, 2, 2, "in"))
  
  ggsave(paste(title,".svg"), p2, device = "svg",
         width = 16, height = 14, units = "in", dpi = 500)
}


#run the code with many input file
csv_files <- list.files(pattern = "_enriched.csv")
for(i in csv_files) {
  tree_string_simple(i,5,3)
  
}
