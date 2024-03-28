library(ape)
library(dendextend)
library(tidyverse)
library(ggtree)
library(gridExtra)
library(ggtreeExtra)
library(svglite)

tree <- function(df, table_name,n = 10) { 
  # Loop through each column starting with 'score'
  for (col_name in names(df)[startsWith(names(df), "score")]) {
    # Select the column and 'other_column'
    sub_df <- df %>% select(col_name, preferredNames,description)
    sub_clean_df <- na.omit(sub_df)
    
    arr.data <- sub_clean_df[order(sub_clean_df[[col_name]]), ]
    
    more_than_0 <- arr.data %>% filter(.data[[col_name]] > 0) %>% tail(n)

    less_than_0 <- arr.data %>% filter(.data[[col_name]] < 0) %>% head(n)
    
    up_down <- rbind(more_than_0, less_than_0)
    
    dff <- up_down %>% mutate(correlated = ifelse(.data[[col_name]] > 0,
                                                  paste("more correlated in", gsub("score|vs.*", "", col_name)),
                                                  ifelse(.data[[col_name]] < 0, paste("more correlated in", sub(".*vs\\s*", "", col_name)),
                                                         "No regulation")))
    dff$abs.score <- abs(round(dff[[col_name]],2))
    all_genes <- unique(unlist(dff$preferredNames))
    
    
    binary_matrix <- matrix(0, nrow = nrow(dff), ncol = length(all_genes), dimnames = list(dff$description, all_genes)) 
    n <- 1
    for (i in dff$preferredNames) {
      binary_matrix[dff$description[n], intersect(all_genes, i)] <- 1 
      n <- n + 1
    }

    binary_position_matrix <- (binary_matrix == 1)

    # Compute Jaccard similarity matrix
    matrix <- dist(binary_position_matrix,  method = "binary")
    hclust_object <- hclust(as.dist(matrix), method = "ward.D2" )
    # Create a dendrogram object
    dendrogram_object <- as.dendrogram(hclust_object)
    newick <- as.phylo(dendrogram_object)
    # make dendrogram
    p1 <- ggtree(newick, branch.length="none", layout="circular",open.angle=15)  + 
      geom_tiplab(align=TRUE, linesize = 0, size = 5, offset = 1.5) + ggtitle(col_name)+
      theme(plot.title = element_text(hjust = 0))+ xlim(0, 50) +
      theme(plot.margin = margin(1, 1, 1, 1, "in"))
    
    p2 <- p1 + 
      geom_fruit(
        data = dff,
        geom = geom_col,
        mapping = aes(y=description, x=abs.score, fill=correlated),
        axis.params=list(
          axis       = "x",
          text.size  = 2,
          hjust      = 0,
          vjust      = 0.70,
          line.size = 0.2, line.alpha = 1.2),
        grid.params=list()for (table_name in table_names) {
  table_summary <- CoPPI_results$resultsCoPPI[[table_name]]$table_summary
      ) 
    
    
    ggsave(paste(col_name,table_name,".svg"), p2, device = "svg",width = 100, height = 98, units = "cm", dpi = 500)

  }
}


tree(dataframe,name.output.svg, 10,2)
