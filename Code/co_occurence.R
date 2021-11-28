# Generate the co-occurence map indicating the number of proteins shared between subcellular locations.

calculate_co_occurence <- function() {
  proteins_with_subcellular_location <- unique(data[c("Protein", "Subcellular_location")])
  
  proteins_with_subcellular_location_no_empty <- proteins_with_subcellular_location[!(proteins_with_subcellular_location$Subcellular_location == ""),]
  
  protein_in_subcellular_location = matrix(ncol = 9, nrow = nrow(proteins_with_subcellular_location_no_empty))
  
  for (j in 1:length(search_terms)) {
    protein_in_subcellular_location[,j] =  grepl(search_terms[[j]], proteins_with_subcellular_location_no_empty$Subcellular_location, fixed=TRUE)
  }
  
  protein_in_subcellular_location <- protein_in_subcellular_location*1
  
  protein_in_subcellular_location_df <- as.data.frame(protein_in_subcellular_location)
  colnames(protein_in_subcellular_location_df) <- c("Cell membrane", "Mitochondrion", "Nucleus", "Endoplasmic Reticulum", "Golgi apparatus", "Lysosome", "Cytoplasm", "Secreted", "Extracellular space")
  rownames(protein_in_subcellular_location_df) <- proteins_with_subcellular_location_no_empty$Protein
  
  protein_in_subcellular_location_df
  
  total_occurences <- colSums(protein_in_subcellular_location_df)
  co_occurence <<- t(protein_in_subcellular_location) %*% protein_in_subcellular_location
}

library(igraph)

plot_co_occurence_absolute <- function() {
  if (! exists("co_occurence")) {
    calculate_co_occurence()
  }
  
  graph <- graph.adjacency(co_occurence,
                           weighted=TRUE,
                           mode="undirected",
                           diag=FALSE) 
  
  co_occurence_plot <- plot(graph,
                            vertex.label=names(protein_in_subcellular_location_df),
                            vertex.size=total_occurences/100,
                            edge.width=E(graph)$weight/100,
                            layout = layout_in_circle(graph, order = V(graph)),
                            margin = -0.4)
}

