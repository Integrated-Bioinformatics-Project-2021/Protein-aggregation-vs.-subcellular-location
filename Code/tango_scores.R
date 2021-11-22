# All functions dealing with the TANGO scores

source("subcellular_location_annotation.R")
# Average TANGO scores --------------------------------------------------------

initialize_by_APR_protein <- function() {
  # Acquire only the vital rows of each protein (APRs)
  only_APR_data <- subset(data, APRdef2_tango > 0)
  
  
  # Get average tango score for each complete protein 
  unique_data <- data[!((data$APRcount_tango == 0) & (data$maxProtscore != 0)),]
  unique_data <- unique_data[!duplicated(unique_data[,c(1,10)]),]
  by_protein <<- aggregate(avgScore ~ Protein, unique_data, mean) # Save by_protein globally
  
  # Get average tango score for each APR protein
  unique_APR_data <- only_APR_data[!duplicated(only_APR_data[,c(1,10)]),]
  by_APR_protein <<- aggregate(avgScore ~ Protein, unique_APR_data, mean) # Save by_APR_protein globally
}
initialize_by_APR_protein()

#get average Tango score for each subcellular locations for complete proteins 
get_average_tango_score_complete_protein <- function(given_protein_list) {
  avgScores = c()
  for (j in 1:nrow(given_protein_list)){
    if (given_protein_list$wanted_proteins[j] %in% by_protein$Protein) {
      avgScores <- append(avgScores,by_protein[which(by_protein$Protein %in% given_protein_list$wanted_proteins[j]), "avgScore"])
    }
    else {avgScores <- append(avgScores,NA)}
  }
  return (avgScores)
}

#get average max Tango score for each subcellular locations for complete proteins 
get_avgmax_tango_score_complete_protein <- function(given_protein_list) {
  avgScores = c()
  for (j in 1:nrow(given_protein_list)){
    if (given_protein_list$wanted_proteins[j] %in% by_protein_max$Protein) {
      avgScores <- append(avgScores,by_protein_max[which(by_protein_max$Protein %in% given_protein_list$wanted_proteins[j]), "maxProtscore"])
    }
    else {avgScores <- append(avgScores,NA)}
  }
  return (avgScores)
}

#get average Tango score for each subcellular locations for APR proteins 
get_average_tango_score_APR_protein <- function(given_protein_list) {
  avgScores = c()
  for (j in 1:nrow(given_protein_list)){
    if (given_protein_list$wanted_proteins[j] %in% by_APR_protein$Protein) {
      avgScores <- append(avgScores,by_APR_protein[which(by_APR_protein$Protein %in% given_protein_list$wanted_proteins[j]), "avgScore"])
    }
    else {avgScores <- append(avgScores,NA)}
  }
  return (avgScores)
}

# Get the average tango score for each subcellular location
calculate_average_tango_scores <- function() {
  tango_scores_complete_protein = data.frame() # For all proteins
  tango_scores_APR_protein = data.frame() # For proteins with an APR region
  for (i in 1:length(search_terms)) {
    given_protein_list = get_proteins_with_given_subcellular_location(search_terms[i])
    tango_scores_complete_protein[nrow(tango_scores_complete_protein)+1:nrow(given_protein_list),1] = given_protein_list
    tango_scores_APR_protein[nrow(tango_scores_APR_protein)+1:nrow(given_protein_list),1] = given_protein_list
    
    newIndex_CP = nrow(tango_scores_complete_protein)+1-nrow(given_protein_list)
    newIndex_AP = nrow(tango_scores_APR_protein)+1-nrow(given_protein_list)
    
    tango_scores_complete_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),2] = rep(search_terms[[i]], nrow(given_protein_list))
    tango_scores_APR_protein[newIndex_AP:(newIndex_AP+nrow(given_protein_list)-1),2] = rep(search_terms[[i]], nrow(given_protein_list))
    #if (! is.na(given_protein_list$wanted_proteins)) {
    avgScores_complete_protein = get_average_tango_score_complete_protein(given_protein_list)
    avgScores_APR_protein = get_average_tango_score_APR_protein(given_protein_list)
    tango_scores_complete_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),3] = avgScores_complete_protein
    tango_scores_APR_protein[newIndex_AP:(newIndex_AP+nrow(given_protein_list)-1),3] = avgScores_APR_protein
    tango_scores_complete_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),4] = rep(check_secretory(search_terms[i]), nrow(given_protein_list))
    tango_scores_APR_protein[newIndex_AP:(newIndex_AP+nrow(given_protein_list)-1),4] = rep(check_secretory(search_terms[i]), nrow(given_protein_list))
    #}
  }
  colnames(tango_scores_complete_protein) <- c("Proteins", "Subcellular_location", "Tango_scores", "Secretory")
  colnames(tango_scores_APR_protein)  <- c("Proteins", "Subcellular_location", "Tango_scores", "Secretory")
  tango_scores_complete_protein <<- tango_scores_complete_protein # Save object globally
  tango_scores_APR_protein <<- tango_scores_APR_protein # Save object globally
}

plots_average_tango_scores <- function() {
  if (! exists("tango_scores_APR_protein")) {
    calculate_average_tango_scores()
  }
  
  #barplot(tango_scores, xlab = "subcellular location", names.arg = search_terms)
  
  #Plot complete proteins
  box_tango_sub_complete_proteins <- ggplot(tango_scores_complete_protein, aes(x=Tango_scores, y = Subcellular_location, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_tango_sub_complete_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
  
  #Plot APR proteins
  box_tango_sub_APR_proteins <- ggplot(tango_scores_APR_protein, aes(x=Tango_scores, y = Subcellular_location, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_tango_sub_APR_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")

  # TODO: join secreted and extracellular
  
  #Plot for difference between secretory and non-secretory 
  #Plot complete proteins
  box_tango_sec_complete_proteins <- ggplot(tango_scores_complete_protein, aes(x=Tango_scores, y = Secretory, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_tango_sec_complete_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
  #Plot APR proteins
  box_tango_sec_APR_proteins <- ggplot(tango_scores_APR_protein, aes(x=Tango_scores, y = Secretory, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_tango_sec_complete_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}

# ----------- Max Tango Scores for each protein and APR -----------------------

plots_max_tango_scores <- function() {
  by_protein_max = aggregate(maxProtscore~Protein, data, max)
  scores_protein = data.frame()
  
  for(i in 1:length(search_terms))
  {
    given_protein_list = get_proteins_with_given_subcellular_location(search_terms[i])
    scores_protein[nrow(scores_protein)+1:nrow(given_protein_list),1] = given_protein_list
    
    newIndex_CP = nrow(scores_protein)+1-nrow(given_protein_list)
    scores_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),2] = rep(search_terms[[i]], nrow(given_protein_list))
    avgMaxScores_protein = get_avgmax_tango_score_complete_protein(given_protein_list)
    
    scores_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),3] = avgMaxScores_protein
    scores_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),4] = rep(check_secretory(search_terms[i]), nrow(given_protein_list))
  }
  colnames(scores_protein) <- c("Proteins", "Subcellular_location", "Tango_scores", "Secretory")
  
  box_avgmax_complete_proteins <- ggplot(scores_protein, aes(x=Tango_scores, y = Subcellular_location, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_avgmax_complete_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}

# Number of APRs per subcellular location -------------------------------------

#get normalized number of APR regions for a protein list
get_normalized_number_APR_regions <- function(given_protein_list) {
  unique_proteins_in_list = unique(given_protein_list$Protein)
  avgNumbers = c()
  for (j in 1:length(unique_proteins_in_list)){ # loop over all unique proteins in the given list
    peptides_for_protein = subset(given_protein_list, given_protein_list$Protein == unique_proteins_in_list[j])
    total_length = nrow(peptides_for_protein) # count every row for the protein representing each residue 
    APR_list <- unique(peptides_for_protein$APRcount_tango) # add all unique APR ideas to the list
    # get the number of APR regions in the current protein
    number_APR_regions = length(APR_list)
    # divide by the total length of the protein
    avgNumber <- number_APR_regions/total_length
    # print(avgNumber)
    # append to the data frame
    avgNumbers <- append(avgNumbers, avgNumber)
  }
  return (avgNumbers)
}

# get normalized number of APR regions for each subcellular location
plots_normalized_number_APR_regions <- function() {
  normalized_number_APR_regions = data.frame()
  for (i in 1:length(search_terms)) {
    given_protein_list = get_proteins_with_given_subcellular_location(search_terms[i])
    given_protein_list_extended = subset(data, data$Protein %in% given_protein_list[[1]])
    
    normalized_number_APR_regions[nrow(normalized_number_APR_regions)+1:nrow(given_protein_list),1] = given_protein_list
    
    newIndex_CP = nrow(normalized_number_APR_regions)+1-nrow(given_protein_list)
    
    normalized_number_APR_regions[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),2] = rep(search_terms[[i]], nrow(given_protein_list))
    
    avg_normalized_number_APR_regions = get_normalized_number_APR_regions(given_protein_list_extended)
    normalized_number_APR_regions[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),3] = avg_normalized_number_APR_regions
    normalized_number_APR_regions[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),4] = rep(check_secretory(search_terms[i]), nrow(given_protein_list))
    
  }
  colnames(normalized_number_APR_regions) <- c("Proteins", "Subcellular_location", "Nb_APRs", "Secretory")
  normalized_number_APR_regions_copy = normalized_number_APR_regions
  normalized_number_APR_regions_copy$Subcellular_location[normalized_number_APR_regions$Subcellular_location == "Secreted" & ! normalized_number_APR_regions$Proteins %in% normalized_number_APR_regions$Proteins[normalized_number_APR_regions$Subcellular_location == "Extracellular space"]] = "Extracellular space" 
  
  #Plot normalized number of APR regions
  box_normalized_number_APR_regions <- ggplot(normalized_number_APR_regions_copy, aes(x=Nb_APRs, y = Subcellular_location, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_normalized_number_APR_regions + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
  
  #Plot for difference between secretory and non-secretory
  box_normalized_number_APR_regions_secretory <- ggplot(normalized_number_APR_regions_copy, aes(x=Nb_APRs, y = Secretory, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_normalized_number_APR_regions_secretory + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}