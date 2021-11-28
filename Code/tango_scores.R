# All functions dealing with the TANGO scores

source("subcellular_location_annotation.R")
source("statistical_analysis.R")

# Statistical significance ---------------------------------------------------

test_significance_secretory <- function(data){
  if (significance_secretory(data) == TRUE){
    print("The observed difference is statistically significant.")
  }
  else{print("The observed difference is not statistically significant.")}
}

# Average TANGO scores --------------------------------------------------------

initialize_by_protein <- function() {
  # Acquire only the vital rows of each protein (APRs)
  only_APR_data <<- subset(data, APRdef2_tango > 0)
  
 
  # Get average tango score for each APR protein
  unique_data <- only_APR_data[!duplicated(only_APR_data[,c(1,10)]),]
  by_protein <<- aggregate(avgScore ~ Protein, unique_data, mean) # Save by_protein globally
}
initialize_by_protein()

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
get_average_tango_score <- function(given_protein_list) {
  avgScores = c()
  for (j in 1:nrow(given_protein_list)){
    if (given_protein_list$wanted_proteins[j] %in% by_protein$Protein) {
      avgScores <- append(avgScores,by_protein[which(by_protein$Protein %in% given_protein_list$wanted_proteins[j]), "avgScore"])
    }
    else {avgScores <- append(avgScores,NA)}
  }
  return (avgScores)
}

# Get the average tango score for each subcellular location
calculate_average_tango_scores <- function() {
    tango_scores = data.frame()
    for (i in 1:length(search_terms)) {
        given_protein_list = get_proteins_with_given_subcellular_location(search_terms[i])
        tango_scores[nrow(tango_scores)+1:nrow(given_protein_list),1] = given_protein_list
        newIndex_AP = nrow(tango_scores)+1-nrow(given_protein_list)
  
        tango_scores[newIndex_AP:(newIndex_AP+nrow(given_protein_list)-1),2] = rep(search_terms[[i]], nrow(given_protein_list))
        #if (! is.na(given_protein_list$wanted_proteins)) 
        avgScores = get_average_tango_score(given_protein_list)
        tango_scores[newIndex_AP:(newIndex_AP+nrow(given_protein_list)-1),3] = avgScores
        tango_scores[newIndex_AP:(newIndex_AP+nrow(given_protein_list)-1),4] = rep(check_secretory(search_terms[i]), nrow(given_protein_list))
        #}
        }
    colnames(tango_scores)  <- c("Proteins", "Subcellular_location", "Tango_scores", "Secretory")
    tango_scores_complete_protein <- tango_scores
    tango_scores_complete_protein[is.na(tango_scores_complete_protein)] <- 0
    tango_scores_APR_protein <- na.omit(tango_scores)
    tango_scores_complete_protein <<- tango_scores_complete_protein # Save object globally
    tango_scores_APR_protein <<- tango_scores_APR_protein # Save object globally
}

plot_average_tango_scores_complete_proteins <- function() {
  if (! exists("tango_scores_complete_protein")) {
    calculate_average_tango_scores()
  }
  #Plot complete proteins
  violin_tango_sub_complete_proteins <- ggplot(tango_scores_complete_protein, aes(x=Tango_scores, y = Subcellular_location, fill = Secretory)) + 
    geom_violin(notch=TRUE) + scale_color_brewer(palette="Dark2")
  violin_tango_sub_complete_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}

plot_average_tango_scores_complete_proteins_joined_secreted <- function() {
  if (! exists("tango_scores_APR_protein")) {
    calculate_average_tango_scores()
  }
  #Plot for difference between secretory and non-secretory 
  #Plot complete proteins
 violin_tango_sec_complete_proteins <- ggplot(tango_scores_complete_protein, aes(x=Tango_scores, y = Secretory, fill = Secretory)) + 
  geom_violin(trim=FALSE) + scale_color_brewer(palette="Dark2")
violin_tango_sec_complete_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}

plot_average_tango_scores_APR_proteins <- function() {
  if (! exists("tango_scores_APR_protein")) {
    calculate_average_tango_scores()
  }
  #Plot APR proteins
violin_tango_sub_APR_proteins <- ggplot(tango_scores_APR_protein, aes(x=Tango_scores, y = Subcellular_location, fill = Secretory)) + 
  geom_violin(trim=TRUE) + scale_color_brewer(palette="Dark2")
violin_tango_sub_APR_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")

}

plot_average_tango_scores_APR_proteins_joined_secreted <- function() {
  if (! exists("tango_scores_APR_protein")) {
    calculate_average_tango_scores()
  }
  #Plot for difference between secretory and non-secretory 
  #Plot APR proteins
 violin_tango_sec_APR_proteins <- ggplot(tango_scores_APR_protein, aes(x=Tango_scores, y = Secretory, fill = Secretory)) + 
  geom_violin(trim=FALSE) + scale_color_brewer(palette="Dark2")
 violin_tango_sec_APR_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}

# ----------- Max Tango Scores for each protein and APR -----------------------

calculate_max_tango_scores <- function() {
  by_protein_max <<- aggregate(maxProtscore~Protein, data, max)
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
  max_tango_scores <<- scores_protein
}

plot_max_tango_scores <- function() {
  if (! exists("max_tango_scores")) {
    calculate_max_tango_scores()
  }
 violin_avgmax_complete_proteins <- ggplot(max_tango_scores, aes(x=Tango_scores, y = Subcellular_location, fill = Secretory)) + 
  geom_violin(trim=FALSE) + scale_color_brewer(palette="Dark2")
 violin_avgmax_complete_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}

#TODO: violin plot max tango score per secretory pathway missing!

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
calculated_normalized_number_APR_regions <- function() {
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
  normalized_number_APR_regions <<- normalized_number_APR_regions
}

plot_normalized_number_APR_regions <- function() {
  if (! exists("normalized_number_APR_regions")) {
    calculated_normalized_number_APR_regions()
  }
  #Plot normalized number of APR regions
  box_normalized_number_APR_regions <- ggplot(normalized_number_APR_regions, aes(x=Nb_APRs, y = Subcellular_location, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_normalized_number_APR_regions + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}

plot_normalized_number_APR_regions_joined_secreted <- function() {
  if (! exists("normalized_number_APR_regions")) {
    calculated_normalized_number_APR_regions()
  }
  #Plot for difference between secretory and non-secretory
  box_normalized_number_APR_regions_secretory <- ggplot(normalized_number_APR_regions, aes(x=Nb_APRs, y = Secretory, fill = Secretory)) + 
    geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
  box_normalized_number_APR_regions_secretory + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")
}