#install.packages("UniprotR")
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("GenomicAlignments")
#install.packages("sjmisc")
#install.packages("hash")
# install.packages("dplyr")
# install.packages("Peptides")
library(UniprotR)
library(sjmisc) # Used for str_contains
library(hash) # Used to make a disctionary
library(dplyr) # Used for aggregate (get avg tango score of protein)
library(Peptides) # Used for charge (get the charge of a peptide sequence) ETC

## Defining the working directory
directory = dirname(rstudioapi::getSourceEditorContext()$path) # Should work when data is placed in same folder
setwd(directory)
source("load_data.R")
source("subcellular_location_annotation.R")

return_value = read_data() # See load_data.R
data = return_value[1]
hashed_proteins = return_value[2]

# ------------------ tango scores ---------------------------

# Acquire only the vital rows of each protein (APRs)
only_APR_data <- subset(data, APRdef2_tango > 0)


# Get average tango score for each complete protein 
unique_data <- data[!((data$APRcount_tango == 0) & (data$maxProtscore != 0)),]
unique_data <- unique_data[!duplicated(unique_data[,c(1,10)]),]
by_protein <- aggregate(avgScore ~ Protein, unique_data, mean)

# Get average tango score for each APR protein
unique_APR_data <- only_APR_data[!duplicated(data[,c(1,10)]),]
by_APR_protein <- aggregate(avgScore ~ Protein, unique_APR_data, mean)

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
tango_scores_complete_protein = data.frame()
tango_scores_APR_protein = data.frame()
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

# normalized_number_APR_regions <- get_normalized_number_APR_regions(data)

#barplot(tango_scores, xlab = "subcellular location", names.arg = search_terms)

library(ggplot2)
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


# ------------------ peptides ---------------------------

# Classify into 3 subset dataframes (peptides with only APR, peptides with APR and GK and peptides with APR, GK and 2 FR )
APR_peptides <- subset(only_APR_data, APRdef2_tango < 2)
GK_peptides <- subset(only_APR_data, APRdef2_tango < 3)
FR_peptides<-  subset(only_APR_data, APRdef2_tango < 5)

# For every APR in the list peptides_name, calculate the net charge
get_charge <- function(peptides_name, ph = 7) {
  protein_charges = hash()
  unique_APRs = unique(peptides_name$APRcount_tango)
  for (i in 1:length(unique_APRs)) {
    peptides_for_APR = subset(peptides_name, peptides_name$APRcount_tango == unique_APRs[i])
    sequence = paste(peptides_for_APR$Residue, collapse = "")
    protein_charges[unique_APRs[i]] <- charge(sequence, pH = ph)
    key = toString(unique_APRs[i])
    print(protein_charges[[key]])
  }
  return (protein_charges)
}

protein_sequences_APR_peptides = get_charge(APR_peptides)


# protein_sequences_GK_peptides = get_charge(GK_peptides)
# protein_sequences_FR_peptides = get_charge(FR_peptides)

# Make a histogram

