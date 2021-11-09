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

## Reading the data
load("Data/human_proteome_df.RData")
data = human_proteome

#Delete all proteins witch were found in transmembrane regions or that are signal peptides (THMM, THMM_domain, SingalIP)
data <- subset(data, TMHMM == "No"& TMHMM_domain == "No"& SignalP == "No")

## Get all unique protein identifiers
proteins <- unique(data$Protein)
proteins



## Get subcellular location of every protein
locations = GetSubcellular_location(proteins[4999:5200], directorypath = NULL)
subcellular_locations = locations[1] # Remove NA's from intramembrane, topological domain and transmembrane information
search_terms <- list("Cell membrane", "Mitochondrion", "Nucleus", "Endoplasmic Reticulum", "Golgi apparatus", "Lysosome", "Cytoplasm", "Secreted", "Extracellular space")

secretory_pathway = matrix(ncol = 1, nrow = 0)

locations_short = data.frame(subcellular_location, secretory_pathway) # creating a dataframe to add the subcellular locations
hashed_proteins = hash() #generate dictionary

secretory <- c("Cell membrane", "Endoplasmatic Reticulum", "Endoplasmic Reticulum", "Golgi apparatus", "Lysosome", "Secreted region", "Secreted", "Extracellular space")
non_secretory <- c("Cytoplasm", "Nucleus", "Mitochondrion")

check_secretory <- function(x) {
  
  locations_str <- as.character(x)
  locations_vector <- strsplit(locations_str, split = ", ")
  result = ""
  
  if(length(locations_vector) != 0)
  {
    for (j in 1:length(locations_vector))
    {
      if (locations_vector[j] %in% secretory)
      {
        result = paste(result, "sec", sep = "")
      }
      else if(locations_vector[j] %in% non_secretory)
      {
        result = paste(result, "nonSec", sep = "")
      }
    }
  }
  
  if (grepl("sec", result, fixed = TRUE) & grepl("nonSec", result, fixed = TRUE)){
    return(NA)
  }
  else if (grepl("sec", result, fixed = TRUE)){
    return(TRUE)
  }
  else if (grepl("nonSec", result, fixed = TRUE)){
    return(FALSE)
  }
  else {
    return(NA)
  }
}

for (i in 1:dim(subcellular_locations)[1]) { #Deleting everything between {} and the "SUBCELLULAR LOCATION"
  string_location <- str_contains(subcellular_locations[i,], search_terms, ignore.case = TRUE)
  x <- search_terms[string_location]
  if (length(x) != 0) {
    if (length(x) > 1) {
      x <- paste(x, collapse = ", ")
    }
    locations_short[nrow(locations_short) + 1,1] = x # adding each value to the dataframe
    locations_short[nrow(locations_short),2] = check_secretory(x) # stating whether the proteins is found in an intracellular or extracellular location
  }
  else {
    locations_short[nrow(locations_short) + 1,] = ""
  }
  current_protein_name = proteins[4999+i]
  hashed_proteins[current_protein_name] <- locations_short[i,]
  current_protein_name = toString(current_protein_name)
  data$Subcellular_location[data$Protein == current_protein_name]  <- hashed_proteins[[current_protein_name]][["subcellular_location"]]
  data$secretory_pathway[data$Protein == current_protein_name]  <- hashed_proteins[[current_protein_name]][["secretory_pathway"]]
}

# save(data, file = "Data/complete_data.RData") # Uncomment to save the results to a file

# Return a list with all the protein IDs of that contain the given subcellular location.
get_proteins_with_given_subcellular_location <- function(given_subcellular_location){
  wanted_proteins <- sapply(keys(hashed_proteins), function(x) grepl(given_subcellular_location, hashed_proteins[[x]][["subcellular_location"]]))
  wanted_proteins <- keys(hashed_proteins)[wanted_proteins]
  return (data.frame(wanted_proteins))
}


# ------------------ tango scores ---------------------------

# Acquire only the vital rows of each protein (APRs)
data <- subset(data, APRdef2_tango > 0)

# Get average tango score for each protein
unique_data <- data[!duplicated(data[,c(1,10)]),]
by_protein <- aggregate(avgScore ~ Protein, unique_data, mean)

#get average Tango score for each subcellular locations 
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
tango_scores = data.frame()
for (i in 1:length(search_terms)) {
  given_protein_list = get_proteins_with_given_subcellular_location(search_terms[i])
  tango_scores[nrow(tango_scores)+1:nrow(given_protein_list),1] = given_protein_list
  newIndex = nrow(tango_scores)+1-nrow(given_protein_list)
  tango_scores[newIndex:(newIndex+nrow(given_protein_list)-1),2] = rep(search_terms[[i]], nrow(given_protein_list))
  #if (! is.na(given_protein_list$wanted_proteins)) {
  avgScores = get_average_tango_score(given_protein_list)
  tango_scores[newIndex:(newIndex+nrow(given_protein_list)-1),3] = avgScores
  tango_scores[newIndex:(newIndex+nrow(given_protein_list)-1),4] = rep(check_secretory(search_terms[i]), nrow(given_protein_list))
    #}
}
colnames(tango_scores) <- c("Proteins", "Subcellular_location", "Tango_scores", "Secretory")

#barplot(tango_scores, xlab = "subcellular location", names.arg = search_terms) # TODO: make it pretty ^_^ library(ggplot2)

library(ggplot2)
box_tango_sub <- ggplot(tango_scores, aes(x=Tango_scores, y = Subcellular_location, fill = Secretory)) + 
  geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
box_tango_sub + theme_minimal() + stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red")

# TODO: join secreted and extracellular

# ------------------ peptides ---------------------------

# Classify into 3 subset dataframes (peptides with only APR, peptides with APR and GK and peptides with APR, GK and 2 FR )
APR_peptides <- subset(data, APRdef2_tango < 2)
GK_peptides <- subset(data, APRdef2_tango < 3)
FR_peptides<-  subset(data, APRdef2_tango < 5)

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
