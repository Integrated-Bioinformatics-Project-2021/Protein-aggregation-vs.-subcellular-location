#install.packages("UniprotR")
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("GenomicAlignments")
#install.packages("sjmisc")
#install.packages("hash")
#install.packages("dplyr")
library(UniprotR)
library(stringr) # Used to get the last word of a string
library(sjmisc) # Used for str_contains
library(hash)
library(dplyr)

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
locations = GetSubcellular_location(proteins[5000:5005], directorypath = NULL)
subcellular_locations = locations[1] # Remove NA's from intramembrane, topological domain and transmembrane information
search_terms <- list("Cell membrane", "Mitochondrion", "Nucleus", "Endoplasmatic Reticulum", "Golgi apparatus", "Lysosomes", "Cytoplasm", "Secreted region", "Extracellular region")

subcellular_location = matrix(ncol = 1, nrow = 0) 
secretory_pathway = matrix(ncol = 1, nrow = 0)

locations_short = data.frame(subcellular_location, secretory_pathway) # creating a dataframe to add the subcellular locations
hashed_proteins = hash() #generate dictionary

secretory <- c("Cell membrane", "Endoplasmatic Reticulum", "Golgi apparatus", "Lysosomes", "Secreted region")
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
  string_location <- str_contains(subcellular_locations[i,], search_terms)
  x <- search_terms[string_location]
  if (length(x) != 0) {
    print(length(x))
    if (length(x) > 1) {
      x <- paste(x, collapse = ", ")
      print(x)
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

proteins_cell_membrane = get_proteins_with_given_subcellular_location("Cell membrane")
proteins_nucleus = get_proteins_with_given_subcellular_location("Nucleus")


# Acquire only the vital rows of each protein (APRs)
data <- subset(data, APRdef2_tango > 0)

# Get average tango score for each protein
unique_data <- data[!duplicated(data[,c(1,10)]),]

by_protein <- aggregate(avgScore ~ Protein, unique_data, mean)

#get average Tango score for each subcellular locations 

get_average_tango_score <- function(given_protein_list) {
  avgTotalScore = 0
  for (j in 1:dim(given_protein_list)[1]){
    if (given_protein_list$wanted_proteins[j] %in% by_protein$Protein) {
      avgTotalScore <- avgTotalScore + by_protein[which(by_protein$Protein %in% given_protein_list$wanted_proteins[j]), "avgScore"]
    }
  }
  return (avgTotalScore)
}




# Classify into 3 subset dataframes (peptides with only APR, peptides with APR and GK and peptides with APR, GK and 2 FR )
APR_peptides <- subset(data, APRdef2_tango < 2)
GK_peptides <- subset(data, APRdef2_tango < 3)
FR_peptides<-  subset(data, APRdef2_tango < 5)



