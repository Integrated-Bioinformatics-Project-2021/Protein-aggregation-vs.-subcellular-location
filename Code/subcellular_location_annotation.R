# All functions dealing with the annotation of subcellular locations

## Get subcellular location of every protein ----------------------------------

# Define all possible subcellular locations
search_terms <- list("Cell membrane", "Mitochondrion", "Nucleus", "Endoplasmic Reticulum", "Golgi apparatus", "Lysosome", "Cytoplasm", "Secreted", "Extracellular space")

# Get the Uniprot information of the subcellular locations of all given proteins
get_subcellular_locations_for_all_proteins <- function(proteins) {
  locations = GetSubcellular_location(proteins, directorypath = NULL)
  subcellular_locations = locations[1] # Remove NA's from intramembrane, topological domain and transmembrane information
  return (subcellular_locations)
}

# Parse the wanted subcellular location information, check if this is secretory
# or non-secretory and add it to the data
add_subcellular_location_and_secretory_information <- function(data, proteins) {
  subcellular_locations = get_subcellular_locations_for_all_proteins(proteins)
  
  subcellular_location = matrix(ncol = 1, nrow = 0)
  secretory_pathway = matrix(ncol = 1, nrow = 0)
  
  locations_short = data.frame(subcellular_location, secretory_pathway) # creating a dataframe to add the subcellular locations
  hashed_proteins = hash() #generate dictionary
  
  for (i in 1:dim(subcellular_locations)[1]) {
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
    current_protein_name = proteins[i]
    hashed_proteins[current_protein_name] <- locations_short[i,]
    current_protein_name = toString(current_protein_name)
    data$Subcellular_location[data$Protein == current_protein_name]  <- hashed_proteins[[current_protein_name]][["subcellular_location"]]
    data$secretory_pathway[data$Protein == current_protein_name]  <- hashed_proteins[[current_protein_name]][["secretory_pathway"]]
  }
  return (c(data, hashed_proteins))
}

# Get hashed_proteins when complete_data.R is loaded ---------------------------
get_hashed_proteins <- function(complete_data) {
  proteins = unique(complete_data$Protein)
  hashed_proteins = hash() #generate dictionary
  subcellular_location = matrix(ncol = 1, nrow = 0)
  secretory_pathway = matrix(ncol = 1, nrow = 0)
  
  locations_short = data.frame(subcellular_location, secretory_pathway) # creating a dataframe to add the subcellular locations
  
  for (i in 1:length(proteins)) {
    current_protein_name = proteins[i]
    locations_short[1,1] = complete_data$Subcellular_location[complete_data$Protein == current_protein_name][1]
    locations_short[1,2] = complete_data$Subcellular_location[complete_data$Protein == secretory_pathway][1]
    hashed_proteins[current_protein_name] <- locations_short[1,]
  }
  return (hashed_proteins)
}

# Check secretory --------------------------------------------------------------

# Define which subcellular locations are secretory and which are non-secretory
secretory <- c("Cell membrane", "Endoplasmic Reticulum", "Golgi apparatus", "Lysosome", "Secreted", "Extracellular space")
non_secretory <- c("Cytoplasm", "Nucleus", "Mitochondrion")

# Check if a given subcellular location group is secretory
# Returns TRUE if the group only contains secretory subcellular locations
# Returns FALSE if the group only contains non-secretory subcellular locations
# Returns NA if the group doesn't contain subcellular locations or if the group
# contains both secretory and non-secretory ones.
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

# Proteins for given subcellular location --------------------------------------

# Return a list with all the protein IDs of that contain the given subcellular location.
get_proteins_with_given_subcellular_location <- function(given_subcellular_location){
  wanted_proteins <- sapply(keys(hashed_proteins), function(x) grepl(given_subcellular_location, hashed_proteins[[x]][["subcellular_location"]]))
  wanted_proteins <- keys(hashed_proteins)[wanted_proteins]
  return (data.frame(wanted_proteins))
}