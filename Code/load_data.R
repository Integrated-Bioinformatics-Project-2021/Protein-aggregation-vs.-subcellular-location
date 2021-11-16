# Load the human_proteome data set, add the subcellular_location information
# retrieved from UniProt, add the secretory information and save it as
# complete_data.RData

## Reading the data
read_data <- function() {
  if (! file.exists("Data/complete_data.RData")) { # Only retrieve the data when it isn't stored yet
    load("Data/human_proteome_df.RData")
    data = human_proteome
    
    #Delete all proteins witch were found in transmembrane regions or that are signal peptides (THMM, THMM_domain, SingalIP)
    data <- subset(data, TMHMM == "No"& TMHMM_domain == "No"& SignalP == "No")
    
    ## Get all unique protein identifiers
    proteins <- unique(data$Protein)
    proteins
    
    
    
    ## Get subcellular location of every protein
    locations = GetSubcellular_location(proteins, directorypath = NULL)
    subcellular_locations = locations[1] # Remove NA's from intramembrane, topological domain and transmembrane information
    search_terms <- list("Cell membrane", "Mitochondrion", "Nucleus", "Endoplasmic Reticulum", "Golgi apparatus", "Lysosome", "Cytoplasm", "Secreted", "Extracellular space")
    
    subcellular_location = matrix(ncol = 1, nrow = 0)
    secretory_pathway = matrix(ncol = 1, nrow = 0)
    
    locations_short = data.frame(subcellular_location, secretory_pathway) # creating a dataframe to add the subcellular locations
    hashed_proteins = hash() #generate dictionary
    
    secretory <- c("Cell membrane", "Endoplasmic Reticulum", "Golgi apparatus", "Lysosome", "Secreted", "Extracellular space")
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
      current_protein_name = proteins[i]
      hashed_proteins[current_protein_name] <- locations_short[i,]
      current_protein_name = toString(current_protein_name)
      data$Subcellular_location[data$Protein == current_protein_name]  <- hashed_proteins[[current_protein_name]][["subcellular_location"]]
      data$secretory_pathway[data$Protein == current_protein_name]  <- hashed_proteins[[current_protein_name]][["secretory_pathway"]]
    }
    
    save(data, file = "Data/complete_data.RData") # Uncomment to save the results to a file
  } else {
    load("Data/complete_data.RData")
  }
  return (data)
}
