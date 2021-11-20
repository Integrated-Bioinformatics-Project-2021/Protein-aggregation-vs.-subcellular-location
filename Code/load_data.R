# Load the human_proteome data set, add the subcellular_location information
# retrieved from UniProt, add the secretory information and save it as
# complete_data.RData

source("subcellular_location_annotation.R")

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
    
    # Retrieve the wanted subcellular location information, check if this is secretory
    # or non-secretory and add it to the data
    return_value = add_subcellular_location_and_secretory_information(data, proteins)
    data = return_value$data
    hashed_proteins = return_value$hashed_proteins
    
    save(data, file = "Data/complete_data.RData") # Uncomment to save the results to a file
    save(hashed_proteins, file = "Data/hashed_proteins.RData")
  } else {
    load("Data/complete_data.RData")
    load("Data/hashed_proteins.RData")
  }
  return (list("data" = data, "hashed_proteins" = hashed_proteins))
}
