# install.packages("UniprotR")
# install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("GenomicAlignments")
# install.packages("sjmisc")
# install.packages("hash")
library(UniprotR)
library(stringr) # Used to get the last word of a string
library(sjmisc) # Used for str_contains
library(hash)

## Defining the working directory
directory = dirname(rstudioapi::getSourceEditorContext()$path) # Should work when data is placed in same folder
setwd(directory)

## Reading the data
load("Data/human_proteome_df.RData")
data = human_proteome

## Get all unique protein identifiers
proteins <- unique(data$Protein)
proteins

## Get subcellular location of every protein
locations = GetSubcellular_location(proteins[5000:5005], directorypath = NULL)
subcellular_locations = locations[1] # Remove NA's from intramembrane, topological domain and transmembrane information
search_terms <- list("Cell membrane", "Mitochondrion", "Nucleus", "Endoplasmatic Reticulum", "Golgi apparatus", "Lysosomes", "Cytoplasm", "Secreted region", "Extracellular region")
subcellular_location = matrix(ncol = 1, nrow = 0) 
locations_short = data.frame(subcellular_location) # creating a dataframe to add the subcellular locations
hashed_proteins = hash() #generate dictionary
for (i in 1:dim(subcellular_locations)[1]) { #Deleting everything between {} and the "SUBCELLULAR LOCATION"
  #print(subcellular_locations[i,])
  string_location <- str_contains(subcellular_locations[i,], search_terms)
  print(string_location)
  x <- search_terms[string_location]
  #print(x)
  if (length(x) != 0) {
    print(length(x))
    if (length(x) > 1) {
      x <- paste(x, collapse = ", ")
      print(x)
    }
    locations_short[nrow(locations_short) + 1,] = x #adding each value to the dataframe
    print(locations_short)
  }
  else {
    locations_short[nrow(locations_short) + 1,] = ""
  }
  current_protein_name = proteins[4999+i]
  hashed_proteins[current_protein_name] <- locations_short[i,]
}

#Adding the subcellular location to the dataset
locations_short<- cbind(locations_short, Protein = row.names(subcellular_locations)) #creating a second column Protein
complete_data = merge(locations_short, human_proteome, by="Protein", all.y = TRUE) #merging the dataframes on column Protein, keeping non matches
save(complete_data, file = "Data/complete_data.RData")

# To check if the merging is working
# complete_data[complete_data$Protein=="P40259",c("Protein","subcellular_location")]

secretory <- c("Cell membrane", "Endoplasmatic Reticulum", "Golgi apparatus", "Lysosomes", "Secreted region")
non_secretory <- c("Cytoplasm", "Nucleus", "Mitochondrion")

for (i in 1:nrow(locations_short)){

  locations_str <- sapply(locations_short[i,]["subcellular_location"], as.character)
  locations_vector <- strsplit(locations_str, split = ", ")
  
  result = ""
  
  if(length(locations_vector$subcellular_location) != 0)
  {
    for (j in 1:length(locations_vector$subcellular_location))
    {
      if (locations_vector$subcellular_location[j] %in% secretory)
      {
        result = paste(result, "a", sep = "")
      }
      else if(locations_vector$subcellular_location[j] %in% non_secretory)
      {
        result = paste(result, "b", sep = "")
      }
    }
  }
  
  if (grepl("a", result, fixed = TRUE) & grepl("b", result, fixed = TRUE)){
    locations_short[i,]["Secretory"] = "3"
  }
  else if (grepl("a", result, fixed = TRUE)){
    locations_short[i,]["Secretory"] = "1"
  }
  else if (grepl("b", result, fixed = TRUE)){
    locations_short[i,]["Secretory"] = "2"
  }
  else{
    locations_short[i,]["Secretory"] = "3"
  }
}





