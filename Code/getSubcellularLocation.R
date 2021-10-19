install.packages("UniprotR")
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("GenomicAlignments")
install.packages("sjmisc") # Used for str_contains
install.packages("hash")
library(UniprotR)
library(stringr) # Used to get the last word of a string
library(sjmisc)
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






