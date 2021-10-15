#install.packages("UniprotR")
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("GenomicAlignments")
library(UniprotR)
library(stringr) # Used to get the last word of a string

## Defining the working directory
directory = dirname(rstudioapi::getSourceEditorContext()$path) # Should work when data is placed in same folder
setwd(directory)

## Reading the data
load("Data/human_proteome_df.RData")
data = human_proteome

## Get all unique protein identifiers
proteins = unique(data$Protein)
proteins

## Get subcellular location of every protein
locations = GetSubcellular_location(proteins[0:5], directorypath = NULL)
subcellular_locations = locations[1] # Remove NA's from intramembrane, topological domain and transmembrane information
for (i in 1:dim(subcellular_locations)[1]) { # Iterate over every protein in the list
  string_location = subcellular_locations[i,]
  while (grepl("\\{", string_location[length(string_location)])) { # As long as there is a "{" left, splitting needs to happen
    string_location[length(string_location)] = strsplit(toString(string_location[length(string_location)]), "{", fixed = TRUE)
  }
  string_location = unlist(string_location)
  n = length(string_location)-1
  string_location = string_location[1:n]
  for (j in 1:n) { # For every subcellular location, remove the unnecessary information
    string_location[j] = sub(".*:", "", string_location[j])
    if (grepl("\\.", string_location[j])) {
      string_location[j] = sub(".*\\.", "", string_location[j])
    }
  }
  print(string_location)
}

# TODO: Add all the subcellular locations to a dictionary mapped to the corresponding Uniprot ID