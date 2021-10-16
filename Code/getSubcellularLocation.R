install.packages("UniprotR")
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("GenomicAlignments")
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
locations = GetSubcellular_location(proteins[5000:5005], directorypath = NULL)
subcellular_locations = locations[1] # Remove NA's from intramembrane, topological domain and transmembrane information
subcellular_location = matrix(ncol = 1, nrow = 0) 
locations_short = data.frame(subcellular_location) # creating a dataframe to add the subcellular locations
for (i in 1:dim(subcellular_locations)[1]) { #Deleting everything between {} and the "SUBCELLULAR LOCATION"
  string_location = subcellular_locations[i,]
  while (grepl("\\{", string_location[length(string_location)])){
    string_location = sub(" \\{.*?\\}.", ",", string_location)
  }
  string_location = gsub("S.*?: ", "", string_location)
  string_location = substr(string_location,1, nchar(string_location)-1)
  print(string_location)
  locations_short[nrow(locations_short) + 1,] = string_location #adding each value to the dataframe
}

#for (i in 1:dim(subcellular_locations)[1]) { # Iterate over every protein in the list
  #string_location = subcellular_locations[i,]
  #while (grepl("\\{", string_location[length(string_location)])) { # As long as there is a "{" left, splitting needs to happen
    #string_location[length(string_location)] = strsplit(toString(string_location[length(string_location)]), "{", perl = TRUE)
  #}
  #string_location = unlist(string_location)
  #n = length(string_location)-1
  #if (n > 0) {
    #string_location = string_location[1:n]
    #for (j in 1:n) { # For every subcellular location, remove the unnecessary information
     #string_location[j] = sub(".*:", "", string_location[j])
      #if (grepl("\\.", string_location[j])) {
        #string_location[j] = sub(".subc*\\.", "", string_location[j])
      #}
    #}
    
  #}
  #print(string_location)
#}

#Adding the subcellular location to the dataset
locations_short<- cbind(locations_short, Protein = row.names(subcellular_locations)) #creating a second column Protein
complete_data = merge(locations_short, human_proteome, by="Protein", all.y = TRUE) #merging the dataframes on column Protein, keeping non matches



