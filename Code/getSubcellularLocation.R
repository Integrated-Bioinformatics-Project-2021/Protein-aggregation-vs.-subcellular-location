#install.packages("UniprotR")
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("GenomicAlignments")
#install.packages("sjmisc")
#install.packages("hash")
#install.packages("dplyr")
#install.packages("Peptides")
library(UniprotR)
library(sjmisc) # Used for str_contains
library(hash) # Used to make a disctionary
library(dplyr) # Used for aggregate (get avg tango score of protein)
library(Peptides) # Used for charge (get the charge of a peptide sequence) ETC
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyverse)

## Defining the working directory
directory = dirname(rstudioapi::getSourceEditorContext()$path) # Should work when data is placed in same folder
setwd(directory)

## Reading the data
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

# Return a list with all the protein IDs of that contain the given subcellular location.
get_proteins_with_given_subcellular_location <- function(given_subcellular_location){
  wanted_proteins <- sapply(keys(hashed_proteins), function(x) grepl(given_subcellular_location, hashed_proteins[[x]][["subcellular_location"]]))
  wanted_proteins <- keys(hashed_proteins)[wanted_proteins]
  return (data.frame(wanted_proteins))
}


# ------------------ tango scores ---------------------------

# Acquire only the vital rows of each protein (APRs)
only_APR_data <- subset(data, APRdef2_tango > 0)


# Get average tango score for each complete protein 
unique_data <- data[!((data$APRcount_tango == 0) & (data$maxProtscore != 0)),]
unique_data <- unique_data[!duplicated(unique_data[,c(1,10)]),]
by_protein <- aggregate(avgScore ~ Protein, unique_data, mean)

# Get average tango score for each APR protein
unique_APR_data <- only_APR_data[!duplicated(only_APR_data[,c(1,10)]),]
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

#get average max Tango score for each subcellular locations for complete proteins 
get_avgmax_tango_score_complete_protein <- function(given_protein_list) {
  avgScores = c()
  for (j in 1:nrow(given_protein_list)){
    if (given_protein_list$wanted_proteins[j] %in% by_protein_max$Protein) {
      avgScores <- append(avgScores,by_protein_max[which(by_protein_max$Protein %in% given_protein_list$wanted_proteins[j]), "maxProtscore"])
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

#barplot(tango_scores, xlab = "subcellular location", names.arg = search_terms)

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

# ----------- Max Tango Scores for each protein and APR -----------------------

by_protein_max = aggregate(maxProtscore~Protein, data, max)
scores_protein = data.frame()

for(i in 1:length(search_terms))
{
  given_protein_list = get_proteins_with_given_subcellular_location(search_terms[i])
  scores_protein[nrow(scores_protein)+1:nrow(given_protein_list),1] = given_protein_list
  
  newIndex_CP = nrow(scores_protein)+1-nrow(given_protein_list)
  scores_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),2] = rep(search_terms[[i]], nrow(given_protein_list))
  avgMaxScores_protein = get_avgmax_tango_score_complete_protein(given_protein_list)
  
  scores_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),3] = avgMaxScores_protein
  scores_protein[newIndex_CP:(newIndex_CP+nrow(given_protein_list)-1),4] = rep(check_secretory(search_terms[i]), nrow(given_protein_list))
}
colnames(scores_protein) <- c("Proteins", "Subcellular_location", "Tango_scores", "Secretory")

box_avgmax_complete_proteins <- ggplot(scores_protein, aes(x=Tango_scores, y = Subcellular_location, fill = Secretory)) + 
  geom_boxplot(notch=TRUE) + scale_color_brewer(palette="Dark2")
box_avgmax_complete_proteins + theme_minimal() + stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red")

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

#--------- Pie Plot for the residues ------------

# Only for the GK regions
GK_residues <- subset(data, APRdef2_tango == 2, select = c("Protein","Residue","APRdef2_tango","Side"))
GK_FL_residues <- subset(data, APRdef2_tango > 1, select = c("Protein","Residue","APRdef2_tango","Side"))
#group <- unique(GK_residues$Residue)
#cts = as.data.frame(table(GK_residues$Residue))
#colnames(cts) <-  c("Residue","Count")
#pie = ggplot(cts, aes(x="", y=Count, fill=Residue)) + geom_bar(stat="identity", width=1)

cts_gk <- GK_residues %>% 
  group_by(Residue) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

cts_gk_side <- GK_residues %>% 
  group_by(Residue, Side) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

cts_gk_fl <- GK_FL_residues %>% 
  group_by(Residue) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

# This one is done for all of the proteins, but now we will be working on all the proteins that
# is belonging to a specific subcellular location
ggplot(cts_gk, aes(x = "", y = perc, fill = Residue)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Residue"))+
  theme_void()
  #scale_fill_brewer(palette="PiYG")+
  #geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE)

# Lysine (Lys) --> K
# Arginine (Arg) --> R
# Aspartic acid (Asp) --> D
# Glutamic acid (Glu) --> E

lys_ct = cts_gk_side[cts_gk_side$Residue == "K",]
arg_ct = cts_gk_side[cts_gk_side$Residue == "R",]
asp_ct = cts_gk_side[cts_gk_side$Residue == "D",]
glu_ct = cts_gk_side[cts_gk_side$Residue == "E",]

lys_ct <- lys_ct %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

arg_ct <- arg_ct %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

asp_ct <- asp_ct %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

glu_ct <- glu_ct %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

ggplot(lys_ct, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()

ggplot(arg_ct, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()

ggplot(asp_ct, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()

ggplot(glu_ct, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()

# For each subcellular location
total_length = length(search_terms)
total_length

#Give i any value to see the statistics for that subcellular location
i = 1

  given_protein_list = get_proteins_with_given_subcellular_location(search_terms[i])
  
  # FOR ONLY GK
  
  wanted <- GK_residues[which(GK_residues$Protein %in% given_protein_list$wanted_proteins),colnames(GK_residues)]
  
  cts_gk_s <- wanted %>% 
    group_by(Residue) %>% # Variable to be transformed
    count() %>% 
    ungroup() %>% 
    mutate(perc = `n` / sum(`n`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  
  ggplot(cts_gk_s, aes(x = "", y = perc, fill = Residue)) +
    geom_col() +
    geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
    coord_polar(theta = "y")+
    guides(fill = guide_legend(title = "Residue"))+
    theme_void()
  
  cts_gk_side <- wanted %>% 
    group_by(Residue, Side) %>% # Variable to be transformed
    count() %>% 
    ungroup() %>% 
    mutate(perc = `n` / sum(`n`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  
    lys_ct = cts_gk_side[cts_gk_side$Residue == "K",]
    arg_ct = cts_gk_side[cts_gk_side$Residue == "R",]
    asp_ct = cts_gk_side[cts_gk_side$Residue == "D",]
    glu_ct = cts_gk_side[cts_gk_side$Residue == "E",]
    
    ggplot(lys_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()
    
    ggplot(arg_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()
    
    ggplot(asp_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()
    
    ggplot(glu_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()
    
    # FOR GK + FR
    
    wanted <- GK_FL_residues[which(GK_FL_residues$Protein %in% given_protein_list$wanted_proteins),colnames(GK_FL_residues)]
    
    cts_gk_fl_s <- wanted %>% 
      group_by(Residue) %>% # Variable to be transformed
      count() %>% 
      ungroup() %>% 
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    ggplot(cts_gk_fl_s, aes(x = "", y = perc, fill = Residue)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Residue"))+
      theme_void()
    
    cts_gk_fl_side <- wanted %>% 
      group_by(Residue, Side) %>% # Variable to be transformed
      count() %>% 
      ungroup() %>% 
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    lys_ct = cts_gk_fl_side[cts_gk_fl_side$Residue == "K",]
    arg_ct = cts_gk_fl_side[cts_gk_fl_side$Residue == "R",]
    asp_ct = cts_gk_fl_side[cts_gk_fl_side$Residue == "D",]
    glu_ct = cts_gk_fl_side[cts_gk_fl_side$Residue == "E",]

    ggplot(lys_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()
    
    ggplot(arg_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()
    
    ggplot(asp_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()
    
    ggplot(glu_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()
    
# Make a histogram

