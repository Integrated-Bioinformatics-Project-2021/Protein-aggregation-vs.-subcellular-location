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
source("load_data.R")
source("subcellular_location_annotation.R")

return_value = read_data() # See load_data.R
data = return_value$data
hashed_proteins = return_value$hashed_proteins

# ------------------ tango scores ---------------------------
source("tango_scores.R")

plot_average_tango_scores_complete_proteins()
plot_average_tango_scores_complete_proteins_joined_secreted()
plot_average_tango_scores_APR_proteins()
plot_average_tango_scores_APR_proteins_joined_secreted()

# TODO for each function: split function calculate_average_tango_scores and plots_average_tango_scores
# this way the plots can be replotted without redoing the calculations
# Check if calc has ran before executing plot
# Plot separate plots
# Improve documentation

plots_max_tango_scores()

plots_normalized_number_APR_regions() # Takes some time to load

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

