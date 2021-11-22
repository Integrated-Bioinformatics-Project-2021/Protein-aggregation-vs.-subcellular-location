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
rm(return_value) # To get rid of the extra memory usage

# ------------------ tango scores ---------------------------
source("tango_scores.R")

plots_average_tango_scores()
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

# Lysine (Lys) --> K
# Arginine (Arg) --> R
# Aspartic acid (Asp) --> D
# Glutamic acid (Glu) --> E

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

cts_gk_fl_side <- GK_FL_residues %>% 
  group_by(Residue, Side) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

#### PERCENTAGE OF ALL RESIDUES FOR ALL OF THE GK REGIONS IN ALL PROTEINS

ggplot(cts_gk, aes(x = "", y = perc, fill = Residue)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Residue"))+
  theme_void()+
  ggtitle("PERCENTAGE OF ALL RESIDUES FOR ALL\n OF THE GK REGIONS IN ALL PROTEINS")
  #scale_fill_brewer(palette="PiYG")+
  #geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE)

#### PERCENTAGE OF ALL RESIDUES FOR ALL OF THE GK + FL REGIONS IN ALL PROTEINS

ggplot(cts_gk_fl, aes(x = "", y = perc, fill = Residue)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Residue"))+
  theme_void()+
  ggtitle("PERCENTAGE OF ALL RESIDUES FOR ALL\n OF THE GK + FL REGIONS IN ALL PROTEINS")
  #scale_fill_brewer(palette="PiYG")+
  #geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE)

#### PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK REGIONS IN ALL PROTEINS

# Lysine (Lys) --> K
# Arginine (Arg) --> R
# Aspartic acid (Asp) --> D
# Glutamic acid (Glu) --> E
# Serine (Ser) --> S
# Proline (Pro) --> P
interested_aa <- c("K","R","D","E","S","P")
cts_interest_gk <- cts_gk[which(cts_gk$Residue %in% interested_aa),colnames(cts_gk)]
cts_interest_gk_fl <- cts_gk_fl[which(cts_gk_fl$Residue %in% interested_aa),colnames(cts_gk_fl)]

ggplot(cts_interest_gk, aes(x = "", y = perc, fill = Residue)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Residue"))+
  theme_void()+
  scale_fill_brewer(palette="Dark2")+
  ggtitle("PERCENTAGE OF INTERESTED RESIDUES FOR ALL\n OF THE GK REGIONS IN ALL PROTEINS")
#geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE)

#### PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK + FL REGIONS IN ALL PROTEINS

ggplot(cts_interest_gk_fl, aes(x = "", y = perc, fill = Residue)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Residue"))+
  theme_void()+
  scale_fill_brewer(palette="Dark2")+
  ggtitle("PERCENTAGE OF INTERESTED RESIDUES FOR ALL\n OF THE GK + FL REGIONS IN ALL PROTEINS")

# This one is done for all of the proteins, but now we will be working on all the proteins that
# is belonging to a specific subcellular location

#### PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK SIDES

lys_ct_gk = cts_gk_side[cts_gk_side$Residue == "K",]
arg_ct_gk = cts_gk_side[cts_gk_side$Residue == "R",]
asp_ct_gk = cts_gk_side[cts_gk_side$Residue == "D",]
glu_ct_gk = cts_gk_side[cts_gk_side$Residue == "E",]
ser_ct_gk = cts_gk_side[cts_gk_side$Residue == "S",]
pro_ct_gk = cts_gk_side[cts_gk_side$Residue == "P",]

lys_ct_gk <- lys_ct_gk %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

arg_ct_gk <- arg_ct_gk %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

asp_ct_gk <- asp_ct_gk %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

glu_ct_gk <- glu_ct_gk %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

ser_ct_gk <- ser_ct_gk %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

pro_ct_gk <- pro_ct_gk %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

ggplot(lys_ct_gk, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF LYS RESIDUES FOR\n ALL OF THE GK SIDES")

ggplot(arg_ct_gk, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF ARG RESIDUES FOR\n ALL OF THE GK SIDES")

ggplot(asp_ct_gk, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF ASP RESIDUES FOR\n ALL OF THE GK SIDES")

ggplot(glu_ct_gk, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF GLU RESIDUES FOR\n ALL OF THE GK SIDES")

ggplot(ser_ct_gk, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF SER RESIDUES FOR\n ALL OF THE GK SIDES")

ggplot(pro_ct_gk, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF PRO RESIDUES FOR\n ALL OF THE GK SIDES")

#### PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK + FL SIDES

lys_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "K",]
arg_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "R",]
asp_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "D",]
glu_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "E",]
ser_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "S",]
pro_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "P",]

lys_ct_gk_fl <- lys_ct_gk_fl %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

arg_ct_gk_fl <- arg_ct_gk_fl %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

asp_ct_gk_fl <- asp_ct_gk_fl %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

glu_ct_gk_fl <- glu_ct_gk_fl %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

ser_ct_gk_fl <- ser_ct_gk_fl %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

pro_ct_gk_fl <- pro_ct_gk_fl %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

ggplot(lys_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF LYS RESIDUES FOR\n ALL OF THE GK + FL SIDES")

ggplot(arg_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF ARG RESIDUES FOR\n ALL OF THE GK + FL SIDES")

ggplot(asp_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF ASP RESIDUES FOR\n ALL OF THE GK + FL SIDES")

ggplot(glu_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF GLU RESIDUES FOR\n ALL OF THE GK + FL SIDES")

ggplot(ser_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF SER RESIDUES FOR\n ALL OF THE GK + FL SIDES")

ggplot(pro_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
  geom_col() +
  geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
  coord_polar(theta = "y")+
  guides(fill = guide_legend(title = "Side"))+
  theme_void()+
  ggtitle("PERCENTAGE OF PRO RESIDUES FOR\n ALL OF THE GK + FL SIDES")

# FOR EACH SUBCELLULAR LOCATION

total_length = length(search_terms)
total_length

#Give i any value to see the statistics for that subcellular location
i = 9

  given_protein_list = get_proteins_with_given_subcellular_location(search_terms[i])
  
  wanted <- GK_residues[which(GK_residues$Protein %in% given_protein_list$wanted_proteins),colnames(GK_residues)]
  wanted2 <- GK_FL_residues[which(GK_FL_residues$Protein %in% given_protein_list$wanted_proteins),colnames(GK_FL_residues)]
  
  # FOR ALL GK IN GIVEN SUBCELLULAR LOCATION
  
  cts_gk_s <- wanted %>% 
    group_by(Residue) %>% # Variable to be transformed
    count() %>% 
    ungroup() %>% 
    mutate(perc = `n` / sum(`n`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  
  cts_gk_fl_s <- wanted2 %>% 
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
    theme_void()+
    ggtitle("Percentage of every residue in GK regions for",search_terms[[i]])
  
  ggplot(cts_gk_fl_s, aes(x = "", y = perc, fill = Residue)) +
    geom_col() +
    geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
    coord_polar(theta = "y")+
    guides(fill = guide_legend(title = "Residue"))+
    theme_void()+
    ggtitle("Percentage of every residue in GK + FL regions for",search_terms[[i]])
  
  # PERCENTAGE OF INTERESTED PROTEINS IN SUBCELLULAR LOCATIONS
  
  cts_interest_gk <- cts_gk_s[which(cts_gk_s$Residue %in% interested_aa),colnames(cts_gk_s)]
  cts_interest_gk_fl <- cts_gk_fl_s[which(cts_gk_fl_s$Residue %in% interested_aa),colnames(cts_gk_fl_s)]
  
  ggplot(cts_interest_gk, aes(x = "", y = perc, fill = Residue)) +
    geom_col() +
    geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
    coord_polar(theta = "y")+
    guides(fill = guide_legend(title = "Residue"))+
    theme_void()+
    scale_fill_brewer(palette="Dark2")+
    ggtitle("Percentage of interested residues in GK regions for",search_terms[[i]])
  
  ggplot(cts_interest_gk_fl, aes(x = "", y = perc, fill = Residue)) +
    geom_col() +
    geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
    coord_polar(theta = "y")+
    guides(fill = guide_legend(title = "Residue"))+
    theme_void()+
    scale_fill_brewer(palette="Dark2")+
    ggtitle("Percentage of interested residues in GK + FL regions for",search_terms[[i]])
  
  # NOW THE SIDES PART
  
  cts_gk_side <- wanted %>% 
    group_by(Residue, Side) %>% # Variable to be transformed
    count() %>% 
    ungroup() %>% 
    mutate(perc = `n` / sum(`n`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  
  cts_gk_fl_side <- wanted2 %>% 
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
    ser_ct = cts_gk_side[cts_gk_side$Residue == "S",]
    pro_ct = cts_gk_side[cts_gk_side$Residue == "P",]
    
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
    
    ser_ct <- ser_ct %>%
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    pro_ct <- pro_ct %>%
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    lys_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "K",]
    arg_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "R",]
    asp_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "D",]
    glu_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "E",]
    ser_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "S",]
    pro_ct_gk_fl = cts_gk_fl_side[cts_gk_fl_side$Residue == "P",]
    
    lys_ct_gk_fl <- lys_ct_gk_fl %>%
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    arg_ct_gk_fl <- arg_ct_gk_fl %>%
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    asp_ct_gk_fl <- asp_ct_gk_fl %>%
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    glu_ct_gk_fl <- glu_ct_gk_fl %>%
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    ser_ct_gk_fl <- ser_ct_gk_fl %>%
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    pro_ct_gk_fl <- pro_ct_gk_fl %>%
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    ggplot(lys_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Lys residue sides in GK regions for",search_terms[[i]])
    
    ggplot(arg_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Arg residue sides in GK regions for",search_terms[[i]])
    
    ggplot(asp_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Asp residue sides in GK regions for",search_terms[[i]])
    
    ggplot(glu_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Glu residue sides in GK regions for",search_terms[[i]])
    
    ggplot(ser_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Ser residue sides in GK regions for",search_terms[[i]])
    
    ggplot(pro_ct, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Pro residue sides in GK regions for",search_terms[[i]])
    
    # FOR GK + FL

    ggplot(lys_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Lys residue sides in GK + FL regions for",search_terms[[i]])
    
    ggplot(arg_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Arg residue sides in GK + FL regions for",search_terms[[i]])
    
    ggplot(asp_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Asp residue sides in GK + FL regions for",search_terms[[i]])
    
    ggplot(glu_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Glu residue sides in GK + FL regions for",search_terms[[i]])
    
    ggplot(ser_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Ser residue sides in GK + FL regions for",search_terms[[i]])
    
    ggplot(pro_ct_gk_fl, aes(x = "", y = perc, fill = Side)) +
      geom_col() +
      geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
      coord_polar(theta = "y")+
      guides(fill = guide_legend(title = "Side"))+
      theme_void()+
      ggtitle("Percentage of Pro residue sides in GK + FL regions for",search_terms[[i]])
    