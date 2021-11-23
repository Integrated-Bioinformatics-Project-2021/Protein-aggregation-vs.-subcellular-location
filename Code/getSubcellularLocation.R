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

# Average tango score distribution in all proteins
plot_average_tango_scores_complete_proteins()
plot_average_tango_scores_complete_proteins_joined_secreted()

# Average tango score distribution in proteins with an APR region
plot_average_tango_scores_APR_proteins()
plot_average_tango_scores_APR_proteins_joined_secreted()

# TODO: Improve documentation

# Distribution of maximal tango scores
plot_max_tango_scores()

# Distribution of number of APR regions per protein normalized for protein length
plot_normalized_number_APR_regions() # Takes some time to load
plot_normalized_number_APR_regions_joined_secreted()

# ------------------ peptides ---------------------------
source("charges.R")
classify_subsets()

# TODO: make plots showing the differences between APR, GK and FR based on charges
# TODO: see if differences when split into subcellular locations
protein_sequences_APR_peptides = get_charge(APR_peptides)
# protein_sequences_GK_peptides = get_charge(GK_peptides)
# protein_sequences_FR_peptides = get_charge(FR_peptides)

#--------- Pie Plot for the residues ------------
# Lysine (Lys) --> K
# Arginine (Arg) --> R
# Aspartic acid (Asp) --> D
# Glutamic acid (Glu) --> E

source("gatekeeper.R")

# Only for the GK regions
GK_analysis = analyse_gate_keeper_regions()

#### PERCENTAGE OF ALL RESIDUES FOR ALL OF THE GK REGIONS IN ALL PROTEINS
pie_plot_percentage_of_all_residues(GK_analysis$cts_gk, "GK REGIONS")

#### PERCENTAGE OF ALL RESIDUES FOR ALL OF THE GK + FL REGIONS IN ALL PROTEINS
pie_plot_percentage_of_all_residues(GK_analysis$cts_gk_fl, "GK + FL REGIONS")

#### PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK REGIONS IN ALL PROTEINS

# Lysine (Lys) --> K
# Arginine (Arg) --> R
# Aspartic acid (Asp) --> D
# Glutamic acid (Glu) --> E
# Serine (Ser) --> S
# Proline (Pro) --> P

GK_analysis_interested_aa = analyse_interested_gate_keeper_regions(GK_analysis)

pie_plot_percentage_of_interested_residues(GK_analysis_interested_aa$cts_interest_gk, "GK REGIONS")

#### PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK + FL REGIONS IN ALL PROTEINS
pie_plot_percentage_of_interested_residues(GK_analysis_interested_aa$cts_interest_gk_fl, "GK + FL REGIONS")

# This one is done for all of the proteins, but now we will be working on all the proteins that
# belong to a specific subcellular location

#### PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK SIDES

GK_analysis_GK_residues = analyse_gate_keeper_residues(GK_analysis$cts_gk_side)
side_string = "GK"

pie_plot_percentage_of_specific_residue(GK_analysis_GK_residues$lys_ct, "LYS", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_residues$arg_ct, "ARG", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_residues$asp_ct, "ASP", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_residues$glu_ct, "GLU", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_residues$ser_ct, "SER", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_residues$pro_ct, "PRO", side_string)

#### PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK + FL SIDES

GK_analysis_GK_and_FL_residues = analyse_gate_keeper_residues(GK_analysis$cts_gk_fl_side)
side_string = "GK + FL"

pie_plot_percentage_of_specific_residue(GK_analysis_GK_and_FL_residues$lys_ct, "LYS", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_and_FL_residues$arg_ct, "ARG", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_and_FL_residues$asp_ct, "ASP", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_and_FL_residues$glu_ct, "GLU", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_and_FL_residues$ser_ct, "SER", side_string)
pie_plot_percentage_of_specific_residue(GK_analysis_GK_and_FL_residues$pro_ct, "PRO", side_string)

# FOR EACH SUBCELLULAR LOCATION

# Statistics for each subcellular location
for (i in 1:length(search_terms)) {
  counts = get_counts_for_subcellular_location(search_terms[i], GK_analysis)
  
  residue_category = "all"
  pie_plot_subcellular_location(counts$cts_gk, residue_category, "GK", search_terms[i])
  pie_plot_subcellular_location(counts$cts_gk_fl, residue_category, "GK + FL", search_terms[i])
  
  # PERCENTAGE OF INTERESTED PROTEINS IN SUBCELLULAR LOCATIONS
  
  counts_interest <- analyse_interested_gate_keeper_regions(counts)
  residue_category = "interested"
  pie_plot_subcellular_location(counts_interest$cts_interest_gk, residue_category, "GK", search_terms[i])
  pie_plot_subcellular_location(counts_interest$cts_interest_gk_fl, residue_category, "GK + FL", search_terms[i])
}
  
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
    