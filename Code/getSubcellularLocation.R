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
  
  # NOW THE SIDES PART
  counts_side_gk <- analyse_sides(counts$cts_gk_side)
  counts_side_gk_fl <- analyse_sides(counts$cts_gk_fl_side)

  region_string = "GK"
  pie_plot_sides(counts_side_gk$lys_ct, "LYS", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk$arg_ct, "ARG", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk$asp_ct, "ASP", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk$glu_ct, "GLU", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk$ser_ct, "SER", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk$pro_ct, "PRO", region_string, search_terms[i])
  
  region_string = "GK + FL"
  pie_plot_sides(counts_side_gk_fl$lys_ct, "LYS", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk_fl$arg_ct, "ARG", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk_fl$asp_ct, "ASP", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk_fl$glu_ct, "GLU", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk_fl$ser_ct, "SER", region_string, search_terms[i])
  pie_plot_sides(counts_side_gk_fl$pro_ct, "PRO", region_string, search_terms[i])
}
    