# All functions dealing with gate keeper regions

analyse_gate_keeper_regions <- function() {
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
  
  return (list("GK_residues" = GK_residues,
               "GK_FL_residues" = GK_FL_residues,
               "cts_gk" = cts_gk,
               "cts_gk_side" = cts_gk_side,
               "cts_gk_fl" = cts_gk_fl,
               "cts_gk_fl_side" = cts_gk_fl_side))
}

pie_plot_percentage_of_all_residues <- function(given_region, region_string) {
  title = paste("PERCENTAGE OF ALL RESIDUES FOR ALL\n OF THE ", region_string,  " IN ALL PROTEINS", sep="", collapse=NULL)
  ggplot(given_region, aes(x = "", y = perc, fill = Residue)) +
    geom_col() +
    geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
    coord_polar(theta = "y")+
    guides(fill = guide_legend(title = "Residue"))+
    theme_void()+
    ggtitle(title)
  #scale_fill_brewer(palette="PiYG")+
  #geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE)
}

# PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE GK REGIONS IN ALL PROTEINS

analyse_interested_gate_keeper_regions <- function(GK_analysis) {
  interested_aa <- c("K","R","D","E","S","P")
  cts_interest_gk <- GK_analysis$cts_gk[which(GK_analysis$cts_gk$Residue %in% interested_aa),colnames(GK_analysis$cts_gk)]
  cts_interest_gk_fl <- GK_analysis$cts_gk_fl[which(GK_analysis$cts_gk_fl$Residue %in% interested_aa),colnames(GK_analysis$cts_gk_fl)]
  return (list("cts_interest_gk" = cts_interest_gk,
               "cts_interest_gk_fl" = cts_interest_gk_fl))
}

pie_plot_percentage_of_interested_residues <- function(given_region, region_string) {
  title = paste("PERCENTAGE OF INTERESTED RESIDUES FOR ALL\n OF THE ", region_string," IN ALL PROTEINS", sep="", collapse=NULL)
  ggplot(given_region, aes(x = "", y = perc, fill = Residue)) +
    geom_col() +
    geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
    coord_polar(theta = "y")+
    guides(fill = guide_legend(title = "Residue"))+
    theme_void()+
    scale_fill_brewer(palette="Dark2")+
    ggtitle(title)
  #geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE)
}

# PERCENTAGE OF INTERESTED RESIDUES FOR ALL OF THE SPECIFIED SIDES
analyse_gate_keeper_residues <- function(sides) {
  lys_ct = sides[sides$Residue == "K",]
  arg_ct = sides[sides$Residue == "R",]
  asp_ct = sides[sides$Residue == "D",]
  glu_ct = sides[sides$Residue == "E",]
  ser_ct = sides[sides$Residue == "S",]
  pro_ct = sides[sides$Residue == "P",]
  
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
  
  return (list("lys_ct" = lys_ct,
               "arg_ct" = arg_ct,
               "asp_ct" = asp_ct,
               "glu_ct" = glu_ct,
               "ser_ct" = ser_ct,
               "pro_ct" = pro_ct))
}

pie_plot_percentage_of_specific_residue <- function(given_region, residue_string, side_string) {
  title = paste("PERCENTAGE OF ", residue_string," RESIDUES FOR\n ALL OF THE ", side_string, " SIDES", sep="", collapse=NULL)
  ggplot(given_region, aes(x = "", y = perc, fill = Side)) +
    geom_col() +
    geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
    coord_polar(theta = "y")+
    guides(fill = guide_legend(title = "Side"))+
    theme_void()+
    ggtitle(title)
}

# FOR EACH SUBCELLULAR LOCATION

get_counts_for_subcellular_location <- function(subcellular_location, GK_analysis) {
  source("subcellular_location_annotation.R")
  
  given_protein_list = get_proteins_with_given_subcellular_location(subcellular_location)
  
  cts_interest_gk <- GK_analysis$GK_residues[which(GK_analysis$GK_residues$Protein %in% given_protein_list$wanted_proteins),colnames(GK_analysis$GK_residues)]
  cts_interest_gk_fl <- GK_analysis$GK_FL_residues[which(GK_analysis$GK_FL_residues$Protein %in% given_protein_list$wanted_proteins),colnames(GK_analysis$GK_FL_residues)]
  
  # FOR ALL GK IN GIVEN SUBCELLULAR LOCATION
  
  cts_gk <- cts_interest_gk %>% 
    group_by(Residue) %>% # Variable to be transformed
    count() %>% 
    ungroup() %>% 
    mutate(perc = `n` / sum(`n`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  
  cts_gk_fl <- cts_interest_gk_fl %>% 
    group_by(Residue) %>% # Variable to be transformed
    count() %>% 
    ungroup() %>% 
    mutate(perc = `n` / sum(`n`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  
  return (list("cts_gk" = cts_gk,
               "cts_gk_fl" = cts_gk_fl))
}

pie_plot_subcellular_location <- function(counts, residue_category, side_string, subcellular_location) {
  title = paste("Percentage of ", residue_category, " residues in ", side_string, " regions\n for ", subcellular_location, sep="", collapse=NULL)
  plot = ggplot(counts, aes(x = "", y = perc, fill = Residue)) +
    geom_col() +
    geom_label_repel(aes(label = labels), max.overlaps = 30,size = 4.5, position = position_stack(vjust = 0.5), show.legend = FALSE)+
    coord_polar(theta = "y")+
    guides(fill = guide_legend(title = "Residue"))+
    theme_void()+
    ggtitle(title)
  print(plot)
}
