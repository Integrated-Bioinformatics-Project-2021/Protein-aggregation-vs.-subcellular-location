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
  
  return (list("cts_gk" = cts_gk,
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