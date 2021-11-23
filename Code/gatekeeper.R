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