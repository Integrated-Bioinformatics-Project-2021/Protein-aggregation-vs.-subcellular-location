# The structural analysis
# Please make sure that Perl is installed on your device:
# https://learn.perl.org/installing/windows.html (for Windows)

# ------------------------------------------------------------------------------
# Steps:
# 1) Download all human AlphaFold structures.
# 2) For each domain in a protein, generate a "sub-structure" based on the Alphafold structures.
# 3) Calculate the contact order of each domain (you can use the perl script from David Bakers lab).
# The R object "domains.RData" contains all domains in human proteins, determined by Gene3D software.
# ------------------------------------------------------------------------------

# 0) Install packages ----------------------------------------------------------
# install.packages("bio3d")
library(bio3d) # Used for read.pdb

# 1) Download all human AlphaFold structures -----------------------------------
# Go to https://www.alphafold.ebi.ac.uk/download and download the entire human dataset.
# Save it to the /Data/AF folder and untar it.
directory = dirname(rstudioapi::getSourceEditorContext()$path) # Place the data in a Data folder on the path of this file
setwd(directory)

# 2) Generate a "sub-structure" ------------------------------------------------
# Check if all proteins in data are present in domains:
load("Data/domains.RData")
proteins_with_subcellular_location <- unique(data$Protein[data$Subcellular_location != ""])
length(proteins_with_subcellular_location)
proteins_with_annotated_domains <- unique(domains$Protein)
length(proteins_with_annotated_domains)
total_list_of_proteins <- Reduce(intersect,list(proteins_with_subcellular_location, proteins_with_annotated_domains))
length(total_list_of_proteins)

# Return a list with for each domain in the given protein the boundaries
get_domain_boundaries <- function(protein) {
  protein_domains = subset(domains, domains$Protein == protein)
  return (protein_domains$boundaries)
}

# Return the structure (pdb file) of the given protein
get_structure_protein <- function(protein) {
  filename = paste(directory,"/Data/AF/AF-", protein, "-F1-model_v1.pdb.gz", sep="")
  file = tryCatch({
    read.pdb(filename)
  }, warning = function(w) {
    warning-handler-code
  }, error = function(e) {
    error_message = paste("There is no AlphaFold structure available for",protein)
    print(error_message)
  }, finally = {
    NULL
  })
  return(file)
}

# Return the sub-structure (pdb file) of the given structure in the given boundaries
get_structure_domain <- function(structure, boundaries) {
  start_position = sub("-.*", "", boundaries)   
  end_position = sub(".*-", "", boundaries)
  all_atoms <- start_position:end_position
  substructure = atom.select(structure, type="ATOM", eleno=all_atoms, value=TRUE)
  filename = paste(directory,"/Output/temp.pdb", sep="")
  write.pdb(pdb=substructure, file=filename)
  return(filename)
}

# For each domain, retrieve the substructure from AlphaFold, calculate its contact and add this to the domains dataframe.
# Save the complete dataframe as "Data/domains_with_contact_order.RData".
save_all_contact_orders <- function() {
  for (i in 1:length(total_list_of_proteins)) {
    print(i)
    current_protein = total_list_of_proteins[i]
    protein_structure = get_structure_protein(current_protein)
    boundaries = get_domain_boundaries(current_protein) # TODO: make more efficient using $begin.domain
    for (j in 1:length(boundaries)) {
      filename = get_structure_domain(protein_structure, boundaries[j])
      # 3) Calculate the contact order of each domain ----------------------------
      # https://depts.washington.edu/bakerpg/contact_order/ (Download Contact Order Program (Perl version))
      # contactOrder.pl <options> <pdbfile>
      # with <options> = -r: returns relative CO
      # with pdbfile the pdb file with the AlphaFold structure of a single domain
      
      cmd = paste("perl contactOrder.pl -r ", filename)
      contact_order_current_structure = system(cmd, intern = TRUE)
      contact_order_current_structure_number = as.numeric(sub(".*: ", "", contact_order_current_structure))
      
      domains$contact_order[domains$Protein == current_protein & domains$boundaries == boundaries[j]] = contact_order_current_structure_number
    }
  }
  save(domains, file = "Data/domains_with_contact_order.RData")
  return (domains)
}

### Mapping: take the maximum tango score for every domain in the dataset ####

# For the given domain, check which APRs from the list APRs_in_protein are within
# the domain. An APR is part of a domain when its middle residue is within the boundaries of the domain.
get_APRs_in_domain <- function(domain, APRs_in_protein) {
  APRs_in_domain = list()
  if (length(APRs_in_protein) == 0) {
    return (APRs_in_domain)
  }
  domain_vector = c(domain$begin.domain:domain$end.domain)
  for (k in 1:length(APRs_in_protein)) {
    residues_in_APR = subset(data, data$APRcount_tango == APRs_in_protein[k])
    min_pos = min(residues_in_APR["Position"])
    max_pos = max(residues_in_APR["Position"])
    # Assume the APR region lies in a domain if its middle residue lies within the domain boundaries
    if (as.integer((max_pos + min_pos)/2) %in% domain_vector) {
      APRs_in_domain = append(APRs_in_domain, APRs_in_protein[k])
    }
  }
  return (APRs_in_domain)
}

# Return the sum of the avgScore attribute of all the APRs in the given list APRs_in_domain.
get_tango_score_domain <- function(APRs_in_domain) {
  tango_score = 0
  if (length(APRs_in_domain) == 0) {
    return (tango_score)
  }
  for (k in 1:length(APRs_in_domain)) {
    current_tango_score = data$avgScore[data$APRcount_tango == APRs_in_domain[k]][1]
    tango_score = tango_score + current_tango_score
  }
  return (tango_score)
}

# For each protein in the list total_list_of_proteins, loop over the domains of the
# protein and add the corresponding tango scores to the domains dataframe.
# Save the complete dataframe as "Data/domains_with_COandTango.RData".
save_tango_per_domain <- function(){
  for (p in 1:length(total_list_of_proteins)) {
    print(p)
    current_protein = total_list_of_proteins[p]
    protein_domains = subset(domains, domains$Protein == current_protein)
    for (d in 1:nrow(protein_domains)){
      APRs_in_protein = unique(data[data$Protein == current_protein & data$APRcount_tango != 0, "APRcount_tango"]) # Get all non-zero APR domains in the protein
      APRs_in_domain = get_APRs_in_domain(protein_domains[d,], APRs_in_protein)
      tango_score_domain = get_tango_score_domain(APRs_in_domain)
      domains[domains$CATH.domain.ID == protein_domains[d,]$CATH.domain.ID, "tango"] = tango_score_domain # Save in tango column of domains dataframe
      }
    }
  save(domains, file = "Data/domains_with_COandTango.RData")
  # return(domains)
}

# MAIN FUNCTION
# Load all the domains information. This function adds a column to the domains data frame
# with the contact order of each domain and the tango score of each domain.
# If this information already has been calculated, load the file "Data/domains_with_COandTango.RData"
# for the full information, or the file "Data/domains_with_contact_order.RData" for the contact order information.
# Calculated the remaining information if needed.
load_full_domains <- function() {
  if (file.exists("Data/domains_with_COandTango.RData")) {
    load("Data/domains_with_COandTango.RData")
  }
  else {
    if (file.exists("Data/domains_with_contact_order.RData")) {
      load("Data/domains_with_contact_order.RData")
    }
    else { # Only retrieve the data when it isn't stored yet
      # load("Data/domains.RData") # Loaded at top of file
      save_all_contact_orders()
    }
    save_tango_per_domain()
  }
  domains <<- domains
}

#General plot 
plot_domains_contact_order <- function() {
  plot(domains$tango, domains$contact_order)
  cor.test(domains$tango, domains$contact_order)
}

#plot grouping the domains for contact_order 
plot_tango_contact_order_grouped <- function(data) {
  labels <- c("<0.05", "0.05-0.1", "0.10-0.15", "0.15-0.20", "0.20-0.25", ">0.25")
  ggplot(data, aes(x= .groups, y = tango)) + 
    geom_boxplot() + scale_color_brewer(palette="Dark2") +
    xlab("contact_order") +
    scale_x_discrete(labels = labels)
}
#creating groupings for contact_order - domains without NA values
get_grouped_domains <- function() {
  domains <- domains[order(domains$contact_order),]
  grouped_domains <- group(domains, n = 0.05, method = "staircase")
  grouped_domains <- grouped_domains[!is.na(grouped_domains$contact_order),]
  return(grouped_domains)
}
#creating groupings for contact_order - NA set to zero
get_grouped_domains_NA_set_zero <- function() {
  domains <- domains[order(domains$contact_order),]
  grouped_domains_NA <- group(domains, n = 0.05, method = "staircase")
  grouped_domains_NA <- grouped_domains_NA[!is.na(grouped_domains_NA$contact_order),]
  grouped_domains_NA$contact_order[is.na(grouped_domains_NA$contact_order)] <- 0
  return(grouped_domains_NA)
}

#density plot for every subcellular location 
plot_2D_contact_order_and_tango_in_subcellular_location <- function() {
  unique_prot_data <- data[!duplicated(data$Protein),]
  domain_and_data <- merge(unique_prot_data, domains)
  for (l in 1: length(search_terms)) {
    domains_in_subcellular_location <- subset(domain_and_data, domain_and_data$Subcellular_location == search_terms[l])
    ggplot(domains_in_subcellular_location, aes(x = tango, y = contact_order)) + 
      geom_density2d()
  }
}

