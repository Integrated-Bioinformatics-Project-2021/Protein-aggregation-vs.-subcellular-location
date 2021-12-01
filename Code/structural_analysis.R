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

save_all_contact_orders <- function() {
  for (i in 1:length(total_list_of_proteins)) {
    print(i)
    current_protein = total_list_of_proteins[i]
    protein_structure = get_structure_protein(current_protein)
    boundaries = get_domain_boundaries(current_protein)
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
      
      # TODO: make a list/hash with every contact_order_current_structure
    }
  }
  save(domains, file = "Data/domains_with_contact_order.RData")
  return (domains)
}

main <- function() {
  if (! file.exists("Data/domains_with_contact_order.RData")) { # Only retrieve the data when it isn't stored yet
    load("Data/domains.RData")
    save_all_contact_orders()
  }
  else {
    load("Data/domains_with_contact_order.RData")
  }
  domains <<- domains
}

main()
