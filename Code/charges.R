# All functions dealing with the charges

# Classify into 3 subset dataframes (peptides with only APR, peptides with APR and GK and peptides with APR, GK and 2 FR )
classify_subsets <- function() {
  if (! exists("only_APR_data")) {
    source("tango_scores.R")
    initialize_by_APR_protein()
  }
  APR_peptides <<- subset(only_APR_data, APRdef2_tango < 2)
  GK_peptides <<- subset(only_APR_data, APRdef2_tango < 3) # Most interesting group + group on intra/extra --> APR part can stay bc prob no charge
  FR_peptides <<-  subset(only_APR_data, APRdef2_tango < 5)
}

# For every APR in the list peptides_name, calculate the net charge
get_charge <- function(peptides_name, ph = 7) {
  protein_charges = hash()
  unique_APRs = unique(peptides_name$APRcount_tango)
  for (i in 1:length(unique_APRs)) {
    peptides_for_APR = subset(peptides_name, peptides_name$APRcount_tango == unique_APRs[i])
    sequence = paste(peptides_for_APR$Residue, collapse = "")
    protein_charges[unique_APRs[i]] <- charge(sequence, pH = ph)
    key = toString(unique_APRs[i])
    #print(protein_charges[[key]])
  }
  return (protein_charges)
}