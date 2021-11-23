# All functions dealing with the charges

# Classify into 3 subset dataframes (peptides with only APR, peptides with APR and GK and peptides with APR, GK and 2 FR )
classify_subsets <- function() {
  if (! exists("only_APR_data")) {
    source("tango_scores.R")
    initialize_by_APR_protein()
  }
  APR_peptides <<- subset(only_APR_data, APRdef2_tango < 2)
  GK_peptides <<- subset(only_APR_data, APRdef2_tango < 3)
  FR_peptides <<-  subset(only_APR_data, APRdef2_tango < 5)
}