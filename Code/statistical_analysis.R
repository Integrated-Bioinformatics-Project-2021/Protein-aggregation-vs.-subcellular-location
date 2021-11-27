# All functions dealing with significance tests

#install.packages("dgof")
library(dgof)
check_significance <- function(x,y){
  res <- ks.test(x, y)
  pvalue <- res$p.value
  if (pvalue < 0.05){
    return(TRUE)
  }
  else{return(FALSE)}
}

## Difference between TANGO scores distributions in subcellular locations belonging to secretory vs non-secretory pathways #### 

# It can take as parameter the following datasets:
# * tango_scores_APR_protein
# * max_tango_scores
# * tango_scores_complete_protein

significance_secretory <- function(data){
  x <- data[data$Secretory==TRUE,]$Tango_scores
  y <- data[data$Secretory==FALSE,]$Tango_scores
  return(check_significance(x,y))
}


