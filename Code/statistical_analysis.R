# All functions dealing with significance tests

library(dgof)
library(ggpubr) # Used for ggline

#SUBCELLULAR LOCATION
#FUNCTION TO CALCULATE SIGNIFICANCE BETWEEN THE SUBCELLULAR LOCATIONS FOR tango_scores_complete_protein, tango_scores_APR_protein, max_tango_scores
#OUTPUT IS THE OUTPUT FROM THE TEST, BUT MORE IMPORTANTLY A TABLE WITH FOR EACH RELATION OF SUBCELLULAR LOCATIONS THE SIGNIFICANCE IN P-VALUE
significance_subcellular <- function(data) {
	sig_subc_APR_prot <- kruskal.test(Tango_scores ~ Subcellular_location, data = data)
	pairwise.wilcox.test(data$Tango_scores, data$Subcellular_location,
                     p.adjust.method = "BH")
}
#FUNCTION TO CALCULATE SIGNIFICANCE BETWEEN THE SUBCELLULAR LOCATIONS FOR normalized_number_APR_regions
#OUTPUT IS THE OUTPUT FROM THE TEST, BUT MORE IMPORTANTLY A TABLE WITH FOR EACH RELATION OF SUBCELLULAR LOCATIONS THE SIGNIFICANCE IN P-VALUE
significance_subcellular_nbAPR <- function(data) {
	sig_subc_APR_prot <- kruskal.test(Nb_APRs ~ Subcellular_location, data = data)
	pairwise.wilcox.test(data$Nb_APRs, data$Subcellular_location,
                     p.adjust.method = "BH")
}

#FUNTION TO CALCULATE THE SIGNIFICANCE IN GENERAL FOR tango_scores_complete_protein, tango_scores_APR_protein, max_tango_scores
get_krustal_test_results <- function(data) {
	sig_subc <- kruskal.test(Tango_scores ~ Subcellular_location, data =  data)
	return (sig_subc)
}

#FUNTION TO CALCULATE THE SIGNIFICANCE IN GENERAL FOR normalized_number_APR_regions
get_krustal_test_results_nbAPR <- function(data) {
	sig_subc <- kruskal.test(Nb_APR ~ Subcellular_location, data =  data)
	return (sig_subc)
}


#FUNTION TO CALCULATE THE SIGNIFICANCE BETWEEN ALL SUBCELLULAR LOCATIONS FOR tango_scores_complete_protein, tango_scores_APR_protein, max_tango_scores
#RETURNS A TABLE WITH ALL THE P-VALUES
get_wilcox_test_table <- function(data) {
	wilcox_sig_subc_com_prot<-pairwise.wilcox.test(data$Tango_scores, data$Subcellular_location,
                     p.adjust.method = "BH")
	table_sig_subc_com_prot<-wilcox_sig_subc_com_prot$p.value

	print(table_sig_subc_com_prot)
	
	file_path = paste(directory, "/Output/table_sig_subc_com_prot.csv", sep="", collapse = NULL)
	write.csv(table_sig_subc_com_prot, file = file_path)
}

library(ggpubr)
library(gtools)

#FUNTION TO CALCULATE THE SIGNIFICANCE BETWEEN ALL SUBCELLULAR LOCATIONS FOR normalized_number_APR_regions
#RETURNS A TABLE WITH ALL THE P-VALUES
get_wilcox_test_table_nbAPR <- function(data) {
	wilcox_sig_subc_com_prot<-pairwise.wilcox.test(data$Nb_APR, data$Subcellular_location,
                     p.adjust.method = "BH")
	table_sig_subc_com_prot<-wilcox_sig_subc_com_prot$p.value

	file_path = paste(directory, "/Output/table_sig_subc_com_prot.csv", sep="", collapse = NULL)
	write.csv(table_sig_subc_com_prot, file = file_path)
}

plot_stat_significance <- function(data){
  ggplot(data, aes(x=Subcellular_location, y = Tango_scores, fill = Secretory),
         color = "supp", palette = "jco",
         add = "jitter") + 
    rotate_x_text(angle = 45)+
    geom_violin() + 
    geom_hline(yintercept = mean(data$Tango_scores), linetype = 2)+
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    stat_compare_means(label = "p.signif", method = "wilcox.test",  ref.group = ".all.",  hide.ns = TRUE)
}

plot_stat_significance_nbAPR <- function(data){
  ggplot(data, aes(x=Subcellular_location, y = Nb_APR, fill = Secretory),
         color = "supp", palette = "jco",
         add = "jitter") + 
    rotate_x_text(angle = 45)+
    geom_violin() + 
    geom_hline(yintercept = mean(data$Nb_APR), linetype = 2)+
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    stat_compare_means(label = "p.signif", method = "wilcox.test",  ref.group = ".all.",  hide.ns = TRUE)
}



group_by(tango_scores_complete_protein, Subcellular_location) %>%
  summarise(
    count = n(),
    mean = mean(Tango_scores, na.rm = TRUE),
    sd = sd(Tango_scores, na.rm = TRUE),
    median = median(Tango_scores, na.rm = TRUE),
    IQR = IQR(Tango_scores, na.rm = TRUE)
  )

ggline(tango_scores_complete_protein, x = "Subcellular_location", y = "Tango_scores", 
       add = c("mean_se", "jitter"), 
       order = c("Cell membrane", "Mitochondrion", "Nucleus", "Endoplasmic Reticulum", "Golgi apparatus", "Lysosome", "Cytoplasm", "Secreted", "Extracellular space"),
       ylab = "Tango_scores", xlab = "Subcellular_location")


#SECRETORY VS NON_SECRETORY
check_significance_sec <- function(x,y){
  res <- ks.test(x, y)
  return (res)
}

## Difference between TANGO scores distributions in subcellular locations belonging to secretory vs non-secretory pathways #### 

# It can take as parameter the following datasets:
# * tango_scores_APR_protein
# * max_tango_scores
# * tango_scores_complete_protein

#FUNCTION TO CALCULATE SIGNIFICANCE, returning D-value and p-value 
#FUNCTION USED FOR tango_scores_complete_protein, tango_scores_APR_protein, max_tango_scores
significance_secretory <- function(data){
  x <- data$Tango_scores[data$Secretory==TRUE]
  y <- data$Tango_scores[data$Secretory==FALSE]
  return(check_significance_sec(x,y))
}
#FUNCTION USED FOR normalized_number_APR_regions
significance_secretory_nbAPR <- function(data){
  x <- data[data$Secretory==TRUE,]$Nb_APRs
  y <- data[data$Secretory==FALSE,]$Nb_APRs
  return(check_significance_sec(x,y))
}

#FUNCTION FOR VISUALISATION OF THE KS-TEST FOR tango_scores_complete_protein, tango_scores_APR_protein, max_tango_scores
plot_ks_test <- function(data) {
	Secretory <- data$Tango_scores[data$Secretory == TRUE]
	Non_secretory <- data$Tango_scores[data$Secretory == FALSE]
	group <- c(rep("Secretory", length(Secretory)), rep("Non_secretory", length(Non_secretory)))
	data_s1_s2<- data.frame(KSD = c(Secretory,Non_secretory), group = group)

	cdf1 <- ecdf(Secretory)
	cdf2 <- ecdf(Non_secretory)

	minMax <- seq(min(Secretory, Non_secretory), max(Secretory, Non_secretory), length.out=length(Secretory)) 
	x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
	y0 <- cdf1(x0) 
	y1 <- cdf2(x0) 

	ggplot(data_s1_s2, aes(x = KSD, group = group, color = group)) + stat_ecdf(size=1) + theme_bw(base_size = 28) + theme(legend.position ="top") + xlab("Group") + ylab("ECDF") + geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]), linetype = "dashed", color = "red") + geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=8) + geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=8) + ggtitle("K-S Test:\nSecretory / Non_secretory") + theme(legend.title=element_blank())
	}

#FUNCTION FOR VISUALISATION OF THE KS-TEST FOR tango_scores_complete_protein, tango_scores_APR_protein, max_tango_scores
plot_ks_test_nbAPR <- function(data) {
	sample1 <- data$Nb_APRs[data$Secretory == TRUE]
	sample2 <- data$Nb_APRs[data$Secretory == FALSE]
	group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
	data_s1_s2<- data.frame(KSD = c(sample1,sample2), group = group)

	cdf1 <- ecdf(sample1)
	cdf2 <- ecdf(sample2)

	minMax <- seq(min(sample1, sample2), max(sample1, sample2), length.out=length(sample1)) 
	x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
	y0 <- cdf1(x0) 
	y1 <- cdf2(x0) 

	ggplot(data_s1_s2, aes(x = KSD, group = group, color = group)) + stat_ecdf(size=1) + theme_bw(base_size = 28) + theme(legend.position ="top") + xlab("Sample") + ylab("ECDF") + geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]), linetype = "dashed", color = "red") + geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=8) + geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=8) + ggtitle("K-S Test: Sample 1 / Sample 2") + theme(legend.title=element_blank())
	}



