#!/usr/bin/env Rscript
library('Quartet')
library('TreeDist')
library(treeio)
library(glue)


# change directory to where the nwk files are
setwd('~/Documents/grad_school/CMSC829A/') 

# Init empty vectors for tree distance metrics to be tested
tree_dis_all_patients <- list()
patient_list<- c('CRC01','CRC02','CRC04','CRC10','CRC11','CRC12','CRC13','CRC14','CRC15')
last_i_list <- c(40,39,40,40,40,24,40,29,40)
#last iteration for each patient is different 
#as some patients are so sparse that they don't have eligible sites left 
#after certain number of iterations

for (p in (1:length(patient_list)))   {
patient <- patient_list[p]
last_i <- last_i_list[p]
#init empty vectors 
RF<-c()
NYE<-c()
MS<-c()
QD<-c()

for (i in (1:last_i)) {
  i_1=i-1
  t_i_file <- glue('{patient}_10pi/t{i}/tree.nwk')
  t_i_1_file <-  glue('{patient}_10pi/t{i_1}/tree.nwk')
  t_i <- read.newick(t_i_file)
  t_i_1 <- read.newick(t_i_1_file)
  
  #Robinson Foulds 
  rf_dis <-RobinsonFoulds(
    tree1 = t_i,
    tree2 = t_i_1,
    similarity = FALSE,
    normalize = FALSE,)
  RF<-append(RF,rf_dis)
  
  #Nye similarity
  nye_dis <-NyeSimilarity(
    tree1 = t_i,
    tree2 = t_i_1,
    similarity = TRUE,
    normalize = FALSE,)
  NYE<-append(NYE,nye_dis)
  
  #Matching Distance
  ms_dis <-MatchingSplitDistance(
    tree1 = t_i,
    tree2 = t_i_1,
    normalize = FALSE,
    reportMatching = FALSE)
  MS<-append(MS, ms_dis)
  
  # Quartet Distance
  qd_dis <-QuartetDistance(
    t_i_1_file,t_i_file)
  QD<-append(QD,qd_dis)
  
}
#store distance results
 tree_dis_all_patients[[patient]] <- list(rf=RF,nye=NYE,ms=MS,qd = QD)
}


# write results to dataframe
for (p in (1:length(patient_list)))   {
  patient <- patient_list[p]
  l <- tree_dis_all_patients[[patient]]
  df <- data.frame(l$rf,l$nye,l$ms,l$qd)
  output_filename = glue('./output/{patient}_tree_dist.tsv')
  write.table(df, file=output_filename, quote=FALSE, sep='\t')
}


