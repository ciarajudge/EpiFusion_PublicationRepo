library(ape)
library(truncnorm)
library(xml2)
library(stringr)
library(castor)
library(magrittr)
library(tidyverse)
library(dplyr)
library(TreeTools)
library(EpiFusionUtilities)
library(HDInterval)

##### THIS SCRIPT AUTOMATES THE CREATION OF MANY REMASTER XMLS, #####
##### RUNS THEM IN BEAST WITH A SUBPROCESS CALL, PARSES THE DATA, #####
##### MAKES EPIFUSION XMLs, THEN RUNS THEM IN EPIFUSION #####


get_tree <- function(treestring, downsamplefactor) {
  tree <- read_tree(treestring)
  ntips <- length(tree$tip.label)
  labels <- tree$tip.label
  for ( i in 1:length(labels)) {
    tmp <- unlist(str_split(str_remove(labels[i], "]"), "_"))
    labels[i] <- paste0("leaf_", tmp[1], "[", tmp[2], "]")
  }
  tree$tip.label <- labels
  
  nodelabels <- tree$node.label
  for (i in 1:length(nodelabels)) {
    tmp <- str_remove(str_remove(nodelabels[i], "_"), "]")
    nodelabels[i] <- paste0("node_", i, "[", tmp, "]")
  }
  tree$node.label <- nodelabels
  
  numsamples <- round(ntips*downsamplefactor)
  
  #Downsample the tree
  sampled_tips <- sample(tree$tip.label, numsamples, replace = F)
  downsampled_tree <- keep.tip(tree, sampled_tips)
  
  return(downsampled_tree)
}

get_incidence <- function(txtfile) {
  table <- read.csv(txtfile, header = F)
  colnames(table) <- c("t", "R", "S", "I", "sample")
  table <- table %>%
    select(t, sample) %>%
    distinct(sample, .keep_all = T) %>%
    filter(sample != 0) %>%
    transmute(t = floor(t)) %>%
    group_by(t) %>%
    mutate(Incidence = n()) %>%
    distinct()
  return(table)
}

get_weekly_incidence <- function(incidence) {
  incidence$t <- as.integer(as.character(incidence$t))
  incidence$Incidence[1] <- 0
  weekly_incidence <- c()
  for (t in seq(0, max(incidence$t+7), 7)) {
    dates <- seq(t,t+6)
    inc <- incidence[which(incidence$t %in% dates),2]
    weekly_incidence <- append(weekly_incidence, sum(inc))
  }
  return(weekly_incidence)
}


numreplicates <- 700

betadraws <- rtruncnorm(numreplicates, a = 0, mean = 0.000035, sd = 0.000003)
gammadraws <- rtruncnorm(numreplicates, a = 0, mean = 0.143, sd = 0.07)
phidraws <- rtruncnorm(numreplicates, a = 0, mean = 0.02, sd = 0.01)
psidraws <- rtruncnorm(numreplicates, a = 0, mean = 0.05, sd = 0.02)

draws <- data.frame(beta = betadraws, gamma = gammadraws, psi = psidraws*phidraws, phi = phidraws)
write.csv(draws, "Simulation_Based_Calibration/paramdraws.csv", row.names = F)

for (i in 1:numreplicates) {
  label <- paste0("rep_", i)
  dir.create(paste0("Simulation_Based_Calibration/", label), showWarnings = F)
  
  #Read template xml and adjust it for drawn parameters
  xmlstem <- as_list(read_xml("Simulation_Based_Calibration/remastertemplate.xml"))
  xml_content <- read_xml("Simulation_Based_Calibration/remastertemplate.xml")
  reaction_node <- xml_find_all(xml_content, "//reaction[@rate='0.000029']")
  xml_attr(reaction_node, "rate") <- as.character(betadraws[i])
  reaction_node2 <- xml_find_all(xml_content, "//reaction[@rate='0.143']")
  xml_attr(reaction_node2, "rate") <- as.character(gammadraws[i])
  reaction_node3 <- xml_find_all(xml_content, "//reaction[@rate='0.02']")
  xml_attr(reaction_node3, "rate") <- as.character(phidraws[i])
  write_xml(xml_content, paste0("Simulation_Based_Calibration/",label,"/simulate.xml"))

  filestem <- paste0("Simulation_Based_Calibration/",label,"/simulate")
  system2("/Applications/BEAST\ 2.7.3/bin/beast", args = c( "-working", '-overwrite', paste0("Simulation_Based_Calibration/",label,"/simulate.xml")))
   
  system(paste0('sed -e "s/\\[\\&type=\\"I\\",time=/_/g" ',filestem,'.trees > ',filestem,'.tree'))
  system(paste0('sed -e "s/;/\\n/g" -e "s/t=//g" -e "s/:R=/,/g" -e "s/:S=/,/g" -e "s/:sample=/,/g" -e "s/:E=/,/g" -e "s/:I=/,/g" ',filestem,'.traj > ',filestem,'placeholder.txt'))

  treefile <- readLines(paste0(filestem, ".tree"), warn = FALSE)
  treestring <- str_replace(treefile[length(treefile)-1], "tree STATE_0 = ", "")
  downsampledtree <- get_tree(treestring, psidraws[i])
  write.tree(downsampledtree, paste0(filestem, "_downsampledtree.tree"))
  
  system(paste0("sed -e '1d' -e '2s/^0\\t//' ",filestem,"placeholder.txt > ",filestem,".txt"))
  incidence <- get_weekly_incidence(get_incidence(paste0(filestem, ".txt")))
  write.table(incidence, paste0(filestem, "_weeklyincidence.txt"), row.names = F, quote = F, col.names = F)
  
  xml_content <- read_xml("Simulation_Based_Calibration/EpiFusion_template.xml")
  reaction_node <- xml_find_all(xml_content, "//fileBase")
  xml_text(reaction_node) <- paste0(filestem)
  reaction_node2 <- xml_find_all(xml_content, "//incidenceVals")
  xml_text(reaction_node2) <- paste0(incidence, collapse = " ")
  reaction_node3 <- xml_find_all(xml_content, "//tree")
  xml_text(reaction_node3) <- readLines(paste0(filestem, "_downsampledtree.tree"))[1]
  write_xml(xml_content, paste0("Simulation_Based_Calibration/",label,"/epifusion.xml"))
}

for (i in 1:numreplicates) {
  label <- paste0("rep_", i)
  filestem <- paste0("Simulation_Based_Calibration/",label,"/")
  epifusionfile <- paste0(filestem,"epifusion.xml")
  system(paste0("java -jar EpiFusion.jar ",epifusionfile))
}



betaovertime <- function(xml_file) {
  xml_doc <- read_xml(xml_file)
  reactions <- xml_find_all(xml_doc, "//reaction")
  reaction <- reactions[[1]]
  rate <- xml_attr(reaction, "rate")
  changeTimes <- xml_attr(reaction, "changeTimes")
  if (is.na(changeTimes)) {
    return(as.numeric(rate))
  } else {
    changeTimes <- as.numeric(unlist(str_split(changeTimes, " ")))
    rate <- as.numeric(unlist(str_split(rate, " ")))
    return(list(rate = rate, changeTimes = changeTimes))
  }
}

getgamma <- function(xml_file, length) {
  xml_doc <- read_xml(xml_file)
  reactions <- xml_find_all(xml_doc, "//reaction")
  reaction <- reactions[[2]]
  rate <- xml_attr(reaction, "rate")
  changeTimes <- xml_attr(reaction, "changeTimes")
  return(rep(as.numeric(rate), length))
}




library(progress)

pb <- progress_bar$new(total = 700)
for (i in 1:numreplicates) {
  pb$tick()
  label <- paste0("rep_", i)
  filestem <- paste0("Simulation_Based_Calibration/",label,"/")
  trajectory <- read.csv(paste0(filestem,"simulate.txt"))
  colnames(trajectory) <- c("t", "R", "S", "I", "psi")
  trajectoryspare <- trajectory %>%
    mutate(t = ceiling(t)) %>%
    group_by(t) %>%
    dplyr::select(t, S, I) %>%
    mutate(S = mean(S), I = mean(I)) %>%
    distinct()
  beta <- betaovertime(paste0(filestem,"simulate.xml"))
  gamma <- getgamma(paste0(filestem,"simulate.xml"), 1)
  
  trajectoryspare <- trajectoryspare %>%
    mutate(Gamma = gamma, Beta = beta) %>%
    mutate(Beta = ((Beta*(S*I))/(Gamma*I))*Gamma)
  write.table(trajectoryspare$Beta, paste0(filestem, "beta.txt"), row.names = F)
}










