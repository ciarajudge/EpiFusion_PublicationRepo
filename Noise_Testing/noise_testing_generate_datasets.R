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
  write.csv(table, txtfile, row.names = F)
  table <- table %>%
    dplyr::select(t, sample) %>%
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


gettruert <- function(truetable, remasterxml) {
  truth <- read.csv(truetable)
  beta <- betaovertime(remasterxml)
  gamma <- getgamma(remasterxml, nrow(truth))
  
  if (length(beta) != 1) {
    time <- 0
    betavec <- c()
    for (i in 1:length(beta$changeTimes)) {
      betavec <- append(betavec, rep(beta$rate[i], length(truth$t[truth$t > time & truth$t <= beta$changeTimes[i]])))
      time <- beta$changeTimes[i]
    }
    l <- length(betavec)
    betavec <- append(betavec, rep(beta$rate[i+1], (nrow(truth)-l)))
  } else {
    betavec <- rep(beta, nrow(truth))
  }
  
  truth <- truth %>%
    mutate(Gamma = gamma, Beta = betavec) %>%
    mutate(Rt = (Beta*(S*I))/(Gamma*I))
  
  return(truth)
  
}


noise_intervals <- seq(0, 1, 0.1)

##### Transmission Noise #####
for (i in 1:length(noise_intervals)) {
  label <- paste0("transmission-noise_", noise_intervals[i])
  dir.create(paste0("Noise_Testing/", label), showWarnings = F)
  
  intervals <- rpois(100, 6)
  intervals <- intervals[intervals!= 0]
  changetimes <- cumsum(intervals)
  changetimes <- changetimes[changetimes<200]
  chtimes <- paste0(changetimes, collapse = " ")
  
  low_noise <- abs(round(rnorm(length(changetimes)+1, mean = 0.000045, sd = 0.000045*noise_intervals[i]),8))
  rates <- paste0(low_noise, collapse = " ")
  
  #Read template xml and adjust it for drawn parameters
  xml_content <- read_xml("Noise_Testing/transmission-noise-remaster-template.xml")
  reaction_node <- xml_find_all(xml_content, "//reaction[@rate='0.000029']")
  xml_attr(reaction_node, "rate") <- rates
  xml_attr(reaction_node, "changeTimes") <- chtimes
  write_xml(xml_content, paste0("Noise_Testing/",label,"/simulate.xml"))

  filestem <- paste0("Noise_Testing/",label,"/simulate")
  system2("/Applications/BEAST\ 2.7.3/bin/beast", args = c( "-working", '-overwrite', paste0("Noise_Testing/",label,"/simulate.xml")))
   
  system(paste0('sed -e "s/\\[\\&type=\\"I\\",time=/_/g" ',filestem,'.trees > ',filestem,'.tree'))
  system(paste0('sed -e "s/;/\\n/g" -e "s/t=//g" -e "s/:R=/,/g" -e "s/:S=/,/g" -e "s/:sample=/,/g" -e "s/:E=/,/g" -e "s/:I=/,/g" ',filestem,'.traj > ',filestem,'placeholder.txt'))

  treefile <- readLines(paste0(filestem, ".tree"), warn = FALSE)
  treestring <- str_replace(treefile[length(treefile)-1], "tree STATE_0 = ", "")
  downsampledtree <- get_tree(treestring, 0.05)
  write.tree(downsampledtree, paste0(filestem, "_downsampledtree.tree"))
  
  system(paste0("sed -e '1d' -e '2s/^0\\t//' ",filestem,"placeholder.txt > ",filestem,".txt"))
  incidence <- get_weekly_incidence(get_incidence(paste0(filestem, ".txt")))
  write.table(incidence, paste0(filestem, "_weeklyincidence.txt"), row.names = F, quote = F, col.names = F)
  
  truth <- gettruert(paste0(filestem, ".txt"), paste0(filestem, ".xml"))
  write.csv(truth, paste0(filestem, "_truth.csv"), row.names = F)
  
  xml_content <- read_xml("Noise_Testing/EpiFusion_combined.xml")
  reaction_node <- xml_find_all(xml_content, "//fileBase")
  xml_text(reaction_node) <- paste0("Noise_Testing/Results/", label, "_combined")
  reaction_node2 <- xml_find_all(xml_content, "//incidenceVals")
  xml_text(reaction_node2) <- paste0(incidence, collapse = " ")
  reaction_node3 <- xml_find_all(xml_content, "//tree")
  xml_text(reaction_node3) <- readLines(paste0(filestem, "_downsampledtree.tree"))[1]
  write_xml(xml_content, paste0("Noise_Testing/EpiFusion_XMLs/", label,"_combined.xml"))
  
  xml_content <- read_xml("Noise_Testing/EpiFusion_epi.xml")
  reaction_node <- xml_find_all(xml_content, "//fileBase")
  xml_text(reaction_node) <- paste0("Noise_Testing/Results/", label, "_epi")
  reaction_node2 <- xml_find_all(xml_content, "//incidenceVals")
  xml_text(reaction_node2) <- paste0(incidence, collapse = " ")
  reaction_node3 <- xml_find_all(xml_content, "//tree")
  xml_text(reaction_node3) <- readLines(paste0(filestem, "_downsampledtree.tree"))[1]
  write_xml(xml_content, paste0("Noise_Testing/EpiFusion_XMLs/", label,"_epi.xml"))
  
  xml_content <- read_xml("Noise_Testing/EpiFusion_phylo.xml")
  reaction_node <- xml_find_all(xml_content, "//fileBase")
  xml_text(reaction_node) <- paste0("Noise_Testing/Results/", label, "_phylo")
  reaction_node2 <- xml_find_all(xml_content, "//incidenceVals")
  xml_text(reaction_node2) <- paste0(incidence, collapse = " ")
  reaction_node3 <- xml_find_all(xml_content, "//tree")
  xml_text(reaction_node3) <- readLines(paste0(filestem, "_downsampledtree.tree"))[1]
  write_xml(xml_content, paste0("Noise_Testing/EpiFusion_XMLs/", label,"_phylo.xml"))
}



##### Observation Noise #####
for (i in 2:length(noise_intervals)) {
  label <- paste0("observation-noise_", noise_intervals[i])
  dir.create(paste0("Noise_Testing/", label), showWarnings = F)
  
  intervals <- rpois(100, 6)
  intervals <- intervals[intervals!= 0]
  changetimes <- cumsum(intervals)
  changetimes <- changetimes[changetimes<200]
  chtimes <- paste0(changetimes, collapse = " ")
  
  low_noise <- abs(round(rnorm(length(changetimes)+1, mean = 0.02, sd = 0.02*noise_intervals[i]),4))
  rates <- paste0(low_noise, collapse = " ")
  
  #Read template xml and adjust it for drawn parameters
  xml_content <- read_xml("Noise_Testing/observation-noise-remaster-template.xml")
  reaction_node <- xml_find_all(xml_content, "//reaction[@rate='0.02']")
  xml_attr(reaction_node, "rate") <- rates
  xml_attr(reaction_node, "changeTimes") <- chtimes
  write_xml(xml_content, paste0("Noise_Testing/",label,"/simulate.xml"))
  
  filestem <- paste0("Noise_Testing/",label,"/simulate")
  system2("/Applications/BEAST\ 2.7.3/bin/beast", args = c( "-working", '-overwrite', paste0("Noise_Testing/",label,"/simulate.xml")))
  
  system(paste0('sed -e "s/\\[\\&type=\\"I\\",time=/_/g" ',filestem,'.trees > ',filestem,'.tree'))
  system(paste0('sed -e "s/;/\\n/g" -e "s/t=//g" -e "s/:R=/,/g" -e "s/:S=/,/g" -e "s/:sample=/,/g" -e "s/:E=/,/g" -e "s/:I=/,/g" ',filestem,'.traj > ',filestem,'placeholder.txt'))
  
  treefile <- readLines(paste0(filestem, ".tree"), warn = FALSE)
  treestring <- str_replace(treefile[length(treefile)-1], "tree STATE_0 = ", "")
  downsampledtree <- get_tree(treestring, 0.05)
  write.tree(downsampledtree, paste0(filestem, "_downsampledtree.tree"))
  
  system(paste0("sed -e '1d' -e '2s/^0\\t//' ",filestem,"placeholder.txt > ",filestem,".txt"))
  incidence <- get_weekly_incidence(get_incidence(paste0(filestem, ".txt")))
  write.table(incidence, paste0(filestem, "_weeklyincidence.txt"), row.names = F, quote = F, col.names = F)
  
  truth <- gettruert(paste0(filestem, ".txt"), paste0(filestem, ".xml"))
  write.csv(truth, paste0(filestem, "_truth.csv"), row.names = F)
  
  xml_content <- read_xml("Noise_Testing/EpiFusion_combined.xml")
  reaction_node <- xml_find_all(xml_content, "//fileBase")
  xml_text(reaction_node) <- paste0("Noise_Testing/Results/", label, "_combined")
  reaction_node2 <- xml_find_all(xml_content, "//incidenceVals")
  xml_text(reaction_node2) <- paste0(incidence, collapse = " ")
  reaction_node3 <- xml_find_all(xml_content, "//tree")
  xml_text(reaction_node3) <- readLines(paste0(filestem, "_downsampledtree.tree"))[1]
  write_xml(xml_content, paste0("Noise_Testing/EpiFusion_XMLs/", label,"_combined.xml"))
  
  xml_content <- read_xml("Noise_Testing/EpiFusion_epi.xml")
  reaction_node <- xml_find_all(xml_content, "//fileBase")
  xml_text(reaction_node) <- paste0("Noise_Testing/Results/", label, "_epi")
  reaction_node2 <- xml_find_all(xml_content, "//incidenceVals")
  xml_text(reaction_node2) <- paste0(incidence, collapse = " ")
  reaction_node3 <- xml_find_all(xml_content, "//tree")
  xml_text(reaction_node3) <- readLines(paste0(filestem, "_downsampledtree.tree"))[1]
  write_xml(xml_content, paste0("Noise_Testing/EpiFusion_XMLs/", label,"_epi.xml"))
  
  xml_content <- read_xml("Noise_Testing/EpiFusion_phylo.xml")
  reaction_node <- xml_find_all(xml_content, "//fileBase")
  xml_text(reaction_node) <- paste0("Noise_Testing/Results/", label, "_phylo")
  reaction_node2 <- xml_find_all(xml_content, "//incidenceVals")
  xml_text(reaction_node2) <- paste0(incidence, collapse = " ")
  reaction_node3 <- xml_find_all(xml_content, "//tree")
  xml_text(reaction_node3) <- readLines(paste0(filestem, "_downsampledtree.tree"))[1]
  write_xml(xml_content, paste0("Noise_Testing/EpiFusion_XMLs/", label,"_phylo.xml"))
}












