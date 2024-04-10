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
library(phytools)

##### Helper Functions #####
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

get_incidence_BEAST273 <- function(txtfile) {
  table <- read.csv(txtfile, header = F)
  colnames(table) <- c("t", "X", "sample")
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

get_incidence_BEAST276 <- function(txtfile) {
  table <- read.table(txtfile, header = T, sep = "\t") %>%
    select(t, population, value) %>%
    pivot_wider(names_from = population, values_from = value)
  write.csv(table, str_replace(txtfile, ".traj", ".txt"), row.names = F)
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


##### Dataset Parameter Combinations #####
betarange <- seq(0.19, 0.27, length.out = 3)
gammarange <- seq(0.08, 0.16, length.out = 3)
psirange <- seq(0.01, 0.05, length.out = 3)

paramcombos <- data.frame(beta = rep(mean(betarange), 9),
                          gamma = rep(mean(gammarange), 9),
                          psi = rep(mean(psirange), 9))

paramcombos$beta[1:3] <- betarange
paramcombos$gamma[4:6] <- gammarange
paramcombos$psi[7:9] <- psirange

write.csv(paramcombos, "Likelihood_Comparison/datasimulationcombos.csv", row.names = F)

##### Make ReMASTER files and run ReMASTER #####
for (i in 1:nrow(paramcombos)) {
  label <- paste0("beta", paramcombos[i,1], "_gamma", paramcombos[i,2], "_psi", paramcombos[i,3])
  dir.create(paste0("Likelihood_Comparison/", label), showWarnings = F)

  #Read template xml and adjust it for drawn parameters
  xmlstem <- as_list(read_xml("Likelihood_Comparison/remastertemplate.xml"))
  xml_content <- read_xml("Likelihood_Comparison/remastertemplate.xml")
  reaction_node <- xml_find_all(xml_content, "//reaction[@rate='0.25']")
  xml_attr(reaction_node, "rate") <- as.character(paramcombos[i, 1])
  reaction_node2 <- xml_find_all(xml_content, "//reaction[@rate='0.143']")
  xml_attr(reaction_node2, "rate") <- as.character(paramcombos[i, 2])
  reaction_node3 <- xml_find_all(xml_content, "//reaction[@rate='0.02']")
  xml_attr(reaction_node3, "rate") <- as.character(paramcombos[i, 3])
  write_xml(xml_content, paste0("Likelihood_Comparison/",label,"/simulate.xml"))

  #Run ReMASTER simulation
  filestem <- paste0("Likelihood_Comparison/",label,"/simulate")
  system2("/Applications/BEAST\ 2.7.6/bin/beast", args = c( "-working", '-overwrite', paste0("Likelihood_Comparison/",label,"/simulate.xml")))
   
   #Parse ReMASTER Output
  system(paste0('sed -e "s/\\[\\&type=\\"X\\",time=/_/g" ',filestem,'.trees > ',filestem,'.tree'))
  #system(paste0('sed -e "s/;/\\n/g" -e "s/t=//g" -e "s/:X=/,/g" -e "s/:S=/,/g" -e "s/:sample=/,/g" -e "s/:E=/,/g" -e "s/:I=/,/g" ',filestem,'.traj > ',filestem,'placeholder.txt'))
   
  treefile <- readLines(paste0(filestem, ".tree"), warn = FALSE)
  treestring <- str_replace(treefile[length(treefile)-1], "tree STATE_0 = ", "")
  downsampledtree <- get_tree(treestring, 1)
  write.tree(downsampledtree, paste0(filestem, "_downsampledtree.tree"))
  
  downsampledtree$root.edge <- min(as.numeric(gsub("\\[|\\]", "", str_extract(downsampledtree$node.label, "\\[([^\\]]+)\\]"))))
  downsampledtree$tip.label <- seq(1, length(downsampledtree$tip.label))
  downsampledtree$node.label <- ""
  write.tree(downsampledtree, paste0(filestem, "_downsampledtree_densitree.tree"))
    
  #system(paste0("sed -e '1d' -e '2s/^0\\t//' ",filestem,"placeholder.txt > ",filestem,".txt"))
  incidence <- get_weekly_incidence(get_incidence_BEAST276(paste0(filestem, ".traj")))
  write.table(incidence, paste0(filestem, "_weeklyincidence.txt"), row.names = F, quote = F, col.names = F)

  #Params to try
  if (i %in% 1:3) {
    pvec <- seq((paramcombos[i, 1] - (paramcombos[i, 1]*0.75)), paramcombos[i, 1] + (paramcombos[i, 1]*0.75), length.out = 50)
  } else if (i %in% 4:6) {
    pvec <- seq((paramcombos[i, 2] - (paramcombos[i, 2]*0.75)), paramcombos[i, 2] + (paramcombos[i, 2]*0.75), length.out = 50)
  } else if (i %in% 7:9) {
    pvec <- seq((paramcombos[i, 3] - (paramcombos[i, 3]*0.75)), paramcombos[i, 3] + (paramcombos[i, 3]*0.75), length.out = 50)
  }


  #Make EpiFusion XMLs
  for (j in 1:length(pvec)) {
    params <- unlist(paramcombos[i,])
    if (i %in% 1:3) {
      params[1] <- pvec[j]
    } else if (i %in% 4:6) {
      params[2] <- pvec[j]
    } else if (i %in% 7:9) {
      params[3] <- pvec[j]
    }
    filebase <- paste0("Likelihood_Comparison/", label, "/beta", params[1], "_gamma", params[2], "_psi", params[3])
    xml_content <- read_xml("Likelihood_Comparison/EpiFusion_template.xml")
    reaction_node <- xml_find_all(xml_content, "//fileBase")
    xml_text(reaction_node) <- filebase
    reaction_node2 <- xml_find_all(xml_content, "//incidenceVals")
    xml_text(reaction_node2) <- paste0(incidence, collapse = " ")
    reaction_node3 <- xml_find_all(xml_content, "//tree")
    xml_text(reaction_node3) <- readLines(paste0(filestem, "_downsampledtree.tree"))[1]
    beta_element <- xml_find_first(xml_content, "//beta/value")
    xml_text(beta_element) <- as.character(params[1])
    gamma_element <- xml_find_first(xml_content, "//gamma/value")
    xml_text(gamma_element) <- as.character(params[2])
    psi_element <- xml_find_first(xml_content, "//psi/value")
    xml_text(psi_element) <- as.character(params[3])

    write_xml(xml_content, paste0("Likelihood_Comparison/",label,"/beta", params[1], "_gamma", params[2], "_psi", params[3],".xml"))

  }
}

##### Run the EpiFusion XMLs #####
folders <- list.dirs("Likelihood_Comparison", recursive = FALSE)

for (i in 1:length(folders)) {
  strings_list <- list.files(folders[i], include.dirs = F, recursive = F)
  epifusionxmls <- strings_list[grepl("beta.*gamma.*psi|gamma.*beta.*psi|beta.*psi.*gamma", strings_list)]
  likelihooddf <- data.frame(matrix(0, ncol = 4, nrow = 0))
  for (j in 1:length(epifusionxmls)) {
    epifusionfile <- paste0(folders[i], "/", epifusionxmls[j])
    xml_content <- read_xml(epifusionfile)
    reaction_node <- xml_find_all(xml_content, "//numSteps")
    xml_text(reaction_node) <- '1000'
    write_xml(xml_content, epifusionfile)
    system(paste0("java -jar Likelihood_Comparison/EpiFusion_SaveAllLikelihoods.jar ",epifusionfile))
    paramvals <- unlist(str_split(epifusionxmls[j], "_"))
    betaval <- as.numeric(str_remove(paramvals[1], "beta"))
    gammaval <- as.numeric(str_remove(paramvals[2], "gamma"))
    psival <- as.numeric(str_remove(str_remove(paramvals[3], "psi"), ".xml"))
    likelihoods <- read.table(paste0(folders[i], "/", str_remove(epifusionxmls[j], ".xml"), "/likelihoods_chain0.txt"))[,1]
    averagelikelihood <- mean(likelihoods[likelihoods!= -Inf])
    likelihooddf <- rbind(likelihooddf, c(betaval, gammaval, psival, averagelikelihood))
  }
  write.csv(likelihooddf, paste0(folders[i], "/likelihooddataframe.csv"), row.names = F)
}

##### Parse EpiFusion Results #####
likelihooddataframe <- data.frame(matrix(0, nrow = 0, ncol = 7))
colnames(likelihooddataframe) <- c("TrueBeta", "TrueGamma", "TruePsi", "ModelBeta", "ModelGamma", "ModelPsi", "AvgLikelihood")
folders <- list.dirs("Likelihood_Comparison", recursive = FALSE)

for (i in 1:length(folders)) {
  trueparams <- unlist(str_split(str_remove(folders[i], "Likelihood_Comparison/"), "_"))
  row <- c()
  row[1] <- as.numeric(str_remove(trueparams[1], "beta"))
  row[2] <- as.numeric(str_remove(trueparams[2], "gamma"))
  row[3] <- as.numeric(str_remove(trueparams[3], "psi"))
  inferred <- list.dirs(folders[i], recursive = FALSE)
  for (j in 1:length(inferred)) {
    modelparams <- unlist(str_split(inferred[j], "/"))[3]
    modelparams <- unlist(str_split(modelparams, "_"))
    row[4] <- as.numeric(str_remove(modelparams[1], "beta"))
    row[5] <- as.numeric(str_remove(modelparams[2], "gamma"))
    row[6] <- as.numeric(str_remove(modelparams[3], "psi"))
    likelihoods <- read.table(paste0(inferred[j], "/likelihoods_chain0.txt"))[,1]
    row[7] <- median(likelihoods[likelihoods!= -Inf])
    names(row) <- c("TrueBeta", "TrueGamma", "TruePsi", "ModelBeta", "ModelGamma", "ModelPsi", "AvgLikelihood")
    likelihooddataframe <- rbind(likelihooddataframe, row)
  }
}

colnames(likelihooddataframe) <- c("TrueBeta", "TrueGamma", "TruePsi", "ModelBeta", "ModelGamma", "ModelPsi", "AvgLikelihood")
write.csv(likelihooddataframe, "EpiFusion_Likelihoods.csv", row.names = F)


##### Make the densitymapper xmls #####
folders <- list.dirs("Likelihood_Comparison",  recursive = FALSE)
folders <- folders[ grepl("beta", folders) ]
densitymapperfiles <- c()
for (i in 1:length(folders)) {
  trueparams <- unlist(str_split(str_remove(folders[i], "Likelihood_Comparison/"), "_"))
  truebeta <- as.numeric(str_remove(trueparams[1], "beta"))
  truegamma <- as.numeric(str_remove(trueparams[2], "gamma"))
  truepsi <- as.numeric(str_remove(trueparams[3], "psi"))
  if (truebeta == 0.23 & truegamma == 0.12 & truepsi == 0.03) {
    # This is the folder where all three things vary
    treeobject <- read.tree(paste0(folders[i], "/simulate_downsampledtree.tree"))
    root <- ceiling(max(nodeHeights(treeobject))) + treeobject$root.edge
    treestring <- readLines(paste0(folders[i], "/simulate_downsampledtree_densitree.tree"))[1]
    
    lower <- truebeta - (0.75*truebeta)
    upper <- truebeta + (0.75*truebeta)
    
    xml_doc <- read_xml("Likelihood_Comparison/densitymapper_birthrate_template.xml")
    tree_template_nodes <- xml_find_all(xml_doc, "//tree")
    xml_attr(tree_template_nodes, "newick") <- treestring
    lower_nodes <- xml_find_all(xml_doc, "//realParam")
    xml_attr(lower_nodes, "lower") <- as.character(lower)
    xml_attr(lower_nodes, "upper") <- as.character(upper)
    dist_nodes <- xml_find_all(xml_doc, "//distribution")
    xml_attr(dist_nodes, "origin") <- as.character(root)
    xml_attr(dist_nodes, "deathRate") <- as.character(truegamma)
    xml_attr(dist_nodes, "samplingRate") <- as.character(truepsi)
    
    write_xml(xml_doc, paste0(folders[i],"/densitymapper_gamma", truegamma, "_psi", truepsi,".xml"))
    densitymapperfiles <- append(densitymapperfiles, paste0(folders[i],"/densitymapper_gamma", truegamma, "_psi", truepsi,".xml"))
    
    lower <- truegamma - (0.75*truegamma)
    upper <- truegamma + (0.75*truegamma)
    
    xml_doc <- read_xml("Likelihood_Comparison/densitymapper_deathrate_template.xml")
    tree_template_nodes <- xml_find_all(xml_doc, "//tree")
    xml_attr(tree_template_nodes, "newick") <- treestring
    lower_nodes <- xml_find_all(xml_doc, "//realParam")
    xml_attr(lower_nodes, "lower") <- as.character(lower)
    xml_attr(lower_nodes, "upper") <- as.character(upper)
    dist_nodes <- xml_find_all(xml_doc, "//distribution")
    xml_attr(dist_nodes, "origin") <- as.character(root)
    xml_attr(dist_nodes, "samplingRate") <- as.character(truepsi)
    xml_attr(dist_nodes, "birthRate") <- as.character(truebeta)
    
    write_xml(xml_doc, paste0(folders[i],"/densitymapper_beta", truebeta, "_psi", truepsi,".xml"))
    densitymapperfiles <- append(densitymapperfiles, paste0(folders[i],"/densitymapper_beta", truebeta, "_psi", truepsi,".xml"))
    
    lower <- truepsi - (0.75*truepsi)
    upper <- truepsi + (0.75*truepsi)
    
    xml_doc <- read_xml("Likelihood_Comparison/densitymapper_samplerate_template.xml")
    tree_template_nodes <- xml_find_all(xml_doc, "//tree")
    xml_attr(tree_template_nodes, "newick") <- treestring
    lower_nodes <- xml_find_all(xml_doc, "//realParam")
    xml_attr(lower_nodes, "lower") <- as.character(lower)
    xml_attr(lower_nodes, "upper") <- as.character(upper)
    dist_nodes <- xml_find_all(xml_doc, "//distribution")
    xml_attr(dist_nodes, "origin") <- as.character(root)
    xml_attr(dist_nodes, "deathRate") <- as.character(truegamma)
    xml_attr(dist_nodes, "birthRate") <- as.character(truebeta)
    
    write_xml(xml_doc, paste0(folders[i],"/densitymapper_beta", truebeta, "_gamma", truegamma,".xml"))
    densitymapperfiles <- append(densitymapperfiles, paste0(folders[i],"/densitymapper_beta", truebeta, "_gamma", truegamma,".xml"))
    
    
  } else if (truegamma == 0.12 & truepsi == 0.03) {
    # betavaries
    treeobject <- read.tree(paste0(folders[i], "/simulate_downsampledtree.tree"))
    root <- ceiling(max(nodeHeights(treeobject))) + treeobject$root.edge
    treestring <- readLines(paste0(folders[i], "/simulate_downsampledtree_densitree.tree"))[1]
    
    lower <- truebeta - (0.75*truebeta)
    upper <- truebeta + (0.75*truebeta)
    
    xml_doc <- read_xml("Likelihood_Comparison/densitymapper_birthrate_template.xml")
    tree_template_nodes <- xml_find_all(xml_doc, "//tree")
    xml_attr(tree_template_nodes, "newick") <- treestring
    lower_nodes <- xml_find_all(xml_doc, "//realParam")
    xml_attr(lower_nodes, "lower") <- as.character(lower)
    xml_attr(lower_nodes, "upper") <- as.character(upper)
    dist_nodes <- xml_find_all(xml_doc, "//distribution")
    xml_attr(dist_nodes, "origin") <- as.character(root)
    xml_attr(dist_nodes, "deathRate") <- as.character(truegamma)
    xml_attr(dist_nodes, "samplingRate") <- as.character(truepsi)
    
    write_xml(xml_doc, paste0(folders[i],"/densitymapper_gamma", truegamma, "_psi", truepsi,".xml"))
    densitymapperfiles <- append(densitymapperfiles, paste0(folders[i],"/densitymapper_gamma", truegamma, "_psi", truepsi,".xml"))
  } else if (truebeta == 0.23 & truepsi == 0.03) {
    # gammavaries
    treeobject <- read.tree(paste0(folders[i], "/simulate_downsampledtree.tree"))
    root <- ceiling(max(nodeHeights(treeobject))) + treeobject$root.edge
    treestring <- readLines(paste0(folders[i], "/simulate_downsampledtree_densitree.tree"))[1]
    lower <- truegamma - (0.75*truegamma)
    upper <- truegamma + (0.75*truegamma)
    
    xml_doc <- read_xml("Likelihood_Comparison/densitymapper_deathrate_template.xml")
    tree_template_nodes <- xml_find_all(xml_doc, "//tree")
    xml_attr(tree_template_nodes, "newick") <- treestring
    lower_nodes <- xml_find_all(xml_doc, "//realParam")
    xml_attr(lower_nodes, "lower") <- as.character(lower)
    xml_attr(lower_nodes, "upper") <- as.character(upper)
    dist_nodes <- xml_find_all(xml_doc, "//distribution")
    xml_attr(dist_nodes, "origin") <- as.character(root)
    xml_attr(dist_nodes, "samplingRate") <- as.character(truepsi)
    xml_attr(dist_nodes, "birthRate") <- as.character(truebeta)
    
    write_xml(xml_doc, paste0(folders[i],"/densitymapper_beta", truebeta, "_psi", truepsi,".xml"))
    densitymapperfiles <- append(densitymapperfiles, paste0(folders[i],"/densitymapper_beta", truebeta, "_psi", truepsi,".xml"))
  
  } else  {
    # psivaries
    treeobject <- read.tree(paste0(folders[i], "/simulate_downsampledtree.tree"))
    root <- ceiling(max(nodeHeights(treeobject))) + treeobject$root.edge
    treestring <- readLines(paste0(folders[i], "/simulate_downsampledtree_densitree.tree"))[1]
    lower <- truepsi - (0.75*truepsi)
    upper <- truepsi + (0.75*truepsi)
    
    xml_doc <- read_xml("Likelihood_Comparison/densitymapper_samplerate_template.xml")
    tree_template_nodes <- xml_find_all(xml_doc, "//tree")
    xml_attr(tree_template_nodes, "newick") <- treestring
    lower_nodes <- xml_find_all(xml_doc, "//realParam")
    xml_attr(lower_nodes, "lower") <- as.character(lower)
    xml_attr(lower_nodes, "upper") <- as.character(upper)
    dist_nodes <- xml_find_all(xml_doc, "//distribution")
    xml_attr(dist_nodes, "origin") <- as.character(root)
    xml_attr(dist_nodes, "deathRate") <- as.character(truegamma)
    xml_attr(dist_nodes, "birthRate") <- as.character(truebeta)
    
    write_xml(xml_doc, paste0(folders[i],"/densitymapper_beta", truebeta, "_gamma", truegamma,".xml"))
    densitymapperfiles <- append(densitymapperfiles, paste0(folders[i],"/densitymapper_beta", truebeta, "_gamma", truegamma,".xml"))
  }
}


##### Run the densitymapper XMLs #####
for (i in 1:length(densitymapperfiles)) {
  system2("/Applications/BEAST\ 2.7.6/bin/beast", args = c( "-working", '-overwrite', densitymapperfiles[i]))
}

##### Parse the DensityMapper Results #####
densitymapperresults <- paste0("Likelihood_Comparison/", list.files("Likelihood_Comparison/", pattern = ".log", recursive = T))
bdskydataframe <- data.frame(matrix(0, nrow = 0, ncol = 7))
colnames(bdskydataframe) <- c("TrueBeta", "TrueGamma", "TruePsi", "ModelBeta", "ModelGamma", "ModelPsi", "BDSkyDensity")


for (i in 1:length(densitymapperresults)) {
  output <- read.table(densitymapperresults[i], sep = "\t", header = T)
  file <- unlist(str_split(densitymapperresults[i], "/"))[3]
  trueparams <- unlist(str_split(str_remove(str_remove(file, ".log"), "densitymapper_"), "_"))
  simulationparams <- unlist(str_split(unlist(str_split(densitymapperresults[i], "/"))[2], "_"))
  if ("birthRate" %in% colnames(output)) {
    output <- output %>%
      transmute(ModelBeta = birthRate, BDSkyDensity = density) %>%
      mutate(TrueGamma = as.numeric(str_remove(trueparams[1], "gamma")),
             TruePsi = as.numeric(str_remove(trueparams[2], "psi")),
             ModelGamma = as.numeric(str_remove(trueparams[1], "gamma")),
             ModelPsi = as.numeric(str_remove(trueparams[2], "psi")),
             TrueBeta = as.numeric(str_remove(simulationparams[1], "beta")))
  } else if ("deathRate" %in% colnames(output)) {
    output <- output %>%
      transmute(ModelGamma = deathRate, BDSkyDensity = density) %>%
      mutate(TrueBeta = as.numeric(str_remove(trueparams[1], "beta")),
             TruePsi = as.numeric(str_remove(trueparams[2], "psi")),
             ModelBeta = as.numeric(str_remove(trueparams[1], "beta")),
             ModelPsi = as.numeric(str_remove(trueparams[2], "psi")),
             TrueGamma = as.numeric(str_remove(simulationparams[2], "gamma")))
  } else {
    output <- output %>%
      transmute(ModelPsi = sampleRate, BDSkyDensity = density) %>%
      mutate(TrueBeta = as.numeric(str_remove(trueparams[1], "beta")),
             TrueGamma = as.numeric(str_remove(trueparams[2], "gamma")),
             ModelBeta = as.numeric(str_remove(trueparams[1], "beta")),
             ModelGamma = as.numeric(str_remove(trueparams[2], "gamma")),
             TruePsi = as.numeric(str_remove(simulationparams[3], "psi")))
  }
  output <- select(output, TrueBeta, TrueGamma, TruePsi, ModelBeta, ModelGamma, ModelPsi, BDSkyDensity)
  bdskydataframe <- rbind(bdskydataframe, output)
}

write.csv(bdskydataframe, "Likelihood_Comparison/BDSKY_Likelihoods.csv", row.names = F)


for (i in 1:ncol(likelihooddataframe)) {
  likelihooddataframe[,i] <- round(as.numeric(likelihooddataframe[,i]), 4)
}
for (i in 1:ncol(bdskydataframe)) {
  bdskydataframe[,i] <- round(as.numeric(bdskydataframe[,i]), 4)
}

phyloonlycomparison <- full_join(bdskydataframe, likelihooddataframe)

write.csv(phyloonlycomparison, "Likelihood_Comparison/EpiFusionvsBDSky_Likelihoods.csv", row.names = F)


