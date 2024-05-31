library(xml2)
library(phytools)



gethpd <- function(str, ind) {
  val <- as.numeric(str_remove(str_remove(unlist(str_split(str, ","))[ind], "\\["), "\\]"))
  return(val)
}

getbdskychangetimes <- function(bdskyxml) {
  xml_doc <- read_xml(bdskyxml)
  tree <- xml_attr(xml_find_first(xml_doc, "//tree"), "newick")
  write.table(tree, str_replace(bdskyxml, ".xml", ".tree"), col.names = F, row.names = F, quote = F)
  tree <- read.tree(str_replace(bdskyxml, ".xml", ".tree"))
  lastleaf <- max(nodeHeights(tree))
  interval_times <- as.numeric(unlist(str_split(xml_text(xml_find_first(xml_doc, "//parameter[@name='intervalTimes']")), " ")))
  changetimes <- lastleaf - interval_times
  return(c(0, round(changetimes)))
}

gethpd <- function(str, ind) {
  val <- as.numeric(str_remove(str_remove(unlist(str_split(str, ","))[ind], "\\["), "\\]"))
  return(val)
}





baseline_bdsky_fulltable <- read.table("Benchmarking/BDSky/baseline_BDSky.tab", sep = "\t", header = T)
baseline_bdsky_means <- as.numeric(baseline_bdsky_fulltable[1,2:(ncol(baseline_bdsky_fulltable))])
baseline_bdsky_hpds <- baseline_bdsky_fulltable[8,2:(ncol(baseline_bdsky_fulltable))]
baseline_bdsky_upperhpds <- sapply(baseline_bdsky_hpds, gethpd, 2)
baseline_bdsky_lowerhpds <- sapply(baseline_bdsky_hpds, gethpd, 1)
baseline_bdsky_changetimes <- getbdskychangetimes("Benchmarking/BDSky/baseline_BDSky_fixedtree_copy.xml")

baseline_bdsky_results <- data.frame(matrix(0, nrow = 0, ncol = 4))

for (i in 1:(length(baseline_bdsky_changetimes)-1)) {
  tmptable <- data.frame(Time = seq(baseline_bdsky_changetimes[i], (baseline_bdsky_changetimes[i+1]-1)),
                         Mean = rep(baseline_bdsky_means[i], length(seq(baseline_bdsky_changetimes[i], (baseline_bdsky_changetimes[i+1]-1)))),
                         Lower95 = rep(baseline_bdsky_lowerhpds[i], length(seq(baseline_bdsky_changetimes[i], (baseline_bdsky_changetimes[i+1]-1)))),
                         Upper95 = rep(baseline_bdsky_upperhpds[i], length(seq(baseline_bdsky_changetimes[i], (baseline_bdsky_changetimes[i+1]-1)))))
  baseline_bdsky_results <- rbind(baseline_bdsky_results, tmptable)
}

baseline_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/baseline/baseline_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()


baseline_results <- full_join(baseline_bdsky_results, baseline_truth)
write.csv(baseline_results, "Benchmarking/BDSky/baseline_BDSKY_fullresults.csv", row.names = F)

sampling_bdsky_fulltable <- read.table("Benchmarking/BDSky/sampling_BDSky.tab", sep = "\t", header = T)
sampling_bdsky_means <- as.numeric(sampling_bdsky_fulltable[1,2:(ncol(sampling_bdsky_fulltable))])
sampling_bdsky_hpds <- sampling_bdsky_fulltable[8,2:(ncol(sampling_bdsky_fulltable))]
sampling_bdsky_upperhpds <- sapply(sampling_bdsky_hpds, gethpd, 2)
sampling_bdsky_lowerhpds <- sapply(sampling_bdsky_hpds, gethpd, 1)
sampling_bdsky_changetimes <- getbdskychangetimes("Benchmarking/BDSky/sampling_BDSky_fixedtree.xml")

sampling_bdsky_results <- data.frame(matrix(0, nrow = 0, ncol = 4))

for (i in 1:(length(sampling_bdsky_changetimes)-1)) {
  tmptable <- data.frame(Time = seq(sampling_bdsky_changetimes[i], (sampling_bdsky_changetimes[i+1]-1)),
                         Mean = rep(sampling_bdsky_means[i], length(seq(sampling_bdsky_changetimes[i], (sampling_bdsky_changetimes[i+1]-1)))),
                         Lower95 = rep(sampling_bdsky_lowerhpds[i], length(seq(sampling_bdsky_changetimes[i], (sampling_bdsky_changetimes[i+1]-1)))),
                         Upper95 = rep(sampling_bdsky_upperhpds[i], length(seq(sampling_bdsky_changetimes[i], (sampling_bdsky_changetimes[i+1]-1)))))
  sampling_bdsky_results <- rbind(sampling_bdsky_results, tmptable)
}

sampling_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/sampling/sampling_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()


sampling_results <- full_join(sampling_bdsky_results, sampling_truth) %>%
  filter(Time > 0)

write.csv(sampling_results, "Benchmarking/BDSky/sampling_BDSKY_fullresults.csv", row.names = F)





transmission_bdsky_fulltable <- read.table("Benchmarking/BDSky/transmission_BDSky.tab", sep = "\t", header = T)
transmission_bdsky_means <- as.numeric(transmission_bdsky_fulltable[1,2:(ncol(transmission_bdsky_fulltable))])
transmission_bdsky_hpds <- transmission_bdsky_fulltable[8,2:(ncol(transmission_bdsky_fulltable))]
transmission_bdsky_upperhpds <- sapply(transmission_bdsky_hpds, gethpd, 2)
transmission_bdsky_lowerhpds <- sapply(transmission_bdsky_hpds, gethpd, 1)
transmission_bdsky_changetimes <- getbdskychangetimes("Benchmarking/BDSky/transmission_BDSky_fixedtree.xml")

transmission_bdsky_results <- data.frame(matrix(0, nrow = 0, ncol = 4))

for (i in 1:(length(transmission_bdsky_changetimes)-1)) {
  tmptable <- data.frame(Time = seq(transmission_bdsky_changetimes[i], (transmission_bdsky_changetimes[i+1]-1)),
                         Mean = rep(transmission_bdsky_means[i], length(seq(transmission_bdsky_changetimes[i], (transmission_bdsky_changetimes[i+1]-1)))),
                         Lower95 = rep(transmission_bdsky_lowerhpds[i], length(seq(transmission_bdsky_changetimes[i], (transmission_bdsky_changetimes[i+1]-1)))),
                         Upper95 = rep(transmission_bdsky_upperhpds[i], length(seq(transmission_bdsky_changetimes[i], (transmission_bdsky_changetimes[i+1]-1)))))
  transmission_bdsky_results <- rbind(transmission_bdsky_results, tmptable)
}

transmission_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/transmission/transmission_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()


transmission_results <- full_join(transmission_bdsky_results, transmission_truth) %>%
  filter(Time > 0)

write.csv(transmission_results, "Benchmarking/BDSky/transmission_BDSKY_fullresults.csv", row.names = F)



baseline_results <- mutate(baseline_results, Scenario = "Baseline")
sampling_results <- mutate(sampling_results, Scenario = "Sampling")
transmission_results <- mutate(transmission_results, Scenario = "Transmission")


full_bdsky_results <- rbind(baseline_results, sampling_results, transmission_results)

write.csv(full_bdsky_results, "Benchmarking/BDSky/full_BDSky_results.csv", row.names = F)


baseline_BDSky_samples <- read.table("Benchmarking/BDSky/baseline_BDSky_samples.tab", sep = "\t", header = T)
baseline_BDSky_samples <- baseline_BDSky_samples[sample(seq(1, nrow(baseline_BDSky_samples)), 10000),] %>%
  dplyr::select(!starts_with("state"))

baseline_bdsky_changetimes <- getbdskychangetimes("Benchmarking/BDSky/baseline_BDSky_fixedtree_copy.xml")

baseline_bdsky_samples_master <- data.frame(matrix(0, nrow = 0, ncol = nrow(baseline_BDSky_samples)+1))

for (i in 1:(length(baseline_bdsky_changetimes)-1)) {
  for (j in seq(baseline_bdsky_changetimes[i], (baseline_bdsky_changetimes[i+1]-1))) {
    baseline_bdsky_samples_master <- rbind(baseline_bdsky_samples_master, c(j, baseline_BDSky_samples[,i]))
  }
}

baseline_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/baseline/baseline_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()

baseline_bdsky_samples_master <- baseline_bdsky_samples_master %>%
  rename(Time = X0) %>%
  left_join(baseline_truth) %>%
  na.omit()

write.csv(baseline_bdsky_samples_master, "Benchmarking/BDSky/baseline_BDSky_samples.csv", row.names = F)




sampling_BDSky_samples <- read.table("Benchmarking/BDSky/sampling_BDSky_samples.tab", sep = "\t", header = T)
sampling_BDSky_samples <- sampling_BDSky_samples[sample(seq(1, nrow(sampling_BDSky_samples)), 10000),] %>%
  dplyr::select(!starts_with("state"))

sampling_bdsky_changetimes <- getbdskychangetimes("Benchmarking/BDSky/sampling_BDSky_fixedtree.xml")

sampling_bdsky_samples_master <- data.frame(matrix(0, nrow = 0, ncol = nrow(sampling_BDSky_samples)+1))

for (i in 1:(length(sampling_bdsky_changetimes)-1)) {
  for (j in seq(sampling_bdsky_changetimes[i], (sampling_bdsky_changetimes[i+1]-1))) {
    sampling_bdsky_samples_master <- rbind(sampling_bdsky_samples_master, c(j, sampling_BDSky_samples[,i]))
  }
}

sampling_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/sampling/sampling_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()

sampling_bdsky_samples_master <- sampling_bdsky_samples_master %>%
  rename(Time = X0) %>%
  left_join(sampling_truth) %>%
  na.omit()

write.csv(sampling_bdsky_samples_master, "Benchmarking/BDSky/sampling_BDSky_samples.csv", row.names = F)





transmission_BDSky_samples <- read.table("Benchmarking/BDSky/transmission_BDSky_samples.tab", sep = "\t", header = T)
transmission_BDSky_samples <- transmission_BDSky_samples[sample(seq(1, nrow(transmission_BDSky_samples)), 10000),] %>%
  dplyr::select(!starts_with("state"))

transmission_bdsky_changetimes <- getbdskychangetimes("Benchmarking/BDSky/transmission_BDSky_fixedtree.xml")

transmission_bdsky_samples_master <- data.frame(matrix(0, nrow = 0, ncol = nrow(transmission_BDSky_samples)+1))

for (i in 1:(length(transmission_bdsky_changetimes)-1)) {
  for (j in seq(transmission_bdsky_changetimes[i], (transmission_bdsky_changetimes[i+1]-1))) {
    transmission_bdsky_samples_master <- rbind(transmission_bdsky_samples_master, c(j, transmission_BDSky_samples[,i]))
  }
}

transmission_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/transmission/transmission_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()

transmission_bdsky_samples_master <- transmission_bdsky_samples_master %>%
  rename(Time = X0) %>%
  left_join(transmission_truth) %>%
  na.omit()

write.csv(transmission_bdsky_samples_master, "Benchmarking/BDSky/transmission_BDSky_samples.csv", row.names = F)












