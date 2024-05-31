library(xml2)
library(phytools)



gethpd <- function(str, ind) {
  val <- as.numeric(str_remove(str_remove(unlist(str_split(str, ","))[ind], "\\["), "\\]"))
  return(val)
}

gettimtamchangetimes <- function(bdskyxml) {
  xml_doc <- read_xml(bdskyxml)
  tree <- xml_attr(xml_find_first(xml_doc, "//tree"), "newick")
  write.table(tree, str_replace(bdskyxml, ".xml", ".tree"), col.names = F, row.names = F, quote = F)
  tree <- read.tree(str_replace(bdskyxml, ".xml", ".tree"))
  lastleaf <- max(nodeHeights(tree))
  print(lastleaf)
  interval_times <- as.numeric(unlist(str_split(xml_text(xml_find_first(xml_doc, "//parameter[@name='r0ChangeTimes']")), " ")))
  changetimes <- lastleaf - interval_times
  return(c(round(changetimes), round(lastleaf)))
}

gethpd <- function(str, ind) {
  val <- as.numeric(str_remove(str_remove(unlist(str_split(str, ","))[ind], "\\["), "\\]"))
  return(val)
}



##### BASELINE #####

baseline_timtam_fulltable <- read.table("Benchmarking/TimTam/baseline_timtam.tab", sep = "\t", header = T)
baseline_timtam_means <- as.numeric(baseline_timtam_fulltable[1,2:(ncol(baseline_timtam_fulltable))])
baseline_timtam_hpds <- baseline_timtam_fulltable[8,2:(ncol(baseline_timtam_fulltable))]
baseline_timtam_upperhpds <- sapply(baseline_timtam_hpds, gethpd, 2)
baseline_timtam_lowerhpds <- sapply(baseline_timtam_hpds, gethpd, 1)
baseline_timtam_changetimes <- gettimtamchangetimes("Benchmarking/timtam/baseline_timtam_with_prevalence_evenR0.xml")

baseline_timtam_results <- data.frame(matrix(0, nrow = 0, ncol = 4))

for (i in 1:(length(baseline_timtam_changetimes)-1)) {
  tmptable <- data.frame(Time = seq(baseline_timtam_changetimes[i], (baseline_timtam_changetimes[i+1]-1)),
                         Mean = rep(baseline_timtam_means[i], length(seq(baseline_timtam_changetimes[i], (baseline_timtam_changetimes[i+1]-1)))),
                         Lower95 = rep(baseline_timtam_lowerhpds[i], length(seq(baseline_timtam_changetimes[i], (baseline_timtam_changetimes[i+1]-1)))),
                         Upper95 = rep(baseline_timtam_upperhpds[i], length(seq(baseline_timtam_changetimes[i], (baseline_timtam_changetimes[i+1]-1)))))
  baseline_timtam_results <- rbind(baseline_timtam_results, tmptable)
}

baseline_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/baseline/baseline_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()


baseline_results <- full_join(baseline_timtam_results, baseline_truth)
write.csv(baseline_results, "Benchmarking/timtam/baseline_timtam_fullresults.csv", row.names = F)



baseline_timtam_samples <- read.table("Benchmarking/timtam/baseline_timtam_samples.tab", sep = "\t", header = T)
baseline_timtam_samples <- baseline_timtam_samples %>%
  dplyr::select(!starts_with("state"))

baseline_timtam_changetimes <- gettimtamchangetimes("Benchmarking/timtam/baseline_timtam_with_prevalence_evenR0.xml")

baseline_timtam_samples_master <- data.frame(matrix(0, nrow = 0, ncol = nrow(baseline_timtam_samples)+1))

for (i in 1:(length(baseline_timtam_changetimes)-1)) {
  for (j in seq(baseline_timtam_changetimes[i], (baseline_timtam_changetimes[i+1]-1))) {
    baseline_timtam_samples_master <- rbind(baseline_timtam_samples_master, c(j, baseline_timtam_samples[,i]))
  }
}

baseline_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/baseline/baseline_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()

baseline_timtam_samples_master <- baseline_timtam_samples_master %>%
  rename(Time = X3) %>%
  left_join(baseline_truth) %>%
  na.omit()

write.csv(baseline_timtam_samples_master, "Benchmarking/timtam/baseline_timtam_samples.csv", row.names = F)



##### SAMPLING #####

sampling_timtam_fulltable <- read.table("Benchmarking/TimTam/sampling_timtam.tab", sep = "\t", header = T)
sampling_timtam_means <- as.numeric(sampling_timtam_fulltable[1,2:(ncol(sampling_timtam_fulltable))])
sampling_timtam_hpds <- sampling_timtam_fulltable[8,2:(ncol(sampling_timtam_fulltable))]
sampling_timtam_upperhpds <- sapply(sampling_timtam_hpds, gethpd, 2)
sampling_timtam_lowerhpds <- sapply(sampling_timtam_hpds, gethpd, 1)
sampling_timtam_changetimes <- gettimtamchangetimes("Benchmarking/timtam/sampling_timtam.xml")

sampling_timtam_results <- data.frame(matrix(0, nrow = 0, ncol = 4))

for (i in 1:(length(sampling_timtam_changetimes)-1)) {
  tmptable <- data.frame(Time = seq(sampling_timtam_changetimes[i], (sampling_timtam_changetimes[i+1]-1)),
                         Mean = rep(sampling_timtam_means[i], length(seq(sampling_timtam_changetimes[i], (sampling_timtam_changetimes[i+1]-1)))),
                         Lower95 = rep(sampling_timtam_lowerhpds[i], length(seq(sampling_timtam_changetimes[i], (sampling_timtam_changetimes[i+1]-1)))),
                         Upper95 = rep(sampling_timtam_upperhpds[i], length(seq(sampling_timtam_changetimes[i], (sampling_timtam_changetimes[i+1]-1)))))
  sampling_timtam_results <- rbind(sampling_timtam_results, tmptable)
}

sampling_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/sampling/sampling_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()


sampling_results <- full_join(sampling_timtam_results, sampling_truth)
write.csv(sampling_results, "Benchmarking/timtam/sampling_timtam_fullresults.csv", row.names = F)



sampling_timtam_samples <- read.table("Benchmarking/timtam/sampling_timtam_samples.tab", sep = "\t", header = T)
sampling_timtam_samples <- sampling_timtam_samples %>%
  dplyr::select(!starts_with("state"))

sampling_timtam_changetimes <- gettimtamchangetimes("Benchmarking/timtam/sampling_timtam.xml")

sampling_timtam_samples_master <- data.frame(matrix(0, nrow = 0, ncol = nrow(sampling_timtam_samples)+1))

for (i in 1:(length(sampling_timtam_changetimes)-1)) {
  for (j in seq(sampling_timtam_changetimes[i], (sampling_timtam_changetimes[i+1]-1))) {
    sampling_timtam_samples_master <- rbind(sampling_timtam_samples_master, c(j, sampling_timtam_samples[,i]))
  }
}

sampling_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/sampling/sampling_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()

sampling_timtam_samples_master <- sampling_timtam_samples_master %>%
  rename(Time = X3) %>%
  left_join(sampling_truth) %>%
  na.omit()

write.csv(sampling_timtam_samples_master, "Benchmarking/timtam/sampling_timtam_samples.csv", row.names = F)




##### TRANSMISSION #####

transmission_timtam_fulltable <- read.table("Benchmarking/TimTam/transmission_timtam.tab", sep = "\t", header = T)
transmission_timtam_means <- as.numeric(transmission_timtam_fulltable[1,2:(ncol(transmission_timtam_fulltable))])
transmission_timtam_hpds <- transmission_timtam_fulltable[8,2:(ncol(transmission_timtam_fulltable))]
transmission_timtam_upperhpds <- sapply(transmission_timtam_hpds, gethpd, 2)
transmission_timtam_lowerhpds <- sapply(transmission_timtam_hpds, gethpd, 1)
transmission_timtam_changetimes <- gettimtamchangetimes("Benchmarking/timtam/transmission_timtam.xml")

transmission_timtam_results <- data.frame(matrix(0, nrow = 0, ncol = 4))

for (i in 1:(length(transmission_timtam_changetimes)-1)) {
  tmptable <- data.frame(Time = seq(transmission_timtam_changetimes[i], (transmission_timtam_changetimes[i+1]-1)),
                         Mean = rep(transmission_timtam_means[i], length(seq(transmission_timtam_changetimes[i], (transmission_timtam_changetimes[i+1]-1)))),
                         Lower95 = rep(transmission_timtam_lowerhpds[i], length(seq(transmission_timtam_changetimes[i], (transmission_timtam_changetimes[i+1]-1)))),
                         Upper95 = rep(transmission_timtam_upperhpds[i], length(seq(transmission_timtam_changetimes[i], (transmission_timtam_changetimes[i+1]-1)))))
  transmission_timtam_results <- rbind(transmission_timtam_results, tmptable)
}

transmission_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/transmission/transmission_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()


transmission_results <- full_join(transmission_timtam_results, transmission_truth)
write.csv(transmission_results, "Benchmarking/timtam/transmission_timtam_fullresults.csv", row.names = F)



transmission_timtam_samples <- read.table("Benchmarking/timtam/transmission_timtam_samples.tab", sep = "\t", header = T)
transmission_timtam_samples <- transmission_timtam_samples %>%
  dplyr::select(!starts_with("state"))

transmission_timtam_changetimes <- gettimtamchangetimes("Benchmarking/timtam/transmission_timtam.xml")

transmission_timtam_samples_master <- data.frame(matrix(0, nrow = 0, ncol = nrow(transmission_timtam_samples)+1))

for (i in 1:(length(transmission_timtam_changetimes)-1)) {
  for (j in seq(transmission_timtam_changetimes[i], (transmission_timtam_changetimes[i+1]-1))) {
    transmission_timtam_samples_master <- rbind(transmission_timtam_samples_master, c(j, transmission_timtam_samples[,i]))
  }
}

transmission_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/transmission/transmission_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()

transmission_timtam_samples_master <- transmission_timtam_samples_master %>%
  rename(Time = X4) %>%
  left_join(transmission_truth) %>%
  na.omit()

write.csv(transmission_timtam_samples_master, "Benchmarking/timtam/transmission_timtam_samples.csv", row.names = F)



baseline_results <- mutate(baseline_results, Scenario = "Baseline")
sampling_results <- mutate(sampling_results, Scenario = "Sampling")
transmission_results <- mutate(transmission_results, Scenario = "Transmission")


full_timtam_results <- rbind(baseline_results, sampling_results, transmission_results)

write.csv(full_timtam_results, "Benchmarking/timtam/full_timtam_results.csv", row.names = F)



