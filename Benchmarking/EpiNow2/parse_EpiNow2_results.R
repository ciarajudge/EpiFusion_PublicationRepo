library(EpiNow2)
library(magrittr)
library(tidyverse)
library(dplyr)

parseoutput <- function(epinow2object) {
  estimates <- epinow2object$estimates
  summarised <- estimates$summarised
  EpiRs <- summarised$mean
  EpiRs <- EpiRs[summarised$variable=="R"]
  EpiRLow90 <- summarised$lower_90[summarised$variable=="R"] 
  EpiRUp90 <- summarised$upper_90[summarised$variable=="R"]
  dates <- summarised$date[summarised$variable=="R"]
  days <- as.numeric(dates - as.Date("01-01-2021", format = "%d-%m-%Y"))
  summarytable <- data.frame(Time = days, Mean = EpiRs, Lower90 = EpiRLow90,
                             Upper90 = EpiRUp90)
  tab <- epinow2object$estimates$samples %>%
    dplyr::filter(variable== "R") %>%
    dplyr::mutate(Time = as.numeric(date - as.Date("01-01-2021", format = "%d-%m-%Y"))) %>%
    dplyr::select(Time, sample, value) %>%
    tidyr::pivot_wider(names_from = sample, values_from = value)
  transformedtab <- t(tab)
  hpd95 <- HDInterval::hdi(transformedtab, 0.95)
  summarytable$Lower95 <- unlist(hpd95[1,])
  summarytable$Upper95 <- unlist(hpd95[2,])
  summarytable <- summarytable %>%
    dplyr::select(Time, Mean, Lower95, Upper95)
  return(list(summary = summarytable, samples = tab))
}



baseline <- readRDS("Benchmarking/EpiNow2/baseline_epinow2.RDS")
parsedbaseline <- parseoutput(baseline)
saveRDS(parsedbaseline, "Benchmarking/EpiNow2/parsed_baseline_epinow2.RDS")

sampling <- readRDS("Benchmarking/EpiNow2/sampling_epinow2.RDS")
parsedsampling <- parseoutput(sampling)
saveRDS(parsedsampling, "Benchmarking/EpiNow2/parsed_sampling_epinow2.RDS")

transmission <- readRDS("Benchmarking/EpiNow2/transmission_epinow2.RDS")
parsedtransmission <- parseoutput(transmission)
saveRDS(parsedtransmission, "Benchmarking/EpiNow2/parsed_transmission_epinow2.RDS")

baseline <- readRDS("Benchmarking/EpiNow2/parsed_baseline_epinow2.RDS")$summary
baseline_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/baseline/baseline_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()

baseline_epinow2_results <- full_join(baseline, baseline_truth) %>%
  mutate(Scenario = "Baseline")


sampling <- readRDS("Benchmarking/EpiNow2/parsed_sampling_epinow2.RDS")$summary
sampling_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/sampling/sampling_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()
sampling_epinow2_results <- full_join(sampling, sampling_truth) %>%
  mutate(Scenario = "Sampling")




transmission <- readRDS("Benchmarking/EpiNow2/parsed_transmission_epinow2.RDS")$summary

transmission_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/transmission/transmission_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()


transmission_epinow2_results <- full_join(transmission, transmission_truth) %>%
  mutate(Scenario = "Transmission")


full_epinow2_results <- rbind(baseline_epinow2_results, sampling_epinow2_results, transmission_epinow2_results)

write.csv(full_epinow2_results, "Benchmarking/EpiNow2/full_epinow2_results.csv", row.names = F)







baseline <- readRDS("Benchmarking/EpiNow2/parsed_baseline_epinow2.RDS")$samples
baseline_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/baseline/baseline_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()

baseline_epinow2_results <- full_join(baseline, baseline_truth) %>%
  mutate(Scenario = "Baseline")


sampling <- readRDS("Benchmarking/EpiNow2/parsed_sampling_epinow2.RDS")$summary
sampling_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/sampling/sampling_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()
sampling_epinow2_results <- full_join(sampling, sampling_truth) %>%
  mutate(Scenario = "Sampling")




transmission <- readRDS("Benchmarking/EpiNow2/parsed_transmission_epinow2.RDS")$summary

transmission_truth <- read.csv("Scenario_Testing/Data_Simulation/Main_Scenarios/transmission/transmission_truth.csv") %>%
  mutate(Time = ceiling(t)) %>%
  dplyr::select(Time, Rt) %>%
  group_by(Time) %>%
  mutate(Rt = mean(Rt)) %>%
  distinct()


transmission_epinow2_results <- full_join(transmission, transmission_truth) %>%
  mutate(Scenario = "Transmission")




