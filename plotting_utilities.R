library(plyr)
library(tidyverse)
library(paletteer)
library(HDInterval)
library(plotly)
library(xml2)
library(truncnorm)
library(RColorBrewer)
library(ape)
library(Metrics)
library(stableGR)
library(dplyr)
library(stringr)
library(fitdistrplus)
library(scoringRules)
library(scoringutils)
library(dplyr)
library(magrittr)
library(stringr)
library(ggh4x)
library(ggpubr)
library(sjmisc)



lshtm_theme <- function() {
  theme(
    # add border 1)
    panel.border = element_rect(colour = "#01454f", fill = NA, size = 0.5),
    # color background 2)
    panel.background = element_rect(fill = "white"),
    # modify grid 3)
    #panel.grid.major.x = element_line(colour = "steelblue", linetype = 3, size = 0.5),
    panel.grid.minor.x = element_line(colour = "aliceblue"),
    #panel.grid.major.y =  element_line(colour = "steelblue", linetype = 3, size = 0.5),
    panel.grid.minor.y = element_line(colour = "aliceblue"),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "#01454f"),
    axis.title = element_text(colour = "#01454f"),
    axis.ticks = element_line(colour = "#01454f"),
    # legend at the bottom 6)
    #legend.position = "bottom"
    strip.text.x = element_text(colour = "white"),
    strip.text.y = element_text(colour = "white"),
    strip.background = element_rect(
      color="#01454f", fill="#01454f", size=1.5, linetype="solid"
    ),
    legend.position = "bottom",
    legend.title = element_text(colour = "#01454f", face = "bold"),
    legend.text = element_text(colour = "#01454f")
  )
}


loadtruth <- function(truthdirectory) {
  filestem <- tail(unlist(str_split(truthdirectory, "/")), 2)[1]
  tree <- read.tree(paste0(truthdirectory, filestem, "_downsampledtree.tree"))
  truthcontinuous <- read.csv(paste0(truthdirectory, filestem, "_truth.csv"))
  weeklyincidence <- read.table(paste0(truthdirectory, filestem, "_weeklyincidence.txt"))[,1]
  beta <- xmlreader(paste0(truthdirectory, filestem, ".xml"), 1)
  gamma <- xmlreader(paste0(truthdirectory, filestem, ".xml"), 2)
  phi <- xmlreader(paste0(truthdirectory, filestem, ".xml"), 3)
  psi <- phi$rate*(length(tree$tip.label)/sum(weeklyincidence))
  truthdiscrete <- truthcontinuous %>%
    dplyr::select(t, I, Rt) %>%
    mutate(t = ceiling(t)) %>%
    group_by(t) %>%
    transmute(Time = t, I = mean(I), Rt = mean(Rt)) %>%
    distinct()
  
  return(list(truthcontinuous = truthcontinuous,
              truthdiscrete = truthdiscrete,
              tree = tree,
              incidence = weeklyincidence,
              beta = beta,
              gamma = gamma,
              phi = phi,
              psi = psi))

  
  
}

xmlreader <- function(xml_file, r) {
  xml_doc <- read_xml(xml_file)
  reactions <- xml_find_all(xml_doc, "//reaction")
  reaction <- reactions[[r]]
  rate <- xml_attr(reaction, "rate")
  changeTimes <- xml_attr(reaction, "changeTimes")
  if (is.na(changeTimes)) {
    return(list(rate = as.numeric(rate)))
  } else {
    changeTimes <- as.numeric(unlist(str_split(changeTimes, " ")))
    rate <- as.numeric(unlist(str_split(rate, " ")))
    return(list(rate = rate, changeTimes = changeTimes))
  }
}

plottruthrow <- function(baselinetruth, title) {
  plot(baselinetruth$truthcontinuous$t, baselinetruth$truthcontinuous$I, type = "l", ylab = "True Number Infected", xlab = "Time", lwd = 1.5)
  plot(seq(7, 7*length(baselinetruth$incidence), 7), baselinetruth$incidence, main = title, pch = 10, col = "orangered3", ylab = "Weekly Case Incidence", xlab = "Time", cex = 1.5)
  plot(baselinetruth$tree, root.edge = TRUE, edge.color = c("darkslateblue"), show.tip.label = F)
  plot(baselinetruth$truthcontinuous$t, baselinetruth$truthcontinuous$Rt, type = "l", ylab = "True Rt", xlab = "Time", lwd = 1.5)
  axis(1)
}

trajectorytable <- function(object, truth, approach, scenario) {
  time <- length(object$infection_trajectories$mean_infection_trajectory)
  table <- data.frame(
    Time = seq(1, time),
    Approach = rep(approach, time),
    Scenario = rep(scenario, time),
    Mean_Infected = object$infection_trajectories$mean_infection_trajectory,
    Lower95_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Lower,
    Upper95_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Upper,
    Lower88_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.88$Lower,
    Upper88_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.88$Upper,
    Lower66_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.66$Lower,
    Upper66_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.66$Upper,
    Mean_Rt = object$rt_trajectories$mean_rt_trajectory,
    Lower95_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.95$Lower,
    Upper95_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.95$Upper,
    Lower88_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.88$Lower,
    Upper88_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.88$Upper,
    Lower66_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.66$Lower,
    Upper66_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.66$Upper
  )
  table <- table %>%
    filter(Time < max(truth$t))
  table <- table %>%
    full_join(truth)
  table$Approach[is.na(table$Approach)] <- approach
  table$Scenario[is.na(table$Scenario)] <- scenario
  return(table)
}


epifusion_rt_crps <- function(object, truth) {
  samples <- object$rt_trajectories$rt_trajectory_samples
  samples <- as.matrix(t(as.matrix(samples)))
  nonnatruth <- truth[!is.na(truth$Rt),]$t
  nonnatruth <- nonnatruth[nonnatruth <= length(object$rt_trajectories$mean_rt_trajectory)]
  true_values <- truth$Rt[!is.na(truth$Rt)]
  samples <- na.omit(samples[nonnatruth,])
  l <- min(length(true_values), nrow(samples))
  crps_value <- scoringutils::crps_sample(true_values[1:l], predictions = samples[1:l,])
  return(crps_value)
}

epifusion_infection_crps <- function(object, truth) {
  samples <- object$infection_trajectories$infection_trajectory_samples
  samples <- as.matrix(t(as.matrix(samples)))
  true_values <- truth$I
  l <- min(length(true_values), nrow(samples))
  crps_value <- scoringutils::crps_sample(true_values[1:l], predictions = samples[1:l,])
  return(crps_value)
}

epifusion_rt_brier <- function(object, truth) {
  samples <- object$rt_trajectories$rt_trajectory_samples
  samples <- as.matrix(t(as.matrix(samples)))
  nonnatruth <- truth[!is.na(truth$Rt),]$t
  nonnatruth <- nonnatruth[nonnatruth <= length(object$rt_trajectories$mean_rt_trajectory)]
  true_values <- truth$Rt[!is.na(truth$Rt)]
  samples <- na.omit(samples[nonnatruth,])
  l <- min(length(true_values), nrow(samples))
  true_values <- true_values[1:l]
  transmission_phase <- as.numeric(true_values >= 1)
  predicted_transmission_phases <- samples[1:l,]
  predicted_transmission_phases[predicted_transmission_phases < 1] <- 0
  predicted_transmission_phases[predicted_transmission_phases >= 1] <- 1
  predictions <- rowMeans(predicted_transmission_phases)
  brier <- mean(scoringutils::brier_score(transmission_phase, predictions = predictions))
  return(brier)
}

performancemetrics <- function(object, truth, approach, scenario) {
  time <- length(object$infection_trajectories$mean_infection_trajectory)
  table <- data.frame(
    Time = seq(1, time),
    Approach = rep(approach, time),
    Scenario = rep(scenario, time),
    Mean_Infected = object$infection_trajectories$mean_infection_trajectory,
    Lower95_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Lower,
    Upper95_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Upper,
    Mean_Rt = object$rt_trajectories$mean_rt_trajectory,
    Lower95_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.95$Lower,
    Upper95_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.95$Upper
  ) 
  
  table <- table %>%
    filter(Time < max(truth$t))
  table <- table %>%
    full_join(truth) %>%
    mutate(Infection_Coverage = ifelse((I >= Lower95_Infected & I <= Upper95_Infected), 1, 0)) %>%
    mutate(Rt_Coverage = ifelse((Rt >= Lower95_Rt & Rt <= Upper95_Rt), 1, 0))
  table$Approach[is.na(table$Approach)] <- approach
  table$Scenario[is.na(table$Scenario)] <- scenario
  
  table <- na.omit(table) %>%
    mutate(Scaled_Inf_HPD_Width  = (Upper95_Infected - Lower95_Infected)/Mean_Infected) %>%
    mutate(Scaled_Rt_HPD_Width = (Upper95_Rt - Lower95_Rt)/Mean_Rt)
  inf_rmse <- round(rmse(table$I, table$Mean_Infected), 2)
  inf_coverage <- round((sum(table$Infection_Coverage)/nrow(table))/0.95, 2)
  scaled_inf_hpd_width <- round(mean(table$Scaled_Inf_HPD_Width), 2)
  rt_rmse <- round(rmse(table$Rt, table$Mean_Rt), 3)
  rt_coverage <- round((sum(table$Rt_Coverage)/nrow(table))/0.95, 2)
  scale_rt_hpd_width <- round(mean(table$Scaled_Rt_HPD_Width), 2)
  infection_crps <- round(mean(epifusion_infection_crps(object, truth)), 2)
  rt_crps <- round(mean(epifusion_rt_crps(object, truth)),3)
  
  brier <- round(epifusion_rt_brier(object, truth), 3)

  return(c(scenario,
           approach,
           inf_rmse,
           inf_coverage,
           scaled_inf_hpd_width,
           rt_rmse,
           rt_coverage,
           scale_rt_hpd_width,
           infection_crps,
           rt_crps,
           brier))
}

performancemetrics_epinow2 <- function(epinow2_fit, truth, approach, scenario) {
  time <- length(object$infection_trajectories$mean_infection_trajectory)
  table <- data.frame(
    Time = seq(1, time),
    Approach = rep(approach, time),
    Scenario = rep(scenario, time),
    Mean_Infected = object$infection_trajectories$mean_infection_trajectory,
    Lower95_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Lower,
    Upper95_Infected = object$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Upper,
    Mean_Rt = object$rt_trajectories$mean_rt_trajectory,
    Lower95_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.95$Lower,
    Upper95_Rt = object$rt_trajectories$rt_trajectory_hpdintervals$HPD0.95$Upper
  ) 
  
  table <- table %>%
    filter(Time < max(truth$t))
  table <- table %>%
    full_join(truth) %>%
    mutate(Infection_Coverage = ifelse((I >= Lower95_Infected & I <= Upper95_Infected), 1, 0)) %>%
    mutate(Rt_Coverage = ifelse((Rt >= Lower95_Rt & Rt <= Upper95_Rt), 1, 0))
  table$Approach[is.na(table$Approach)] <- approach
  table$Scenario[is.na(table$Scenario)] <- scenario
  
  table <- na.omit(table) %>%
    mutate(Scaled_Inf_HPD_Width  = (Upper95_Infected - Lower95_Infected)/Mean_Infected) %>%
    mutate(Scaled_Rt_HPD_Width = (Upper95_Rt - Lower95_Rt)/Mean_Rt)
  inf_rmse <- round(rmse(table$I, table$Mean_Infected), 2)
  inf_coverage <- round((sum(table$Infection_Coverage)/nrow(table))/0.95, 2)
  scaled_inf_hpd_width <- round(mean(table$Scaled_Inf_HPD_Width), 2)
  rt_rmse <- round(rmse(table$Rt, table$Mean_Rt), 3)
  rt_coverage <- round((sum(table$Rt_Coverage)/nrow(table))/0.95, 2)
  scale_rt_hpd_width <- round(mean(table$Scaled_Rt_HPD_Width), 2)
  infection_crps <- round(mean(epifusion_infection_crps(object, truth)), 2)
  rt_crps <- round(mean(epifusion_rt_crps(object, truth)),3)
  
  brier <- round(epifusion_rt_brier(object, truth), 3)
  
  return(c(scenario,
           approach,
           inf_rmse,
           inf_coverage,
           scaled_inf_hpd_width,
           rt_rmse,
           rt_coverage,
           scale_rt_hpd_width,
           infection_crps,
           rt_crps,
           brier))
}

epinow2_crps <- function(samples) {
  predictions <- as.matrix(dplyr::select(samples, !c(Time, Rt)))
  truth <- samples$Rt
  crps_value <- mean(scoringutils::crps_sample(truth, predictions = predictions))
  return(round(crps_value, 3))
}

epinow2_brier <- function(samples) {
  predictions <- as.matrix(dplyr::select(samples, !c(Time, Rt)))
  predictions[predictions < 0] <- 0
  predictions[predictions >= 1] <- 1
  truth <- as.numeric(samples$Rt >= 1)
  crps_value <- mean(scoringutils::brier_score(truth, predictions = rowMeans(predictions)))
  return(round(crps_value, 3))
}

casesfit <- function(object, truth, approach, scenario) {
  time <- seq(7, 7*(length(object$fitted_epi_cases$mean_fitted_epi_cases)), 7)
  table <- data.frame(
    Time = time,
    Approach = rep(approach, length(time)),
    Scenario = rep(scenario, length(time)),
    Mean_Cases_Fit = object$fitted_epi_cases$mean_fitted_epi_cases,
    Lower95_Cases_Fit = object$fitted_epi_cases$fitted_epi_cases_hpdintervals$HPD0.95$Lower,
    Upper95_Cases_Fit = object$fitted_epi_cases$fitted_epi_cases_hpdintervals$HPD0.95$Upper,
    Lower88_Cases_Fit = object$fitted_epi_cases$fitted_epi_cases_hpdintervals$HPD0.88$Lower,
    Upper88_Cases_Fit = object$fitted_epi_cases$fitted_epi_cases_hpdintervals$HPD0.88$Upper,
    Lower66_Cases_Fit = object$fitted_epi_cases$fitted_epi_cases_hpdintervals$HPD0.66$Lower,
    Upper66_Cases_Fit = object$fitted_epi_cases$fitted_epi_cases_hpdintervals$HPD0.66$Upper,
    Observed_Cases = truth$incidence[1:length(time)]
  )
  return(table)
}

modelresultstable <- function(object, approach, scenario) {
  params <- names(object$parameters)
  df <- data.frame(matrix(NA, ncol = 7, nrow = length(params)))
  for (i in 1:length(params)) {
    df[i,1] <- scenario
    df[i,2] <- approach
    df[i,3] <- params[i]
    df[i,4] <- round(mean(object$parameters[[i]]$samples), 5)
    df[i,5] <- paste0(hdi(round(object$parameters[[i]]$samples, 5), 0.95), collapse = ",")
    df[i,6] <- round(object$parameters[[i]]$rhat, 4)
    df[i,7] <- object$parameters[[i]]$ess
  }
  colnames(df) <- c("Approach", "Scenario", "Parameter", "Mean", "HPD Interval", "Rhat", "ESS")
  df <- df %>%
    filter(!grepl("Refactor_distribs_0", Parameter))
  if (approach == "Epi") {
    df <- df %>%
      filter(!grepl("psi", Parameter))
  } else if (approach == "Phylo") {
    df <- df %>%
      filter(!grepl("phi", Parameter))
  }
  
  if (scenario == "Sampling") {
    df <- df %>%
      filter(!grepl("changetime", Parameter))
  } 
  
  return(df)
}




