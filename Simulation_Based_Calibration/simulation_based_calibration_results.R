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
library(ggplot2)

##### THIS SCRIPT PARSES THE RESULTS OF THE SIMULATION BASED CALIBRATION #####

numreplicates <- 700

draws <- read.csv("Simulation_Based_Calibration/paramdraws.csv", header = T)

for (i in 1:numreplicates) {
  label <- paste0("rep_", i)
  filestem <- paste0("Simulation_Based_Calibration/",label,"/simulate")
  if (!file.exists(paste0(filestem, "/likelihoods_chain1.txt"))) {
    print(label)
  }
}

simulationsummarytable <- data.frame(matrix(0, nrow = 700, ncol = 28))
gammacalibrationtable1 <- data.frame(matrix(0, nrow = 700, ncol = 19))
gammacalibrationtable2 <- data.frame(matrix(0, nrow = 700, ncol = 19))
phicalibrationtable1 <- data.frame(matrix(0, nrow = 700, ncol = 19))
phicalibrationtable2 <- data.frame(matrix(0, nrow = 700, ncol = 19))
psicalibrationtable1 <- data.frame(matrix(0, nrow = 700, ncol = 19))
psicalibrationtable2 <- data.frame(matrix(0, nrow = 700, ncol = 19))
betacalibrationtable <- data.frame(matrix(0, nrow = 700, ncol = 19))


alpha <- seq(0.05, 0.95, 0.05)
z <- c(0.065, 0.125, 0.19, 0.254, 0.319, 0.385, 0.455, 0.525, 0.6, 0.675, 0.755, 0.845, 0.935, 1.138, 1.15, 1.28, 1.44, 1.645, 1.96)


for (i in 1:numreplicates) {
  label <- paste0("rep_", i)
  print(label)
  epifusionoutputfolder <- paste0("Simulation_Based_Calibration/",label,"/simulate/")
  simulationparams <- paste0("Simulation_Based_Calibration/",label,"/simulate.xml")
  truthtrajectories <- paste0("Simulation_Based_Calibration/",label,"/simulate.txt")
  treefile <- paste0("Simulation_Based_Calibration/",label,"/simulate_downsampledtree.tree")
  caseincidence <- paste0("Simulation_Based_Calibration/",label,"/simulate_weeklyincidence.txt")

  #Rep
  simulationsummarytable[i, 1] <- i
  #Acceptance Rate
  simulationsummarytable[i, 2] <- mean(c(read.table(paste0(epifusionoutputfolder, "acceptance_chain0.txt"))[,1],
                                         read.table(paste0(epifusionoutputfolder, "acceptance_chain1.txt"))[,1]))
  #Likelihood Chain 1
  chain0likelihood <- read.table(paste0(epifusionoutputfolder, "likelihoods_chain0.txt"))[,1]
  simulationsummarytable[i, 3] <- mean(chain0likelihood[100:1001][chain0likelihood[100:1001]!= -Inf])
  #Likelihood Chain 2
  chain1likelihood <- read.table(paste0(epifusionoutputfolder, "likelihoods_chain1.txt"))[,1]
  simulationsummarytable[i, 4] <- mean(chain1likelihood[100:1001][chain1likelihood[100:1001]!= -Inf])
  #Delta
  simulationsummarytable[i, 5] <- abs(simulationsummarytable[i, 4] - simulationsummarytable[i, 3])
  #True peak
  simulationsummarytable[i, 6] <- max(read.csv(truthtrajectories, header = F)[,4])
  #True gamma
  simulationsummarytable[i, 7] <- draws$gamma[i]
  #True phi
  simulationsummarytable[i, 8] <- draws$phi[i]
  #True psi
  simulationsummarytable[i, 9] <- draws$psi[i]
  #Num cases
  simulationsummarytable[i, 10] <- sum(read.table(caseincidence, header = F)[,1])
  #Tree size
  simulationsummarytable[i, 11] <- length(read.tree(treefile)$tip.label)

  epifusionoutput <- extract_posterior_epifusion(load_raw_epifusion(epifusionoutputfolder), 0.1)
  #% traj in hpd
  truth <- read.csv(truthtrajectories, header = F)
  colnames(truth) <- c("Time", "R", "S", "I", "sample")
  truth <- truth %>%
    mutate(Time = ceiling(Time)) %>%
    group_by(Time) %>%
    mutate(Truth = mean(I)) %>%
    select(Time, Truth) %>%
    distinct() %>%
    left_join(data.frame(Time = seq(1, nrow(truth))))
  l <- min(nrow(truth), length(epifusionoutput$infection_trajectories$mean_infection_trajectory))
  simulationsummarytable[i, 12] <- sum(truth > epifusionoutput$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Lower &
        truth < epifusionoutput$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Upper)/l

  #gamma mean 
  simulationsummarytable[i, 13] <- mean(epifusionoutput$parameters$gamma$samples)
  confints <- hdi(epifusionoutput$parameters$gamma$samples, 0.95)
  #gamma lower 95
  simulationsummarytable[i, 14] <- confints[1]
  #gamma upper 95
  simulationsummarytable[i, 15] <- confints[2]
  #phi mean 
  simulationsummarytable[i, 16] <- mean(epifusionoutput$parameters$phi$samples)
  confints <- hdi(epifusionoutput$parameters$phi$samples, 0.95)
  #phi lower 95
  simulationsummarytable[i, 17] <- confints[1]
  #phi upper 95
  simulationsummarytable[i, 18] <- confints[2]
  #psi mean 
  simulationsummarytable[i, 19] <- mean(epifusionoutput$parameters$psi$samples)
  confints <- hdi(epifusionoutput$parameters$psi$samples, 0.95)
  #psi lower 95
  simulationsummarytable[i, 20] <- confints[1]
  #psi upper 95
  simulationsummarytable[i, 21] <- confints[2]
  
  
  for (j in 1:19) {
    gammahpd <- hdi(epifusionoutput$parameters$gamma$samples, alpha[j])
    if (draws$gamma[i] >= gammahpd[1] & draws$gamma[i] <= gammahpd[2]) {
      gammacalibrationtable1[i,j] <- 1
    }
    phihpd <- hdi(epifusionoutput$parameters$phi$samples, alpha[j])
    if (draws$phi[i] >= phihpd[1] & draws$phi[i] <= phihpd[2]) {
      phicalibrationtable1[i,j] <- 1
    }
    psihpd <- hdi(epifusionoutput$parameters$psi$samples, alpha[j])
    if (draws$psi[i] >= psihpd[1] & draws$psi[i] <= psihpd[2]) {
      psicalibrationtable1[i,j] <- 1
    }
  }
  
  for (j in 1:19) {
    gammaupper <- mean(epifusionoutput$parameters$gamma$samples) + (z[j]*sd(epifusionoutput$parameters$gamma$samples))
    gammalower <- mean(epifusionoutput$parameters$gamma$samples) - (z[j]*sd(epifusionoutput$parameters$gamma$samples))
    if (draws$gamma[i] >= gammalower & draws$gamma[i] <= gammaupper) {
      gammacalibrationtable2[i,j] <- 1
    }
    phiupper <- mean(epifusionoutput$parameters$phi$samples) + (z[j]*sd(epifusionoutput$parameters$phi$samples))
    philower <- mean(epifusionoutput$parameters$phi$samples) - (z[j]*sd(epifusionoutput$parameters$phi$samples))
    if (draws$phi[i] >= philower & draws$phi[i] <= phiupper) {
      phicalibrationtable2[i,j] <- 1
    }
    psiupper <- mean(epifusionoutput$parameters$psi$samples) + (z[j]*sd(epifusionoutput$parameters$psi$samples))
    psilower <- mean(epifusionoutput$parameters$psi$samples) - (z[j]*sd(epifusionoutput$parameters$psi$samples))
    if (draws$psi[i] >= psilower & draws$psi[i] <= psiupper) {
      psicalibrationtable2[i,j] <- 1
    }
  }
  
  simulationsummarytable[i, 22] <- sum(gammacalibrationtable1[i,]/19)
  simulationsummarytable[i, 23] <- sum(phicalibrationtable1[i,]/19)
  simulationsummarytable[i, 24] <- sum(psicalibrationtable1[i,]/19)
  simulationsummarytable[i, 25] <- sum(gammacalibrationtable2[i,]/19)
  simulationsummarytable[i, 26] <- sum(phicalibrationtable2[i,]/19)
  simulationsummarytable[i, 27] <- sum(psicalibrationtable2[i,]/19)
  
  simulationsummarytable[i, 28] <- suppressWarnings(read.table(paste0(epifusionoutputfolder, "timings.txt"), header = F)[1,1])
}

##### Assign Colnames to our newly made tables #####
colnames(simulationsummarytable) <- c('Rep', 'Acceptance', 'Likelihood1', 'Likelihood2', 'LikelihoodDelta', 'TruePeak', 'TrueGamma',
                                      'TruePhi', 'TruePsi', 'Cases', 'TreeSize', 'PercentTrajCaptured', 'GammaMean', 'GammaLower95',
                                      'GappaUpper95', 'PhiMean', 'PhiLower95', 'PhiUpper95', 'PsiMean', 'PsiLower95', 
                                      'PsiUpper95', 'GammaCalibration1', 'PhiCalibration1', 'PsiCalibration1', 
                                      'GammaCalibration2', 'PhiCalibration2', 'PsiCalibration2', 'RunTime')
colnames(gammacalibrationtable1) <- paste0("alpha", seq(0.05, 0.95, 0.05))
colnames(gammacalibrationtable2) <- paste0("alpha", seq(0.05, 0.95, 0.05))
colnames(phicalibrationtable1) <- paste0("alpha", seq(0.05, 0.95, 0.05))
colnames(phicalibrationtable2) <- paste0("alpha", seq(0.05, 0.95, 0.05))
colnames(psicalibrationtable1) <- paste0("alpha", seq(0.05, 0.95, 0.05))
colnames(psicalibrationtable2) <- paste0("alpha", seq(0.05, 0.95, 0.05))

##### For the actual paper, the runs were done in batches on difference machines (which have different speeds) so useful to assign these #####
batches <- c(rep(1, 50), rep(2, 50), rep(3, 100), rep(4, 100), rep(5, 100), rep(4, 100), rep(6, 200))
simulationsummarytable <- simulationsummarytable %>%
  mutate(Batch = as.factor(batches))

##### Save complete csvs #####
write.csv(simulationsummarytable, "Simulation_Based_Calibration/simulationsummarytable.csv", row.names = F)
write.csv(gammacalibrationtable1, "Simulation_Based_Calibration/gammacalibrationtableHDI.csv", row.names = F)
write.csv(gammacalibrationtable1, "Simulation_Based_Calibration/gammacalibrationtableZ.csv", row.names = F)
write.csv(phicalibrationtable1, "Simulation_Based_Calibration/phicalibrationtableHDI.csv", row.names = F)
write.csv(phicalibrationtable2, "Simulation_Based_Calibration/phicalibrationtableZ.csv", row.names = F)
write.csv(psicalibrationtable1, "Simulation_Based_Calibration/psicalibrationtableHDI.csv", row.names = F)
write.csv(psicalibrationtable2, "Simulation_Based_Calibration/psicalibrationtableZ.csv", row.names = F)

##### Plot the raw results #####
layout(matrix(seq(1, 9), nrow = 3, byrow = T))
plot(simulationsummarytable$TrueGamma, simulationsummarytable$GammaMean, main = "True Gamma x Mean Gamma", ylim = c(0,0.3), xlim = c(0, 0.3))
abline(lm(simulationsummarytable$GammaMean ~ simulationsummarytable$TrueGamma), col = "red")
lines(c(0, 1), c(0, 1),lty = 2)
plot(simulationsummarytable$TruePhi, simulationsummarytable$PhiMean, main = "True Phi x Mean Phi", ylim = c(0,0.06), xlim = c(0, 0.06))
lines(c(0, 1), c(0, 1),lty = 2)
abline(lm(simulationsummarytable$PhiMean ~ simulationsummarytable$TruePhi), col = "red")
plot(simulationsummarytable$TruePsi, simulationsummarytable$PsiMean, main = "True Psi x Mean Psi", ylim = c(0,0.005), xlim = c(0, 0.005))
lines(c(0, 1), c(0, 1),lty = 2)
abline(lm(simulationsummarytable$PsiMean ~ simulationsummarytable$TruePsi), col = "red")
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(gammacalibrationtable1)), main = "Gamma Calibration (HDI)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(phicalibrationtable1)), main = "Phi Calibration (HDI)",  ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(psicalibrationtable1)), main = "Psi Calibration (HDI)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(gammacalibrationtable2)), main = "Gamma Calibration (Z)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(phicalibrationtable2)), main = "Phi Calibration (Z)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(psicalibrationtable2)), main = "Psi Calibration (Z)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)

##### Discard non convergent runs or runs with too small ESS' (low acceptance rate) #####
simulationsummarytable_acceptancefilter <- simulationsummarytable %>%
  mutate(LikelihoodDelta_Percent = LikelihoodDelta/abs(Likelihood2)) %>%
  filter(Acceptance > 0.25 & LikelihoodDelta_Percent < 0.02)
gammacalibrationtable1_acceptancefilter <- gammacalibrationtable1[simulationsummarytable_acceptancefilter$Rep,]
phicalibrationtable1_acceptancefilter <- phicalibrationtable1[simulationsummarytable_acceptancefilter$Rep,]
psicalibrationtable1_acceptancefilter <- psicalibrationtable1[simulationsummarytable_acceptancefilter$Rep,]
gammacalibrationtable2_acceptancefilter <- gammacalibrationtable2[simulationsummarytable_acceptancefilter$Rep,]
phicalibrationtable2_acceptancefilter <- phicalibrationtable2[simulationsummarytable_acceptancefilter$Rep,]
psicalibrationtable2_acceptancefilter <- psicalibrationtable2[simulationsummarytable_acceptancefilter$Rep,]


layout(matrix(seq(1, 9), nrow = 3, byrow = T))
plot(simulationsummarytable_acceptancefilter$TrueGamma, simulationsummarytable_acceptancefilter$GammaMean, main = "True Gamma x Mean Gamma", ylim = c(0,0.3), xlim = c(0, 0.3))
abline(lm(simulationsummarytable_acceptancefilter$GammaMean ~ simulationsummarytable_acceptancefilter$TrueGamma), col = "red")
lines(c(0, 1), c(0, 1),lty = 2)
plot(simulationsummarytable_acceptancefilter$TruePhi, simulationsummarytable_acceptancefilter$PhiMean, main = "True Phi x Mean Phi", ylim = c(0,0.06), xlim = c(0, 0.06))
lines(c(0, 1), c(0, 1),lty = 2)
abline(lm(simulationsummarytable_acceptancefilter$PhiMean ~ simulationsummarytable_acceptancefilter$TruePhi), col = "red")
plot(simulationsummarytable_acceptancefilter$TruePsi, simulationsummarytable_acceptancefilter$PsiMean, main = "True Psi x Mean Psi", ylim = c(0,0.005), xlim = c(0, 0.005))
lines(c(0, 1), c(0, 1),lty = 2)
abline(lm(simulationsummarytable_acceptancefilter$PsiMean ~ simulationsummarytable_acceptancefilter$TruePsi), col = "red")
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(gammacalibrationtable1_acceptancefilter)), main = "Gamma Calibration (HDI)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(phicalibrationtable1_acceptancefilter)), main = "Phi Calibration (HDI)",  ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(psicalibrationtable1_acceptancefilter)), main = "Psi Calibration (HDI)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(gammacalibrationtable2_acceptancefilter)), main = "Gamma Calibration (Z)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(phicalibrationtable2_acceptancefilter)), main = "Phi Calibration (Z)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)
plot(seq(0.05, 0.95, 0.05), colMeans(na.omit(psicalibrationtable2_acceptancefilter)), main = "Psi Calibration (Z)", ylim = c(0, 1), xlim = c(0, 1))
lines(c(0, 1), c(0, 1), lty = 2)


write.csv(simulationsummarytable_acceptancefilter, "Simulation_Based_Calibration/simulationsummarytable_acceptanceandconvergencefilter.csv", row.names = F)
write.csv(gammacalibrationtable2_acceptancefilter, "Simulation_Based_Calibration/gammacalibrationtable_acceptanceandconvergencefilter.csv", row.names = F)
write.csv(phicalibrationtable2_acceptancefilter, "Simulation_Based_Calibration/phicalibrationtable_acceptanceandconvergencefilter.csv", row.names = F)
write.csv(psicalibrationtable2_acceptancefilter, "Simulation_Based_Calibration/psicalibrationtable_acceptanceandconvergencefilter.csv", row.names = F)


##### Make two new tables, one with all the model fits, one with the trajectory calibration #####
model_fits <- data.frame(matrix(0, nrow = 0, ncol = 10))
colnames(model_fits) <- c("Rep", "Time", "Mean", "lower95", "upper95", "lower88", "upper88", "lower66", "upper66", "Truth")
trajcalibrationtable <- data.frame(matrix(0, nrow = 700, ncol = 19))


for (i in 1:numreplicates) {
  label <- paste0("rep_", i)
  print(label)
  filestem <- paste0("Simulation_Based_Calibration/",label,"/simulate/")
  epifusionoutput <- extract_posterior_epifusion(load_raw_epifusion(filestem), 0.1)
  tmpdf <- data.frame(Rep = rep(i, length(epifusionoutput$infection_trajectories$mean_infection_trajectory)),
                      Time = seq(0, length(epifusionoutput$infection_trajectories$mean_infection_trajectory)-1),
                      Mean = epifusionoutput$infection_trajectories$mean_infection_trajectory,
                      lower95 = epifusionoutput$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Lower,
                      upper95 = epifusionoutput$infection_trajectories$infection_trajectory_hpdintervals$HPD0.95$Upper,
                      lower88 = epifusionoutput$infection_trajectories$infection_trajectory_hpdintervals$HPD0.88$Lower,
                      upper88 = epifusionoutput$infection_trajectories$infection_trajectory_hpdintervals$HPD0.88$Upper,
                      lower66 = epifusionoutput$infection_trajectories$infection_trajectory_hpdintervals$HPD0.66$Lower,
                      upper66 = epifusionoutput$infection_trajectories$infection_trajectory_hpdintervals$HPD0.6$Upper)
  truth <- read.csv(paste0("Simulation_Based_Calibration/",label,"/discretetruth.csv"), header = T) 
  tmpdf <- tmpdf %>%
    left_join(truth)
  model_fits <- rbind(model_fits, tmpdf)
  l <- min(nrow(truth), length(epifusionoutput$infection_trajectories$mean_infection_trajectory))
  
  for (j in 1:length(alpha)) {
    trajcalibrationtable[i, j] <- sum(truth$Truth[1:l] >= hdi(epifusionoutput$infection_trajectories$infection_trajectory_samples, alpha[j])[1,1:l] & 
                                           truth$Truth[1:l] <= hdi(epifusionoutput$infection_trajectories$infection_trajectory_samples, alpha[j])[2,1:l])/l
  }
  
}

##### Save these tables and plot them #####

write.csv(model_fits, "Simulation_Based_Calibration/trajectoryfits.csv", row.names = F)

master_model_fits <- model_fits %>%
  left_join(simulationsummarytable)
write.csv(model_fits, "Simulation_Based_Calibration/trajectoryfitsandrepdata.csv", row.names = F)

model_fits_100 <- filter(master_model_fits, Rep %in% seq(1, 100))
model_fits_200 <- filter(master_model_fits, Rep %in% seq(101, 200))
model_fits_300 <- filter(master_model_fits, Rep %in% seq(201, 300))
model_fits_400 <- filter(master_model_fits, Rep %in% seq(301, 400))
model_fits_500 <- filter(master_model_fits, Rep %in% seq(401, 500))
model_fits_600 <- filter(master_model_fits, Rep %in% seq(501, 600))
model_fits_700 <- filter(master_model_fits, Rep %in% seq(601, 700))


pdf("Simulation_Based_Calibration/trajectory_fits.pdf", height = 20, width = 20, onefile = TRUE)
ggplot(model_fits_100, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower88, ymax = upper88, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower66, ymax = upper66, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_line(aes(y = Mean, col = as.factor(Rep)), show.legend = F) +
  geom_line(aes(y = Truth), col = "black") +
  facet_wrap(~Rep, scales = "free")
ggplot(model_fits_200, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower88, ymax = upper88, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower66, ymax = upper66, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_line(aes(y = Mean, col = as.factor(Rep)), show.legend = F) +
  geom_line(aes(y = Truth), col = "black") +
  facet_wrap(~Rep, scales = "free")
ggplot(model_fits_300, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower88, ymax = upper88, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower66, ymax = upper66, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_line(aes(y = Mean, col = as.factor(Rep)), show.legend = F) +
  geom_line(aes(y = Truth), col = "black") +
  facet_wrap(~Rep, scales = "free")
ggplot(model_fits_400, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower88, ymax = upper88, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower66, ymax = upper66, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_line(aes(y = Mean, col = as.factor(Rep)), show.legend = F) +
  geom_line(aes(y = Truth), col = "black") +
  facet_wrap(~Rep, scales = "free")
ggplot(model_fits_500, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower88, ymax = upper88, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower66, ymax = upper66, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_line(aes(y = Mean, col = as.factor(Rep)), show.legend = F) +
  geom_line(aes(y = Truth), col = "black") +
  facet_wrap(~Rep, scales = "free")
ggplot(model_fits_600, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower88, ymax = upper88, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower66, ymax = upper66, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_line(aes(y = Mean, col = as.factor(Rep)), show.legend = F) +
  geom_line(aes(y = Truth), col = "black") +
  facet_wrap(~Rep, scales = "free")
ggplot(model_fits_700, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower88, ymax = upper88, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower66, ymax = upper66, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_line(aes(y = Mean, col = as.factor(Rep)), show.legend = F) +
  geom_line(aes(y = Truth), col = "black") +
  facet_wrap(~Rep, scales = "free")
dev.off()


##### Filter to the same reps that were filtered out earlier #####
filtered_model_fits <- master_model_fits %>%
  filter(Rep %in% simulationsummarytable_acceptancefilter$Rep)
write.csv(model_fits, "Simulation_Based_Calibration/trajectoryfitsandrepdata_acceptanceandconvergencefilter.csv", row.names = F)


pdf("Simulation_Based_Calibration/sampledfits.pdf", height = 20, width = 20, onefile = TRUE)
ggplot(random_sample_fits, aes(x = Time)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower88, ymax = upper88, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_ribbon(aes(ymin = lower66, ymax = upper66, fill = as.factor(Rep)), alpha = 0.2, show.legend = F) +
  geom_line(aes(y = Mean, col = as.factor(Rep)), show.legend = F) +
  geom_line(aes(y = Truth), col = "black") +
  facet_wrap(~Rep, scales = "free")
dev.off()


##### Plot the Trajectory Convergence for the full runs #####
layout(matrix(c(1)))
plot(seq(0.05, 0.95, 0.05), colMeans(trajcalibrationtable), ylim = c(0, 1), xlim = c(0, 1))
lines(c(0,1), c(0,1), lty = 2)

write.csv(trajcalibrationtable, "Simulation_Based_Calibration/trajcalibrationtable.csv", row.names = F)

trajcalibrationtable_acceptanceandconvergencefilter <- trajcalibrationtable[simulationsummarytable_acceptancefilter$Rep,]
plot(seq(0.05, 0.95, 0.05), colMeans(trajcalibrationtable_acceptanceandconvergencefilter), ylim = c(0, 1), xlim = c(0, 1))
lines(c(0,1), c(0,1), lty = 2)
write.csv(trajcalibrationtable_acceptanceandconvergencefilter, "Simulation_Based_Calibration/trajcalibrationtable_acceptanceandconvergencefilter.csv", row.names = F)


##### Get the real vs inferred values of beta #####
betavstruth <- matrix(nrow = 0, ncol = 5)
for (i in 1:699) {
  label <- paste0("rep_", i)
  print(label)
  epifusionoutputfolder <- paste0("Simulation_Based_Calibration/",label,"/simulate/")
  truth <- read.csv(paste0("Simulation_Based_Calibration/",label,"/beta.txt"))[,1]
  beta1 <- read.csv(paste0(epifusionoutputfolder, "betas_chain0.txt"))[100:1001,]
  beta2 <- read.csv(paste0(epifusionoutputfolder, "betas_chain1.txt"))[100:1001,]
  beta <- rbind(beta1, beta2)
  coverage <- c()
  len <- min((ncol(beta)-1), length(truth))
  truth <- truth[1:len]
  HPD <- HDInterval::hdi(beta, 0.95)
  lower <- HPD[1, 1:len]
  upper <- HPD[2, 1:len]
  for (k in 1:length(lower)) {
    betavstruth <- rbind(betavstruth, c(i, truth[k], mean(na.omit(beta[,k])), lower[k], upper[k]))
  }
}

write.csv(betavstruth, "Simulation_Based_Calibration/betavstruth_raw.csv")

repsused <- read.csv("Simulation_Based_Calibration/simulationsummarytable_acceptanceandconvergencefilter.csv")$Rep

colnames(betavstruth) <- c("Rep", "True_Beta", "Mean_Beta", "Lower95_Beta", "Upper95_Beta")

beta_vs_truth_filtered <- data.frame(betavstruth) %>%
  filter(Rep %in% repsused)

#Subsample as otherwise there's a lot of data points:
beta_vs_truth_filtered <- beta_vs_truth_filtered[sample(nrow(beta_vs_truth_filtered), 10000),]
write.csv(beta_vs_truth_filtered, "Simulation_Based_Calibration/betavstruth_acceptanceandconvergencefilter.csv")


betacoverage <- matrix(nrow = 700, ncol = 19)
for (i in 1:698) {
  label <- paste0("rep_", i)
  print(label)
  epifusionoutputfolder <- paste0("Simulation_Based_Calibration/",label,"/simulate/")
  truth <- read.csv(paste0("Simulation_Based_Calibration/",label,"/beta.txt"))[,1]
  beta1 <- read.csv(paste0(epifusionoutputfolder, "betas_chain0.txt"))[100:1001,]
  beta2 <- read.csv(paste0(epifusionoutputfolder, "betas_chain1.txt"))[100:1001,]
  beta <- rbind(beta1, beta2)
  len <- min((ncol(beta)-1), length(truth))
  truth <- truth[1:len]
  iter <- 0
  for (j in seq(0.05, 0.95, 0.05)) {
    iter <- iter + 1
    coverage <- c()
    HPD <- HDInterval::hdi(beta, j)
    lower <- HPD[1, 1:len]
    upper <- HPD[2, 1:len]
    for (k in 1:length(lower)) {
      coverage[k] <- ifelse(lower[k] < truth[k] & upper[k] > truth[k], 1, 0)
    }
    cov <- mean(coverage)
    betacoverage[i, iter] <- cov
  }
}

betacoverage_filtered <- betacoverage[repsused,]
betacoverage_filtered <- na.omit(betacoverage_filtered)

write.csv(betacoverage_filtered, "Simulation_Based_Calibration/betacalibrationtable_acceptanceandconvergencefilter.csv", row.names = F)



##### Get ESS x runtime and epidemic length x runtime #####
numreplicates <- 700

more_runtime_metrics <- data.frame(matrix(0, nrow = 700, ncol = 7))

for (i in 1:numreplicates) {
  label <- paste0("rep_", i)
  print(label)
  epifusionoutputfolder <- paste0("Simulation_Based_Calibration/",label,"/simulate/")
  
  #Rep
  more_runtime_metrics[i, 1] <- i
  
  epifusionoutput <- extract_posterior_epifusion(load_raw_epifusion(epifusionoutputfolder), 0.1)
  
  more_runtime_metrics[i, 2] <- epifusionoutput$parameters$gamma$ess
  more_runtime_metrics[i, 3] <- epifusionoutput$parameters$psi$ess
  more_runtime_metrics[i, 4] <- epifusionoutput$parameters$phi$ess
  more_runtime_metrics[i, 5] <- epifusionoutput$parameters$initialBeta$ess
  more_runtime_metrics[i, 6] <- length(epifusionoutput$infection_trajectories$mean_infection_trajectory)
  
  more_runtime_metrics[i, 7] <- suppressWarnings(read.table(paste0(epifusionoutputfolder, "timings.txt"), header = F)[1,1])
}

summarytable <- read.csv("Simulation_Based_Calibration/simulationsummarytable_acceptanceandconvergencefilter.csv")
colnames(more_runtime_metrics) <- c("Rep", "Gamma_ESS", "Psi_ESS", "Phi_ESS", "InitialBeta_ESS", "Epidemic_Length", "RunTime")

runtimexlength <- dplyr::select(more_runtime_metrics, Rep, Epidemic_Length)
summarytable <- left_join(summarytable, runtimexlength)

write.csv(summarytable, "Simulation_Based_Calibration/simulationsummarytable_acceptanceandconvergencefilter.csv", row.names = F)

ESSpermin <- more_runtime_metrics %>%
  dplyr::select(Rep, Gamma_ESS, Psi_ESS, Phi_ESS, InitialBeta_ESS, RunTime) %>%
  mutate(RunTime = RunTime / 6e10) %>%
  mutate(Gamma_ESS = Gamma_ESS/RunTime,
         Psi_ESS = Psi_ESS/RunTime,
         Phi_ESS = Phi_ESS/RunTime,
         InitialBeta_ESS = InitialBeta_ESS/RunTime)

write.csv(ESSpermin, "Simulation_Based_Calibration/ESS_per_min.csv", row.names = F)


