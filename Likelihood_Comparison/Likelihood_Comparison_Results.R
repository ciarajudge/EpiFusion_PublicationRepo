
##### Plots #####
#####Separated out by param #####
variedbeta <- phyloonlycomparison %>%
  filter(TrueGamma == 0.12 & TruePsi == 0.03 & ModelGamma == 0.12 & ModelPsi == 0.03 & TrueBeta == 0.19)  %>%
  filter(AvgLikelihood > -600) %>%
  mutate(MaxEpiFusion = max(AvgLikelihood), MaxBDSky = max(BDSkyDensity)) %>%
  mutate(MaxEpiFusion_x = ModelBeta[which(AvgLikelihood == max(AvgLikelihood))]) %>%
  mutate(MaxBDSky_x = ModelBeta[which(BDSkyDensity == max(BDSkyDensity))]) %>%
  mutate(TrueBeta_EpiFusion_y = AvgLikelihood[which.min(abs(ModelBeta-TrueBeta))]) %>%
  mutate(TrueBeta_BDSky_y = BDSkyDensity[which.min(abs(ModelBeta-TrueBeta))])

ggplot(variedbeta, aes(x = ModelBeta)) +
  geom_line(aes(y = AvgLikelihood), col = "darkolivegreen3", linewidth = 1) +
  geom_point(aes(y = AvgLikelihood), col = "darkolivegreen3", shape = 17, size = 2) +
  geom_line(aes(y = BDSkyDensity), col = "red", linewidth = 1) +
  geom_point(aes(y = BDSkyDensity), col = "red") +
  geom_vline(xintercept = 0.19)

variedgamma <- phyloonlycomparison %>%
  filter(TrueBeta == 0.23 & TruePsi == 0.03 & ModelBeta == 0.23 & ModelPsi == 0.03 & TrueGamma == 0.16)  %>%
  filter(AvgLikelihood > -650) 

ggplot(variedgamma, aes(x = ModelGamma)) +
  geom_line(aes(y = AvgLikelihood), col = "darkolivegreen3", linewidth = 1) +
  geom_point(aes(y = AvgLikelihood), col = "darkolivegreen3", shape = 17, size = 2) +
  geom_line(aes(y = BDSkyDensity), col = "red", linewidth = 1) +
  geom_point(aes(y = BDSkyDensity), col = "red") +
  geom_vline(xintercept = 0.16)

variedpsi <- phyloonlycomparison %>%
  filter(TrueBeta == 0.23 & TrueGamma == 0.12 & ModelBeta == 0.23 & ModelGamma == 0.12 & TruePsi == 0.03)  %>%
  filter(AvgLikelihood > -650) 

ggplot(variedpsi, aes(x = ModelPsi)) +
  geom_line(aes(y = AvgLikelihood), col = "darkolivegreen3", linewidth = 1) +
  geom_point(aes(y = AvgLikelihood), col = "darkolivegreen3", shape = 17, size = 2) +
  geom_line(aes(y = BDSkyDensity), col = "red", linewidth = 1) +
  geom_point(aes(y = BDSkyDensity), col = "red") +
  geom_vline(xintercept = 0.03)




##### Plot all together #####
phyloonlycomparisonmaster <- phyloonlycomparison %>%
  filter(AvgLikelihood > -700) %>%
  mutate(TrueParams = paste0("Beta", TrueBeta, "_Gamma", TrueGamma, "_Psi", TruePsi))

ggplot(phyloonlycomparisonmaster, aes(x = BDSkyDensity, y = AvgLikelihood)) +
  geom_smooth(method = "lm", col = "black") +
  geom_point(aes(col = TrueParams)) 

summary(lm(phyloonlycomparisonmaster$BDSkyDensity ~ phyloonlycomparisonmaster$AvgLikelihood))












