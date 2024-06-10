library(EpiNow2)

simdubcensored_SIR <- function(max, rate) {
  primary <- runif(1e6)
  secondary <- primary + rexp(1e6, rate)
  delay <- floor(secondary) - floor(primary)
  cdf <- ecdf(delay)(1:max)
  pmf <- c(cdf[1], diff(cdf))
  return(pmf)
}

simdubcensored_sampledelays_SIR <- function(max, rate1, samplerate) {
  sampledelays <- rexp(1e7, samplerate)
  recovdelays <- rexp(1e7, rate1)
  propobserved <- sum(sampledelays < recovdelays) / 1e7
  secondary <- sampledelays[sampledelays < recovdelays]
  print(length(secondary))
  primary <- runif(length(secondary))
  secondary <- primary + secondary
  delay <- floor(secondary) - floor(primary)
  cdf <- ecdf(delay)(1:max)
  pmf <- c(cdf[1], diff(cdf))
  return(list(pmf, propobserved))
}

simdubcensored_sampledelays_samplingscenario <- function(max, rate1, samplerate1, samplerate2, propsplit) {
  sampledelays1 <- rexp(1e7*propsplit, samplerate1)
  sampledelays2 <- rexp(1e7*(1-propsplit), samplerate2)
  sampledelays <- c(sampledelays1, sampledelays2)
  recovdelays <- rexp(1e7, rate1)
  propobserved <- sum(sampledelays < recovdelays) / 1e7
  secondary <- sampledelays[sampledelays < recovdelays]
  print(length(secondary))
  primary <- runif(length(secondary))
  secondary <- primary + secondary
  delay <- floor(secondary) - floor(primary)
  cdf <- ecdf(delay)(1:max)
  pmf <- c(cdf[1], diff(cdf))
  return(list(pmf, propobserved))
}

#####Introduction Scenario#####
#Read in the weekly incidence data, give it dates and put it into a compatible data frame
incidence <- c(0, unlist(read.table("Scenario_Testing/Data_Simulation/Main_Scenarios/baseline/baseline_weeklyincidence.txt")[,1]))
dates <- seq(as.Date(paste0(c("2021-01-01"), collapse = "")), 
             as.Date(paste0(c("2021-04-20"), collapse = "")), by="weeks")
reported_cases <- data.frame(date = dates, confirm = incidence)
reported_cases <- reported_cases[2:nrow(reported_cases),]

#Generation time
gen_time_epinow2 <- generation_time_opts(
  dist = Gamma(shape = 1, rate = 0.45*0.143, max = 14)
)


samplingstats <- simdubcensored_sampledelays_SIR(21, 1/7, 0.02)
sampling_delay <- NonParametric(samplingstats[[1]])
sampling_delay_epinow2 <- delay_opts(sampling_delay)


prop_observed <- samplingstats[[2]]
obs_epinow2 <- obs_opts(scale = list(mean = prop_observed, sd = 0.02), na = "accumulate")

options(mc.cores = 8)
out <- epinow(reported_cases = reported_cases, 
              generation_time = gen_time_epinow2, #From above
              delays = sampling_delay_epinow2, 
              obs = obs_epinow2,
              rt = rt_opts(prior = list(mean = 3, sd = 1.0)), 
              gp = gp_opts(basis_prop = 0.2), 
              #gp = NULL,
              stan = stan_opts(samples = 4000),
              horizon = 0, 
              target_folder = "results",
              logs = file.path("logs", Sys.Date()),
              return_output = TRUE, 
              verbose = TRUE)
plot(out)
saveRDS(out, "Benchmarking/EpiNow2/baseline_epinow2.RDS")


#####Transmission step-change scenario#####
#Read in the weekly incidence data, give it dates and put it into a compatible data frame
incidence <- c(0, unlist(read.table("Scenario_Testing/Data_Simulation/Main_Scenarios/transmission/transmission_weeklyincidence.txt")[,1]))
dates <- seq(as.Date(paste0(c("2021-01-01"), collapse = "")), 
             as.Date(paste0(c("2021-08-01"), collapse = "")), by="weeks")
reported_cases <- data.frame(date = dates, confirm = incidence)
reported_cases <- reported_cases[2:nrow(reported_cases),]

#Generation time
gen_time_epinow2 <- generation_time_opts(
  dist = Gamma(shape = 1, rate = 0.45*0.143, max = 14)
)


samplingstats <- simdubcensored_sampledelays_SIR(21, 1/7, 0.02)
sampling_delay <- NonParametric(samplingstats[[1]])
sampling_delay_epinow2 <- delay_opts(sampling_delay)


prop_observed <- samplingstats[[2]]
obs_epinow2 <- obs_opts(scale = list(mean = prop_observed, sd = 0.02), na = "accumulate")

options(mc.cores = 8)
out <- epinow(reported_cases = reported_cases, 
              generation_time = gen_time_epinow2, #From above
              delays = sampling_delay_epinow2,
              obs = obs_epinow2,
              rt = rt_opts(prior = list(mean = 3, sd = 1.0)), #Rt is 1 for a lot of this scenario
              gp = gp_opts(basis_prop = 0.2), 
              stan = stan_opts(samples = 4000),
              horizon = 0, 
              target_folder = "results",
              logs = file.path("logs", Sys.Date()),
              return_output = TRUE, 
              verbose = TRUE)
plot(out)
saveRDS(out, "Benchmarking/EpiNow2/transmission_epinow2.RDS")


#####Sampling stepchange scenario#####
#Read in the weekly incidence data, give it dates and put it into a compatible data frame
incidence <- c(0, unlist(read.table("Scenario_Testing/Data_Simulation/Main_Scenarios/sampling/sampling_weeklyincidence.txt")[,1]))
dates <- seq(as.Date(paste0(c("2021-01-01"), collapse = "")), 
             as.Date(paste0(c("2021-04-15"), collapse = "")), by="weeks")
reported_cases <- data.frame(date = dates, confirm = incidence)
reported_cases <- reported_cases[2:nrow(reported_cases),]

#Get generation time
gen_time_epinow2 <- generation_time_opts(
  dist = Gamma(shape = 1, rate = 0.45*0.143, max = 14)
)

#Get sampling delays, this one is more complicated due to the stepchange 
samplingstats <- simdubcensored_sampledelays_samplingscenario(21, 1/7, 0.005, 0.05, 0.3)

sampling_delay <- NonParametric(samplingstats[[1]])
sampling_delay_epinow2 <- delay_opts(sampling_delay)

prop_observed <- samplingstats[[2]]
obs_epinow2 <- obs_opts(scale = list(mean = prop_observed, sd = 0.05), na = "accumulate")





options(mc.cores = 8)
out <- epinow(reported_cases = reported_cases, 
              generation_time = gen_time_epinow2, #From above
              delays = sampling_delay_epinow2,
              obs = obs_epinow2,
              rt = rt_opts(prior = list(mean = 3, sd = 1.0)), #I think this makes sense, Rt starts high and gets low, mean over the whole time is 1.5ish
              gp = gp_opts(basis_prop = 0.2), 
              stan = stan_opts(samples = 4000),
              horizon = 0, 
              target_folder = "results",
              logs = file.path("logs", Sys.Date()),
              return_output = TRUE, 
              verbose = TRUE)
plot(out)
saveRDS(out, "Benchmarking/EpiNow2/sampling_EpiNow2.RDS")
