library(ape)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtree)
library(phytools)
library(sf)
library(RColorBrewer)
library(castor)


##### Load and clean tree #####
#Load tree, and data, get factor levels
big_mcc_tree <- read.nexus("EBOV_SierraLeone/Raw_Data/Makona_1610_cds_ig.MCC.tree")
region_lookup <- read.csv("EBOV_SierraLeone/Raw_Data/regional_lookup_adm1.csv")
SL_shape <- read_sf("EBOV_SierraLeone/Raw_Data/geoBoundaries-SLE-ADM2-all/", layer = "geoBoundaries-SLE-ADM2")

regionlabels_map <- unique(SL_shape$shapeName)
regionlabels_tree <- regionlabels_map
regionlabels_tree[2] <- "PortLoko"
regionlabels_tree[13] <- "WesternRural"
regionlabels_tree[14] <- "WesternUrban"
regionlabels_tree[15] <- "WesternArea"

SL_shape <- mutate(SL_shape, Region = factor(shapeName, levels = regionlabels_map))


#Turn tree labels into metadata table
alllabels <- big_mcc_tree$tip.label
sequence_metadata <- data.frame(matrix(0, ncol = 7, nrow = length(alllabels)))
for (i in 1:length(alllabels)) {
  label <- str_remove_all(alllabels[i], "'")
  row <- c(label, unlist(str_split(label, "\\|")))
  sequence_metadata[i,] <- row
}
colnames(sequence_metadata) <- c("tip", "Pathogen", "ID", "Accession", "Country", "Region", "Date")
big_mcc_tree$tip.label <- str_remove_all(big_mcc_tree$tip.label, "'")

sequence_metadata_final <- sequence_metadata %>%
  mutate(Region = ifelse(Country == "SLE", Region, NA)) %>%
  mutate(Region = ifelse(Region == "?", NA, Region)) %>%
  mutate(Region = ifelse(Region == "", NA, Region)) %>%
  filter(Country != "?") %>%
  left_join(region_lookup) %>%
  mutate(Region = factor(Region, levels = regionlabels_tree))

sierra_leone_sequence_metadata <- sequence_metadata %>%
  dplyr::filter(Country == "SLE") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))


# Plot coloured by country
colourednational <- ggtree(big_mcc_tree) %<+% sequence_metadata_final +
  geom_tippoint(aes(color = Country)) +
  #scale_color_brewer("Region", palette="Spectral") +
  theme_tree2()
colourednational
ggsave("EBOV_SierraLeone/Intermediate_Figures/MCC_colourcountry_tree.pdf", plot = colourednational, width = 10, height = 20, units = "in" )


# Plot with shape by country, colour by region (within sierra leone)
map <- ggplotGrob(ggplot(SL_shape) +
  geom_sf(aes(fill = Region), show.legend = F) +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill = NA)))
mapgrob <- list("1373" = map)

colouredmccregional <- ggtree(big_mcc_tree) %<+% sequence_metadata_final +
  geom_tippoint(aes(color = Region, shape = Country)) +
  #scale_color_brewer("Region", palette="Spectral") +
  theme_tree2()
colouredmccregional <- inset(colouredmccregional, mapgrob, width = 0.4, height = 0.4, vjust = -1, hjust = -0.6)
colouredmccregional

ggsave("EBOV_SierraLeone/Intermediate_Figures/MCC_shapecountry_colourregion_tree.pdf", plot = colouredmccregional, width = 10, height = 20, units = "in" )

# Pick clade in mcc tree and highlight
big_mcc_tree$node.label <- seq(1, big_mcc_tree$Nnode)
ggtree(big_mcc_tree) %<+% sequence_metadata_final +
  geom_tippoint(aes(color = Country)) +
  geom_nodelab(size = 3) +
  theme_tree2()

highlightedclade <- ggtree(big_mcc_tree) %<+% sequence_metadata_final +
  geom_hilight(node = "219", fill = "yellow", alpha = 0.2, type = "rect") +
  geom_tippoint(aes(color = Country)) +
  theme_tree2()
highlightedclade

ggsave("EBOV_SierraLeone/Intermediate_Figures/MCC_highlightedclade.pdf", plot = highlightedclade, width = 10, height = 20, units = "in" )


# Subset to the clade and filter out the GIN clusters
subset_sle_tree <- extract.clade(big_mcc_tree, "219", root.edge = T)
sierra_leone_sequence_metadata <- sequence_metadata_final %>%
  filter(tip %in% subset_sle_tree$tip.label) %>%
  filter(Country == "SLE")
write.csv(sierra_leone_sequence_metadata, "EBOV_SierraLeone/Clean_Data/SLE_tree_metadata.csv", row.names = F)
subset_sle_tree <- keep.tip(subset_sle_tree, sierra_leone_sequence_metadata$tip)
write.tree(subset_sle_tree, "EBOV_SierraLeone/Clean_Data/SLE_tree.tree")

subset_sle_tree$node.label <- seq(1, subset_sle_tree$Nnode)
ggtree(subset_sle_tree, mrsd = max(sierra_leone_sequence_metadata$Date)) %<+% sequence_metadata_final +
  geom_tippoint(aes(color = ADM1, shape = Country)) +
  geom_nodelab(size = 3) +
  #geom_tiplab(size = 3) +
  theme_tree2()

sierraleoneregional <- ggtree(subset_sle_tree, mrsd = max(sierra_leone_sequence_metadata$Date)) %<+% sequence_metadata_final +
  geom_tippoint(aes(color = Region, shape = Country)) +
  theme_tree2()
mapgrob <- list("266" = map)
sierraleoneregional <- inset(sierraleoneregional, mapgrob, width = 0.4, height = 0.4, vjust = 0, hjust = -1)
sierraleoneregional
ggsave("EBOV_SierraLeone/Intermediate_Figures/SLE_regionalandmap.pdf", plot = sierraleoneregional, width = 10, height = 20, units = "in" )


sierraleoneadmin1 <- ggtree(subset_sle_tree, mrsd = max(sierra_leone_sequence_metadata$Date)) %<+% sequence_metadata_final +
  geom_tippoint(aes(color = ADM1, shape = Country)) +
  scale_color_brewer("Region", palette="Spectral") +
  theme_tree2()
sierraleoneadmin1
ggsave("EBOV_SierraLeone/Intermediate_Figures/SLE_admin1.pdf", plot = sierraleoneadmin1, width = 10, height = 20, units = "in" )


# Make into EpiFusion friendly tree
root <- as.Date(max(sierra_leone_sequence_metadata$Date), format = "%Y-%m-%d") - (365*(max(nodeHeights(subset_sle_tree)) + subset_sle_tree$root.edge))
fulltreelength <- (365*(max(nodeHeights(subset_sle_tree)) + subset_sle_tree$root.edge))


SLE_epifusion_tree <- subset_sle_tree
distances <- get_all_distances_to_root(SLE_epifusion_tree) + rnorm(length(get_all_distances_to_root(SLE_epifusion_tree)), 0.001, 0.001) 
node_distances <- distances[(length(SLE_epifusion_tree$tip.label)+1):(length(SLE_epifusion_tree$tip.label)+SLE_epifusion_tree$Nnode)]
tip_distances <- distances[1:length(SLE_epifusion_tree$tip.label)] 
SLE_epifusion_tree$node.label <- paste0("X[", (node_distances + subset_sle_tree$root.edge)*365, "]")
SLE_epifusion_tree$tip.label <- paste0("X[", (tip_distances + subset_sle_tree$root.edge)*365, "]")


ggtree(SLE_epifusion_tree) +
  geom_tiplab(size = 1) +
  geom_nodelab(size = 1)

test_downsampled_tree <- keep.tip(SLE_epifusion_tree, sample(SLE_epifusion_tree$tip.label, round(0.5*length(SLE_epifusion_tree$tip.label))))

write.tree(SLE_epifusion_tree, "EBOV_SierraLeone/Clean_Data/SLE_epifusiontree.tree")
write.tree(test_downsampled_tree, "EBOV_SierraLeone/Clean_Data/downsampled_SLE_epifusiontree.tree")

# Get some tree analytics 
distances <- get_all_distances_to_root(SLE_epifusion_tree)
births <- data.frame(event = "birth", time = (distances[(length(SLE_epifusion_tree$tip.label)+1):(length(SLE_epifusion_tree$tip.label)+SLE_epifusion_tree$Nnode)]+ subset_sle_tree$root.edge)*365)
samplings <- data.frame(event = "sampling", time = (distances[1:length(SLE_epifusion_tree$tip.label)] + rnorm(length(SLE_epifusion_tree$tip.label), 0.001, 0.001) + subset_sle_tree$root.edge)*365)
events <- rbind(births, samplings) %>%
  mutate(day = ceiling(time)) %>%
  group_by(day, event) %>%
  mutate(count = n()) %>%
  select(day, event, count) %>%
  distinct()

treeevents <- ggplot(events, aes(x = day)) +
  geom_bar(aes(y = count, col = event, fill = event), stat = "identity", position = "stack")
ggsave("EBOV_SierraLeone/Intermediate_Figures/SLE_tree_eventbarplot.pdf", treeevents, width = 8, height = 6, units = c("in"))

weeklysequences <- events %>%
  ungroup() %>%
  #select(day, event, count) %>%
  mutate(Week = ceiling(day/7)) %>%
  filter(event == "sampling") %>%
  group_by(Week) %>%
  mutate(total_seqs = sum(count)) %>%
  select(Week,  total_seqs) %>%
  distinct()


##### Incidence Data #####
caselinelist <- read.csv("EBOV_SierraLeone/Raw_Data/pnas.1518587113.sd02.csv") %>%
  mutate(Date = as.Date(Date.of.symptom.onset, format = "%d-%b-%y")) %>%
  filter(Date < max(sierra_leone_sequence_metadata$Date)) %>%
  mutate(Week = ceiling(as.numeric(Date - root)/7)) %>%
  group_by(District, Week) %>%
  mutate(Cases = n()) %>%
  dplyr::select(Week, District, Cases) %>%
  mutate(District = str_remove(District, " ")) %>%
  mutate(District = factor(District, levels = regionlabels_tree)) %>%
  distinct() %>%
  left_join(mutate(region_lookup, District = Region)) %>%
  ungroup()
write.csv(caselinelist, "EBOV_SierraLeone/Clean_Data/tidy_regional_incidence.csv", row.names = F)


regionalcases <- ggplot(caselinelist, aes(x = Week)) +
  geom_bar(aes(y = Cases, col = District, fill = District), position = "stack", stat = "identity")
ggsave("EBOV_SierraLeone/Intermediate_Figures/SLE_regional_weekly_cases.pdf", regionalcases, width = 8, height = 6, units = c("in"))


suspected_caselinelist <- read.csv("EBOV_SierraLeone/Raw_Data/pnas.1518587113.sd02_suspected.csv") %>%
  mutate(Date = as.Date(Date.of.symptom.onset, format = "%d-%b-%y")) %>%
  filter(Date < max(sierra_leone_sequence_metadata$Date)) %>%
  mutate(Week = ceiling(as.numeric(Date - root)/7)) %>%
  group_by(District, Week) %>%
  mutate(Cases = n()) %>%
  dplyr::select(Week, District, Cases) %>%
  rename(Suspected_Cases = Cases) %>%
  mutate(District = str_remove(District, " ")) %>%
  mutate(District = factor(District, levels = regionlabels_tree)) %>%
  distinct() %>%
  ungroup()

suspectedregionalcases <- ggplot(suspected_caselinelist, aes(x = Week)) +
  geom_bar(aes(y = Suspected_Cases, col = District, fill = District), position = "stack", stat = "identity")
ggsave("EBOV_SierraLeone/Intermediate_Figures/SLE_suspected_regional_weekly_cases.pdf", suspectedregionalcases, width = 8, height = 6, units = c("in"))


combined <- full_join(caselinelist, suspected_caselinelist) %>%
  mutate(Suspected_Cases = ifelse(is.na(Suspected_Cases), 0, Suspected_Cases)) %>%
  mutate(Cases = ifelse(is.na(Cases), 0, Cases)) %>%
  group_by(Week) %>%
  mutate(Cases = sum(Cases), Suspected_Cases = sum(Suspected_Cases)) %>%
  dplyr::select(Week, Cases, Suspected_Cases) %>%
  distinct() %>%
  ungroup() %>%
  mutate(Total = Cases + Suspected_Cases) %>%
  mutate(Proportion = Cases/Total)

proportion_confirmed_vs_suspected <- ggplot(combined, aes(x = Week, y = Proportion)) +
  geom_point() +
  geom_smooth()


total_incidence <- caselinelist %>%
  group_by(Week) %>%
  mutate(Total = sum(Cases)) %>%
  select(Week, Total) %>%
  ungroup() %>%
  distinct()

write.csv(total_incidence, "EBOV_SierraLeone/Clean_Data/tidy_national_incidence.csv", row.names = F)

casesandsequences <- total_incidence %>%
  full_join(weeklysequences) %>%
  mutate(total_seqs = ifelse(is.na(total_seqs), 0, total_seqs)) %>%
  mutate(sample_prop = total_seqs/Total)

ggplot(casesandsequences, aes(x = Week, y = sample_prop)) +
  geom_point() +
  geom_smooth()



paste0(total_incidence$Week*7, collapse = " ")
paste0(total_incidence$Total, collapse = " ")

sum(total_incidence$Total)
length(SLE_epifusion_tree$tip.label)


plot(total_incidence$Week, total_incidence$Total)


dailyconfirmedcases <- read.csv("EBOV_SierraLeone/Raw_Data/pnas.1518587113.sd02.csv") %>%
  mutate(Date = as.Date(Date.of.symptom.onset, format = "%d-%b-%y")) %>%
  filter(Date < max(sierra_leone_sequence_metadata$Date)) %>%
  mutate(Day = ceiling(as.numeric(Date - root))) %>%
  group_by(Day) %>%
  mutate(Cases = n()) %>%
  dplyr::select(Day, Cases) %>%
  distinct() %>%
  ungroup()

dailysuspectedcases <- read.csv("EBOV_SierraLeone/Raw_Data/pnas.1518587113.sd02_suspected.csv") %>%
  mutate(Date = as.Date(Date.of.symptom.onset, format = "%d-%b-%y")) %>%
  filter(Date < max(sierra_leone_sequence_metadata$Date)) %>%
  mutate(Day = ceiling(as.numeric(Date - root))) %>%
  group_by(Day) %>%
  mutate(SuspectedCases = n()) %>%
  dplyr::select(Day, SuspectedCases) %>%
  distinct() %>%
  ungroup()

tree <- read.tree("EBOV_SierraLeone/Clean_Data/SLE_epifusiontree.tree")
tree$edge.length <- tree$edge.length*365
tree$root.edge <- tree$root.edge*365
lineagesthroughtime <- count_lineages_through_time(tree, times = seq(0, 550, 1))
lineages <- data.frame(Day = lineagesthroughtime$times+round(tree$root.edge), Lineages = lineagesthroughtime$lineages)

masterdailyevents <- events %>%
  pivot_wider(names_from = event, values_from = count) %>%
  full_join(data.frame(day = seq(1, max(events$day)))) %>%
  arrange(day) %>%
  mutate(Date = root + day) %>%
  transmute(Date = Date, Day = day, BirthEvents = birth, Samplings = sampling) %>%
  left_join(dailyconfirmedcases) %>%
  left_join(dailysuspectedcases) %>%
  left_join(lineages) %>%
  replace_na(list(BirthEvents = 0, Samplings = 0, Cases = 0, SuspectedCases = 0)) %>%
  mutate(TotalCases = Cases + SuspectedCases) %>%
  mutate(PropSequencedTotal = Samplings / TotalCases)

write.csv(masterdailyevents, "EBOV_SierraLeone/Clean_Data/master_daily_events.csv", row.names = F)

masterweeklyevents <- masterdailyevents %>%
  ungroup() %>%
  mutate(Week = ceiling(Day/7)) %>%
  dplyr::select(!Day) %>%
  summarise(.by = Week, BirthEvents = sum(BirthEvents), 
            Samplings = sum(Samplings), Cases = sum(Cases), 
            SuspectedCases = sum(SuspectedCases),
            TotalCases = sum(TotalCases),
            Lineages = mean(Lineages))
  
ggplot(masterdailyevents, aes(x = Day)) +
  geom_point(aes(y = PropSequencedTotal)) 




