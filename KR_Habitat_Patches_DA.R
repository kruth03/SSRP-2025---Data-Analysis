# Habitat Patches Productivity and Diversity Data Analysis 
# Created by Kaylee Ruth 

# **************************************************************************
# R FILE PREPARATION ####
# **************************************************************************

# Load packages
library(tidyr)
library(ggplot2)
library(ggstatsplot)
library(vegan)
library(FD)
library(dplyr)
library(myClim)
library(plyr)
library(sjstats)
library(afex)
library(ggfortify)
library(grid)


# Clear R's memory
rm(list = ls())

# Set working directory (remember to update filepath)
setwd("/Users/kayleeruth/Desktop/R Bootcamp 2025/Habitat Patches Data Analysis")

# Load data
habitat_patches <- read.csv("KR-D-027 KR_data.csv", na.strings = c("", "NA"), header = TRUE)

# **************************************************************************
# INSPECT AND CLEAN THE DATA SET ####
# **************************************************************************

str(habitat_patches)
View(habitat_patches)

habitat_patches = data.frame(habitat_patches)


# ***************************************************************************
# CALCULATE MESOCOSM PRODUCTIVITY
# ***************************************************************************

# In order to calculate mesocosm productivity, need to sum the biomass of all species
# this requires skipping columns 1-5 within the dataset 

habitat_patches$Productivity = rowSums(habitat_patches[, -c(1,2,3,4,5)]) # this skips the irrelevant columns in the sum
View(habitat_patches) # view the data frame to check if the productivity column was generated correctly.


# Plot Total Mesocosm productivity vs mesocosm number  -- let's put 
ggplot(data = habitat_patches, aes(x = Code, y = Productivity)) +
  geom_point(alpha = 1) + 
      labs(x = "Block", y = "Productivity (g)") +
      theme_bw(base_size = 16)
    
# -- Raphanus raphanistrum sample has a biomass of 49 g
# I need to actually measure productivity vs deltaT -- I have temp and precip normals for each site but
# need to calculate the difference in field site temp vs experimental soil temp


# ****************************************************************************
# Average Temperature and precipitation of Carnegie Science hot/dry conditions
# ****************************************************************************
tms <- mc_read_files(c("data_95322909_2025_07_24_0.csv", "data_95322908_2025_07_24_0.csv", "data_95322902_2025_07_24_0.csv"), dataformat_name = "TOMST", silent = TRUE)

# crop the time-series to dates of interest (some of these were sitting in the lab, we don't care about those)
start <- as.POSIXct("2024-10-01", tz = "UTC") # as.POSIXct is a standard date format.
end <- as.POSIXct("2025-07-24", tz = "UTC")

tms <- mc_prep_crop(tms, start, end) # this crops tms data by start/end dates we defined above
tms.calc <- mc_calc_vwc(tms, soiltype = "loam") # calculate soil moisture
head(tms.calc)

tms.long <- mc_reshape_long(tms.calc, sensors = c("TMS_T3"))
View(tms.long)
# now want to extract the temperature of the air -- this should be from the sensor furthest from the soil.
# just checked the TOMST website, TMS_T3 = air temperature -- need to extract this

# Main issue: trying to extract T3 and VWC data into a new data frame -- it worked, tms.long = df
# now, I can calculate the average temperature over all sensors (all T3)

avg.temp = mean(tms.long$value)
# we have an average temp for the entire time period now = 16.66479 C (approx 62.6 F)
# can calculate ΔT values with this

# I received precipitation data from Trevor Herbert (thank you!) at Jasper Ridge:
# jasper_ridge <- read.csv("JR_Weather_Oct2024-July2025.csv", na.strings = c("", "NA"), header = TRUE) # import weather data
# View(jasper_ridge) # column X.3 is the rainfall amount -- so take the average of this
# avg.precip = sum(jasper_ridge$X.3[-1]) # error -- saying that this column is not numeric, so convert X.3 to numeric

# jasper_ridge$X.3 = as.numeric(jasper_ridge$X.3) # convert precip column to all numeric values (1st row = NA)
# jasper_ridge_cleaned = jasper_ridge[!is.na(jasper_ridge$X.3), ] # filter out the rows with NA (this is 1st row)

# total.precip = sum(jasper_ridge_cleaned$X.3) # take sum of precip column = total experimental precip.

# total.precip = total.precip / 10 # convert to cm (currently in mm)
# we now have total precip for the time period in cm = 51.7 cm rainfall
# this is very odd, about 20 cm more rain than expected. Andrea gave me 
# Stanford's weather data: 34.8488 cm to use for our total.precip instead

total.precip = 34.8488

# ******************************************************************************
# Calculate Δprecip and ΔT values using avg experimental precip, T, and field site 
# T and precip values
# ******************************************************************************

# first, turn our different field sites from the Habitat Patches ds into factors 
habitat_patches$Source <- as.factor(habitat_patches$Source)

# now, let's calculate Δprecip and ΔT values

# ΔT = T(experimental site) - T(field site)
habitat_patches$delta.T = avg.temp - habitat_patches$temp_C

# Δprecip = precip(experimental site) - precip(field site)
habitat_patches$delta.precip = total.precip - habitat_patches$ppt_cm 


# ******************************************************************************
# Calculate Shannon and Simpson Diversity Indices for the Habitat Patches ds
# NOTE: The Shannon index is more sensitive to species richness (number of species) and the evenness of their distribution, 
# while the Simpson index is more influenced by the dominant species
# ******************************************************************************
habitat_patches$Shannon <- diversity(habitat_patches[,-c(1,2,3,4,5,77,78,79)], index = "shannon")

# Calculate Simpson Diversity Index
habitat_patches$Simpson <- diversity(habitat_patches[,-c(1,2,3,4,5,77,78,79)], index = "simpson")
View(habitat_patches) # it worked... asked Andrea for help/checking my data.



# ******************************************************************************
# Statistical Analyses -- GLMM of Productivity
# ******************************************************************************

# First, ensure that Source and Block (Code) are factors 
habitat_patches$Code = as.factor(habitat_patches$Code)
habitat_patches$delta.precip = as.numeric(habitat_patches$delta.precip)
habitat_patches$delta.T = as.numeric(habitat_patches$delta.T)

# try to make delta.T and delta.precip numeric... what happens?


# 1) Linear Mixed Model Analysis of Productivity

  # I got an error saying: number of levels of each grouping factor < number of obs
    # in Code, we are sorting by mesocosm. Let's sort by Block: 
habitat_patches$Block <- sub("_.*", "", habitat_patches$Code)
View(habitat_patches) # now we have the rows sorted by block, try LMM again:
# also got a singularity error (coming from Source being treated as a )

# also tried using block, source as random effects. drop block completely, only use
# source as a random effects variable:
prod_model_1 <- mixed(Productivity ~ delta.T + delta.precip + (1 | Source),
                      data = habitat_patches, method = "KR",
                      control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

prod_model_assump <- lme4::lmer(Productivity ~ delta.T + delta.precip + (1 | Source),
                           data = habitat_patches)

plot(prod_model_assump)

# Check the normality of the residuals assumption 
qqnorm(residuals(prod_model_assump)) 
# Normal Q-Q Plot is not very linear, log transform data
anova(prod_model_1)

### LOG transform analysis with same LMM: ###

prod_model_log <- mixed(log(Productivity) ~ delta.T + delta.precip  +  (1 | Source),
                      data = habitat_patches, method = "KR",
                      control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

prod_model_log_assump <- lme4::lmer(log(Productivity) ~ delta.T + delta.precip + (1 | Source),
                                data = habitat_patches)

# Check the normality of the residuals assumption 
qqnorm(residuals(prod_model_log_assump))  # slightly more linear? 

hb_biomass_R2 <- piecewiseSEM::rsquared(prod_model_log_assump, method = nagelkerke)
hb_biomass_R2

# Calculate the F, df, and P value for each independent factor
anova(prod_model_log) # seems that p value of log(delta.precip) is significant?

# I reran the code and I am not seeing significance anymore-- I must have 
# overwritten values while I was coding this past week.

# ******************************************************************************
# Statistical Analyses -- GLMM of Shannon and Simpson Diversity Indices
# ******************************************************************************

# GLMM for Shannon Diversity Index 

shannon_model_1 = mixed(Shannon ~ delta.T + delta.precip + (1 | Source),
                                        data = habitat_patches, method = "KR",
                                        control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

shannon_model_assump <- lme4::lmer(Shannon ~ delta.T + delta.precip + (1 | Source),
                                data = habitat_patches)

# Q-Q Plot of Residuals to check for Normality
qqnorm(residuals(shannon_model_assump)) # slope is nearly unity, this looks good

plot(shannon_model_assump)

hp_shannon_R2 <- piecewiseSEM::rsquared(shannon_model_assump, method = nagelkerke)
hp_shannon_R2 

# Calculate the F, df, and P value for each independent factor using anova
anova(shannon_model_1)

# GLMM for Simpson Diversity Index
simpson_model_1 = mixed(Simpson ~ delta.T + delta.precip + (1 | Source),
                        data = habitat_patches, method = "KR",
                        control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

simpson_model_assump <- lme4::lmer(Simpson ~ delta.T + delta.precip + (1 | Source),
                                   data = habitat_patches)

hp_simpson_R2 <- piecewiseSEM::rsquared(simpson_model_assump, method = nagelkerke)
hp_simpson_R2

# Q-Q Plot of Residuals to check for Normality
qqnorm(residuals(simpson_model_assump))

# Calculate the F, df, and P value for each independent factor
anova(simpson_model_1)


# ******************************************************************************
# Statistical Analyses -- GLMM of Functional Traits Data (Height & SLA)
# ******************************************************************************

# import trait data sent from Andrea Nebhut: 
hp_traits <- read.csv("trait_cwm_KR.csv", na.strings = c("", "NA"), header = TRUE)
View(hp_traits)

# add a column for Block
hp_traits$Block <- sub("_.*", "", hp_traits$Code)
View(hp_traits)

# add a column for temperature and precipitation changes experienced by the block
# Δprecipitation
hp_traits <- hp_traits %>%
  left_join(habitat_patches %>% select(Code, delta.precip), by = "Code")
View(hp_traits)

# ΔT
hp_traits <- hp_traits %>%
  left_join(habitat_patches %>% select(Code, delta.T), by = "Code")
View(hp_traits)

# add a column for Source of community
hp_traits <- hp_traits %>%
  left_join(habitat_patches %>% select(Code, Source), by = "Code")
View(hp_traits)

# Set Source as a factor for GLMM analysis later
hp_traits$Source <- as.factor(hp_traits$Source)

###### GLMMs #######

# GLMM for CWM Plant Height 

height_model_1 = mixed(log(heightveg_m_CWM) ~ delta.T + delta.precip + (1 | Source),
                        data = hp_traits, method = "KR",
                        control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

height_model_assump <- lme4::lmer(log(heightveg_m_CWM) ~ delta.T + delta.precip + (1 | Source),
                                   data = hp_traits)


# Q-Q Plot of Residuals to check for Normality
qqnorm(residuals(height_model_assump)) # slope is nearly unity, this looks good

plot(height_model_assump)

hp_height_R2 <- piecewiseSEM::rsquared(height_model_assump, method = nagelkerke)
hp_height_R2 


# Calculate the F, df, and P value for each independent factor using anova
anova(height_model_1)


# GLMM for CWM SLA 
SLA_model_1 = mixed(log(SLA_mm2.mg_CWM) ~ delta.T + delta.precip + (1 | Source),
                       data = hp_traits, method = "KR",
                       control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

SLA_model_assump <- lme4::lmer(log(SLA_mm2.mg_CWM) ~ delta.T + delta.precip + (1 | Source),
                                  data = hp_traits)


# Q-Q Plot of Residuals to check for Normality
qqnorm(residuals(SLA_model_assump)) # slope is nearly unity, this looks good

plot(SLA_model_assump)

hp_SLA_R2 <- piecewiseSEM::rsquared(SLA_model_assump, method = nagelkerke)
hp_SLA_R2 

# Calculate the F, df, and P value for each independent factor using anova
anova(SLA_model_1)


### GLMM for CWM seed mass 
seedmass_model_1 = mixed(log(seedmass_mg_CWM) ~ delta.T + delta.precip + (1 | Source),
                    data = hp_traits, method = "KR",
                    control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)

seedmass_model_assump <- lme4::lmer(log(seedmass_mg_CWM) ~ delta.T + delta.precip + (1 | Source),
                               data = hp_traits)


# Q-Q Plot of Residuals to check for Normality
qqnorm(residuals(seedmass_model_assump)) # slope is nearly unity, this looks good

plot(seedmass_model_assump)

hp_seedmass_R2 <- piecewiseSEM::rsquared(seedmass_model_assump, method = nagelkerke)
hp_seedmass_R2 

# Calculate the F, df, and P value for each independent factor using anova
anova(seedmass_model_1)


# ******************************************************************************
# Statistical Analyses -- Community Functional Traits NMDS
# ******************************************************************************

# Import ggrepel() for arrow labels 

library(ggrepel)

# get rid of NA entries: 
traits_nmds <- hp_traits[c(2,4,6)] # only look at height, SLA, seed mass
hp_traits_complete <- hp_traits[complete.cases(traits_nmds), ]

# Perform a NMDS on the trait data
hp.mds <- metaMDS(hp_traits_complete[c(2,4,6)])
hp.data.scores <- as.data.frame(scores(hp.mds)$sites)
hp.data.scores$Label <- hp_traits_complete$Source
# adonis2(traits[c(2:13)] ~ traits$FctGrp_Origin, data = traits)$aov.tab # adonis test -- are the shapes we see significant 
# figure out how to do adonis test on google/chatgpt to get p value

hp.data.scores_traits <- as.data.frame(scores(hp.mds)$species)
hp.data.scores_traits <- cbind(hp.data.scores_traits, trait = rownames(hp.data.scores_traits))
hp.data.scores <- cbind(hp.data.scores, Source = hp_traits_complete$Source)

# Function to sort points by convex hull # founds outside points and draws lines appropriately
sort_by_hull <- function(df) {
  hull_indices <- chull(df$NMDS1, df$NMDS2)
  df[hull_indices, ]
}

# Apply the sorting function to each group
sorted_data <- hp.data.scores %>%
  group_by(Label) %>%
  do(sort_by_hull(.))

# add columns for delta.T and delta.precip that correspond to Source to sorted_data
hp_traits_unique <- hp_traits %>%
  distinct(Source, .keep_all = TRUE) # prevent duplicate row entries being made

# delta.precip column creation in sorted_data 
sorted_data <- sorted_data %>%
  left_join(hp_traits_unique %>% select(Source, delta.precip), by = "Source" )
View(sorted_data)

# delta.T column creation in sorted_data
sorted_data <- sorted_data %>%
  left_join(hp_traits_unique %>% select(Source, delta.T), by = "Source")
View(sorted_data)

# Run adonis2() on the NMDS data to get a p-value for functional traits: 
adonis_result_traits <- adonis2(hp_traits_complete[c(2,4,6)] ~ Source, data = hp_traits_complete, method = "bray")
print(adonis_result_traits)
# Hey! The p-value is 0.02 * 

adonis_result_traits_temp <- adonis2(hp_traits_complete[c(2,4,6)] ~ delta.T, data = hp_traits_complete, method = "bray")
print(adonis_result_traits_temp)

adonis_result_traits_precip <- adonis2(hp_traits_complete[c(2,4,6)] ~ delta.precip, data = hp_traits_complete, method = "bray")
print(adonis_result_traits_precip)

# rename hp.data.scores_traits trait column to give sensical names to the functional traits 
hp.data.scores_traits$trait <- recode(hp.data.scores_traits$trait,
                              "heightveg_m_CWM" = "CWM Height (m)",
                              "SLA_mm2.mg_CWM" = "CWM SLA (mm²/mg)",
                              "seedmass_mg_CWM" = "CWM Seed Mass (g)"
)

# Plot results on a color gradient based upon delta.T
temp_func <- ggplot() +
  geom_point(data = sorted_data, aes(x = NMDS1, y = NMDS2, color = delta.T), size = 3) +
  geom_segment(
    data = hp.data.scores_traits,
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    arrow = arrow(length = unit(0.25, "cm")), colour = "grey30"
  ) +
  geom_text_repel(
    data = hp.data.scores_traits,
    aes(x = NMDS1, y = NMDS2, label = trait),
    size = 15, color = "black",
    max.overlaps = Inf
  ) + 
  scale_color_gradientn(colors = c("#ff0000" ,"#ff4d00","#ff7400", "#ff9a00", "#ffc100")) +
  theme_bw(base_size = 60) + 
  theme(
    axis.text = element_text(size = 40)  # changes both x and y axis numbers
  ) + 
  theme(legend.position = "bottom") +
  annotation_custom(grobTree(textGrob("p = 0.09",
    gp = gpar(fontsize = 60),
    x = 0.15,
    y = 0.90
)))

  ggsave(
    filename = "NMDS_FunctionalTraits_Temperature_final.png",  # Specify the output filename and desired format
    plot = temp_func,       # Save the last displayed ggplot object (or specify a plot object)
    width = 10,                # Set the desired width
    height = 10,               # Set the desired height
    units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
    dpi = 300                 # Set the desired resolution in dots per inch
  )   
  
# Why low = orange, high = blue? A larger delta.T = cooler source community,
# smaller delta.T = hotter source community

# Plot results on a color gradient based upon delta.precip
precip_func <- ggplot() +
  geom_point(data = sorted_data, aes(x = NMDS1, y = NMDS2, color = delta.precip), size = 3) +
  geom_segment(
    data = hp.data.scores_traits,
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    arrow = arrow(length = unit(0.25, "cm")), colour = "grey30"
  ) +
  geom_text_repel(
    data = hp.data.scores_traits,
    aes(x = NMDS1, y = NMDS2, label = trait),
    size = 15, color = "black",
    max.overlaps = Inf
  ) + 
  scale_color_gradientn(colors = c("#1d20d0","#1348a3", "#2361a7", "#3279ae", "lightblue")) +
  theme_bw(base_size = 60) +
  theme(
    axis.text = element_text(size = 40)  # changes both x and y axis numbers
  ) + 
  theme(legend.position = "bottom") +
  annotation_custom(grobTree(textGrob("p = 0.78",
    gp = gpar(fontsize = 60),
    x = 0.15,
    y = 0.90
  )))

  
ggsave(
  filename = "NMDS_FunctionalTraits_Precipitation_final.png",  # Specify the output filename and desired format
  plot = precip_func,       # Save the last displayed ggplot object (or specify a plot object)
  width = 10,                # Set the desired width
  height = 10,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
  )  


# Plot results on a color blind-friendly color scheme based upon Source.
func_traits = ggplot() +
  geom_point(data = sorted_data, aes(x = NMDS1, y = NMDS2, color = Source), size = 3) +
  geom_segment(
    data = hp.data.scores_traits,
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    arrow = arrow(length = unit(0.25, "cm")), colour = "grey30"
  ) +
  geom_text_repel(
    data = hp.data.scores_traits,
    aes(x = NMDS1, y = NMDS2, label = trait),
    size = 15, color = "black",
    max.overlaps = Inf
  ) + 
  scale_color_manual(values = c("#0072b2","#e69f00", "#56b4e9", "#009e73", "#9C9C9C", "#cc79a7","#000000")) +
  theme_bw(base_size = 60) +
  theme(
    axis.text = element_text(size = 40)  # changes both x and y axis numbers
  ) + 
  theme(legend.position = "bottom") +
  annotation_custom(grobTree(textGrob("p = 0.02 *",
    gp = gpar(fontsize = 60),
    x = 0.10,
    y = 0.90
   
  ))) 
  

ggsave(
  filename = "NMDS_FunctionalTraits_Source_final.png",  # Specify the output filename and desired format
  plot = func_traits,       # Save the last displayed ggplot object (or specify a plot object)
  width = 12,                # Set the desired width
  height = 10,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  

# Why low = blue, high = orange? A larger delta.precip = drier source community,
# smaller delta.precip = wetter source community 
  
# ******************************************************************************
# Community Composition NMDS
# ******************************************************************************

# using same code as above for the functional traits, run an NMDS on the 
# community composition data (which is using the habitat_patches dataset)

# filter out columns of zeros 
habitat_patches_filtered <- habitat_patches[, colSums(habitat_patches != 0) > 0]
View(habitat_patches_filtered) # check manually to ensure code executed correctly 


# Perform a NMDS on the composition data -- only select species present in the blocks
comp.mds <- metaMDS(habitat_patches_filtered[-c(1,2,3,4,5,56,57,58,59,60,61)]) 
comp.data.scores <- as.data.frame(scores(comp.mds)$sites)
comp.data.scores$Label <- habitat_patches_filtered$Source
  # adonis2(traits[c(2:13)] ~ traits$FctGrp_Origin, data = traits)$aov.tab # adonis test -- are the shapes we see significant 
  # figure out how to do adonis test on google/chatgpt to get p value
  
  comp.data.scores_traits <- as.data.frame(scores(comp.mds)$species)
  comp.data.scores_traits <- cbind(comp.data.scores_traits, trait = rownames(comp.data.scores_traits))
  comp.data.scores <- cbind(comp.data.scores, Source = habitat_patches_filtered$Source)
  
  # Function to sort points by convex hull # founds outside points and draws lines appropriately
sort_by_hull <- function(df) {
  hull_indices <- chull(df$NMDS1, df$NMDS2)
  df[hull_indices, ]
}
  
  # Apply the sorting function to each group
comp_sorted_data <- comp.data.scores %>%
  group_by(Label) %>%
  do(sort_by_hull(.))
  
# add columns for delta.T and delta.precip that correspond to Source to sorted_data
habitat_patches_unique <- habitat_patches %>%
  distinct(Source, .keep_all = TRUE) # prevent duplicate row entries being made

# delta.precip column creation in sorted_data 
comp_sorted_data <- comp_sorted_data %>%
  left_join(habitat_patches_unique %>% select(Source, delta.precip), by = "Source" )
View(comp_sorted_data)

# delta.T column creation in sorted_data
comp_sorted_data <- comp_sorted_data %>%
  left_join(habitat_patches_unique %>% select(Source, delta.T), by = "Source")
View(sorted_data)

# Run adonis2() on the NMDS data to get a p-value for community composition: 
adonis_result_comp <- adonis2(habitat_patches_filtered[-c(1,2,3,4,5,56,57,58,59,60,61)] 
                      ~ Source, data = habitat_patches_filtered, method = "bray")
print(adonis_result_comp)
# Hey! The p-value is 0.001 *** 

# Plot results on a color gradient based upon by delta.T
ggplot() +
  geom_point(data = comp_sorted_data, aes(x = NMDS1, y = NMDS2, color = delta.T), size = 3) +
  # geom_segment(
  #   data = comp.data.scores_traits,
  #   aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
  #   arrow = arrow(length = unit(0.25, "cm")), colour = "grey30"
  # ) +
  # geom_text_repel(
  #   data = comp.data.scores_traits,
  #   aes(x = NMDS1, y = NMDS2, label = trait),
  #   size = 3, color = "black",
  #   max.overlaps = Inf
  # ) + 
  scale_color_gradientn(colors = c("#ff0000" ,"#ff4d00","#ff7400", "#ff9a00", "#ffc100")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

  # annotation_custom(grobTree(textGrob("p = 0.001 ***",
    # gp = gpar(fontsize = 16),
    # x = 0.90,
    # y = 0.95
    # )))

# Why low = darker orange, high = lighter orange? A larger delta.T = cooler source community,
# smaller delta.T = hotter source community

# Plot results on a color gradient based upon delta.precip
ggplot() +
  geom_point(data = comp_sorted_data, aes(x = NMDS1, y = NMDS2, color = delta.precip), size = 3) +
  # geom_segment(
  #   data = comp.data.scores_traits,
  #   aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
  #   arrow = arrow(length = unit(0.25, "cm")), colour = "grey30"
  # ) +
  # geom_text_repel(
  #   data = comp.data.scores_traits,
  #   aes(x = NMDS1, y = NMDS2, label = trait),
  #   size = 3, color = "black",
  #   max.overlaps = Inf
  # ) + 
  scale_color_gradientn(colors = c("#1d20d0","#1348a3", "#2361a7", "#3279ae", "lightblue")) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom")
# annotation_custom(grobTree(textGrob("p = 0.001 ***",
  # gp = gpar(fontsize = 16),
  # x = 0.90,
  # y = 0.95
# )))

# Why low = blue, high = orange? A larger delta.precip = drier source community,
# smaller delta.precip = wetter source community 

# Plot results using a color palette (colorblind friendly) based upon Source.
nmds_comp = ggplot() +
  geom_point(data = comp_sorted_data, aes(x = NMDS1, y = NMDS2, color = Source), size = 3) +
  scale_color_manual(values = c("#0072b2","#e69f00", "#56b4e9", "#009e73", "#9C9C9C", "#cc79a7","#000000")) +
  theme_bw(base_size = 60) +
  theme(legend.position = "bottom") + 
  annotation_custom(grobTree(textGrob("p = 0.001 **",
    gp = gpar(fontsize = 60),
    x = 0.85,
    y = 0.92
  ))) + 
  theme(
    axis.text = element_text(size = 40)  # changes both x and y axis numbers
  )
  

ggsave(
  filename = "NMDS_Composition_Source_final.png",  # Specify the output filename and desired format
  plot = nmds_comp,       # Save the last displayed ggplot object (or specify a plot object)
  width = 12,                # Set the desired width
  height = 10,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  
###################################### Notes ###################################

# planning to use linear mixed model for statistical analysis

#Serial numbers: M = 95322903, Q = 95322907, O = 95322902 (hot block), I = 95322908 (hot block), G = 95322906, A = 95322909 (hot block)
# Average temperatures from A, I, O and that is our experimental temperature
# need to ensure that the dates match up with the temp normals for field sites = Oct 1 - July 25th
# current dates for TMS4 data: 2024-03-18 to 2025-07-24


# Question: Should I make the code(block), source(field site), factors? 
# could need to make source fixed factor, but for now it is random effects 
# start with maximal number of terms (max model) then drop terms if it does not 
# converge
# first response variable -- productivity 
# productivity ~ deltaT + deltaP + Random effects = block, field site 
# diversity (same as productivity)
# composition = NMDS 

# ******************************************************************************
# Plots
# ******************************************************************************
# Productivity vs ΔT
p1 = ggplot(data = habitat_patches, aes(x = delta.T, y = log(Productivity), 
                                   color = delta.T)) +
  geom_point(alpha = 1,, size = 3) + 
  scale_color_gradientn(colors = c("#ff0000" ,"#ff4d00","#ff7400", "#ff9a00", "#ffc100")) +
  theme_bw(base_size = 55) +
  # theme(legend.position = "bottom") +
  labs(x = "ΔTemperature (experimental - source; °C)", y = "Productivity (log(g))") + 
  theme(legend.position = "none") + 
  theme(
  axis.text = element_text(size = 35)  # changes both x and y axis numbers
)
  

ggsave(
  filename = "Productivity_vs_deltaT_final.png",  # Specify the output filename and desired format
  plot = p1,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  




# Productivity vs Δprecip
p2 = ggplot(data = habitat_patches, aes(x = delta.precip, y = log(Productivity), color = delta.precip)) +
  geom_point(alpha = 1, size = 3) +
  scale_color_gradientn(colors = c("#1d20d0","#1348a3", "#2361a7", "#3279ae", "lightblue")) +
  theme_bw(base_size = 55) +
  labs(x = "ΔPrecipitation (experimental - source; cm)", y = "Productivity (log(g))") +
  theme(legend.position = "none") + 
  theme(
  axis.text = element_text(size = 35)  # changes both x and y axis numbers
)



ggsave(
  filename = "Productivity_vs_deltaPrecip_final.png",  # Specify the output filename and desired format
  plot = p2,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  


# Shannon Diversity Index vs ΔT
p3 = ggplot(data = habitat_patches, aes(x = delta.T, y = Shannon, 
                                  color = delta.T)) +
  geom_point(alpha = 1,, size = 3) + 
  scale_color_gradientn(colors = c("#ff0000" ,"#ff4d00","#ff7400", "#ff9a00", "#ffc100")) +
  labs(x = "ΔTemperature (experimental - source; °C)", y = "Shannon Index") + 
  theme_bw(base_size = 55) + 
  theme(legend.position = "none") + 
  theme(
    axis.text = element_text(size = 35)  # changes both x and y axis numbers
  )

ggsave(
  filename = "Shannon_vs_deltaT_final.png",  # Specify the output filename and desired format
  plot = p3,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  
# Shannon Diversity index vs Δprecip
p4 = ggplot(data = habitat_patches, aes(x = delta.precip, y = Shannon, color = delta.precip)) +
  geom_point(alpha = 1, size = 3) + 
  scale_color_gradientn(colors = c("#1d20d0","#1348a3", "#2361a7", "#3279ae", "lightblue")) +
  labs(x = "ΔPrecipitation (experimental - source; cm)", y = "Shannon Index") +
  theme_bw(base_size = 55)  + 
  theme(legend.position = "none") + 
  theme(
    axis.text = element_text(size = 35)  # changes both x and y axis numbers
  )

ggsave(
  filename = "Shannon_vs_deltaPrecip_final.png",  # Specify the output filename and desired format
  plot = p4,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  

# Simpson Diversity index vs ΔT
p5 = ggplot(data = habitat_patches, aes(x = delta.T, y = Simpson, 
                                    color = delta.T)) +
  geom_point(alpha = 1, size = 3) + 
  scale_color_gradientn(colors = c("#ff0000" ,"#ff4d00","#ff7400", "#ff9a00", "#ffc100")) +
  theme_bw(base_size = 55) +
  labs(x = "ΔTemperature (experimental - source; °C)", y = "Simpson Index")  + 
  theme(legend.position = "none") + 
  theme(
    axis.text = element_text(size = 35)  # changes both x and y axis numbers
  )


ggsave(
  filename = "Simpson_vs_deltaT_final.png",  # Specify the output filename and desired format
  plot = p5,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  


##### Simpson Diversity vs Δprecip #####
p6 = ggplot(data = habitat_patches, aes(x = delta.precip, y = Simpson, 
                                   color = delta.precip)) +
  geom_point(alpha = 1, size = 3) + 
  scale_color_gradientn(colors = c("#1d20d0","#1348a3", "#2361a7", "#3279ae", "lightblue")) +
  labs(x = "ΔPrecipitation (experimental - source; cm)", y = "Simpson Index") +
  theme_bw(base_size = 55) + 
  theme(legend.position = "none") + 
  theme(
    axis.text = element_text(size = 35)  # changes both x and y axis numbers
  )

ggsave(
  filename = "Simpson_vs_deltaPrecip_final.png",  # Specify the output filename and desired format
  plot = p6,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  

##### Plant height vs ΔT #####
p7 <- ggplot(data = hp_traits, aes(x = delta.T, y = heightveg_m_CWM, 
                                  color = delta.T)) +
  geom_point(alpha = 1, size = 3) + 
  scale_color_gradientn(colors = c("#ff0000" ,"#ff4d00","#ff7400", "#ff9a00", "#ffc100")) +
  theme_bw(base_size = 55) +
  theme(legend.position = "none") +
  geom_smooth(method = "lm") +
  labs(x = "ΔTemperature (experimental - source; °C)", y = "CWM Plant Height (m)") +
  theme(axis.text = element_text(size = 35)) 

ggsave(
  filename = "Height_vs_deltaT_final.png",  # Specify the output filename and desired format
  plot = p7,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  

##### Plant height vs Δprecip #####
p8 <- ggplot(data = hp_traits, aes(x = delta.precip, y = heightveg_m_CWM, 
                             color = delta.precip)) +
  geom_point(alpha = 1, size = 3) + 
  scale_color_gradientn(colors = c("#1d20d0","#1348a3", "#2361a7", "#3279ae", "lightblue")) +
  theme_bw(base_size = 55) +
  theme(legend.position = "none") +
  labs(x = "ΔPrecipitation (experimental - source; cm)", y = "CWM Plant Height (m)") +
  theme(axis.text = element_text(size = 35))

ggsave(
  filename = "Height_vs_deltaPrecip_final.png",  # Specify the output filename and desired format
  plot = p8,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  

###### Plant SLA vs ΔT #####
p9 <- ggplot(data = hp_traits, aes(x = delta.T, y = SLA_mm2.mg_CWM, 
                             color = delta.T)) +
  geom_point(alpha = 1, size = 3) + 
  scale_color_gradientn(colors = c("#ff0000" ,"#ff4d00","#ff7400", "#ff9a00", "#ffc100")) +
  theme_bw(base_size = 55) +
  theme(legend.position = "none") +
  labs(x = "ΔTemperature (experimental - source; °C)", y = expression(CWM~Plant~SLA~(mm^2*mg^-1))) +
  theme(axis.text = element_text(size = 35))

ggsave(
  filename = "SLA_vs_deltaT_final.png",  # Specify the output filename and desired format
  plot = p9,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  

##### Plant SLA vs Δprecip #####
p10 <- ggplot(data = hp_traits, aes(x = delta.precip, y = SLA_mm2.mg_CWM, 
                             color = delta.precip)) +
  geom_point(alpha = 1, size = 3) + 
  scale_color_gradientn(colors = c("#1d20d0","#1348a3", "#2361a7", "#3279ae", "lightblue")) +
  theme_bw(base_size = 55) +
  theme(legend.position = "none") +
  labs(x = "ΔPrecipitation (experimental - source; cm)", y = expression(CWM~Plant~SLA~(mm^2*mg^-1))) +
  theme(axis.text = element_text(size = 35))

ggsave(
  filename = "SLA_vs_deltaPrecip_final.png",  # Specify the output filename and desired format
  plot = p10,       # Save the last displayed ggplot object (or specify a plot object)
  width = 6,                # Set the desired width
  height = 6,               # Set the desired height
  units = "in",             # Specify the units (e.g., "in", "cm", "mm", "px")
  dpi = 300                 # Set the desired resolution in dots per inch
)  





