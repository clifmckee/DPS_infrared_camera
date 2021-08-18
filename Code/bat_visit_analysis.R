# Code for "Seasonality of bats' date palm sap feeding behavior in an area of Bangladesh associated with Nipah virus spillover"
# Code author: Clif McKee
# Works as of: 18 August 2021
# R version: 4.1.0 "Camp Pontanezen"

#####################
### Configuration ###
#####################

# Load required packages
library(boot)
library(broom)
library(cowplot)
library(dplyr)
library(emmeans)
library(ggeffects)
library(ggplot2)
library(lme4)
library(lubridate)
library(lunar)
library(MASS)
library(multcomp)
library(multcompView)
library(MuMIn)
library(RColorBrewer)
library(reshape2)
library(tidyr)
library(tidyverse)
library(PupillometryR)

## Set appropriate working directory
# Mac
setwd("~/Dropbox/JHSPH postdoc/Projects/Date palm infrared camera analysis")

# Change when scientific notation is used
options(scipen = 12)

# Set dodge position for plotting
set_dodge <- position_dodge(width = 0.5)

#################
### Functions ###
#################

## Function for extracting legend from ggplot object
g_legend <- function(a.gplot) {
  tmp <-
    ggplot_gtable(ggplot_build(a.gplot))
  leg <-
    which(sapply(tmp$grobs,
                 function(x)
                   x$name) == "guide-box")
  legend <-
    tmp$grobs[[leg]]
  
  return(legend)
}

# Create function for bootstrap resampler
days_sampler <- function(df, samples){
  # Make empty file for storing results of loop
  df_day <- NULL
  # Outer loop for iterating through year-month index
  for(i in 1:length(unique(df$index1))){
    # Resampling step: from a particular year-month index,
    # sample with replacement from the available dtobs-treeNo indices
    choose <- base::sample(x = unique(df[df$index1 == unique(df$index1)[i],]$index2),
                           size = samples,
                           replace = TRUE)
    # Inner loop for iterating through sampled dtobs-treeNo indices
    for(j in 1:length(choose)){
      # Summarize the number of bat visits (i.e., those that actually land) per tree per day
      df_filter <- df %>%
        filter(index2 == choose[j]) %>%
        group_by(year = year(as.Date(dtobs, format = "%m/%d/%Y")),
                 month = month(as.Date(dtobs, format = "%m/%d/%Y")),
                 date = as.Date(dtobs, format = "%m/%d/%Y"),
                 treeNo,
                 index1,
                 index2,
                 pnp,
                 Month_Season) %>%
        summarize(visits = sum(nstay_fixed, na.rm = TRUE))
      # Add summarized data to larger file
      df_day <- rbind(df_day, df_filter) %>%
        # Add another index to label dtobs-treeNo indices that are chosen twice in the sampler
        mutate(index3 = row_number())
    }
  }
  # Return the summarized data as output
  return(df_day)
}

## Summarize bootstrap observations by month or season
boots_summary <- function(bootstraps, variables){
  # Make empty file for storing results
  boots_summary <- NULL
  # Loop through the list of bootstrap resamples
  for(i in 1:length(lengths(bootstraps))){
    # Spread the bat species column into separate columns and replace missing values with zeroes
    infra_boot_spread <- spread(data = bootstraps[[i]], key = pnp, value = visits, fill = 0)
    # Rename bat species columns
    colnames(infra_boot_spread)[9:12] <- c("no_bat", "non_Pteropus", "Pteropus", "unidentified")
    # Group bootstrap resampled data by season and month
    infra_boot_summary <- infra_boot_spread %>%
      group_by_(.dots = variables) %>%
      # Calculate mean and sum of bat visits by species
      summarize(mean_non_Pteropus = mean(non_Pteropus, na.rm = TRUE),
                sum_non_Pteropus = sum(non_Pteropus, na.rm = TRUE),
                mean_Pteropus = mean(Pteropus, na.rm = TRUE),
                sum_Pteropus = sum(Pteropus, na.rm = TRUE))
    # Add summarized data to larger file
    boots_summary <- rbind(boots_summary, infra_boot_summary)
  }
  return(boots_summary)
}

# Summarize bootstrap observation summary over months or seasons
boots_grand_summary <- function(boots_summary, variables){
  boots_grand_summary <- boots_summary %>%
    group_by_(.dots = variables) %>%
    # Calculate summary statistics of MEAN bat visits by species (bootstrap estimator of mean):
    # mean, median, minimum, maximum, lower 95% confidence interval, and upper 95% confidence interval
    summarize(Mean_mean_non_Pteropus = mean(mean_non_Pteropus, na.rm = TRUE),
              Med_mean_non_Pteropus = median(mean_non_Pteropus, na.rm = TRUE),
              Min_mean_non_Pteropus = min(mean_non_Pteropus, na.rm = TRUE),
              Max_mean_non_Pteropus = max(mean_non_Pteropus, na.rm = TRUE),
              Lower_mean_non_Pteropus = quantile(x = mean_non_Pteropus, probs = 0.025, na.rm = TRUE),
              Upper_mean_non_Pteropus = quantile(x = mean_non_Pteropus, probs = 0.975, na.rm = TRUE),
              
              Mean_sum_non_Pteropus = mean(sum_non_Pteropus, na.rm = TRUE),
              Med_sum_non_Pteropus = median(sum_non_Pteropus, na.rm = TRUE),
              Min_sum_non_Pteropus = min(sum_non_Pteropus, na.rm = TRUE),
              Max_sum_non_Pteropus = max(sum_non_Pteropus, na.rm = TRUE),
              Lower_sum_non_Pteropus = quantile(x = sum_non_Pteropus, probs = 0.025, na.rm = TRUE),
              Upper_sum_non_Pteropus = quantile(x = sum_non_Pteropus, probs = 0.975, na.rm = TRUE),
              
              Mean_mean_Pteropus = mean(mean_Pteropus, na.rm = TRUE),
              Med_mean_Pteropus = median(mean_Pteropus, na.rm = TRUE),
              Min_mean_Pteropus = min(mean_Pteropus, na.rm = TRUE),
              Max_mean_Pteropus = max(mean_Pteropus, na.rm = TRUE),
              Lower_mean_Pteropus = quantile(x = mean_Pteropus, probs = 0.025, na.rm = TRUE),
              Upper_mean_Pteropus = quantile(x = mean_Pteropus, probs = 0.975, na.rm = TRUE),
              
              Mean_sum_Pteropus = mean(sum_Pteropus, na.rm = TRUE),
              Med_sum_Pteropus = median(sum_Pteropus, na.rm = TRUE),
              Min_sum_Pteropus = min(sum_Pteropus, na.rm = TRUE),
              Max_sum_Pteropus = max(sum_Pteropus, na.rm = TRUE),
              Lower_sum_Pteropus = quantile(x = sum_Pteropus, probs = 0.025, na.rm = TRUE),
              Upper_sum_Pteropus = quantile(x = sum_Pteropus, probs = 0.975, na.rm = TRUE))
  return(boots_grand_summary)
}

###################################
### Loading and formatting data ###
###################################

# Read in bat visit data from file
infra <- read.csv("./Data/13. Bat_image_and_temperature_data_with_treeinfo_27_Feb_2020.csv",
                  header = TRUE)

# Define approach methods
infra <- mutate(infra, way = case_when(way == 0 ~ "no contamination",
                                way == 1 ~ "left/right",
                                way == 2 ~ "direct landing",
                                way == 3 ~ "landing on branches and climbing down"))
# Define contamination methods
infra <- mutate(infra, conPlace = case_when(conPlace == 0 ~ "no contamination",
                                            conPlace == 1 ~ "shaved part",
                                            conPlace == 2 ~ "sap stream",
                                            conPlace == 3 ~ "tap",
                                            conPlace == 4 ~ "collection pot"))
# Define species
infra <- mutate(infra, pnp = case_when(pnp == 0 ~ "no bat",
                                       pnp == 1 ~ "Pteropus",
                                       pnp == 2 ~ "Non-pteropus",
                                       pnp == 3 ~ "Unidentified"))
# Define months
infra <- mutate(infra, mon = case_when(mon == 1 ~ "January",
                                       mon == 2 ~ "February",
                                       mon == 3 ~ "March",
                                       mon == 4 ~ "April",
                                       mon == 5 ~ "May",
                                       mon == 6 ~ "June",
                                       mon == 7 ~ "July",
                                       mon == 8 ~ "August",
                                       mon == 9 ~ "September",
                                       mon == 10 ~ "October",
                                       mon == 11 ~ "November",
                                       mon == 12 ~ "December"))
# Define seasons
infra <- mutate(infra, Month_Season = case_when(mon %in% c("December", "January", "February") ~ "Winter",
                                       mon %in% c("March", "April", "May") ~ "Spring",
                                       mon %in% c("June", "July", "August", "September") ~ "Monsoon",
                                       mon %in% c("October", "November") ~ "Postmonsoon"))
# Calculate shaved area of tree
infra <- mutate(infra, shavearea = treeLength * treeWidth)

## Create formatted "intime" and "outtime" columns
# Unite columns for hours, minutes, and seconds into new "intime" column
infra <- unite(data = infra, col = "intime", c("inH", "inM", "inS"), sep = ":", remove = FALSE)
# Add AM/PM to "intime"
infra <- unite(data = infra, col = "intime_ampm", c("intime", "ampm"), sep = " ", remove = FALSE)
# Change "intime" column to 24-hour POSIXct format
infra$intime_POSIXct <- as.POSIXct(infra$intime_ampm, format = "%I:%M:%S %p")
# Unite columns for hours, minutes, and seconds into new "outtime" column
infra <- unite(data = infra, col = "outtime", c("outH", "outM", "outS"), sep = ":", remove = FALSE)
# Add AM/PM to "outtime"
infra <- unite(data = infra, col = "outtime_ampm", c("outtime", "ampm"), sep = " ", remove = FALSE)
# Change "outtime" column to 24-hour POSIXct format
infra$outtime_POSIXct <- as.POSIXct(infra$outtime_ampm, format = "%I:%M:%S %p")

# Make new column for hour of observation based on formatted "intime"
infra$hour <- hour(infra$intime_POSIXct)
# Rename durCount column as durCont
infra <- rename(infra, durCont = durCount)

# Fix some erroneously large values in the stay and contact columns
infra$nstay_fixed <- as.numeric(infra$nStay>0)
infra$nconta_fixed <- as.numeric(infra$nConta>0)

# Replace bad values in shaving date
infra <- mutate(infra, dtshaving = replace(dtshaving, dtobs == "5/5/2013" & dtshaving == "5/5/2014", "5/5/2013"))

##################################
### Observation summary output ###
##################################

# Quick summary for abstract
infra_summary <- data.frame(stat = c(length(unique(infra$dtobs)),
                                     sum(infra$nVisit, na.rm = TRUE),
                                     sum(infra[infra$pnp == "Non-pteropus",]$nVisit, na.rm = TRUE),
                                     sum(infra[infra$pnp == "Non-pteropus",]$nVisit, na.rm = TRUE)/sum(infra$nVisit, na.rm = TRUE)*100,
                                     sum(infra[infra$pnp == "Pteropus",]$nVisit, na.rm = TRUE),
                                     sum(infra[infra$pnp == "Pteropus",]$nVisit, na.rm = TRUE)/sum(infra$nVisit, na.rm = TRUE)*100,
                                     sum(infra[infra$pnp == "Unidentified",]$nVisit, na.rm = TRUE),
                                     sum(infra[infra$pnp == "Unidentified",]$nVisit, na.rm = TRUE)/sum(infra$nVisit, na.rm = TRUE)*100,
                                     sum(infra$nstay_fixed, na.rm = TRUE),
                                     sum(infra[infra$pnp == "Non-pteropus",]$nstay_fixed, na.rm = TRUE),
                                     sum(infra[infra$pnp == "Non-pteropus",]$nstay_fixed, na.rm = TRUE)/sum(infra$nstay_fixed, na.rm = TRUE)*100,
                                     sum(infra[infra$pnp == "Pteropus",]$nstay_fixed, na.rm = TRUE),
                                     sum(infra[infra$pnp == "Pteropus",]$nstay_fixed, na.rm = TRUE)/sum(infra$nstay_fixed, na.rm = TRUE)*100,
                                     sum(infra[infra$pnp == "Unidentified",]$nstay_fixed, na.rm = TRUE),
                                     sum(infra[infra$pnp == "Unidentified",]$nstay_fixed, na.rm = TRUE)/sum(infra$nstay_fixed, na.rm = TRUE)*100
                                     ))
rownames(infra_summary) <- c("Total observation nights",
                             "Total bat visits",
                             "Total non-Pteropus visits",
                             "Percentage of total visits, non-Pteropus",
                             "Total Pteropus visits",
                             "Percentage of total visits, Pteropus",
                             "Total unidentified visits",
                             "Percentage of total visits, unidentified",
                             "Total bat stays",
                             "Total non-Pteropus stays",
                             "Percentage of total stays, non-Pteropus",
                             "Total Pteropus stays",
                             "Percentage of total stays, Pteropus",
                             "Total unidentified stays",
                             "Percentage of total stays, unidentified"
                             )
write.csv(infra_summary, "./Results/abstract_summary_statistics.csv")

# Summarize the number of unique trees and observation days per year and month
infra_days_trees <- infra %>%
  group_by(yr, mon) %>%
  summarize(unique_trees = n_distinct(treeNo, na.rm = TRUE),
            unique_days = n_distinct(dtobs, na.rm = TRUE))

# Summarize the number of unique observation days per tree, year, and month
infra_count <- infra %>%
  group_by(yr, mon, treeNo) %>%
  summarize(count = n_distinct(dtobs, na.rm = TRUE))

# Summarize the number of observation days per year and month
infra_totals <- infra_count %>%
  group_by(mon) %>%
  summarize(total_count = sum(count, na.rm = TRUE))

#######################################
### Raw data analysis of bat visits ###
#######################################

# Summarize data about observation nights and bat contacts
raw_summary <- data.frame(count = c(n_distinct(infra$dtobs, na.rm = TRUE),
                                   n_distinct(infra %>% filter(nVisit>0) %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Pteropus", nVisit>0) %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Non-pteropus", nVisit>0) %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Unidentified", nVisit>0) %>% dplyr::select(dtobs), na.rm = TRUE),
                                   
                                   n_distinct(infra %>% filter(nstay_fixed>0) %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Pteropus", nstay_fixed>0) %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Non-pteropus", nstay_fixed>0) %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Unidentified", nstay_fixed>0) %>% dplyr::select(dtobs), na.rm = TRUE),
                                   
                                   n_distinct(infra %>% filter(nconta_fixed>0, conPlace!="no contamination") %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Pteropus", nconta_fixed>0, conPlace!="no contamination") %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Non-pteropus", nconta_fixed>0, conPlace!="no contamination") %>% dplyr::select(dtobs), na.rm = TRUE),
                                   n_distinct(infra %>% filter(pnp == "Unidentified", nconta_fixed>0, conPlace!="no contamination") %>% dplyr::select(dtobs), na.rm = TRUE),
                                   
                                   sum(infra$nVisit, na.rm = TRUE),
                                   sum(infra$nstay_fixed, na.rm = TRUE),
                                   sum(infra %>% filter(conPlace!="no contamination") %>% dplyr::select(nconta_fixed), na.rm = TRUE),
                                   
                                   sum(infra %>% filter(pnp == "Pteropus") %>% dplyr::select(nVisit), na.rm = TRUE),
                                   sum(infra %>% filter(pnp == "Pteropus") %>% dplyr::select(nstay_fixed), na.rm = TRUE),
                                   sum(infra %>% filter(pnp == "Pteropus", conPlace!="no contamination") %>% dplyr::select(nconta_fixed), na.rm = TRUE),
                                   
                                   sum(infra %>% filter(pnp == "Non-pteropus") %>% dplyr::select(nVisit), na.rm = TRUE),
                                   sum(infra %>% filter(pnp == "Non-pteropus") %>% dplyr::select(nstay_fixed), na.rm = TRUE),
                                   sum(infra %>% filter(pnp == "Non-pteropus", conPlace!="no contamination") %>% dplyr::select(nconta_fixed), na.rm = TRUE),
                                   
                                   sum(infra %>% filter(pnp == "Unidentified") %>% dplyr::select(nVisit), na.rm = TRUE),
                                   sum(infra %>% filter(pnp == "Unidentified") %>% dplyr::select(nstay_fixed), na.rm = TRUE),
                                   sum(infra %>% filter(pnp == "Unidentified", conPlace!="no contamination") %>% dplyr::select(nconta_fixed), na.rm = TRUE)))
# Rename rows
row.names(raw_summary) <- c("Total nights of observation",
                            "Nights with 1+ bats present at night",
                            "Nights with 1+ Pteropus bats present at night",
                            "Nights with 1+ Non-Pteropus bats present at night",
                            "Nights with 1+ unidentified bats present at night",
                            
                            "Nights with 1+ bats staying at night",
                            "Nights with 1+ Pteropus bats staying at night",
                            "Nights with 1+ Non-Pteropus bats staying at night",
                            "Nights with 1+ Unidentified bats staying at night",
                            
                            "Nights with 1+ bats contacting sap at night",
                            "Nights with 1+ Pteropus bats contacting sap at night",
                            "Nights with 1+ Non-Pteropus bats contacting sap at night",
                            "Nights with 1+ Unidentified bats contacting sap at night",
                            
                            "Total bats present at night",
                            "Total bats staying at night",
                            "Total bats contacting sap",
                            
                            "Pteropus bats present at night",
                            "Pteropus bats staying at night",
                            "Pteropus bats contacting sap",
                            
                            "Non-Pteropus bats present at night",
                            "Non-Pteropus bats staying at night",
                            "Non-Pteropus bats contacting sap",
                            
                            "Unidentified bats present at night",
                            "Unidentified bats staying at night",
                            "Unidentified bats contacting sap")
# Calculate percentages
raw_summary$percent = c(raw_summary$count[1]/raw_summary$count[1]*100,
         raw_summary$count[2]/raw_summary$count[1]*100,
         raw_summary$count[3]/raw_summary$count[1]*100,
         raw_summary$count[4]/raw_summary$count[1]*100,
         raw_summary$count[5]/raw_summary$count[1]*100,
         
         raw_summary$count[6]/raw_summary$count[1]*100,
         raw_summary$count[7]/raw_summary$count[1]*100,
         raw_summary$count[8]/raw_summary$count[1]*100,
         raw_summary$count[9]/raw_summary$count[1]*100,
         
         raw_summary$count[10]/raw_summary$count[1]*100,
         raw_summary$count[11]/raw_summary$count[1]*100,
         raw_summary$count[12]/raw_summary$count[1]*100,
         raw_summary$count[13]/raw_summary$count[1]*100,
         
         raw_summary$count[14]/raw_summary$count[14]*100,
         raw_summary$count[15]/raw_summary$count[14]*100,
         raw_summary$count[16]/raw_summary$count[14]*100,
         
         raw_summary$count[17]/raw_summary$count[17]*100,
         raw_summary$count[18]/raw_summary$count[17]*100,
         raw_summary$count[19]/raw_summary$count[17]*100,
         
         raw_summary$count[20]/raw_summary$count[20]*100,
         raw_summary$count[21]/raw_summary$count[20]*100,
         raw_summary$count[22]/raw_summary$count[20]*100,
         
         raw_summary$count[23]/raw_summary$count[23]*100,
         raw_summary$count[24]/raw_summary$count[23]*100,
         raw_summary$count[25]/raw_summary$count[23]*100)
# Output observation night summary to file
write.csv(raw_summary, "./Results/observation_night_summary.csv")

# Make pivot table for methods bats use to approach date palm sap before contamination
approach_summary <- infra %>% filter(nconta_fixed>0) %>% group_by(pnp, way) %>% summarize(count = n()) %>%
  spread(data = ., key = pnp, value = count, fill = 0)
# Rename columns
colnames(approach_summary) <- c("way", "non_Pteropus", "Pteropus", "unidentified")
# Calculate percentages by species
approach_summary$non_Pteropus_percent <- approach_summary$non_Pteropus/sum(approach_summary$non_Pteropus, na.rm = TRUE)*100
approach_summary$Pteropus_percent <- approach_summary$Pteropus/sum(approach_summary$Pteropus, na.rm = TRUE)*100
approach_summary$unidentified_percent <- approach_summary$unidentified/sum(approach_summary$unidentified, na.rm = TRUE)*100
# Reorder columns
approach_summary <- approach_summary[, c(1, 2, 5, 3, 6, 4, 7)]
# Output approach method summary to file
write.csv(approach_summary, "./Results/approach_method_summary.csv")

# Make pivot table for methods bats contaminate date palm sap
contam_summary <- infra %>% filter(nconta_fixed>0) %>% group_by(pnp, conPlace) %>% summarize(count = n()) %>%
  spread(data = ., key = pnp, value = count, fill = 0)
# Rename columns
colnames(contam_summary) <- c("conplace", "non_Pteropus", "Pteropus", "unidentified")
# Calculate percentages by species
contam_summary$non_Pteropus_percent <- contam_summary$non_Pteropus/sum(contam_summary$non_Pteropus, na.rm = TRUE)*100
contam_summary$Pteropus_percent <- contam_summary$Pteropus/sum(contam_summary$Pteropus, na.rm = TRUE)*100
contam_summary$unidentified_percent <- contam_summary$unidentified/sum(contam_summary$unidentified, na.rm = TRUE)*100
# Reorder columns
contam_summary <- contam_summary[, c(1, 2, 5, 3, 6, 4, 7)]
# Output approach method summary to file
write.csv(contam_summary, "./Results/contamination_method_summary.csv")

# Count the number of visits per tree per day
infra_day <- infra %>%
  mutate(year = year(as.Date(dtobs, format = "%m/%d/%Y")),
         month = month(as.Date(dtobs, format = "%m/%d/%Y")),
         date = as.Date(dtobs, format = "%m/%d/%Y"),
         date_shaved = as.Date(dtshaving, format = "%m/%d/%Y"),
         days_since_shaving = date - date_shaved,
         moon_phase = lunar.phase(date, shift = 12),
         moon_illum = lunar.illumination(date, shift = 12)) %>%
  # Retain information about year, month, tree number, season, temperature, and tree measurements
  group_by(year,
           month,
           date,
           date_shaved,
           days_since_shaving,
           moon_phase,
           moon_illum,
           treeNo,
           pnp,
           Month_Season,
           min_temp,
           mean_temp,
           median_temp,
           sd_temp,
           max_temp,
           treeHeight,
           treeAge,
           shavearea) %>%
  # Count only those bats that stayed on the tree
  summarize(visits = sum(nstay_fixed, na.rm = TRUE))

# Spread the bat species column into separate columns and replace missing values with zeroes
infra_day_spread <- spread(data = infra_day, key = pnp, value = visits, fill = 0) %>%
  # Calculate total bat visits across all species
  rowwise(.) %>%
  mutate(total = sum(`Non-pteropus` + Pteropus + Unidentified, na.rm = TRUE))
# Rename columns
colnames(infra_day_spread)[18:22] <- c("no_bat", "non_Pteropus", "Pteropus", "unidentified", "total")

# Calculate number of tree nights with observations of Pteropus and non-Pteropus bats
tree_nights <- data.frame(tree_nights_w_both = nrow(infra_day_spread %>%
                                                      filter(non_Pteropus > 0, Pteropus > 0)),
                          tree_nights_w_none = nrow(infra_day_spread %>%
                                                      filter(non_Pteropus == 0, Pteropus == 0)),
                          tree_nights_w_onlyP = nrow(infra_day_spread %>%
                                                       filter(non_Pteropus == 0, Pteropus > 0)),
                          tree_nights_w_onlyNP = nrow(infra_day_spread %>%
                                                        filter(non_Pteropus > 0, Pteropus == 0)),
                          tree_nights_w_multiP = nrow(infra_day_spread %>%
                                                        filter(Pteropus > 1)),
                          tree_nights_w_multiNP = nrow(infra_day_spread %>%
                                                         filter(non_Pteropus > 1)))
round(tree_nights/nrow(infra_day_spread)*100, 1)

# Calculate mean number and total number of visits for each month and species
infra_month_obs_summary <- infra_day_spread %>%
  group_by(month) %>%
  summarize(mean_non_Pteropus = mean(non_Pteropus, na.rm = TRUE),
            sum_non_Pteropus = sum(non_Pteropus, na.rm = TRUE),
            mean_Pteropus = mean(Pteropus, na.rm = TRUE),
            sum_Pteropus = sum(Pteropus, na.rm = TRUE))
# Save output to file
write.csv(infra_month_obs_summary, "./Results/visits_summary_by_month.csv")

# Calculate mean number and total number of visits for each season and species
infra_season_obs_summary <- infra_day_spread %>%
  group_by(Month_Season) %>%
  summarize(mean_no_bat = mean(no_bat, na.rm = TRUE),
            sum_no_bat = sum(no_bat, na.rm = TRUE),
            mean_non_Pteropus = mean(non_Pteropus, na.rm = TRUE),
            sum_non_Pteropus = sum(non_Pteropus, na.rm = TRUE),
            mean_Pteropus = mean(Pteropus, na.rm = TRUE),
            sum_Pteropus = sum(Pteropus, na.rm = TRUE),
            mean_unidentified = mean(unidentified, na.rm = TRUE),
            sum_unidentified = sum(unidentified, na.rm = TRUE)) %>%
  rowwise() %>%
  mutate(season_total = sum(sum_non_Pteropus, sum_Pteropus, sum_unidentified, na.rm = TRUE))
infra_season_obs_summary$non_Pteropus_percent = infra_season_obs_summary$sum_non_Pteropus / sum(infra_season_obs_summary$sum_non_Pteropus, na.rm = TRUE) * 100
infra_season_obs_summary$Pteropus_percent = infra_season_obs_summary$sum_Pteropus / sum(infra_season_obs_summary$sum_Pteropus, na.rm = TRUE) * 100
infra_season_obs_summary$unidentified_percent = infra_season_obs_summary$sum_unidentified / sum(infra_season_obs_summary$sum_unidentified, na.rm = TRUE) * 100
infra_season_obs_summary$total_percent = infra_season_obs_summary$season_total / sum(infra_season_obs_summary$season_total, na.rm = TRUE) * 100
# Save output to file
write.csv(infra_season_obs_summary, "./Results/visits_summary_by_season.csv")

## Plots for raw data analysis: mean of bat visits by month

# Subplot for Pteropus mean visits
obs_month_mean_visitsA <- ggplot(infra_day_spread, aes(x = as.factor(month), y = Pteropus,
                                                       color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Number of bat visits", 
                     limits = c(0, 40),
                     breaks = seq(0, 40, 10)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Subplot for non-Pteropus mean visits
obs_month_mean_visitsB <- ggplot(infra_day_spread, aes(x = as.factor(month), y = non_Pteropus,
                                                       color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Number of bat visits", 
                     limits = c(0, 1500),
                     breaks = seq(0, 1500, 250)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Combine mean visits subplots vertically (removing legends from both plots)
obs_month_mean_visits <- plot_grid(obs_month_mean_visitsA + theme(legend.position = "none"),
                                   obs_month_mean_visitsB + theme(legend.position = "none"),
                                    nrow = 2, align = "v")
# Extract legend from first plot
obs_month_mean_visits_legend <- g_legend(obs_month_mean_visitsA)
# Arrange combined subplots and legend
plot_grid(obs_month_mean_visits, obs_month_mean_visits_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/observed_month_mean_bat_visits.png",
       height = 6, width = 8, dpi = 300, bg = "white")
ggsave("./Results/observed_month_mean_bat_visits.tiff",
       height = 6, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: mean of bat visits by month, but separated by years

# Plot for Pteropus mean visits over all years
ggplot(infra_day_spread, aes(x = as.factor(month), y = Pteropus,
                             color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Number of bat visits", 
                     limits = c(0, 40),
                     breaks = seq(0, 40, 10)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        strip.background = element_blank()) +
  # Separate years across panels in figure
  facet_wrap(~year)
# Save plot directly to file
ggsave("./Results/observed_month_mean_bat_visits_by_year_Pteropus.png",
       height = 6, width = 10, dpi = 300, bg = "white")
ggsave("./Results/observed_month_mean_bat_visits_by_year_Pteropus.tiff",
       height = 6, width = 10, dpi = 300, bg = "white", compression = "lzw")

# Plot for non-Pteropus mean visits over all years
ggplot(infra_day_spread, aes(x = as.factor(month), y = non_Pteropus,
                             color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  ## Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Number of bat visits", 
                     limits = c(0, 1500),
                     breaks = seq(0, 1500, 250)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        strip.background = element_blank()) +
  # Separate years across panels in figure
  facet_wrap(~year)
# Save plot directly to file
ggsave("./Results/observed_month_mean_bat_visits_by_year_nonPteropus.png",
       height = 6, width = 10, dpi = 300, bg = "white")
ggsave("./Results/observed_month_mean_bat_visits_by_year_nonPteropus.tiff",
       height = 6, width = 10, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: mean of bat visits by season

# Subplot for Pteropus mean visits
obs_season_mean_visitsA <- ggplot(infra_day_spread, aes(x = as.factor(Month_Season), y = Pteropus,
                                                        color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Number of bat visits",
                     limits = c(0, 40),
                     breaks = seq(0, 40, 10)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus mean visits
obs_season_mean_visitsB <- ggplot(infra_day_spread, aes(x = as.factor(Month_Season), y = non_Pteropus,
                                                        color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Number of bat visits",
                     limits = c(0, 1500),
                     breaks = seq(0, 1500, 250)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine mean visits subplots vertically (removing legends from both plots)
plot_grid(obs_season_mean_visitsA, obs_season_mean_visitsB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/observed_season_mean_bat_visits.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/observed_season_mean_bat_visits.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: sum of bat visits by month

# Subplot for Pteropus sum of visits
obs_month_sum_visitsA <- infra_day_spread %>%
  group_by(month, Month_Season) %>%
  summarize(sum_visits = sum(Pteropus)) %>%
  # Bar plot for sum of visits by month, with labels
  ggplot() +
  geom_col(aes(x = as.factor(month), y = sum_visits, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = as.factor(month), y = sum_visits + (300/25), label = sum_visits), size = 3) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month", labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Total number of bat visits",
                     limits = c(0, 300),
                     breaks = seq(0, 300, 50)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Subplot for non-Pteropus sum of visits
obs_month_sum_visitsB <- infra_day_spread %>%
  group_by(month, Month_Season) %>%
  summarize(sum_visits = sum(non_Pteropus)) %>%
  # Bar plot for sum of visits by month, with labels
  ggplot() +
  geom_col(aes(x = as.factor(month), y = sum_visits, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = as.factor(month), y = sum_visits + (15000/25), label = sum_visits), size = 3) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month", labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Total number of bat visits",
                     limits = c(0, 15000),
                     breaks = seq(0, 15000, 2500)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Combine sum of visits subplots vertically (removing legends from both plots)
obs_month_sum_visits <- plot_grid(obs_month_sum_visitsA + theme(legend.position = "none"),
                                   obs_month_sum_visitsB + theme(legend.position = "none"),
                                   nrow = 2, align = "v")
# Extract legend from first plot
obs_month_sum_visits_legend <- g_legend(obs_month_sum_visitsA)
# Arrange combined subplots and legend
plot_grid(obs_month_sum_visits, obs_month_sum_visits_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/observed_month_sum_bat_visits.png",
       height = 6, width = 8, dpi = 300, bg = "white")
ggsave("./Results/observed_month_sum_bat_visits.tiff",
       height = 6, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: sum of bat visits by season

# Subplot for Pteropus sum of visits
obs_season_sum_visitsA <- infra_day_spread %>%
  group_by(Month_Season) %>%
  summarize(sum_visits = sum(Pteropus)) %>%
  # Bar plot for sum of visits by season, with labels
  ggplot() +
  geom_col(aes(x = as.factor(Month_Season), y = sum_visits, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = as.factor(Month_Season), y = sum_visits + (500/25), label = sum_visits), size = 3) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Total number of bat visits",
                     limits = c(0, 500),
                     breaks = seq(0, 500, 100)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus sum of visits
obs_season_sum_visitsB <- infra_day_spread %>%
  group_by(Month_Season) %>%
  summarize(sum_visits = sum(non_Pteropus)) %>%
  # Bar plot for sum of visits by season, with labels
  ggplot() +
  geom_col(aes(x = as.factor(Month_Season), y = sum_visits, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = as.factor(Month_Season), y = sum_visits + (25000/25), label = sum_visits), size = 3) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Total number of bat visits",
                     limits = c(0, 25000),
                     breaks = seq(0, 25000, 5000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine sum of visits subplots vertically (removing legends from both plots)
plot_grid(obs_season_sum_visitsA, obs_season_sum_visitsB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/observed_season_sum_bat_visits.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/observed_season_sum_bat_visits.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

########################################################
### Raw data analysis of bat visits by time of night ###
########################################################

# Summarize Pteropus bat visits by hour of night
infra_prop_Pteropus <- infra %>%
  # Filter for nights and trees with visits by Pteropus bats
  filter(nVisit == 1, pnp == "Pteropus") %>%
  # Group data by hour, season, and species
  group_by(hour = hour(intime_POSIXct),
           season = Month_Season,
           species = pnp) %>%
  # Count the total number of visits in a given hour
  summarize(visits = n()) %>%
  # Create new variable to order the hours for plotting
  mutate(hour_trans = case_when(hour == 18 ~ 1,
                                hour == 19 ~ 2,
                                hour == 20 ~ 3,
                                hour == 21 ~ 4,
                                hour == 22 ~ 5,
                                hour == 23 ~ 6,
                                hour == 0 ~ 7,
                                hour == 1 ~ 8,
                                hour == 2 ~ 9,
                                hour == 3 ~ 10,
                                hour == 4 ~ 11,
                                hour == 5 ~ 12))
# Add in zeroes for missing rows for particular hours and seasons
infra_prop_Pteropus[47,] <- list(5, "Postmonsoon", "Pteropus", 0, 12)
infra_prop_Pteropus[48,] <- list(6, "Winter", "Pteropus", 0, 13)
infra_prop_Pteropus[49,] <- list(6, "Spring", "Pteropus", 0, 13)
infra_prop_Pteropus[50,] <- list(6, "Monsoon", "Pteropus", 0, 13)
infra_prop_Pteropus[51,] <- list(6, "Postmonsoon", "Pteropus", 0, 13)
infra_prop_Pteropus[52,] <- list(17, "Winter", "Pteropus", 0, 0)
infra_prop_Pteropus[53,] <- list(17, "Spring", "Pteropus", 0, 0)
infra_prop_Pteropus[54,] <- list(17, "Monsoon", "Pteropus", 0, 0)
infra_prop_Pteropus[55,] <- list(17, "Postmonsoon", "Pteropus", 0, 0)
infra_prop_Pteropus[56,] <- list(18, "Spring", "Pteropus", 0, 1)

# Calculate total and proportion of total visits by hour for Pteropus bats
infra_prop_Pteropus <- merge(infra_prop_Pteropus, infra_prop_Pteropus %>%
                               group_by(season, species) %>%
                               summarize(total_visits = sum(visits, na.rm = TRUE))) %>%
  mutate(prop_visits = visits / total_visits)
# Check the data so that each hour has 4 points (for each season)
infra_prop_Pteropus %>% group_by(hour) %>% summarize(n_distinct(season, na.rm = TRUE))

# Summarize non-Pteropus bat visits by hour of night
infra_prop_nonPteropus <- infra %>%
  # Filter for nights and trees with visits by non-Pteropus bats
  filter(nVisit == 1, pnp == "Non-pteropus") %>%
  # Group data by hour, season, and species
  group_by(hour = hour(intime_POSIXct),
           season = Month_Season,
           species = pnp) %>%
  # Count the total number of visits in a given hour
  summarize(visits = n()) %>%
  # Create new variable to order the hours for plotting
  mutate(hour_trans = case_when(hour == 17 ~ 0,
                                hour == 18 ~ 1,
                                hour == 19 ~ 2,
                                hour == 20 ~ 3,
                                hour == 21 ~ 4,
                                hour == 22 ~ 5,
                                hour == 23 ~ 6,
                                hour == 0 ~ 7,
                                hour == 1 ~ 8,
                                hour == 2 ~ 9,
                                hour == 3 ~ 10,
                                hour == 4 ~ 11,
                                hour == 5 ~ 12,
                                hour == 6 ~ 13))
# Add in zeroes for missing rows for particular hours and seasons
infra_prop_nonPteropus[50,] <- list(5, "Monsoon", "Non-pteropus", 0, 12)
infra_prop_nonPteropus[51,] <- list(6, "Spring", "Non-pteropus", 0, 13)
infra_prop_nonPteropus[52,] <- list(6, "Monsoon", "Non-pteropus", 0, 13)
infra_prop_nonPteropus[53,] <- list(6, "Postmonsoon", "Non-pteropus", 0, 13)
infra_prop_nonPteropus[54,] <- list(17, "Spring", "Non-pteropus", 0, 0)
infra_prop_nonPteropus[55,] <- list(17, "Monsoon", "Non-pteropus", 0, 0)
infra_prop_nonPteropus[56,] <- list(17, "Postmonsoon", "Non-pteropus", 0, 0)

# Calculate total and proportion of total visits by hour for Pteropus bats
infra_prop_nonPteropus <- merge(infra_prop_nonPteropus, infra_prop_nonPteropus %>%
                               group_by(season, species) %>%
                               summarize(total_visits = sum(visits, na.rm = TRUE))) %>%
  mutate(prop_visits = visits / total_visits)
# Check the data so that each hour has 4 points (for each season)
infra_prop_nonPteropus %>% group_by(hour) %>% summarize(n_distinct(season, na.rm = TRUE))

## Plots for raw data analysis: bat visits by time of night

# Subplot for Pteropus visits by time of night
obs_hour_visitsA <- ggplot(infra_prop_Pteropus, aes(x = as.factor(hour_trans), y = prop_visits,
                                                    color = as.factor(season),
                                                    linetype = as.factor(season),
                                                    group = as.factor(season))) +
  # Plot separate lines for each season
  geom_line() +
  # Set x axis (translating back to original hours)
  scale_x_discrete(name = "Hour", labels = c(17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6)) +
  # Plot the mean proportion of visits across all seasons
  stat_summary_bin(aes(x = as.factor(hour_trans), y = prop_visits, group = 1),
                   fun = "mean", na.rm = TRUE,
                   geom = "point", size = 1) +
  # Set y axis
  scale_y_continuous(name = "Proportion of bat visits",
                     limits = c(0, 0.25),
                     breaks = seq(0, 0.25, 0.05)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Set the line types and arrangement of seasons
  scale_linetype_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                        labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                 values = c("solid", "dashed", "dotdash", "dotted")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.key.width = unit(3, "line"))

# Subplot for non-Pteropus visits by time of night
obs_hour_visitsB <- ggplot(infra_prop_nonPteropus, aes(x = as.factor(hour_trans), y = prop_visits,
                                                       color = as.factor(season),
                                                       linetype = as.factor(season),
                                                       group = as.factor(season))) +
  # Plot separate lines for each season
  geom_line() +
  # Set x axis (translating back to original hours)
  scale_x_discrete(name = "Hour", labels = c(17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6)) +
  # Plot the mean proportion of visits across all seasons
  stat_summary_bin(aes(x = as.factor(hour_trans), y = prop_visits, group = 1), fun = "mean", na.rm = TRUE, geom = "point", size = 1) +
  # Set y axis
  scale_y_continuous(name = "Proportion of bat visits",
                     limits = c(0, 0.25),
                     breaks = seq(0, 0.25, 0.05)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the line types and arrangement of seasons
  scale_linetype_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                        labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                        values = c("solid", "dashed", "dotdash", "dotted")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.key.width = unit(3, "line"))

# Combine visit time of night subplots vertically (removing legends from both plots)
obs_hour_visits <- plot_grid(obs_hour_visitsA + theme(legend.position = "none"),
                             obs_hour_visitsB + theme(legend.position = "none"),
                                   nrow = 2, align = "v")
# Extract legend from first plot
obs_hour_visits_legend <- g_legend(obs_hour_visitsA)
# Arrange combined subplots and legend
plot_grid(obs_hour_visits, obs_hour_visits_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/observed_bat_visits_time_of_night.png",
       height = 6, width = 9, dpi = 300, bg = "white")
ggsave("./Results/observed_bat_visits_time_of_night.tiff",
       height = 6, width = 9, dpi = 300, bg = "white", compression = "lzw")

###############################################################
### Raw data analysis of bat visit duration and sap contact ###
###############################################################

# Filter visit data for visits where bats stayed
infra_stay <- infra %>%
  filter(nstay_fixed == 1, durStay > 0)

# Filter visit data for visits where bats contacted sap
infra_cont <- infra %>%
  filter(nconta_fixed == 1, durCont > 0, conPlace!="no contamination")

## Plots for raw data analysis: visit duration by month

# Subplot for Pteropus visit duration
obs_month_durstayA <- ggplot(infra_stay %>% filter(pnp == "Pteropus"),
                             aes(x = month(as.Date(dtobs, format = "%m/%d/%Y")), y = durStay,
                                 color = Month_Season, fill = Month_Season, 
                                 group = month(as.Date(dtobs, format = "%m/%d/%Y")))) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename numbered months as abbreviated names
  scale_x_continuous(name = "Month", 
                     breaks = seq(1, 12, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Duration of visit (sec)", 
                     limits = c(0, 15000),
                     breaks = seq(0, 15000, 2500)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Subplot for non-Pteropus visit duration
obs_month_durstayB <- ggplot(infra_stay %>% filter(pnp == "Non-pteropus"),
                             aes(x = month(as.Date(dtobs, format = "%m/%d/%Y")), y = durStay,
                                 color = Month_Season, fill = Month_Season,
                                 group = month(as.Date(dtobs, format = "%m/%d/%Y")))) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename numbered months as abbreviated names
  scale_x_continuous(name = "Month", 
                   breaks = seq(1, 12, 1),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Duration of visit (sec)", 
                     limits = c(0, 800),
                     breaks = seq(0, 800, 200)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Combine duration of bat visit subplots vertically (removing legends from both plots)
obs_month_durstay <- plot_grid(obs_month_durstayA + theme(legend.position = "none"),
                               obs_month_durstayB + theme(legend.position = "none"),
                               nrow = 2, align = "v")
# Extract legend from first plot
obs_month_durstay_legend <- g_legend(obs_month_durstayA)
# Arrange combined subplots and legend
plot_grid(obs_month_durstay, obs_month_durstay_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/observed_month_bat_visit_duration.png",
       height = 6, width = 8, dpi = 300, bg = "white")
ggsave("./Results/observed_month_bat_visit_duration.tiff",
       height = 6, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: visit duration by season

# Subplot for Pteropus visit duration
obs_season_durstayA <- ggplot(infra_stay %>% filter(pnp == "Pteropus"),
                              aes(x = as.factor(Month_Season), y = durStay,
                                  color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Duration of visit (sec)",
                     limits = c(0, 15000),
                     breaks = seq(0, 15000, 2500)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus visit duration
obs_season_durstayB <- ggplot(infra_stay %>% filter(pnp == "Non-pteropus"),
                              aes(x = as.factor(Month_Season), y = durStay,
                                  color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Duration of visit (sec)",
                     limits = c(0, 800),
                     breaks = seq(0, 800, 100)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine duration of bat visit subplots vertically (removing legends from both plots)
plot_grid(obs_season_durstayA, obs_season_durstayB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/observed_season_bat_visit_duration.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/observed_season_bat_visit_duration.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: cumulative visit durations by month

# Subplot for Pteropus cumulative duration
obs_month_sum_durstayA <- infra_stay %>%
  filter(pnp == "Pteropus") %>%
  group_by(month = month(as.Date(dtobs, format = "%m/%d/%Y")), Month_Season) %>%
  summarize(cum_duration = sum(durStay)) %>%
  # Bar plot for cumulative duration of visits by month, with labels
  ggplot() +
  geom_col(aes(x = month, y = cum_duration, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = month, y = cum_duration + (300000/25), label = cum_duration), size = 3) +
  # Rename numbered months as abbreviated names
  scale_x_continuous(name = "Month",
                   breaks = seq(1, 12, 1),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Cumulative visit duration (sec)",
                     limits = c(0, 300000),
                     breaks = seq(0, 300000, 100000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Subplot for non-Pteropus cumulative duration
obs_month_sum_durstayB <- infra_stay %>%
  filter(pnp == "Non-pteropus") %>%
  group_by(month = month(as.Date(dtobs, format = "%m/%d/%Y")), Month_Season) %>%
  summarize(cum_duration = sum(durStay)) %>%
  ungroup() %>%
  add_row(month = 6, Month_Season = "Monsoon", cum_duration = 0) %>%
  # Bar plot for cumulative duration of visits by month, with labels
  ggplot() +
  geom_col(aes(x = month, y = cum_duration, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = month, y = cum_duration + (500000/25), label = cum_duration), size = 3) +
  # Rename numbered months as abbreviated names
  scale_x_continuous(name = "Month",
                     breaks = seq(1, 12, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Cumulative visit duration (sec)",
                     limits = c(0, 500000),
                     breaks = seq(0, 500000, 100000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Combine cumulative duration of bat visits subplots vertically (removing legends from both plots)
obs_month_sum_durstay <- plot_grid(obs_month_sum_durstayA + theme(legend.position = "none"),
                               obs_month_sum_durstayB + theme(legend.position = "none"),
                               nrow = 2, align = "v")
# Extract legend from first plot
obs_month_sum_durstay_legend <- g_legend(obs_month_sum_durstayA)
# Arrange combined subplots and legend
plot_grid(obs_month_sum_durstay, obs_month_sum_durstay_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/observed_month_cumulative_bat_visit_duration.png",
       height = 6, width = 8, dpi = 300, bg = "white")
ggsave("./Results/observed_month_cumulative_bat_visit_duration.tiff",
       height = 6, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: cumulative visit durations by season

# Subplot for Pteropus cumulative duration
obs_season_sum_durstayA <- infra_stay %>%
  filter(pnp == "Pteropus") %>%
  group_by(Month_Season) %>%
  summarize(cum_duration = sum(durStay)) %>%
  # Bar plot for cumulative duration of visits by season, with labels
  ggplot() +
  geom_col(aes(x = as.factor(Month_Season), y = cum_duration, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = as.factor(Month_Season), y = cum_duration + (400000/25), label = cum_duration), size = 3) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Cumulative visit duration (sec)",
                     limits = c(0, 400000),
                     breaks = seq(0, 400000, 100000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus cumulative duration
obs_season_sum_durstayB <- infra_stay %>%
  filter(pnp == "Non-pteropus") %>%
  group_by(Month_Season) %>%
  summarize(cum_duration = sum(durStay)) %>%
  # Bar plot for cumulative duration of visits by season, with labels
  ggplot() +
  geom_col(aes(x = as.factor(Month_Season), y = cum_duration, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = as.factor(Month_Season), y = cum_duration + (800000/25), label = cum_duration), size = 3) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Cumulative visit duration (sec)",
                     limits = c(0, 800000),
                     breaks = seq(0, 800000, 200000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine cumulative duration of bat visits subplots vertically (removing legends from both plots)
plot_grid(obs_season_sum_durstayA, obs_season_sum_durstayB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/observed_season_cumulative_bat_visit_duration.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/observed_season_cumulative_bat_visit_duration.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: sap contact duration by month

# Subplot for Pteropus sap contact duration
obs_month_durcontA <- ggplot(infra_cont %>% filter(pnp == "Pteropus"),
                             aes(x = month(as.Date(dtobs, format = "%m/%d/%Y")), y = durCont,
                                 color = Month_Season, fill = Month_Season,
                                 group = month(as.Date(dtobs, format = "%m/%d/%Y")))) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename numbered months as abbreviated names
  scale_x_continuous(name = "Month", 
                     breaks = seq(1, 12, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Duration of sap contact (sec)", 
                     limits = c(0, 5000),
                     breaks = seq(0, 5000, 1000)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Subplot for non-Pteropus sap contact duration
obs_month_durcontB <- ggplot(infra_cont %>% filter(pnp == "Non-pteropus"),
                             aes(x = month(as.Date(dtobs, format = "%m/%d/%Y")), y = durCont,
                                 color = Month_Season, fill = Month_Season,
                                 group = month(as.Date(dtobs, format = "%m/%d/%Y")))) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename numbered months as abbreviated names
  scale_x_continuous(name = "Month", 
                     breaks = seq(1, 12, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Duration of sap contact (sec)", 
                     limits = c(0, 800),
                     breaks = seq(0, 800, 200)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Combine duration of bat sap contact subplots vertically (removing legends from both plots)
obs_month_durcont <- plot_grid(obs_month_durcontA + theme(legend.position = "none"),
                               obs_month_durcontB + theme(legend.position = "none"),
                               nrow = 2, align = "v")
# Extract legend from first plot
obs_month_durcont_legend <- g_legend(obs_month_durcontA)
# Arrange combined subplots and legend
plot_grid(obs_month_durcont, obs_month_durcont_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/observed_month_bat_sap_contact_duration.png",
       height = 6, width = 8, dpi = 300, bg = "white")
ggsave("./Results/observed_month_bat_sap_contact_duration.tiff",
       height = 6, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: sap contact duration by season

# Subplot for Pteropus sap contact duration
obs_season_durcontA <- ggplot(infra_cont %>% filter(pnp == "Pteropus"),
                              aes(x = as.factor(Month_Season), y = durCont,
                                  color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Duration of sap contact (sec)",
                     limits = c(0, 5000),
                     breaks = seq(0, 5000, 1000)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus sap contact duration
obs_season_durcontB <- ggplot(infra_cont %>% filter(pnp == "Non-pteropus"),
                              aes(x = as.factor(Month_Season), y = durCont,
                                  color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Duration of sap contact (sec)",
                     limits = c(0, 800),
                     breaks = seq(0, 800, 200)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine duration of bat sap contact subplots vertically (removing legends from both plots)
plot_grid(obs_season_durcontA, obs_season_durcontB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/observed_season_bat_sap_contact_duration.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/observed_season_bat_sap_contact_duration.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: cumulative sap contact durations by month

# Subplot for Pteropus cumulative duration
obs_month_sum_durcontA <- infra_stay %>%
  filter(pnp == "Pteropus") %>%
  group_by(month = month(as.Date(dtobs, format = "%m/%d/%Y")), Month_Season) %>%
  summarize(cum_contact = sum(durCont)) %>%
  # Bar plot for cumulative duration of sap contacts by month, with labels
  ggplot() +
  geom_col(aes(x = month, y = cum_contact, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = month, y = cum_contact + (30000/25), label = cum_contact), size = 3) +
  # Rename numbered months as abbreviated names
  scale_x_continuous(name = "Month",
                     breaks = seq(1, 12, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Cumulative sap contact duration (sec)",
                     limits = c(0, 30000),
                     breaks = seq(0, 30000, 10000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Subplot for non-Pteropus cumulative duration
obs_month_sum_durcontB <- infra_stay %>%
  filter(pnp == "Non-pteropus") %>%
  group_by(month = month(as.Date(dtobs, format = "%m/%d/%Y")), Month_Season) %>%
  summarize(cum_contact = sum(durCont)) %>%
  ungroup() %>%
  add_row(month = 6, Month_Season = "Monsoon", cum_contact = 0) %>%
  # Bar plot for cumulative duration of sap contacts by month, with labels
  ggplot() +
  geom_col(aes(x = month, y = cum_contact, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = month, y = cum_contact + (500000/25), label = cum_contact), size = 3) +
  # Rename numbered months as abbreviated names
  scale_x_continuous(name = "Month",
                     breaks = seq(1, 12, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Cumulative sap contact duration (sec)",
                     limits = c(0, 500000),
                     breaks = seq(0, 500000, 100000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Combine cumulative duration of bat sap contacts subplots vertically (removing legends from both plots)
obs_month_sum_durcont <- plot_grid(obs_month_sum_durcontA + theme(legend.position = "none"),
                                   obs_month_sum_durcontB + theme(legend.position = "none"),
                                   nrow = 2, align = "v")
# Extract legend from first plot
obs_month_sum_durcont_legend <- g_legend(obs_month_sum_durcontA)
# Arrange combined subplots and legend
plot_grid(obs_month_sum_durcont, obs_month_sum_durcont_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/observed_month_cumulative_bat_sap_contact_duration.png",
       height = 7, width = 8, dpi = 300, bg = "white")
ggsave("./Results/observed_month_cumulative_bat_sap_contact_duration.tiff",
       height = 7, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for raw data analysis: cumulative sap contact durations by season

# Subplot for Pteropus cumulative duration
obs_season_sum_durcontA <- infra_stay %>%
  filter(pnp == "Pteropus") %>%
  group_by(Month_Season) %>%
  summarize(cum_contact = sum(durCont)) %>%
  # Bar plot for cumulative duration of sap contacts by season, with labels
  ggplot() +
  geom_col(aes(x = as.factor(Month_Season), y = cum_contact, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = as.factor(Month_Season), y = cum_contact + (60000/25), label = cum_contact), size = 3) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Cumulative sap contact duration (sec)",
                     limits = c(0, 60000),
                     breaks = seq(0, 60000, 10000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus cumulative duration
obs_season_sum_durcontB <- infra_stay %>%
  filter(pnp == "Non-pteropus") %>%
  group_by(Month_Season) %>%
  summarize(cum_contact = sum(durCont)) %>%
  # Bar plot for cumulative duration of sap contacts by season, with labels
  ggplot() +
  geom_col(aes(x = as.factor(Month_Season), y = cum_contact, fill = Month_Season), width = 0.5) +
  geom_text(aes(x = as.factor(Month_Season), y = cum_contact + (800000/25), label = cum_contact), size = 3) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Cumulative sap contact duration (sec)",
                     limits = c(0, 800000),
                     breaks = seq(0, 800000, 200000)) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine cumulative duration of bat sap contacts subplots vertically (removing legends from both plots)
plot_grid(obs_season_sum_durcontA, obs_season_sum_durcontB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/observed_season_cumulative_bat_sap_contact_duration.png",
       height = 7, width = 6, dpi = 300, bg = "white")
ggsave("./Results/observed_season_cumulative_bat_sap_contact_duration.tiff",
       height = 7, width = 6, dpi = 300, bg = "white", compression = "lzw")

##################################
### Modeling of bat visit data ###
##################################

# Gather bat species visits together into one column with species names in a second column
infra_day_gather <- melt(infra_day_spread, measure.vars = c("Pteropus", "non_Pteropus", "unidentified"))

# Add a column for ordering seasons so that winter is the reference
infra_stay <- mutate(infra_stay, season_order = case_when(Month_Season == "Winter" ~ 1,
                                                                Month_Season == "Spring" ~ 2,
                                                                Month_Season == "Monsoon" ~ 3,
                                                                Month_Season == "Postmonsoon" ~ 4))
infra_cont <- mutate(infra_cont, season_order = case_when(Month_Season == "Winter" ~ 1,
                                                          Month_Season == "Spring" ~ 2,
                                                          Month_Season == "Monsoon" ~ 3,
                                                          Month_Season == "Postmonsoon" ~ 4))
infra_day_spread <- mutate(infra_day_spread, season_order = case_when(Month_Season == "Winter" ~ 1,
                                                                      Month_Season == "Spring" ~ 2,
                                                                      Month_Season == "Monsoon" ~ 3,
                                                                      Month_Season == "Postmonsoon" ~ 4))

## Modeling the difference between bat species in average number of visits

# Model the number of visits per tree per day as a negative binomial variable
species_visits_mod <- glm(formula = value ~ variable, data = infra_day_gather, na.action = "na.fail",
                          family = negative.binomial(1))
# Model summary, coefficients, and predictions
species_visits_mod_coef <- tidy(species_visits_mod, conf.int = TRUE)
species_visits_mod_coef$term <- c("Pteropus", "non_Pteropus", "unidentified")
species_visits_mod_pred <- ggpredict(model = species_visits_mod, terms = "variable")
species_visits_mod_pred$x <- c("Pteropus", "non_Pteropus", "unidentified")
species_visits_mod_out <- merge(species_visits_mod_coef, species_visits_mod_pred,
                                by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(species_visits_mod_out, "./Results/species_visits_model_output.csv")
# Post-hoc hypothesis test for effect of species
species_visits_mod_pairwise <- emmeans(species_visits_mod, ~variable)
pairs(species_visits_mod_pairwise)
# Output results to file
write.csv(pairs(species_visits_mod_pairwise), "./Results/species_visits_model_pairwise.csv")

## Modeling the difference between bat species in duration of visits

# Model the duration of visits as an inverse Gaussian variable
species_visit_duration_mod <- glm(formula = durStay ~ pnp, data = infra_stay, na.action = "na.fail",
                                  family = inverse.gaussian(link = "identity"))
# Model summary, coefficients, and predictions
species_visit_duration_mod_coef <- tidy(species_visit_duration_mod, conf.int = TRUE)
species_visit_duration_mod_coef$term <- c("non_Pteropus", "Pteropus", "unidentified")
species_visit_duration_mod_pred <- ggpredict(model = species_visit_duration_mod, terms = "pnp")
species_visit_duration_mod_pred$x <- c("non_Pteropus", "Pteropus", "unidentified")
species_visit_duration_mod_out <- merge(species_visit_duration_mod_coef, species_visit_duration_mod_pred,
                                        by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(species_visit_duration_mod_out, "./Results/species_visit_duration_model_output.csv")
# Post-hoc hypothesis test for effect of species
species_visit_duration_mod_pairwise <- emmeans(species_visit_duration_mod, ~pnp)
pairs(species_visit_duration_mod_pairwise)
# Output results to file
write.csv(pairs(species_visit_duration_mod_pairwise), "./Results/species_visit_duration_model_pairwise.csv")

## Modeling the difference between bat species in duration of sap contact

# Model the duration of sap contact as an inverse Gaussian variable
species_contact_duration_mod <- glm(formula = durCont ~ pnp, data = infra_cont, na.action = "na.fail",
                                    family = inverse.gaussian(link = "identity"))
# Model summary, coefficients, and predictions
species_contact_duration_mod_coef <- tidy(species_contact_duration_mod, conf.int = TRUE)
species_contact_duration_mod_coef$term <- c("non_Pteropus", "Pteropus", "unidentified")
species_contact_duration_mod_pred <- ggpredict(model = species_contact_duration_mod, terms = "pnp")
species_contact_duration_mod_pred$x <- c("non_Pteropus", "Pteropus", "unidentified")
species_contact_duration_mod_out <- merge(species_contact_duration_mod_coef, species_contact_duration_mod_pred,
                                          by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(species_contact_duration_mod_out, "./Results/species_contact_duration_model_output.csv")
# Post-hoc hypothesis test for effect of species
species_contact_duration_mod_pairwise <- emmeans(species_contact_duration_mod, ~pnp)
pairs(species_contact_duration_mod_pairwise)
# Output results to file
write.csv(pairs(species_contact_duration_mod_pairwise), "./Results/species_contact_duration_model_pairwise.csv")

## Modeling the difference between seasons in average number of visits

# Model the number of Pteropus visits per tree per day as a negative binomial variable
Pteropus_season_visits_mod <- glm(formula = Pteropus ~ as.factor(season_order), data = infra_day_spread, na.action = "na.fail",
                                  family = negative.binomial(1))
# Model summary, coefficients, and predictions
Pteropus_season_visits_mod_coef <- tidy(Pteropus_season_visits_mod, conf.int = TRUE)
Pteropus_season_visits_mod_coef$term <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
Pteropus_season_visits_mod_pred <- ggpredict(model = Pteropus_season_visits_mod, terms = "season_order")
Pteropus_season_visits_mod_pred$x <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
Pteropus_season_visits_mod_out <- merge(Pteropus_season_visits_mod_coef, Pteropus_season_visits_mod_pred,
                                        by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(Pteropus_season_visits_mod_out, "./Results/Pteropus_season_visits_model_output.csv")
# Post-hoc hypothesis test for effect of seasons for Pteropus
Pteropus_season_visits_mod_pairwise <- emmeans(Pteropus_season_visits_mod, ~as.factor(season_order))
pairs(Pteropus_season_visits_mod_pairwise)
# Output results to file
write.csv(pairs(Pteropus_season_visits_mod_pairwise), "./Results/Pteropus_season_visits_model_pairwise.csv")

# Model the number of non-Pteropus visits per tree per day as a negative binomial variable
non_Pteropus_season_visits_mod <- glm(formula = non_Pteropus ~ as.factor(season_order), data = infra_day_spread, na.action = "na.fail",
                                         family = negative.binomial(1))
# Model summary, coefficients, and predictions
non_Pteropus_season_visits_mod_coef <- tidy(non_Pteropus_season_visits_mod, conf.int = TRUE)
non_Pteropus_season_visits_mod_coef$term <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
non_Pteropus_season_visits_mod_pred <- ggpredict(model = non_Pteropus_season_visits_mod, terms = "season_order")
non_Pteropus_season_visits_mod_pred$x <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
non_Pteropus_season_visits_mod_out <- merge(non_Pteropus_season_visits_mod_coef, non_Pteropus_season_visits_mod_pred,
                                            by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(non_Pteropus_season_visits_mod_out, "./Results/non_Pteropus_season_visits_model_output.csv")
# Post-hoc hypothesis test for effect of seasons for non-Pteropus
non_Pteropus_season_visits_mod_pairwise <- emmeans(non_Pteropus_season_visits_mod, ~as.factor(season_order))
pairs(non_Pteropus_season_visits_mod_pairwise)
# Output results to file
write.csv(pairs(non_Pteropus_season_visits_mod_pairwise), "./Results/non_Pteropus_season_visits_model_pairwise.csv")

# Subplot for Pteropus mean visits
mod_season_mean_visitsA <- ggplot(Pteropus_season_visits_mod_out, aes(x = term, y = predicted, color = term)) +
  # Points are the model estimates of the mean of visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = term, ymin = conf.low.y, ymax = conf.high.y),
                 size = 1) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean number of bat visits",
                     limits = c(0, 3),
                     breaks = seq(0, 3, 0.5)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus mean visits
mod_season_mean_visitsB <- ggplot(non_Pteropus_season_visits_mod_out, aes(x = term, y = predicted, color = term)) +
  # Points are the model estimates of the mean of visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = term, ymin = conf.low.y, ymax = conf.high.y),
                 size = 1) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean number of bat visits",
                     limits = c(0, 150),
                     breaks = seq(0, 150, 25)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine mean visits subplots vertically
plot_grid(mod_season_mean_visitsA, mod_season_mean_visitsB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/NBmodel_season_mean_bat_visits.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/NBmodel_season_mean_bat_visits.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

## Modeling the difference between seasons in average visit duration

# Model the duration of Pteropus visits as an inverse Gaussian variable
Pteropus_visit_duration_mod <- glm(formula = durStay ~ as.factor(season_order), data = infra_stay %>% filter(pnp == "Pteropus"), na.action = "na.fail",
                                   family = inverse.gaussian(link = "identity"))
# Model summary, coefficients, and predictions
Pteropus_visit_duration_mod_coef <- tidy(Pteropus_visit_duration_mod, conf.int = TRUE)
Pteropus_visit_duration_mod_coef$term <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
Pteropus_visit_duration_mod_pred <- ggpredict(model = Pteropus_visit_duration_mod, terms = "season_order")
Pteropus_visit_duration_mod_pred$x <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
Pteropus_visit_duration_mod_out <- merge(Pteropus_visit_duration_mod_coef, Pteropus_visit_duration_mod_pred,
                                         by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(Pteropus_visit_duration_mod_out, "./Results/Pteropus_season_visit_duration_model_output.csv")
# Post-hoc hypothesis test for effect of seasons for Pteropus
Pteropus_visit_duration_mod_pairwise <- emmeans(Pteropus_visit_duration_mod, ~as.factor(season_order))
pairs(Pteropus_visit_duration_mod_pairwise)
# Output results to file
write.csv(pairs(Pteropus_visit_duration_mod_pairwise), "./Results/Pteropus_season_visit_duration_model_pairwise.csv")

# Model the duration of non-Pteropus visits as an inverse Gaussian variable
non_Pteropus_visit_duration_mod <- glm(formula = durStay ~ as.factor(season_order), data = infra_stay %>% filter(pnp == "Non-pteropus"), na.action = "na.fail",
                                       family = inverse.gaussian(link = "identity"))
# Model summary, coefficients, and predictions
non_Pteropus_visit_duration_mod_coef <- tidy(non_Pteropus_visit_duration_mod, conf.int = TRUE)
non_Pteropus_visit_duration_mod_coef$term <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
non_Pteropus_visit_duration_mod_pred <- ggpredict(model = non_Pteropus_visit_duration_mod, terms = "season_order")
non_Pteropus_visit_duration_mod_pred$x <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
non_Pteropus_visit_duration_mod_out <- merge(non_Pteropus_visit_duration_mod_coef, non_Pteropus_visit_duration_mod_pred,
                                             by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(non_Pteropus_visit_duration_mod_out, "./Results/non_Pteropus_season_visit_duration_model_output.csv")
# Post-hoc hypothesis test for effect of seasons for non-Pteropus
non_Pteropus_visit_duration_mod_pairwise <- emmeans(non_Pteropus_visit_duration_mod, ~as.factor(season_order))
pairs(non_Pteropus_visit_duration_mod_pairwise)
# Output results to file
write.csv(pairs(non_Pteropus_visit_duration_mod_pairwise), "./Results/non_Pteropus_season_visit_duration_model_pairwise.csv")

# Subplot for Pteropus mean visit duration
mod_season_mean_visit_durationA <- ggplot(Pteropus_visit_duration_mod_out, aes(x = term, y = predicted, color = term)) +
  # Points are the model estimates of the mean duration of visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = term, ymin = conf.low.y, ymax = conf.high.y),
                 size = 1) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean visit duration (sec)",
                     limits = c(0, 1600),
                     breaks = seq(0, 1600, 400)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus mean visit duration
mod_season_mean_visit_durationB <- ggplot(non_Pteropus_visit_duration_mod_out, aes(x = term, y = predicted, color = term)) +
  # Points are the model estimates of the mean duration of visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = term, ymin = conf.low.y, ymax = conf.high.y),
                 size = 1) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean visit duration (sec)",
                     limits = c(0, 60),
                     breaks = seq(0, 60, 20)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine mean visit duration subplots vertically
plot_grid(mod_season_mean_visit_durationA, mod_season_mean_visit_durationB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/IGmodel_season_mean_bat_visit_duration.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/IGmodel_season_mean_bat_visit_duration.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

## Modeling the difference between seasons in average sap contact duration

# Model the duration of Pteropus sap contact as an inverse Gaussian variable
Pteropus_contact_duration_mod <- glm(formula = durCont ~ as.factor(season_order), data = infra_cont %>% filter(pnp == "Pteropus"), na.action = "na.fail",
                                     family = inverse.gaussian(link = "identity"))
# Model summary, coefficients, and predictions
Pteropus_contact_duration_mod_coef <- tidy(Pteropus_contact_duration_mod, conf.int = TRUE)
Pteropus_contact_duration_mod_coef$term <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
Pteropus_contact_duration_mod_pred <- ggpredict(model = Pteropus_contact_duration_mod, terms = "season_order")
Pteropus_contact_duration_mod_pred$x <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
Pteropus_contact_duration_mod_out <- merge(Pteropus_contact_duration_mod_coef, Pteropus_contact_duration_mod_pred,
                                           by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(Pteropus_contact_duration_mod_out, "./Results/Pteropus_season_sap_contact_duration_model_output.csv")
# Post-hoc hypothesis test for effect of seasons for Pteropus
Pteropus_contact_duration_mod_pairwise <- emmeans(Pteropus_contact_duration_mod, ~as.factor(season_order))
pairs(Pteropus_contact_duration_mod_pairwise)
# Output results to file
write.csv(pairs(Pteropus_contact_duration_mod_pairwise), "./Results/Pteropus_season_contact_duration_model_pairwise.csv")

# Model the duration of non-Pteropus sap contact as an inverse Gaussian variable
non_Pteropus_contact_duration_mod <- glm(formula = durCont ~ as.factor(season_order), data = infra_cont %>% filter(pnp == "Non-pteropus"), na.action = "na.fail",
                                         family = inverse.gaussian(link = "identity"))
# Model summary, coefficients, and predictions
non_Pteropus_contact_duration_mod_coef <- tidy(non_Pteropus_contact_duration_mod, conf.int = TRUE)
non_Pteropus_contact_duration_mod_coef$term <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
non_Pteropus_contact_duration_mod_pred <- ggpredict(model = non_Pteropus_contact_duration_mod, terms = "season_order")
non_Pteropus_contact_duration_mod_pred$x <- c("Winter", "Spring", "Monsoon", "Postmonsoon")
non_Pteropus_contact_duration_mod_out <- merge(non_Pteropus_contact_duration_mod_coef, non_Pteropus_contact_duration_mod_pred,
                                               by.x = "term", by.y = "x", all = TRUE)
# Output results to file
write.csv(non_Pteropus_contact_duration_mod_out, "./Results/non_Pteropus_season_sap_contact_duration_model_output.csv")
# Post-hoc hypothesis test for effect of seasons for non-Pteropus
non_Pteropus_contact_duration_mod_pairwise <- emmeans(non_Pteropus_contact_duration_mod, ~as.factor(season_order))
pairs(non_Pteropus_contact_duration_mod_pairwise)
# Output results to file
write.csv(pairs(non_Pteropus_contact_duration_mod_pairwise), "./Results/non_Pteropus_season_contact_duration_mod_pairwise.csv")

# Subplot for Pteropus mean sap contact duration
mod_season_mean_contact_durationA <- ggplot(Pteropus_contact_duration_mod_out, aes(x = term, y = predicted, color = term)) +
  # Points are the model estimates of the mean duration of sap contact
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = term, ymin = conf.low.y, ymax = conf.high.y),
                 size = 1) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean sap contact duration (sec)",
                     limits = c(0, 400),
                     breaks = seq(0, 400, 100)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus mean sap contact duration
mod_season_mean_contact_durationB <- ggplot(non_Pteropus_contact_duration_mod_out, aes(x = term, y = predicted, color = term)) +
  # Points are the model estimates of the mean duration of sap contact
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = term, ymin = conf.low.y, ymax = conf.high.y),
                 size = 1) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean sap contact duration (sec)",
                     limits = c(0, 60),
                     breaks = seq(0, 60, 20)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine mean contact duration subplots vertically
plot_grid(mod_season_mean_contact_durationA, mod_season_mean_contact_durationB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/IGmodel_season_mean_bat_contact_duration.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/IGmodel_season_mean_bat_contact_duration.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

# Create a reduced data frame that drops rows with NAs for predictor variables
infra_day_spread_noNA <- drop_na(infra_day_spread, c("treeHeight", "shavearea"))

## Demonstration of strong correlation between tree height and tree age

# Collect only unique tree heights and tree ages
infra_day_spread_filtertrees <- infra_day_spread %>%
  ungroup() %>%
  dplyr::select(treeHeight, treeAge) %>%
  distinct(treeHeight, treeAge)

# Make linear model for tree height and tree age
treeage_lm <- lm(treeAge ~ treeHeight, infra_day_spread_filtertrees, na.action = "na.fail")

# Plot correlation between tree height and tree age
ggplot(infra_day_spread_filtertrees, aes(x = treeHeight, y = treeAge)) +
  # Fit a linear model to the data
  stat_smooth(method = "lm", color = "#808080", alpha = 0.5, se = TRUE) +
  # Plot points over fit
  geom_point(shape = 1) +
  # Add R2 value
  annotate(geom = "text", x = 500, y = 35, label = paste("R squared =", round(glance(treeage_lm)$adj.r.squared, 2))) +
  # Add p-value
  annotate(geom = "text", x = 500, y = 33, label = paste("P <", min(0.001, formatC(tidy(treeage_lm)$p.value[2], digits = 2)))) +
  # Set x axis
  scale_x_continuous(name = "Tree height (cm)", breaks = seq(250, 1250, 250), limits = c(250, 1250)) +
  # Set y axis
  scale_y_continuous(name = "Tree age (yr)", breaks = seq(10, 40, 10), limits = c(10, 40)) +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12)
# Output plot to file
ggsave("./Results/tree_height+age_correlation.png",
       height = 4, width = 4, dpi = 300, bg = "white")
ggsave("./Results/tree_height+age_correlation.tiff",
       height = 4, width = 4, dpi = 300, bg = "white", compression = "lzw")

## Demonstration of yearly temperature trend

# Tests to show significant variation in minimum temperature by month and season
KW_month <- tidy(kruskal.test(min_temp ~ as.factor(month), infra_day_spread_noNA))
KW_season <- tidy(kruskal.test(min_temp ~ Month_Season, infra_day_spread_noNA))

# Plot distribution of minimum temperatures across months
min_temp_trendA <- ggplot(infra_day_spread_noNA, aes(x = as.factor(month), y = min_temp,
                                                     color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Add Kruskal-Wallis statistic
  annotate(geom = "text", x = 8, y = 13, label = paste("Kruskal-Wallis H =", round(KW_month$statistic, 1))) +
  # Add df
  annotate(geom = "text", x = 8, y = 10, label = paste("df =", KW_month$parameter)) +
  # Add p-value
  annotate(geom = "text", x = 8, y = 7, label = paste("P <", min(0.001, formatC(KW_month$p.value, digits = 2)))) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Minimum nightly temperature (ºC)",
                     limits = c(0, 40),
                     breaks = seq(0, 40, 20)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Plot distribution of minimum temperatures across seasons
min_temp_trendB <- ggplot(infra_day_spread_noNA, aes(x = Month_Season, y = min_temp,
                                                     color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .2, y = 0)) +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Add Kruskal-Wallis statistic
  annotate(geom = "text", x = 3, y = 13, label = paste("Kruskal-Wallis H =", round(KW_season$statistic, 1))) +
  # Add df
  annotate(geom = "text", x = 3, y = 10, label = paste("df =", KW_season$parameter)) +
  # Add p-value
  annotate(geom = "text", x = 3, y = 7, label = paste("P <", min(0.001, formatC(KW_season$p.value, digits = 2)))) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)")) +
  # Set y axis
  scale_y_continuous(name = "Minimum nightly temperature (ºC)",
                     limits = c(0, 40),
                     breaks = seq(0, 40, 20)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine minimum temperature subplots vertically
plot_grid(min_temp_trendA, min_temp_trendB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/minimum_temperature_trend.png",
       height = 7, width = 8, dpi = 300, bg = "white")
ggsave("./Results/minimum_temperature_trend.tiff",
       height = 7, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Model selection for Pteropus visits to trees per night, all seasons

# Create global model with all potential explanatory variables
Pteropus_all_mod_nb <- glm(formula = Pteropus ~ as.factor(season_order)+
                             treeHeight+
                             shavearea+
                             as.numeric(days_since_shaving) - 1,
                           data = infra_day_spread_noNA,
                           family = negative.binomial(1), na.action = "na.fail", maxit = 1000)
# Run model selection by AICc
Pteropus_all_mod_nb_select <- dredge(Pteropus_all_mod_nb, rank = "AICc")
# Select top models by delta AICc
Pteropus_all_mod_nb_top <- get.models(Pteropus_all_mod_nb_select, subset = delta < 2)

## Run repeated K-fold and leave-one-out cross-validation and output results within model selection table
# Create empty columns for results
Pteropus_all_mod_nb_select$cv.mean.rmse <- NA
Pteropus_all_mod_nb_select$cv.mean.adj.rmse <- NA
Pteropus_all_mod_nb_select$loo.rmse <- NA
# Set seed
set.seed(1234)
# Loop through top models
for(i in 1:length(Pteropus_all_mod_nb_top)){
  # Run repeated K-fold cross-validation and calculated average root mean square error across repeats
  out <- sqrt(cv.glm(data = infra_day_spread_noNA, glmfit = Pteropus_all_mod_nb_top[[i]], K=10)$delta)
  # Output average prediction error
  Pteropus_all_mod_nb_select$cv.mean.rmse[i] <- out[1]
  Pteropus_all_mod_nb_select$cv.mean.adj.rmse[i] <- out[2]
  # Calculate root mean square error from leave-one-out cross-validation
  Pteropus_all_mod_nb_select$loo.rmse[i] <- loo(Pteropus_all_mod_nb_top[[i]], type = "rmse")
}
# Write output to file
write.csv(as.data.frame(Pteropus_all_mod_nb_select), "./Results/Pteropus_visits_NBmodel_selection.csv")

# Output coefficients and predictions from best model by AICc and cross-validation
Pteropus_all_mod_nb_top_coef <- tidy(Pteropus_all_mod_nb_top[[2]], conf.int = TRUE)
Pteropus_all_mod_nb_top_coef$term <- c("winter", "spring", "monsoon", "post-monsoon", "shavearea", "treeHeight")
# Write output to file
write.csv(Pteropus_all_mod_nb_top_coef, "./Results/Pteropus_visits_NBmodel_best_coef.csv")
# Post-hoc hypothesis test for effect of seasons for Pteropus
Pteropus_all_mod_nb_top_pairwise <- emmeans(Pteropus_all_mod_nb_top[[2]], ~as.factor(season_order))
pairs(Pteropus_all_mod_nb_top_pairwise)
# Output results to file
write.csv(pairs(Pteropus_all_mod_nb_top_pairwise), "./Results/Pteropus_visits_NBmodel_best_pairwise.csv")

# Make predictions for top model
Pteropus_all_mod_nb_top_pred1 <- as.data.frame(ggpredict(Pteropus_all_mod_nb_top[[2]],
                                                        c("season_order",
                                                          "treeHeight [457.20, 591.82, 703.58]",
                                                          "shavearea [1131]")))
Pteropus_all_mod_nb_top_pred2 <- as.data.frame(ggpredict(Pteropus_all_mod_nb_top[[2]],
                                                         c("season_order",
                                                           "shavearea [960, 1131, 1320]",
                                                           "treeHeight [591.82]")))
# Write output to file
write.csv(Pteropus_all_mod_nb_top_pred1, "./Results/Pteropus_visits_NBmodel_pred_season+treeheight.csv")
write.csv(Pteropus_all_mod_nb_top_pred2, "./Results/Pteropus_visits_NBmodel_pred_season+shavearea.csv")

# Model predictions
infra_day_spread_noNA$Pteropus_pred <- predict(Pteropus_all_mod_nb_top[[2]], type = "response")
m_Pteropus <- infra_day_spread_noNA %>%
  dplyr::select(Pteropus, Pteropus_pred) %>%
  gather()

# Model performance (prediction interval)
(Pteropus_mod_perform <- ggplot(m_Pteropus, aes(x = key, y = value)) +
    geom_point(position = position_jitter(width = 0.15, height = 0), shape = 21, size = 0.1) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), color = "black", fill = "black") +
    stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
                 size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
    scale_x_discrete(name = "Value", labels = c("observed", "predicted")) +
    scale_y_continuous(name = "Number of bat visits", breaks = seq(0, 40, 10), limits = c(0, 40)) +
    ggtitle("Pteropus") +
    theme_cowplot(font_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)))

# Plot predictions for season + tree height
(mod_Pteropus1 <- ggplot(Pteropus_all_mod_nb_top_pred1, aes(x = as.factor(x), y = predicted, group = group, color = group)) +
  # Estimated means for each level of tree height
  geom_point(size = 2, position = set_dodge) +
  # Estimated confidence intervals for each level of tree height
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), size = 1, position = set_dodge) +
  # Set the color scale for levels of tree height
  scale_color_manual(name = "Tree height (cm)",
                     values = c("grey60", "grey40", "black"),
                     labels = c("25th %ile", "50th %ile", "75th %ile")) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("1", "2", "3", "4"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean Pteropus visits", breaks = seq(0, 5, 1), limits = c(0, 5.5)) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Center horizontal legend at bottom
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5)))

# Plot predictions for season + shaved area
(mod_Pteropus2 <- ggplot(Pteropus_all_mod_nb_top_pred2, aes(x = as.factor(x), y = predicted, group = group, color = group)) +
    # Estimated means for each level of tree height
    geom_point(size = 2, position = set_dodge) +
    # Estimated confidence intervals for each level of tree height
    geom_linerange(aes(ymin = conf.low, ymax = conf.high), size = 1, position = set_dodge) +
    # Set the color scale for levels of shaved area
    scale_color_manual(name = "Shaved area (sq cm)",
                       values = c("grey60", "grey40", "black"),
                       labels = c("25th %ile", "50th %ile", "75th %ile")) +
    # Rename and order seasons
    scale_x_discrete(name = "Season", limits = c("1", "2", "3", "4"),
                     labels = c("winter", "spring", "monsoon", "post-monsoon")) +
    # Set y axis
    scale_y_continuous(name = "Mean Pteropus visits", breaks = seq(0, 4, 1), limits = c(0, 4)) +
    # Plot theme (simplified axes in all black and white)
    theme_cowplot(font_size = 12) +
    # Center horizontal legend at bottom
    theme(legend.position = "bottom",
          legend.justification = "center",
          legend.direction = "horizontal"))

# Combine Pteropus subplots
(mod_Pteropus <- plot_grid(mod_Pteropus1, mod_Pteropus2, nrow = 2, align = "v", rel_heights = c(0.525, 0.475)))
# Output plot to file
ggsave("./Results/Pteropus_visits_NBmodel_pred_season+treeheight+shavearea.png",
       height = 5, width = 5, dpi = 300, bg = "white")
ggsave("./Results/Pteropus_visits_NBmodel_pred_season+treeheight+shavearea.tiff",
       height = 5, width = 5, dpi = 300, bg = "white", compression = "lzw")

## Model selection for non-Pteropus visits to trees per night, all seasons

# Create global model with all potential explanatory variables
non_Pteropus_all_mod_nb <- glm(formula = non_Pteropus ~ as.factor(season_order)+
                                 treeHeight+
                                 shavearea+
                                 as.numeric(days_since_shaving) - 1,
                               data = infra_day_spread_noNA,
                               family = negative.binomial(1), na.action = "na.fail", maxit = 1000)
summary(non_Pteropus_all_mod_nb)
# Run model selection by AICc
non_Pteropus_all_mod_nb_select <- dredge(non_Pteropus_all_mod_nb, rank = "AICc")
# Select top models by delta AICc
non_Pteropus_all_mod_nb_top <- get.models(non_Pteropus_all_mod_nb_select, subset = delta < 10)

## Run repeated K-fold and leave-one-out cross-validation and output results within model selection table
# Create empty columns for results
non_Pteropus_all_mod_nb_select$cv.mean.rmse <- NA
non_Pteropus_all_mod_nb_select$cv.mean.adj.rmse <- NA
non_Pteropus_all_mod_nb_select$loo.rmse <- NA
# Set seed
set.seed(1234)
# Loop through top models
for(i in 1:length(non_Pteropus_all_mod_nb_top)){
  # Run repeated K-fold cross-validation and calculated average root mean square error across repeats
  out <- sqrt(cv.glm(data = infra_day_spread_noNA, glmfit = non_Pteropus_all_mod_nb_top[[i]], K=10)$delta)
  # Output average prediction error
  non_Pteropus_all_mod_nb_select$cv.mean.rmse[i] <- out[1]
  non_Pteropus_all_mod_nb_select$cv.mean.adj.rmse[i] <- out[2]
  # Calculate root mean square error from leave-one-out cross-validation
  non_Pteropus_all_mod_nb_select$loo.rmse[i] <- loo(non_Pteropus_all_mod_nb_top[[i]], type = "rmse")
}
# Write output to file
write.csv(as.data.frame(non_Pteropus_all_mod_nb_select), "./Results/non_Pteropus_visits_NBmodel_selection.csv")

# Output coefficients and predictions from best model by AICc and cross-validation
non_Pteropus_all_mod_nb_top_coef <- tidy(non_Pteropus_all_mod_nb_top[[1]], conf.int = TRUE)
non_Pteropus_all_mod_nb_top_coef$term <- c("winter", "spring", "monsoon", "post-monsoon", "days_since_shaving", "shavearea")
# Write output to file
write.csv(non_Pteropus_all_mod_nb_top_coef, "./Results/non_Pteropus_visits_NBmodel_best_coef.csv")
# Post-hoc hypothesis test for effect of seasons for non-Pteropus
non_Pteropus_all_mod_nb_top_pairwise <- emmeans(non_Pteropus_all_mod_nb_top[[1]], ~as.factor(season_order))
pairs(non_Pteropus_all_mod_nb_top_pairwise)
# Output results to file
write.csv(pairs(non_Pteropus_all_mod_nb_top_pairwise), "./Results/non_Pteropus_visits_NBmodel_best_pairwise.csv")
# Make predictions for top model
non_Pteropus_all_mod_nb_top_pred1 <- as.data.frame(ggpredict(non_Pteropus_all_mod_nb_top[[1]],
                                                             c("season_order",
                                                               "shavearea [960, 1131, 1320]",
                                                               "days_since_shaving [3]")))
non_Pteropus_all_mod_nb_top_pred2 <- as.data.frame(ggpredict(non_Pteropus_all_mod_nb_top[[1]],
                                                             c("season_order",
                                                               "days_since_shaving [0, 3, 6]",
                                                               "shavearea [1131]")))
# Write output to file
write.csv(non_Pteropus_all_mod_nb_top_pred1, "./Results/non_Pteropus_visits_NBmodel_pred_season+shavearea.csv")
write.csv(non_Pteropus_all_mod_nb_top_pred2, "./Results/non_Pteropus_visits_NBmodel_pred_season+days_since_shaving.csv")

# Model predictions
infra_day_spread_noNA$non_Pteropus_pred <- predict(non_Pteropus_all_mod_nb_top[[1]], type = "response")
m_non_Pteropus <- infra_day_spread_noNA %>%
  dplyr::select(non_Pteropus, non_Pteropus_pred) %>%
  gather()

# Model performance (prediction interval)
(non_Pteropus_mod_perform <- ggplot(m_non_Pteropus, aes(x = key, y = value)) +
    geom_point(position = position_jitter(width = 0.15, height = 0), shape = 21, size = 0.1) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), color = "black", fill = "black") +
    stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
                 size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
    scale_x_discrete(name = "Value", labels = c("observed", "predicted")) +
    scale_y_continuous(name = "Number of bat visits", breaks = seq(0, 1500, 250), limits = c(0, 1500)) +
    ggtitle("Non-Pteropus") +
    theme_cowplot(font_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)))

# Plot predictions for season + shaved area
(mod_non_Pteropus1 <- ggplot(non_Pteropus_all_mod_nb_top_pred1, aes(x = as.factor(x), y = predicted, group = group, color = group)) +
  # Estimated means for each level of shaved area
  geom_point(size = 2, position = set_dodge) +
  # Estimated confidence intervals for each level of shaved area
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), size = 1, position = set_dodge) +
  # Set the color scale for levels of shaved area
  scale_color_manual(name = "Shaved area (sq cm)",
                     values = c("grey60", "grey40", "black"),
                     labels = c("25th %ile", "50th %ile", "75th %ile")) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("1", "2", "3", "4"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean non-Pteropus visits", breaks = seq(0, 150, 50), limits = c(0, 150)) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Center horizontal legend at bottom
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5)))

# Plot predictions for season + days since shaving
(mod_non_Pteropus2 <- ggplot(non_Pteropus_all_mod_nb_top_pred2, aes(x = as.factor(x), y = predicted, group = group, color = group)) +
  # Estimated means for each level of days since shaving
  geom_point(size = 2, position = set_dodge) +
  # Estimated confidence intervals for each level of days since shaving
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), size = 1, position = set_dodge) +
  # Set the color scale for levels of days since shaving
  scale_color_manual(name = "Days since shaving",
                     values = c("grey60", "grey40", "black")) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("1", "2", "3", "4"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean non-Pteropus visits", breaks = seq(0, 200, 50), limits = c(0, 200)) +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Center horizontal legend at bottom
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal"))

# Combine non-Pteropus subplots
(mod_non_Pteropus <- plot_grid(mod_non_Pteropus1, mod_non_Pteropus2, nrow = 2, align = "v", rel_heights = c(0.525, 0.475)))
# Output plot to file
ggsave("./Results/non_Pteropus_visits_NBmodel_pred_season+shavearea+days_since_shaving.png",
       height = 5, width = 5, dpi = 300, bg = "white")
ggsave("./Results/non_Pteropus_visits_NBmodel_pred_season+shavearea+days_since_shaving.tiff",
       height = 5, width = 5, dpi = 300, bg = "white", compression = "lzw")

# Combine Pteropus and non-Pteropus subplots
plot_grid(mod_Pteropus, mod_non_Pteropus, ncol = 2)
# Output plot to file
ggsave("./Results/all_bats_NBmodels.png",
       height = 5, width = 10, dpi = 300, bg = "white")
ggsave("./Results/all_bats_NBmodels.tiff",
       height = 5, width = 10, dpi = 300, bg = "white", compression = "lzw")

# Within winter, look at the additional effect of minimum temperature

# Create a reduced data frame that drops rows with NAs for predictor variables
infra_day_spread_noNA2 <- drop_na(infra_day_spread, c("treeHeight", "shavearea", "min_temp")) %>%
  filter(Month_Season == "Winter")

## Model selection for Pteropus visits to trees per night, winter only

# Create global model with all potential explanatory variables
Pteropus_all_mod_nb_winter <- glm(Pteropus ~ min_temp+
                                    treeHeight+
                                    shavearea+
                                    as.numeric(days_since_shaving),
                                  infra_day_spread_noNA2,
                                  family = negative.binomial(1), na.action = "na.fail", maxit = 1000)
# Run model selection by AICc
Pteropus_all_mod_nb_winter_select <- dredge(Pteropus_all_mod_nb_winter, rank = "AICc")
# Select top models by delta AICc
Pteropus_all_mod_nb_winter_top <- get.models(Pteropus_all_mod_nb_winter_select, subset = delta < 2)

## Run repeated K-fold and leave-one-out cross-validation and output results within model selection table
# Create empty columns for results
Pteropus_all_mod_nb_winter_select$cv.mean.rmse <- NA
Pteropus_all_mod_nb_winter_select$cv.mean.adj.rmse <- NA
Pteropus_all_mod_nb_winter_select$loo.rmse <- NA
# Set seed
set.seed(1234)
# Loop through top models
for(i in 1:length(Pteropus_all_mod_nb_winter_top)){
  # Run repeated K-fold cross-validation and calculated average root mean square error across repeats
  out <- sqrt(cv.glm(data = infra_day_spread_noNA2, glmfit = Pteropus_all_mod_nb_winter_top[[i]], K=10)$delta)
  # Output average prediction error
  Pteropus_all_mod_nb_winter_select$cv.mean.rmse[i] <- out[1]
  Pteropus_all_mod_nb_winter_select$cv.mean.adj.rmse[i] <- out[2]
  # Calculate root mean square error from leave-one-out cross-validation
  Pteropus_all_mod_nb_winter_select$loo.rmse[i] <- loo(Pteropus_all_mod_nb_winter_top[[i]], type = "rmse")
}
# Write output to file
write.csv(as.data.frame(Pteropus_all_mod_nb_winter_select), "./Results/Pteropus_visits_winter_NBmodel_selection.csv")

# Output coefficients and predictions from best model by AICc and cross-validation
Pteropus_all_mod_nb_top_coef_winter <- tidy(Pteropus_all_mod_nb_winter_top[[2]], conf.int = TRUE)
Pteropus_all_mod_nb_top_coef_winter$term <- c("intercept", "min_temp", "treeheight")
# Write output to file
write.csv(Pteropus_all_mod_nb_top_coef_winter, "./Results/Pteropus_visits_winter_NBmodel_best_coef.csv")
# Make predictions for top model
Pteropus_all_mod_nb_top_pred_winter <- as.data.frame(ggpredict(Pteropus_all_mod_nb_winter_top[[2]],
                                                               c("min_temp",
                                                                 "treeHeight [513.08]")))
# Write output to file
write.csv(Pteropus_all_mod_nb_top_pred_winter, "./Results/Pteropus_visits_winter_NBmodel_pred_mintemp.csv")

# Model predictions
infra_day_spread_noNA2$winter_Pteropus_pred <- predict(Pteropus_all_mod_nb_winter_top[[2]], type = "response")
m_winter_Pteropus <- infra_day_spread_noNA2 %>%
  dplyr::select(Pteropus, winter_Pteropus_pred) %>%
  gather()

# Model performance (prediction interval)
(winter_Pteropus_mod_perform <- ggplot(m_winter_Pteropus, aes(x = key, y = value)) +
    geom_point(position = position_jitter(width = 0.15, height = 0), shape = 21, size = 0.1) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), color = "black", fill = "black") +
    stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
                 size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
    scale_x_discrete(name = "Value", labels = c("observed", "predicted")) +
    scale_y_continuous(name = "Number of winter bat visits", breaks = seq(0, 40, 10), limits = c(0, 40)) +
    ggtitle("Pteropus (winter)") +
    theme_cowplot(font_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)))

# Plot predictions
(mod_winter_Pteropus <- ggplot(Pteropus_all_mod_nb_top_pred_winter, aes(x = x, y = predicted, group = group)) +
  # Estimated means for single level of tree height
  geom_line(color = "#808080", size = 1) +
  # Estimated confidence intervals for single level of tree height
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#808080", alpha = 0.5) +
  # Set x axis
  scale_x_continuous(name = "Minimum nightly temperature (ºC)", breaks = seq(8, 30, 2), limits = c(8, 30)) +
  # Set y axis
  scale_y_continuous(name = "Mean winter Pteropus visits", breaks = seq(0, 5.1, 1), limits = c(0, 5.1)) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none"))

# Output plot to file
ggsave("./Results/Pteropus_visits_winter_NBmodel_pred_mintemp.png",
       height = 4, width = 4, dpi = 300, bg = "white")
ggsave("./Results/Pteropus_visits_winter_NBmodel_pred_mintemp.tiff",
       height = 4, width = 4, dpi = 300, bg = "white", compression = "lzw")

## Model selection for non-Pteropus visits to trees per night, winter only

# Create global model with all potential explanatory variables
non_Pteropus_all_mod_nb_winter <- glm(non_Pteropus ~ min_temp+
                                        treeHeight+
                                        shavearea+
                                        as.numeric(days_since_shaving),
                                      infra_day_spread_noNA2,
                                      family = negative.binomial(1), na.action = "na.fail", maxit = 1000)
# Run model selection by AICc
non_Pteropus_all_mod_nb_winter_select <- dredge(non_Pteropus_all_mod_nb_winter, rank = "AICc")
# Select top models by delta AICc
non_Pteropus_all_mod_nb_winter_top <- get.models(non_Pteropus_all_mod_nb_winter_select, subset = delta < 2)

## Run repeated K-fold and leave-one-out cross-validation and output results within model selection table
# Create empty columns for results
non_Pteropus_all_mod_nb_winter_select$cv.mean.rmse <- NA
non_Pteropus_all_mod_nb_winter_select$cv.mean.adj.rmse <- NA
non_Pteropus_all_mod_nb_winter_select$loo.rmse <- NA
# Set seed
set.seed(1234)
# Loop through top models
for(i in 1:length(non_Pteropus_all_mod_nb_winter_top)){
  # Run repeated K-fold cross-validation and calculated average root mean square error across repeats
  out <- sqrt(cv.glm(data = infra_day_spread_noNA2, glmfit = non_Pteropus_all_mod_nb_winter_top[[i]], K=10)$delta)
  # Output average prediction error
  non_Pteropus_all_mod_nb_winter_select$cv.mean.rmse[i] <- out[1]
  non_Pteropus_all_mod_nb_winter_select$cv.mean.adj.rmse[i] <- out[2]
  # Calculate root mean square error from leave-one-out cross-validation
  non_Pteropus_all_mod_nb_winter_select$loo.rmse[i] <- loo(non_Pteropus_all_mod_nb_winter_top[[i]], type = "rmse")
}
# Write output to file
write.csv(as.data.frame(non_Pteropus_all_mod_nb_winter_select), "./Results/non_Pteropus_visits_winter_NBmodel_selection.csv")

# Output coefficients and predictions from best model by AICc and cross-validation
non_Pteropus_all_mod_nb_top_coef_winter <- cbind(as.data.frame(tidy(non_Pteropus_all_mod_nb_winter_top[[2]])), confint(non_Pteropus_all_mod_nb_winter_top[[2]]))
non_Pteropus_all_mod_nb_top_coef_winter$term <- c("intercept", "days_since_shaving", "min_temp")
rownames(non_Pteropus_all_mod_nb_top_coef_winter) <- c()
colnames(non_Pteropus_all_mod_nb_top_coef_winter)[6:7] <- c("conf.lower", "conf.upper")
# Write output to file
write.csv(non_Pteropus_all_mod_nb_top_coef_winter, "./Results/non_Pteropus_visits_winter_NBmodel_best_coef.csv")
# Make predictions for top model
non_Pteropus_all_mod_nb_top_pred_winter <- as.data.frame(ggpredict(non_Pteropus_all_mod_nb_winter_top[[2]],
                                                                   c("min_temp",
                                                                     "days_since_shaving [3]")))
# Write output to file
write.csv(non_Pteropus_all_mod_nb_top_pred_winter, "./Results/non_Pteropus_visits_winter_NBmodel_pred_mintemp.csv")

# Model predictions
infra_day_spread_noNA2$winter_non_Pteropus_pred <- predict(non_Pteropus_all_mod_nb_winter_top[[2]], type = "response")
m_winter_non_Pteropus <- infra_day_spread_noNA2 %>%
  dplyr::select(non_Pteropus, winter_non_Pteropus_pred) %>%
  gather()

# Model performance (prediction interval)
(winter_non_Pteropus_mod_perform <- ggplot(m_winter_non_Pteropus, aes(x = key, y = value)) +
    geom_point(position = position_jitter(width = 0.15, height = 0), shape = 21, size = 0.1) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), color = "black", fill = "black") +
    stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .2, y = 0),
                 size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
    scale_x_discrete(name = "Value", labels = c("observed", "predicted")) +
    scale_y_continuous(name = "Number of winter bat visits", breaks = seq(0, 600, 100), limits = c(0, 600)) +
    ggtitle("Non-Pteropus (winter)") +
    theme_cowplot(font_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)))

# Plot predictions
(mod_winter_non_Pteropus <- ggplot(non_Pteropus_all_mod_nb_top_pred_winter, aes(x = x, y = predicted, group = group)) +
  # Estimated means for single level of tree height
  geom_line(color = "#808080", size = 1) +
  # Estimated confidence intervals for single level of tree height
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#808080", alpha = 0.5) +
  # Set x axis
  scale_x_continuous(name = "Minimum nightly temperature (ºC)", breaks = seq(8, 30, 2), limits = c(8, 30)) +
  # Set y axis
  scale_y_continuous(name = "Mean winter non-Pteropus visits", breaks = seq(0, 150, 50), limits = c(0, 150)) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none"))

# Output plot to file
ggsave("./Results/non_Pteropus_visits_winter_NBmodel_pred_mintemp.png",
       height = 4, width = 4, dpi = 300, bg = "white")
ggsave("./Results/non_Pteropus_visits_winter_NBmodel_pred_mintemp.tiff",
       height = 4, width = 4, dpi = 300, bg = "white", compression = "lzw")

# Combine Pteropus and non-Pteropus subplots
plot_grid(mod_winter_Pteropus, mod_winter_non_Pteropus, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/all_bats_winter_NBmodels.png",
       height = 6, width = 6, dpi = 300, bg = "white")
ggsave("./Results/all_bats_winter_NBmodels.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

# Combine model performance plots
plot_grid(Pteropus_mod_perform, non_Pteropus_mod_perform,
          winter_Pteropus_mod_perform, winter_non_Pteropus_mod_perform, nrow = 2, ncol = 2, align = "hv")
# Output plot to file
ggsave("./Results/all_bats_winter_NBmodels_performance.png",
       height = 6, width = 6, dpi = 300, bg = "white",)
ggsave("./Results/all_bats_winter_NBmodels_performance.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

########################################################
### Looking at distributions of covariates by season ###
########################################################

# Plot of average tree height by season
tree_height <- ggplot(infra_day_spread, aes(x = as.factor(Month_Season), y = treeHeight,
                             color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .125, y = 0)) +
  # Line to show median
  geom_hline(yintercept = mean(infra_day_spread$treeHeight, na.rm = TRUE), color = "grey50") +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .125, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .125, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Tree height (cm)",
                     limits = c(0, 1500),
                     breaks = seq(0, 1500, 300)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Plot of average shaved area by season
shaved_area <- ggplot(infra_day_spread, aes(x = as.factor(Month_Season), y = shavearea,
                                            color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .125, y = 0)) +
  # Line to show median
  geom_hline(yintercept = mean(infra_day_spread$shavearea, na.rm = TRUE), color = "grey50") +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .125, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .125, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Shaved area (sq cm)",
                     limits = c(0, 1500),
                     breaks = seq(0, 1500, 300)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Plot of average shaved area by season
days_since_shaving <- ggplot(infra_day_spread, aes(x = as.factor(Month_Season), y = days_since_shaving,
                                            color = Month_Season, fill = Month_Season)) +
  # Flattened violin plot for distribution of points
  geom_flat_violin(scale = "width", position = position_nudge(x = .125, y = 0)) +
  # Line to show median
  geom_hline(yintercept = mean(infra_day_spread$days_since_shaving, na.rm = TRUE), color = "grey50") +
  # Individual data points
  geom_point(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 0.5) +
  # Plot calculated mean number of visits
  stat_summary(fun = "mean", na.rm = TRUE, geom = "point", position = position_nudge(x = .125, y = 0),
               size = 2, shape = 21, color = "white", fill = "black", alpha = 0.8) +
  # Plot calculated median number of visits
  stat_summary(fun = "median", na.rm = TRUE, geom = "point", position = position_nudge(x = .125, y = 0),
               size = 2, shape = 23, color = "white", fill = "black", alpha = 0.8) +
  # Rename seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Days since shaving",
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Set the fill palette and arrangement of seasons
  scale_fill_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                    labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                    values = c("grey80", "grey60", "grey40", "black")) +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine covariate plots
plot_grid(tree_height, shaved_area, days_since_shaving, nrow = 3, align = "hv")
# Output plot to file
ggsave("./Results/visit_covariates.png",
       height = 6, width = 6, dpi = 300, bg = "white",)
ggsave("./Results/visit_covariates.tiff",
       height = 6, width = 6, dpi = 300, bg = "white", compression = "lzw")

##########################################
### Observation bootstrapping analysis ###
##########################################

# Formatting data frame with indices necessary for filtering within resampling loop
infra <- infra %>%
  mutate(year = year(as.Date(dtobs, format = "%m/%d/%Y")),
         month = month(as.Date(dtobs, format = "%m/%d/%Y"))) %>%
  # Creates "index1" for year-month
  unite(data = ., col = "index1", c("year", "month"), sep = "-", remove = TRUE)
# Creates "index2" for dtobs-treeNo (day and tree number)
infra <- unite(data = infra, col = "index2", c("dtobs", "treeNo"), sep = "-", remove = FALSE)

## Run bootstrapping
# Check if the file already already exists and read it in
if(file.exists("./Data/infra_bootstrap_output.rds")){
  # Read in bootstrap data from file
  boots <- readRDS(file = "./Data/infra_bootstrap_output.rds")
  # Otherwise, run the bootsrapping code
} else {
  # Make empty file for storing results
  boots <- NULL
  # Loop through number of bootstrap iterations
  for(i in 1:100){
    # Record time at start of iteration
    start.time <- Sys.time()
    # Bootrap resampling on infrared camera data, sampling 28 observation days
    boots[[i]] <- days_sampler(df = infra, samples = 28)
    # Print the iteration
    print(i)
    # Print time at end of iteration
    print(Sys.time() - start.time)
  }
  # Save bootstrap data to file for later analysis
  saveRDS(boots, file = "./Data/infra_bootstrap_output.rds")
}

## Summarize bootstrap observations by month
boots_month_summary <- boots_summary(boots, c("Month_Season", "month"))
# Summarize bootstrap observations by month (averaging across all bootstrap resamples)
boots_month_grand_summary <- boots_grand_summary(boots_month_summary, c("Month_Season", "month"))
# Output month summary to file
write.csv(boots_month_grand_summary, "./Results/bootstrap_month_summary_statistics.csv")

## Summarize bootstrap observations by season
boots_season_summary <- boots_summary(boots, "Month_Season")
# Summarize bootstrap observations by season (averaging across all bootstrap resamples)
boots_season_grand_summary <- boots_grand_summary(boots_season_summary, "Month_Season")
# Output season summary to file
write.csv(boots_season_grand_summary, "./Results/bootstrap_season_summary_statistics.csv")

# Pteropus pairwise tests of means by season
Pteropus_season_means_test <-
  array(0, c(
    nrow(boots_season_grand_summary),
    nrow(boots_season_grand_summary)
  ))
for (i in 1:nrow(boots_season_grand_summary)) {
  for (j in 1:nrow(boots_season_grand_summary)) {
    Pteropus_season_means_test[i, j] <-
      as.numeric((
        boots_season_grand_summary$Upper_mean_Pteropus[i] <= boots_season_grand_summary$Upper_mean_Pteropus[j] &
          boots_season_grand_summary$Upper_mean_Pteropus[i] >= boots_season_grand_summary$Lower_mean_Pteropus[j]
      ) |
        (
          boots_season_grand_summary$Lower_mean_Pteropus[i] <= boots_season_grand_summary$Upper_mean_Pteropus[j] &
            boots_season_grand_summary$Lower_mean_Pteropus[i] >= boots_season_grand_summary$Lower_mean_Pteropus[j]
        )
      )
  }
}
colnames(Pteropus_season_means_test) <-
  boots_season_grand_summary$Month_Season
rownames(Pteropus_season_means_test) <-
  boots_season_grand_summary$Month_Season
# Output test results to file
write.csv(Pteropus_season_means_test, "./Results/bootstrap_season_means_test_Pteropus.csv")

# Non-Pteropus pairwise tests of means by season
nonPteropus_season_means_test <- array(0, c(nrow(boots_season_grand_summary), nrow(boots_season_grand_summary)))
for(i in 1:nrow(boots_season_grand_summary)){
  for(j in 1:nrow(boots_season_grand_summary)){
    nonPteropus_season_means_test[i, j] <- as.numeric((boots_season_grand_summary$Upper_mean_non_Pteropus[i] <= boots_season_grand_summary$Upper_mean_non_Pteropus[j] & boots_season_grand_summary$Upper_mean_non_Pteropus[i] >= boots_season_grand_summary$Lower_mean_non_Pteropus[j]) |
                                                        (boots_season_grand_summary$Lower_mean_non_Pteropus[i] <= boots_season_grand_summary$Upper_mean_non_Pteropus[j] & boots_season_grand_summary$Lower_mean_non_Pteropus[i] >= boots_season_grand_summary$Lower_mean_non_Pteropus[j]))
  }
}
colnames(nonPteropus_season_means_test) <- boots_season_grand_summary$Month_Season
rownames(nonPteropus_season_means_test) <- boots_season_grand_summary$Month_Season
# Output test results to file
write.csv(nonPteropus_season_means_test, "./Results/bootstrap_season_means_test_nonPteropus.csv")

## Plots for bootstrapping analysis: mean bat visits by month

# Subplot for Pteropus mean visits
boot_month_mean_visitsA <- ggplot(boots_month_grand_summary, aes(x = as.factor(month), y = Mean_mean_Pteropus, color = as.factor(Month_Season))) +
  # Points are the bootstrap estimates of mean visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = as.factor(month), ymin = Lower_mean_Pteropus, ymax = Upper_mean_Pteropus),
                 size = 1) +
  # Thin bars are the range (min, max)
  geom_linerange(aes(x = as.factor(month), ymin = Min_mean_Pteropus, ymax = Max_mean_Pteropus)) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Mean number of bat visits", 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Subplot for non-Pteropus mean visits
boot_month_mean_visitsB <- ggplot(boots_month_grand_summary, aes(x = as.factor(month), y = Mean_mean_non_Pteropus, color = as.factor(Month_Season))) +
  # Points are the bootstrap estimates of mean visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = as.factor(month), ymin = Lower_mean_non_Pteropus, ymax = Upper_mean_non_Pteropus),
                 size = 1) +
  # Thin bars are the range (min, max)
  geom_linerange(aes(x = as.factor(month), ymin = Min_mean_non_Pteropus, ymax = Max_mean_non_Pteropus)) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Mean number of bat visits",
                     limits = c(0, 250),
                     breaks = seq(0, 250, 50)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Combine mean visits subplots vertically (removing legends from both plots)
boot_month_mean_visits <- plot_grid(boot_month_mean_visitsA + theme(legend.position = "none"),
                                    boot_month_mean_visitsB + theme(legend.position = "none"),
                                    nrow = 2, align = "v")
# Extract legend from first plot
boot_month_mean_visits_legend <- g_legend(boot_month_mean_visitsA)
# Arrange combined subplots and legend
plot_grid(boot_month_mean_visits, boot_month_mean_visits_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/bootstrap_month_mean_bat_visits.png",
       height = 8, width = 8, dpi = 300, bg = "white")
ggsave("./Results/bootstrap_month_mean_bat_visits.tiff",
       height = 8, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for bootstrapping analysis: sum of bat visits by month

# Subplot for sum of Pteropus visits
boot_month_sum_visitsA <- ggplot(boots_month_grand_summary, aes(x = as.factor(month), y = Mean_sum_Pteropus, color = as.factor(Month_Season))) +
  # Points are the bootstrap estimates of mean visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = as.factor(month), ymin = Lower_sum_Pteropus, ymax = Upper_sum_Pteropus),
                 size = 1) +
  # Thin bars are the range (min, max)
  geom_linerange(aes(x = as.factor(month), ymin = Min_sum_Pteropus, ymax = Max_sum_Pteropus)) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Total number of bat visits",
                     limits = c(0, 500),
                     breaks = seq(0, 500, 100)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Subplot for sum of non-Pteropus visits
boot_month_sum_visitsB <- ggplot(boots_month_grand_summary, aes(x = as.factor(month), y = Mean_sum_non_Pteropus, color = as.factor(Month_Season))) +
  # Points are the bootstrap estimates of mean visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = as.factor(month), ymin = Lower_sum_non_Pteropus, ymax = Upper_sum_non_Pteropus),
                 size = 1) +
  # Thin bars are the range (min, max)
  geom_linerange(aes(x = as.factor(month), ymin = Min_sum_non_Pteropus, ymax = Max_sum_non_Pteropus)) +
  # Rename numbered months as abbreviated names
  scale_x_discrete(name = "Month",
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # Set y axis
  scale_y_continuous(name = "Total number of bat visits",
                     limits = c(0, 20000),
                     breaks = seq(0, 20000, 4000)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and center horizontal legend at bottom
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal")

# Combine sum of visit subplots vertically (removing legends from both plots)
boot_month_sum_visits <- plot_grid(boot_month_sum_visitsA + theme(legend.position = "none"),
                                   boot_month_sum_visitsB + theme(legend.position = "none"),
                                   nrow = 2, align = "v")
# Extract legend from first plot
boot_month_sum_visits_legend <- g_legend(boot_month_sum_visitsA)
# Arrange combined subplots and legend
plot_grid(boot_month_sum_visits, boot_month_sum_visits_legend, nrow = 2, rel_heights = c(0.95, 0.05))
# Output plot to file
ggsave("./Results/bootstrap_month_sum_bat_visits.png",
       height = 8, width = 8, dpi = 300, bg = "white")
ggsave("./Results/bootstrap_month_sum_bat_visits.tiff",
       height = 8, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for bootstrapping analysis: mean bat visits by season

# Subplot for Pteropus mean visits
boot_season_mean_visitsA <- ggplot(boots_season_grand_summary, aes(x = as.factor(Month_Season), y = Mean_mean_Pteropus, color = as.factor(Month_Season))) +
  # Points are the bootstrap estimates of the sum of visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = as.factor(Month_Season), ymin = Lower_mean_Pteropus, ymax = Upper_mean_Pteropus),
                 size = 1) +
  # Thin bars are the range (min, max)
  geom_linerange(aes(x = as.factor(Month_Season), ymin = Min_mean_Pteropus, ymax = Max_mean_Pteropus)) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean number of bat visits",
                     limits = c(0, 3),
                     breaks = seq(0, 3, 0.5)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for non-Pteropus mean visits
boot_season_mean_visitsB <- ggplot(boots_season_grand_summary, aes(x = as.factor(Month_Season), y = Mean_mean_non_Pteropus, color = as.factor(Month_Season))) +
  # Points are the bootstrap estimates of the sum of visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = as.factor(Month_Season), ymin = Lower_mean_non_Pteropus, ymax = Upper_mean_non_Pteropus),
                 size = 1) +
  # Thin bars are the range (min, max)
  geom_linerange(aes(x = as.factor(Month_Season), ymin = Min_mean_non_Pteropus, ymax = Max_mean_non_Pteropus)) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Mean number of bat visits",
                     limits = c(0, 150),
                     breaks = seq(0, 150, 25)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine mean visits subplots vertically
plot_grid(boot_season_mean_visitsA, boot_season_mean_visitsB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/bootstrap_season_mean_bat_visits.png",
       height = 8, width = 8, dpi = 300, bg = "white")
ggsave("./Results/bootstrap_season_mean_bat_visits.tiff",
       height = 8, width = 8, dpi = 300, bg = "white", compression = "lzw")

## Plots for bootstrapping analysis: sum of bat visits by season

# Subplot for sum of Pteropus visits
boot_season_sum_visitsA <- ggplot(boots_season_grand_summary, aes(x = as.factor(Month_Season), y = Mean_sum_Pteropus, color = as.factor(Month_Season))) +
  # Points are the bootstrap estimates of the sum of visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = as.factor(Month_Season), ymin = Lower_sum_Pteropus, ymax = Upper_sum_Pteropus),
                 size = 1) +
  # Thin bars are the range (min, max)
  geom_linerange(aes(x = as.factor(Month_Season), ymin = Min_sum_Pteropus, ymax = Max_sum_Pteropus)) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Total number of bat visits", 
                     limits = c(0, 600),
                     breaks = seq(0, 600, 100)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Subplot for sum of non-Pteropus visits
boot_season_sum_visitsB <- ggplot(boots_season_grand_summary, aes(x = as.factor(Month_Season), y = Mean_sum_non_Pteropus, color = as.factor(Month_Season))) +
  # Points are the bootstrap estimates of the sum of visits
  geom_point(size = 2) +
  # Wide bars are the 95% confidence intervals
  geom_linerange(aes(x = as.factor(Month_Season), ymin = Lower_sum_non_Pteropus, ymax = Upper_sum_non_Pteropus),
                 size = 1) +
  # Thin bars are the range (min, max)
  geom_linerange(aes(x = as.factor(Month_Season), ymin = Min_sum_non_Pteropus, ymax = Max_sum_non_Pteropus)) +
  # Rename and order seasons
  scale_x_discrete(name = "Season", limits = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                   labels = c("winter", "spring", "monsoon", "post-monsoon")) +
  # Set y axis
  scale_y_continuous(name = "Total number of bat visits",
                     limits = c(0, 40000),
                     breaks = seq(0, 40000, 10000)) +
  # Set the color palette and arrangement of seasons
  scale_color_manual(name = "Season", breaks = c("Winter", "Spring", "Monsoon", "Postmonsoon"),
                     labels = c("winter (Dec-Feb)", "spring (Mar-May)", "monsoon (Jun-Sep)", "post-monsoon (Oct-Nov)"),
                     values = c("grey80", "grey60", "grey40", "black")) +
  # Subplot title
  ggtitle("Non-Pteropus") +
  # Plot theme (simplified axes in all black and white)
  theme_cowplot(font_size = 12) +
  # Justify subplot title to the center and remove legend
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Combine sum of visits subplots vertically
plot_grid(boot_season_sum_visitsA, boot_season_sum_visitsB, nrow = 2, align = "v")
# Output plot to file
ggsave("./Results/bootstrap_season_sum_bat_visits.png",
       height = 8, width = 8, dpi = 300, bg = "white")
ggsave("./Results/bootstrap_season_sum_bat_visits.tiff",
       height = 8, width = 8, dpi = 300, bg = "white", compression = "lzw")

###################
### End of code ###
###################
