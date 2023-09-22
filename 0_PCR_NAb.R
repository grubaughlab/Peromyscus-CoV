
# Code Description --------------------------------------------------------

# 1) Create Fig 1B (sVNT results)
# 2) Create Figure S4 (human SC2 transmission context for swab sampling)
# Note: raw data available upon request

# Load libraries + data + set paths ---------------------------------------------------

library(readxl)
library(dplyr)
library(plyr)
library(flextable)
library(lubridate)
library(ggplot2)
library(readr)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(zoo)

file_loc <- getwd()
sheet_date <- "053023" # date of downloaded metadata sheet

# dat_2020 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                        sheet = "Mice_Guilford_2020")
# dat_2021 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                        sheet = "Mice_Guilford_2021")
# dat_2022_gu <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                           sheet = "Mice_Guilford_2022") # differentiate 2022 naming b/c 2 locations
# dat_2022_nb <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                           sheet = "Mice_North_Branford_2022")
# dat_2022_deer <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                             sheet = "Deer_2021_2022")

path_figures <- paste0(file_loc, "/Figures/")

# Format Data -------------------------------------------------------------

# Standardize tag number label 
dat_2022_deer <- dat_2022_deer %>% dplyr::rename(`Tag #` = `Ear Tag #`)

# Add empty column for age to mice so can standardize column selection later (did not collect age data for mice)
dat_2020$Age <- rep(NA, nrow(dat_2020))
dat_2021$Age <- rep(NA, nrow(dat_2021))
dat_2022_gu$Age <- rep(NA, nrow(dat_2022_gu))
dat_2022_nb$Age <- rep(NA, nrow(dat_2022_nb))

# Add dataset name
dat_2020$Dataset <- "Mice_2020"
dat_2021$Dataset <- "Mice_2021"
dat_2022_gu$Dataset <- "Mice_2022_Guilford"
dat_2022_nb$Dataset <- "Mice_2022_NB"

dat_2022_deer[which(dat_2022_deer$Location == "Norwalk" & (year(dat_2022_deer$Date) == "2021")), "Dataset"] <- "Deer_Norwalk_2021"
dat_2022_deer[which(dat_2022_deer$Location == "Bridgeport" & (dat_2022_deer$Date <= "2022-02-01")), "Dataset"] <- "Deer_Bridgeport_2021" # includes some January 2022
dat_2022_deer[which(dat_2022_deer$Location == "Norwalk" & (year(dat_2022_deer$Date) == "2022")), "Dataset"] <- "Deer_Norwalk_2022"
dat_2022_deer[which(dat_2022_deer$Location == "Bridgeport" & (dat_2022_deer$Date > "2022-02-01")), "Dataset"] <- "Deer_Bridgeport_2022"

# Break out deer groups
dat_2021_nor <- dat_2022_deer %>% dplyr::filter(Dataset == "Deer_Norwalk_2021")
dat_2021_bri <- dat_2022_deer %>% dplyr::filter(Dataset == "Deer_Bridgeport_2021")
dat_2022_nor <- dat_2022_deer %>% dplyr::filter(Dataset == "Deer_Norwalk_2022")
dat_2022_bri <- dat_2022_deer %>% dplyr::filter(Dataset == "Deer_Bridgeport_2022")

# Select needed columns 
keep_selected_col <- function(dat) {
  dat_new <- dat %>% dplyr::select(Date,
                                   `Tag #`,
                                   Sex,
                                   Age,
                                   `Recapture?`,
                                   ELISA_Plate_No, # so can separate ELISA tested vs totally untested samples in case of recaptures
                                   cPass_Result_Run1, # sVNT initial run
                                   cPass_Result_Run2, # sVNT confirmatory duplicate
                                   Dataset)
  dat_new
}

## mice
dat_2020_sub <- keep_selected_col(dat_2020)
dat_2021_sub <- keep_selected_col(dat_2021)
dat_2022_gu_sub <- keep_selected_col(dat_2022_gu)
dat_2022_nb_sub <- keep_selected_col(dat_2022_nb)

## deer
dat_2021_nor_sub <- keep_selected_col(dat_2021_nor)
dat_2021_bri_sub <- keep_selected_col(dat_2021_bri)
dat_2022_nor_sub <- keep_selected_col(dat_2022_nor)
dat_2022_bri_sub <- keep_selected_col(dat_2022_bri)

# Combine into single dataframe
dat_all <- rbind.data.frame(dat_2020_sub, # mice
                            dat_2021_sub,
                            dat_2022_gu_sub,
                            dat_2022_nb_sub,
                            dat_2021_nor_sub, # deer
                            dat_2021_bri_sub,
                            dat_2022_nor_sub,
                            dat_2022_bri_sub)

# Add descriptive group by dataset & make a factor so can order in plots
dat_all[which(dat_all$Dataset %in% c("Mice_2020", "Mice_2021", "Mice_2022_Guilford")), "Group"] <- "Mice_Residential"
dat_all[which(dat_all$Dataset %in% c("Mice_2022_NB")), "Group"] <- "Mice_Forested"
dat_all[which(dat_all$Dataset %in% c("Deer_Norwalk_2021", "Deer_Norwalk_2022")), "Group"] <- "Deer_Norwalk"
dat_all[which(dat_all$Dataset %in% c("Deer_Bridgeport_2021", "Deer_Bridgeport_2022")), "Group"] <- "Deer_Bridgeport"

dat_all$Group <- factor(dat_all$Group, levels = c("Mice_Residential", 
                                                  "Mice_Forested",
                                                  "Deer_Norwalk",
                                                  "Deer_Bridgeport"))

# Create new column to delineate confirmed positives
dat_all$Year <- as.character(year(dat_all$Date))
dat_all$NAb_result <- rep(NA, nrow(dat_all))

dat_all[which((dat_all$cPass_Result_Run1 == "Pos") & (dat_all$cPass_Result_Run2 == "Pos")), "NAb_result"] <- "Pos duplicate" # positive on 2 subsequent runs
dat_all[which((dat_all$cPass_Result_Run1 == "Pos") & (dat_all$cPass_Result_Run2 == "Elevated")), "NAb_result"] <- "Pos elevated" # positive on first run; elevated on second but does not cross positivity threshold
dat_all[which((dat_all$cPass_Result_Run1 == "Pos") & (dat_all$cPass_Result_Run2 == "Insufficient sample")), "NAb_result"] <- "Pos singlicate" # not enough sample left to confirm
dat_all[which((dat_all$cPass_Result_Run1 == "Pos") & (dat_all$cPass_Result_Run2 == "Neg")), "NAb_result"] <- "Neg" # false positive

dat_all[which((is.na(dat_all$ELISA_Plate_No) == FALSE) & (is.na(dat_all$cPass_Result_Run1) == TRUE) & (dat_all$Group %in% c("Mice_Residential", "Mice_Forested"))), "NAb_result"] <- "Did not pass pre-screening" 
dat_all[which((is.na(dat_all$ELISA_Plate_No) == TRUE) & (dat_all$Year != "2022") & (dat_all$Group %in% c("Mice_Residential", "Mice_Forested"))), "NAb_result"] <- "Untested" # in case of recaptures (only applicable to 2020/21)
dat_all[which((is.na(dat_all$NAb_result) == TRUE) & (dat_all$cPass_Result_Run1 == "Neg")), "NAb_result"] <- "Neg"

# Drop untested samples
dat_all <- dat_all %>% dplyr::filter(NAb_result != "Untested")

# Run Fishers Exact Test for Independence --------------------------------
# Note: Chi-squared doesn't work well w/ <5 obs in a category
# No. positive INDIVIDUALS vs. no. negative INDIVIDUALS full time period 

fisher_by_ind <- function(dat, group_1, group_2) {
  
  dat$Year <- year(dat$Date)
  
  print(unique(dat$Group))
  print(unique(dat$NAb_result))
  
  dat_fish <- dat %>% dplyr::filter(Group %in% c(group_1, group_2), 
                                    NAb_result != "Untested",
                                    NAb_result != "Did not pass pre-screening")
  
  ## group 1 - pos and neg counts
  dat_pos_g1 <- dat_fish %>% dplyr::filter(Group == group_1, 
                                           NAb_result != "Neg")
  print(ddply(dat_pos_g1, .(Year), summarize, Unique_No = length(unique(`Tag #`))))
  no_pos_g1 <- length(unique(dat_pos_g1$`Tag #`))
  
  dat_neg_g1 <- dat_fish %>% dplyr::filter(Group == group_1, 
                                           NAb_result == "Neg") 
  print(ddply(dat_neg_g1, .(Year), summarize, Unique_No = length(unique(`Tag #`))))
  no_neg_g1 <-length(unique(dat_neg_g1$`Tag #`))
  
  ## group 2 - pos and neg counts
  dat_pos_g2 <- dat_fish %>% dplyr::filter(Group == group_2, 
                                           NAb_result != "Neg")
  print(ddply(dat_pos_g2, .(Year), summarize, Unique_No = length(unique(`Tag #`))))
  no_pos_g2 <- length(unique(dat_pos_g2$`Tag #`))
  
  dat_neg_g2 <- dat_fish %>% dplyr::filter(Group == group_2, 
                                           NAb_result == "Neg") 
  print(ddply(dat_neg_g2, .(Year), summarize, Unique_No = length(unique(`Tag #`))))
  no_neg_g2 <-length(unique(dat_neg_g2$`Tag #`))
  
  ## test
  dat_fish_test <- cbind.data.frame(group_1 = c(no_neg_g1, no_pos_g1), 
                                    group_2 = c(no_neg_g2, no_pos_g2))
  rownames(dat_fish_test) <- c("Neg", "Pos")
  colnames(dat_fish_test) <- c(group_1, group_2)
  print(dat_fish_test)
  dat_result <- fisher.test(dat_fish_test) 
  dat_result
}

dat_fish_result_mice_IND <- fisher_by_ind(dat = dat_all, group_1 = "Mice_Residential", group_2 = "Mice_Forested")
dat_fish_result_deer_IND <- fisher_by_ind(dat = dat_all, group_1 = "Deer_Norwalk", group_2 = "Deer_Bridgeport")

# Count Individual Mice Tested by Dataset (i.e. does not include recaptures) ---------------------------------

count_by_ind <- function(dat, dataset_name) {
  
  dat$Year <- year(dat$Date)
  dat_count <- dat %>% dplyr::filter(Dataset %in% dataset_name, 
                                     NAb_result != "Untested",
                                     NAb_result != "Did not pass pre-screening")
  
  ## dataset - pos and neg counts
  dat_pos_g1 <- dat_count %>% dplyr::filter(NAb_result != "Neg")
  no_pos_g1 <- length(unique(dat_pos_g1$`Tag #`))
  
  dat_neg_g1 <- dat_count %>% dplyr::filter(NAb_result == "Neg") 
  no_neg_g1 <-length(unique(dat_neg_g1$`Tag #`))
  
  # number of negative samples that were later positive (so don't double count)
  no_pos_orig_neg <- length(which(unique(dat_pos_g1$`Tag #`) %in% unique(dat_neg_g1$`Tag #`)) == TRUE)
  no_neg_g1 <- no_neg_g1 - no_pos_orig_neg
  
  ## test
  dat_count_total <- cbind.data.frame(no_neg_g1, no_pos_g1)
  colnames(dat_count_total) <- c("Neg", "Pos")
  dat_count_total$Total <- rowSums(dat_count_total)
  print(dat_count_total)
}

unique(dat_all$Dataset)
unique(dat_all$NAb_result)

# residential mice

## across all years to prevent double-counting recaptures
count_by_ind_g1 <- count_by_ind(dat = dat_all, dataset_name = c("Mice_2020", "Mice_2021", "Mice_2022_Guilford"))
all_years <- "All"
count_by_ind_mice_res <- cbind.data.frame(all_years, count_by_ind_g1)
per_pos <- round(((count_by_ind_mice_res[which(count_by_ind_mice_res$all_years == "All"), "Pos"] / 
                     count_by_ind_mice_res[which(count_by_ind_mice_res$all_years == "All"), "Total"]) * 100), digits = 0)

## by year - if totals per year add up to overall totals above, no recaptures across years
count_by_ind_g1 <- count_by_ind(dat = dat_all, dataset_name = "Mice_2020")
count_by_ind_g2 <- count_by_ind(dat = dat_all, dataset_name = "Mice_2021")
count_by_ind_g3 <- count_by_ind(dat = dat_all, dataset_name = "Mice_2022_Guilford")
count_by_ind_combo <- rbind.data.frame(count_by_ind_g1, count_by_ind_g2, count_by_ind_g3)
totals <- colSums(count_by_ind_combo)
count_by_ind_combo <- rbind.data.frame(count_by_ind_combo, totals)
all_years <- c("2020", "2021", '2022', "All")
count_by_ind_mice_res_year <- cbind.data.frame(all_years, count_by_ind_combo)
per_pos <- round(((count_by_ind_mice_res_year[which(count_by_ind_mice_res_year$all_years == "All"), "Pos"] / 
                     count_by_ind_mice_res_year[which(count_by_ind_mice_res_year$all_years == "All"), "Total"]) * 100), digits = 0)

# forested mice
count_by_ind_g1 <- count_by_ind(dat = dat_all, dataset_name = "Mice_2022_NB")
all_years <- "2022"
count_by_ind_mice_for <- cbind.data.frame(all_years, count_by_ind_g1)

# deer

## across all years to prevent double-counting recaptures
count_by_ind_g1 <- count_by_ind(dat = dat_all, dataset_name = c("Deer_Norwalk_2021", "Deer_Norwalk_2022"))
count_by_ind_g2 <- count_by_ind(dat = dat_all, dataset_name = c("Deer_Bridgeport_2021", "Deer_Bridgeport_2022"))
count_by_ind_combo <- rbind.data.frame(count_by_ind_g1, count_by_ind_g2)

totals <- colSums(count_by_ind_combo)
count_by_ind_combo <- rbind.data.frame(count_by_ind_combo, totals)
all_years <- c(rep("2021_2022", 2), "All")
all_locations <- c("Norwalk", "Bridgeport", "All")
count_by_ind_deer <- cbind.data.frame(all_years, all_locations, count_by_ind_combo)
per_pos <- round(((count_by_ind_deer[which(count_by_ind_deer$all_years == "All"), "Pos"] / 
                     count_by_ind_deer[which(count_by_ind_deer$all_years == "All"), "Total"]) * 100), digits = 0)

## by year - if totals per year add up to overall totals above, no recaptures across years
count_by_ind_g1 <- count_by_ind(dat = dat_all, dataset_name = "Deer_Norwalk_2021")
count_by_ind_g2 <- count_by_ind(dat = dat_all, dataset_name = "Deer_Norwalk_2022")
count_by_ind_g3 <- count_by_ind(dat = dat_all, dataset_name = "Deer_Bridgeport_2021")
count_by_ind_g4 <- count_by_ind(dat = dat_all, dataset_name = "Deer_Bridgeport_2022")

count_by_ind_combo <- rbind.data.frame(count_by_ind_g1, count_by_ind_g2, count_by_ind_g3, count_by_ind_g4)
totals <- colSums(count_by_ind_combo)
count_by_ind_combo <- rbind.data.frame(count_by_ind_combo, totals)
all_years <- c(rep(c("2021", "2022"), 2), "All")
all_locations <- c(rep(c("Norwalk", "Bridgeport"), each = 2), "All")

count_by_ind_deer_year <- cbind.data.frame(all_years, all_locations, count_by_ind_combo)

# Count Individual Mice Tested by Dataset & Describe Demographics (does not include recaptures) ---------------------------------

dem_ind <- function(dat, dataset_name, species) {
  
  dat$Year <- year(dat$Date)
  
  # print(unique(dat$NAb_result))
  
  # only include tested individuals
  dat <- dat %>% dplyr::filter(Dataset %in% dataset_name,
                               NAb_result != "Untested",
                               NAb_result != "Did not pass pre-screening")
  
  # sex stats
  dat_sex <- dat %>% dplyr::select(`Tag #`, Sex)
  dat_sex_check_dupes <- dat_sex[which(duplicated(dat_sex) == FALSE), ] # checks duplicated rows (ID + sex)
  
  num_ids <- length(dat_sex_check_dupes$`Tag #`)
  num_ids_unique <- length(unique(dat_sex_check_dupes$`Tag #`))
  
  if (num_ids != num_ids_unique) {
    
    dat_check_dupes <- dat_sex_check_dupes
    
    dupes <- dat_check_dupes[which(duplicated(dat_check_dupes$`Tag #`) == TRUE), ] # checks just duplicated IDs
    
    i <- 1
    for (each_id in dupes$`Tag #`) {
      
      single_id <- dupes[which(dupes$`Tag #` == each_id), ]
      num_dupes_each <- nrow(single_id)
      
      if (i == 1) {
        
        df_dupes <- data.frame("ID" = each_id, "Num_Dupes" = num_dupes_each)
        
      } else {
        df_dupes_single <- cbind.data.frame("ID" = each_id, "Num_Dupes" = num_dupes_each)
        df_dupes <- rbind.data.frame(df_dupes, df_dupes_single)
      }
      
      i <- i + 1
      
    }
    print(df_dupes)
    print(sum(df_dupes$Num_Dupes))
    
    dat_check_dupes[which(dat_check_dupes$`Tag #` %in% df_dupes$ID), "Sex"] <- "Mismatched"
    
    ## remove dupes for mismatched category now that ID + mismatched are same
    no_dupes <- dat_check_dupes[which(duplicated(dat_check_dupes$`Tag #`) == FALSE), ]
    
    dat_ct_sex <- table(no_dupes$Sex)
    dat_ct_sex_sum <- sum(table(no_dupes$Sex))
    dat_pr_sex <- round(prop.table(table(no_dupes$Sex)), digits = 2)
    
  } else {
    
    dat_ct_sex <- table(dat_sex_check_dupes$Sex)
    dat_ct_sex_sum <- sum(table(dat_sex_check_dupes$Sex))
    dat_pr_sex <- round(prop.table(table(dat_sex_check_dupes$Sex)), digits = 2)
  }
  
  print(dat_ct_sex)
  print(dat_ct_sex_sum)
  print(dat_pr_sex)
  
  if (species == "deer") {
    
    ## age stats
    dat_age <- dat %>% dplyr::select(`Tag #`, Age)
    dat_age_check_dupes <- dat_age[which(duplicated(dat_age) == FALSE), ] # checks duplicated rows (ID + age)
    
    num_ids <- length(dat_age_check_dupes$`Tag #`)
    num_ids_unique <- length(unique(dat_age_check_dupes$`Tag #`))
    
    if (num_ids != num_ids_unique) {
      
      dat_check_dupes <- dat_age_check_dupes
      
      dupes <- dat_check_dupes[which(duplicated(dat_check_dupes$`Tag #`) == TRUE), ] # checks just duplicated IDs
      
      i <- 1
      for (each_id in dupes$`Tag #`) {
        
        single_id <- dupes[which(dupes$`Tag #` == each_id), ]
        num_dupes_each <- nrow(single_id)
        
        if (i == 1) {
          
          df_dupes <- data.frame("ID" = each_id, "Num_Dupes" = num_dupes_each)
          
        } else {
          df_dupes_single <- cbind.data.frame("ID" = each_id, "Num_Dupes" = num_dupes_each)
          df_dupes <- rbind.data.frame(df_dupes, df_dupes_single)
        }
        
        i <- i + 1
        
      }
      print(df_dupes)
      print(sum(df_dupes$Num_Dupes))
      
      dat_check_dupes[which(dat_check_dupes$`Tag #` %in% df_dupes$ID), "Age"] <- "Age_Progression" # e.g., yearling to adult
      
      ## remove dupes for mismatched category now that ID + age_progression are same
      no_dupes <- dat_check_dupes[which(duplicated(dat_check_dupes$`Tag #`) == FALSE), ]
      
      dat_ct_age <- table(no_dupes$Age)
      dat_ct_age_sum <- sum(table(no_dupes$Age))
      dat_pr_age <- round(prop.table(table(no_dupes$Age)), digits = 2)
      
    } else {
      
      dat_ct_age <- table(dat_age_check_dupes$Age)
      dat_ct_age_sum <- sum(table(dat_age_check_dupes$Age))
      dat_pr_age <- round(prop.table(table(dat_age_check_dupes$Age)), digits = 2)
    }
    
    print(dat_ct_age)
    print(dat_ct_age_sum)
    print(dat_pr_age)
  }
  
}

## across all years to prevent double-counting recaptures
dem_ind(dat = dat_all, dataset_name = c("Mice_2020", "Mice_2021", "Mice_2022_Guilford"), species = "mice") # residential mice
dem_ind(dat = dat_all, dataset_name = "Mice_2022_NB", species = "mice") # forested mice
dem_ind(dat = dat_all, dataset_name = c("Deer_Norwalk_2021", "Deer_Norwalk_2022"), species = "deer") # Norwalk deer - note 1 yearling became and adult
dem_ind(dat = dat_all, dataset_name = c("Deer_Bridgeport_2021", "Deer_Bridgeport_2022"), species = "deer") # Bridgeport deer

# Count Mice & Deer Sera (or Pooled Swab) Samples Tested by Dataset ---------------------------------

## 2021 samples
nrow(dat_all[which(dat_all$Dataset == "Deer_Norwalk_2021"), ])
nrow(dat_all[which(dat_all$Dataset == "Deer_Bridgeport_2021"), ])

## 2022 samples
nrow(dat_all[which(dat_all$Dataset == "Deer_Norwalk_2022"), ])
nrow(dat_all[which(dat_all$Dataset == "Deer_Bridgeport_2022"), ])

## All deer swabs (2022 only)
nrow(dat_all[which(dat_all$Dataset == "Deer_Norwalk_2022"), ]) +
  nrow(dat_all[which(dat_all$Dataset == "Deer_Bridgeport_2022"), ])

# Plot NAb results by group -----------------------------------------------

# Format fields
dat_all$Date <- as.Date(dat_all$Date)
dat_all$Month_Yr <- as.yearmon(dat_all$Date)
dat_all$Month_Yr <- as.Date(dat_all$Month_Yr, format = "%Y-%m")

# Make factor to label in order of result & combine positives & drop 2022 untested
unique(dat_all$NAb_result)
dat_all[which(dat_all$NAb_result %in% c("Pos duplicate", "Pos elevated", "Pos singlicate")), "NAb_result"] <- "Pos" # combine into all positive
unique(dat_all$NAb_result)

dat_all$NAb_result <- factor(dat_all$NAb_result,
                             levels = c("Pos", "Neg", "Did not pass pre-screening"))

# Calculate sample size for each group
ss_mice_res <- nrow(dat_all %>% dplyr::filter(Group == "Mice_Residential"))
ss_mice_for <- nrow(dat_all %>% dplyr::filter(Group == "Mice_Forested"))
ss_deer_nor <- nrow(dat_all %>% dplyr::filter(Group == "Deer_Norwalk"))
ss_deer_bri <- nrow(dat_all %>% dplyr::filter(Group == "Deer_Bridgeport"))

# Create labels for each plot w/ sample size number
group_labels <- c("Mice_Residential" = paste0("Mice 1-Guilford (Residential) (N=", ss_mice_res, ")"), 
                  "Mice_Forested" = paste0("Mice 2-North Branford (Forested) (N=", ss_mice_for, ")"), 
                  "Deer_Norwalk" = paste0("Deer 1-Norwalk (Residential/Forested) (N=", ss_deer_nor, ")"), 
                  "Deer_Bridgeport" = paste0("Deer 2-Bridgeport (Residential/Forested) (N=", ss_deer_bri, ")"))

# Split dataset by group
mice1 <- dat_all %>% dplyr::filter(Group == "Mice_Residential")
mice2 <- dat_all %>% dplyr::filter(Group == "Mice_Forested")
deer1 <- dat_all %>% dplyr::filter(Group == "Deer_Norwalk")
deer2 <- dat_all %>% dplyr::filter(Group == "Deer_Bridgeport")

# Select colors
viridis_names <-c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
n <- 6
par(mfrow=c(4,2), mar = c(1,1,1,1))
f <- sapply(1:8, function(x) image(matrix(1:n, n, 1), col = viridis(n=n, option = LETTERS[x]), axes =FALSE, main = viridis_names[x]))

NAb_result_colors <- c(viridis::rocket(6)[4],
                       "dark grey",
                       "light grey")

# Create vectors of labels / categories to assign colors
NAb_result_cats <- c("Pos", "Neg", "Did not pass pre-screening")
NAb_result_labels <- c("Positive", "Negative", "Did not pass pre-screening")

# Specify axis start & end dates to give 1 month buffer on each side of date range
start_date <- floor_date(min(dat_all$Month_Yr), unit = "month") %m-% months(1)
end_date <- floor_date(max(dat_all$Month_Yr), unit = "month") %m+% months(1)

# Function to plot NAb results for each group
plot_nab_results <- function(dat, 
                             species, 
                             group_labels, 
                             NAb_result_cats, 
                             NAb_result_labels,
                             NAb_result_colors) {
  
  unique_results <- unique(dat$NAb_result) # which results included for this group
  keep_cats <- which(NAb_result_cats %in% unique_results)
  selected_colors <- NAb_result_colors[keep_cats] # pull colors responding to included sVNT result categories 
  selected_labels <- NAb_result_labels[keep_cats] ## "" for labels
  
  p <- ggplot(dat) + 
    geom_bar(stat = "count", width = 25, aes(x=Month_Yr, fill = NAb_result)) +
    facet_wrap(~Group, nrow = 3, scales = "free_y", 
               labeller = labeller(Group = group_labels)) +
    labs(fill = "sVNT Result") +
    xlab("Month-Year") +
    ylab("No. Samples") +
    scale_fill_manual(labels = selected_labels,
                      values = selected_colors) +
    theme_bw() +
    scale_x_date(date_labels = "%b-%Y", date_breaks  ="3 month", limits = c(as.Date(start_date), as.Date(end_date))) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.text = element_text(size = 11),
          strip.text.x = element_text(size = 14),
          legend.title = element_text(size=12), 
          legend.text = element_text(size=11),
          axis.title = element_text(size = 12), 
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow=1, byrow=TRUE))
  p
}

# Set max y axis bounds for mice and deer results
ylim_mice <- 315
ylim_deer <- 15

# Plot each group and annotate with positive counts for each NAb result category

## mice1
count_pos <- mice1 %>% dplyr::filter(NAb_result == "Pos")
ddply(count_pos, .(Month_Yr, NAb_result), nrow)

p_mice1 <- plot_nab_results(dat = mice1, 
                            species = "mice",
                            group_labels = group_labels, 
                            NAb_result_cats = NAb_result_cats, 
                            NAb_result_labels = NAb_result_labels,
                            NAb_result_colors = NAb_result_colors) + 
  ylim(0, ylim_mice) +
  annotate("text", x = as.Date("2021-07-01"), y = 105, label = "N=1", size = 4, color = NAb_result_colors[1]) +
  annotate("text", x = as.Date("2022-06-15"), y = 160, label = "N=4", size = 4, color = NAb_result_colors[1]) +
  annotate("text", x = as.Date("2022-08-01"), y = 260, label = "N=1", size = 4, color = NAb_result_colors[1])

## mice 2
count_pos <- mice2 %>% dplyr::filter(NAb_result == "Pos")
ddply(count_pos, .(Month_Yr, NAb_result), nrow) # no positives

p_mice2 <- plot_nab_results(mice2, 
                            species = "mice",
                            group_labels = group_labels, 
                            NAb_result_cats = NAb_result_cats, 
                            NAb_result_labels = NAb_result_labels,
                            NAb_result_colors = NAb_result_colors) + 
  ylim(0, ylim_mice)

## deer 1
count_pos <- deer1 %>% dplyr::filter(NAb_result == "Pos")
ddply(count_pos, .(Month_Yr, NAb_result), nrow)

p_deer1 <- plot_nab_results(deer1, 
                            species = "deer",
                            group_labels = group_labels, 
                            NAb_result_cats = NAb_result_cats, 
                            NAb_result_labels = NAb_result_labels,
                            NAb_result_colors = NAb_result_colors) + 
  ylim(0, ylim_deer) +
  annotate("text", x = as.Date("2021-05-01"), y = 10.5, label = "N=2", size = 4, color = NAb_result_colors[1]) +
  annotate("text", x = as.Date("2021-08-01"), y = 10.5, label = "N=1*", size = 4, color = NAb_result_colors[1]) +
  annotate("text", x = as.Date("2022-07-01"), y = 9.5, label = "N=2*", size = 4, color = NAb_result_colors[1])

## deer 2
count_pos <- deer2 %>% dplyr::filter(NAb_result == "Pos")
ddply(count_pos, .(Month_Yr, NAb_result), nrow)

p_deer2 <- plot_nab_results(deer2, 
                            species = "deer",
                            group_labels = group_labels, 
                            NAb_result_cats = NAb_result_cats, 
                            NAb_result_labels = NAb_result_labels,
                            NAb_result_colors = NAb_result_colors) + 
  ylim(0, ylim_deer) +
  annotate("text", x = as.Date("2022-01-20"), y = 4.5, label = "N=1", size = 4, color = NAb_result_colors[1])

# Function for grabbing legend for use in a combination plot
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(p_mice1) # includes all sVNT result categories

# Combine all group plots into one single plot
p <- grid.arrange(ggarrange(
  p_mice1 + theme(legend.position="none"),
  p_mice2 + theme(legend.position="none"),
  p_deer1 + theme(legend.position="none"),
  p_deer2 + theme(legend.position="none")),
  mylegend,
  heights = c(9, 1))

ggsave(paste(path_figures, "Fig_NAb_pos.jpg", sep=""), p, height = 10, width =10) 

# Plot sampling context: estimated infections and swab collection ---------

# Count number of swabs
dat_all_swabs_mice <- dat_all %>% dplyr::filter(Dataset %in% c("Mice_2022_Guilford", "Mice_2022_NB"))
dat_all_swabs_mice$No_swabs <- rep(2, nrow(dat_all_swabs_mice)) # to count individual swabs instead of pooled; 2 swabs per mouse per collection
dat_all_swabs_deer <- dat_all %>% dplyr::filter(Dataset %in% c("Deer_Norwalk_2022", "Deer_Bridgeport_2022"))
dat_all_swabs_deer$No_swabs <- rep(3, nrow(dat_all_swabs_deer)) # 3 swabs per mouse per collection

# Combine swab data for both groups
dat_all_swabs <- rbind.data.frame(dat_all_swabs_mice, dat_all_swabs_deer)

# Calculate sample size / total swabs for each group
ss_mice_res <- nrow(dat_all_swabs %>% dplyr::filter(Dataset == "Mice_2022_Guilford")) * 2
ss_mice_for <- nrow(dat_all_swabs %>% dplyr::filter(Dataset == "Mice_2022_NB")) * 2
ss_deer_nor <- nrow(dat_all_swabs %>% dplyr::filter(Dataset == "Deer_Norwalk_2022")) * 3
ss_deer_bri <- nrow(dat_all_swabs %>% dplyr::filter(Dataset == "Deer_Bridgeport_2022")) * 3

# Create labels for plotting
group_labels <- c("Mice_Residential" = paste0("Mice-Guilford (Residential) (N=", ss_mice_res, ")"), 
                  "Mice_Forested" = paste0("Mice-North Branford (Forested) (N=", ss_mice_for, ")"), 
                  "Deer_Norwalk" = paste0("Deer-Norwalk (Residential/Forested) (N=", ss_deer_nor, ")"), 
                  "Deer_Bridgeport" = paste0("Deer-Bridgeport (Residential/Forested) (N=", ss_deer_bri, ")"))

# Break into separate datasets
mice1 <- dat_all_swabs %>% dplyr::filter(Group == "Mice_Residential")
mice2 <- dat_all_swabs %>% dplyr::filter(Group == "Mice_Forested")
deer1 <- dat_all_swabs %>% dplyr::filter(Group == "Deer_Norwalk")
deer2 <- dat_all_swabs %>% dplyr::filter(Group == "Deer_Bridgeport")

# Visually inspect counts
ddply(dat_all_swabs, .(Group), summarize, Total_Swabs = sum(No_swabs))

# Read in weekly estimated infections (Covidestim) in CT
dat_ce <- read_csv("covidestim_estimates_012623.csv")
dat_ce <- dat_ce %>% dplyr::filter(state == "Connecticut") %>% 
  dplyr::select(state, date, infections) %>%
  dplyr::rename(Date = date)

# Read in CT state population
dat_pop <- read_csv("NST-EST2022-ALLDATA.csv")
dat_pop <- dat_pop %>% 
  dplyr::filter(STATE == "09") %>% 
  dplyr::select(POPESTIMATE2022)

# Calculate infections per 100K pop in CT
dat_ce$infections_per_100K <- (dat_ce$infections / dat_pop$POPESTIMATE2022) * 100000

# Create full range of dates
dat_ce <- dat_ce[order(dat_ce$Date, decreasing = FALSE), ]
dat_ce$Week <- seq(1, nrow(dat_ce)) # not calendar weeks
daily_dates <- seq.Date(min(dat_ce$Date), max(dat_ce$Date), by = "day")

# Summarize number swabs by date and group
dat_swab_count <- ddply(dat_all_swabs, .(Date, Group), summarize, No_Samp = sum(No_swabs)) 
dat_swab_count_cast <- reshape2::dcast(dat_swab_count, Date ~ Group, value.var = "No_Samp") 

# Join daily dates, weekly estimated infections, and daily number of swabs
dat_daily_inf <- data.frame(Date = daily_dates) %>% dplyr::left_join(dat_ce, by = "Date")
dat_daily_inf_swab <- dat_daily_inf %>% dplyr::left_join(dat_swab_count_cast, by = "Date")

# Select columns of interest and format
dat_daily_inf <- dat_daily_inf_swab %>% dplyr::select(Date, 
                                                      Week,
                                                      infections_per_100K) 
dat_daily_swab <- dat_daily_inf_swab %>% dplyr::select(Date, 
                                                       Mice_Residential, 
                                                       Mice_Forested,
                                                       Deer_Norwalk,
                                                       Deer_Bridgeport)

dat_daily_swab_melt <- reshape2::melt(dat_daily_swab, id.vars = "Date")
dat_daily_both <- dat_daily_swab_melt %>% dplyr::left_join(dat_daily_inf, by = "Date")
all_weeks <- rep(seq(1, 63, by = 1), each = 7) # fill in missing week info
all_weeks <- c(all_weeks, 64)
dat_daily_both$Week_Rep <- all_weeks

# Summarize swab counts by group and week to match estimated infections
dat_daily_both$value <- as.numeric(dat_daily_both$value)
dat_weekly_both <- ddply(dat_daily_both, .(variable, Week_Rep), summarize, No_Swabs = sum(value, na.rm=TRUE))
dat_weekly_both <- dat_weekly_both %>% dplyr::rename(Week = Week_Rep)
dat_weekly_both <- dat_weekly_both %>% dplyr::left_join(dat_daily_inf, by = "Week") # add est infections back
dat_weekly_both$Month_Yr <- as.yearmon(dat_weekly_both$Date)

# Determine coefficient of monthly max infections per 100K vs. max no. of swabs to set secondary y axis
swabs_per_mth <- ddply(dat_weekly_both, .(Month_Yr), summarize, No_Swabs_Month = sum(No_Swabs, na.rm = TRUE))
swabs_per_mth_max <- max(swabs_per_mth$No_Swabs_Month, na.rm = TRUE)
coeff <- max(dat_weekly_both$infections_per_100K, na.rm = TRUE)/swabs_per_mth_max

# Format variables
dat_weekly_both$Month_Yr <- as.Date(dat_weekly_both$Month_Yr, format = "%Y-%m")
dat_weekly_both$Date <- as.Date(dat_weekly_both$Date)

# Select colors
viridis_names <-c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
n <- 10
par(mfrow=c(4,2), mar = c(1,1,1,1))
f <- sapply(1:8, function(x) image(matrix(1:n, n, 1), col = viridis(n=n, option = LETTERS[x]), axes =FALSE, main = viridis_names[x]))

swab_group_colors <- c(viridis::viridis(4)[1],
                       viridis::viridis(4)[2],
                       viridis::viridis(4)[3],
                       viridis::viridis(4)[4])

p <- ggplot(dat_weekly_both) + 
  theme_bw() + 
  geom_line(aes(x = Date, y = infections_per_100K), color = "black") +
  geom_bar(stat = "identity", aes(x = Month_Yr, y = No_Swabs*coeff, fill = variable)) +
  scale_y_continuous(
    name = "Est. Weekly Human Inf. per 100K Pop. in CT",
    sec.axis = sec_axis(~./coeff, name="Monthly No. Swab Samples")) +
  xlab("Month-Year") +
  labs(fill = "Sampled Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_date(date_labels = "%b-%Y", date_breaks  ="1 month") +
  scale_fill_manual(labels = c("Mice 1-Guilford (Residential)", 
                               "Mice 2-North Branford (Forested)", 
                               "Deer 1-Norwalk (Residential/Forested)",
                               "Deer 2-Bridgeport (Residential/Forested)"),
                    values = swab_group_colors) +
  theme(strip.text.x = element_text(size = 14),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11), 
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow=2, byrow=FALSE))

ggsave(paste(path_figures, "Fig_S4.jpg", sep=""), p, height = 5, width =10)
