
# Code Description --------------------------------------------------------

# 1) Create Fig 1 (Sampling Map + sVNT results)
# 2) Create Fig 3 (PanCoV PCR Map + Phylogenetic Tree)
# 3) Create Supp Fig 2 (sVNT Map)
# Note: Raw data available upon request

# Load libraries + data + set paths ---------------------------------------------------

library(readxl)
library(ggmap)
library(plyr)
library(dplyr)
library(magick)
library(cowplot)

file_loc <- getwd()
sheet_date <- "053023"
# register_google(key = "", write = TRUE) # insert key for geocoding


# dat_2020 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                        sheet = "Mice_Guilford_2020")
# dat_2021 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                        sheet = "Mice_Guilford_2021")
# dat_2022 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                        sheet = "Mice_Guilford_2022")
# dat_2022_nb <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                           sheet = "Mice_North_Branford_2022")
# dat_2021_2022 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                             sheet = "Deer_2021_2022")

path_figures <- paste0(file_loc, "/Figures/")

# Summarize Pan-CoV PCR Results -------------------------------------------------

# Select columns
select_geocode_col <- function(dat) {
  dat_new <- dat %>% dplyr::select(Date,
                                   Address,
                                   `Tag #`,
                                   `Recapture?`,
                                   `PanCov_Gel_Result_Final`)
  dat_new
}

dat_2022 <- select_geocode_col(dat_2022)

# Geocode addresses
dat_2022$Address_Full <- paste0(dat_2022$Address, ", Guilford", ", Connecticut") # create full address
# dat_gc <- mutate_geocode(dat_2022, location = Address_Full, output = "latlona")
# save(dat_gc,file="dat_pancov_geocoded.Rda") # save so don't have to repeat
load(paste0(file_loc, "/dat_pancov_geocoded.Rda"))

# Change maybe results to negative (metagenomic sequencing revealed they were indeed negative)
dat_gc[which(dat_gc$PanCov_Gel_Result_Final == "Maybe"), "PanCov_Gel_Result_Final"] <- "Neg"

# Set labeling
dat_gc$PanCov_Gel_Result_Final <- factor(dat_gc$PanCov_Gel_Result_Final, 
                                         levels = c("Pos", "Neg"),
                                         labels = c("Positive", "Negative"))

dat_gc$PanCoV_Lineage <- "Negative"
dat_gc[which(dat_gc$`Tag #` %in% c(216, 219)), "PanCoV_Lineage"] <- "Positive (PCoV1)"
dat_gc[which((dat_gc$`Tag #` == 225) & (dat_gc$`Recapture?` == "No")), "PanCoV_Lineage"] <- "Positive (PCoV1)"
dat_gc[which(dat_gc$`Tag #` %in% c(445, 446)), "PanCoV_Lineage"] <- "Positive (PCoV2)"
dat_gc$PanCoV_Lineage <- as.factor(dat_gc$PanCoV_Lineage)

# Fix column names so successfully import into ArcGIS Online
dat_gc <- dat_gc %>% dplyr::rename("Tag_No" = "Tag #",
                                   "Recapture" = "Recapture?")
dat_gc <- dat_gc %>% dplyr::select(!(Address))

# Jitter coordinates to improve visibility
dat_gc$lat <- jitter(dat_gc$lat, amount = 0.0005)
dat_gc$lon <- jitter(dat_gc$lon, amount = 0.0005)

# Save for import into ArcGIS online to create map
write.csv(dat_gc, file = paste0("dat_pancov", ".csv"))

# Create figure for pancov results (Fig 3) --------------------------------

img_pancov <-  image_read(paste0(path_figures, "pancov_map.jpg")) # load pancov map created in ArcGIS
img_pancov <- image_trim(img_pancov) # trim white space around image
p_pancov2 <- image_ggplot(img_pancov, interpolate = TRUE) # create ggplot object
p_pancov2 <- p_pancov2 + theme(plot.margin = unit(c(0,0,10,0), "pt")) # add small amnt of white space around margin

img_tree <- image_read(paste0(path_figures, "tree_entire.jpg")) # load phylogenetic tree image
img_tree <- image_trim(img_tree) 
p_tree <- image_ggplot(img_tree, interpolate = TRUE)

p <- plot_grid(p_pancov2, p_tree, # create combined image
               axis = "tl",
               nrow = 2, 
               ncol = 1,
               labels = c("A", "B"),
               hjust = -2)
ggsave(paste(path_figures, "Fig_3.jpg", sep=""), p)

# Summarize sVNT Results -------------------------------------------------

# Select columns 
select_geocode_col <- function(dat) {
  dat_new <- dat %>% dplyr::select(Date,
                                   Address,
                                   `Tag #`,
                                   `Recapture?`,
                                   `cPass_Result_Run1`,
                                   `cPass_Result_Run2`)
  dat_new
}

dat_2020 <- select_geocode_col(dat_2020)
dat_2021 <- select_geocode_col(dat_2021)

dat_2022 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
                       sheet = "Mice_Guilford_2022")
dat_2022 <- select_geocode_col(dat_2022)

# Combine all datasets across years
dat_all <- rbind.data.frame(dat_2020, dat_2021, dat_2022)

# Geocode addresses
dat_all$Address_Full <- paste0(dat_all$Address, ", Guilford", ", Connecticut") # create full address for geocoding
# dat_gc <- mutate_geocode(dat_all, location = Address_Full, output = "latlona")
# save(dat_gc, file = "dat_nab_tested_geocoded.Rda")
load(paste0(file_loc, "/dat_nab_tested_geocoded.Rda"))

# Format result categories
dat_all[which(is.na(dat_all$cPass_Result_Run1) == TRUE), "cPass_Result_Run1"] <- "Untested" # will only plot tested samples

dat_gc[which(dat_gc$cPass_Result_Run1 == "Untested"), "cPass_Result_Run2"] <- "Untested" 
dat_gc[which(dat_gc$cPass_Result_Run1 == "Neg"), "cPass_Result_Run2"] <- "Negative" # these are samples negative on 1st run
dat_gc[which(dat_gc$cPass_Result_Run2 == "Neg"), "cPass_Result_Run2"] <- "Negative" # standardizing term (both neg and negative)

dat_gc[which(dat_gc$cPass_Result_Run2 == "Pos"), "cPass_Result_Run2"] <- "Positive" # combining positive categories
dat_gc[which(dat_gc$cPass_Result_Run2 == "Elevated"), "cPass_Result_Run2"] <- "Positive" 
dat_gc[which(dat_gc$cPass_Result_Run2 == "Insufficient sample"), "cPass_Result_Run2"] <- "Positive" 

dat_gc <- dat_gc %>% dplyr::filter(dat_gc$cPass_Result_Run2 != "Untested") # remove untested

dat_gc$cPass_Result_Run2 <- factor(dat_gc$cPass_Result_Run2, 
                                   levels = c("Positive", "Negative")) # will use final run 2 result

# Fix column names so successfully import into ArcGIS Online
dat_gc <- dat_gc %>% dplyr::rename("Tag_No" = "Tag #",
                                   "Recapture" = "Recapture?")
dat_gc <- dat_gc %>% dplyr::select(!(Address))

# Jitter coordinates for visibility
dat_gc$lat <- jitter(dat_gc$lat, amount = 0.0005)
dat_gc$lon <- jitter(dat_gc$lon, amount = 0.0005)

# Label group
dat_gc[which((dat_gc$Tag_No %in% c(312, 320)) & (dat_gc$cPass_Result_Run2 == "Positive")), "Group"] <- "Group 1"

# Save for import into ArcGIS online to create map
write.csv(dat_gc, file = paste0("dat_svnt", ".csv"))

# Create figure for svnt map (Supp Fig 2) --------------------------------

img_svnt <-  image_read(paste0(path_figures, "svnt_map.jpg"))
img_svnt <- image_trim(img_svnt)
p_svnt <- image_ggplot(img_svnt, interpolate = TRUE)

p <- plot_grid(p_svnt,
               axis = "tl",
               nrow = 1, 
               ncol = 1)

ggsave(paste(path_figures, "Fig_S2.jpg", sep=""), p)

# Summarize sampling locations --------------------------------------------

# Load geocoded mice data
load(paste0(file_loc, "/dat_nab_tested_geocoded.Rda"))
dat_gc_mice_gu <- dat_gc
dat_gc_mice_gu$Group <- "Mice 1-Guilford (Residential)"
dat_gc_mice_gu <- dat_gc_mice_gu %>% dplyr::select(lon, lat, Group)

# Combine w/ forested mice & deer datasets
## do not have specific addresses; pick centroid of sampling area
dat_2022_nb$lon <- -72.775
dat_2022_nb$lat <- 41.350
dat_2022_nb$Group <- "Mice 2-North Branford (Forested)"
dat_2022_nb <- dat_2022_nb %>% dplyr::select(lon, lat, Group)

dat_2021_2022[which(dat_2021_2022$Location == "Norwalk"), "lon"] <- -73.410400
dat_2021_2022[which(dat_2021_2022$Location == "Norwalk"), "lat"] <- 41.076163
dat_2021_2022[which(dat_2021_2022$Location == "Norwalk"), "Group"] <- "Deer 1-Norwalk (Residential/Forested)"

dat_2021_2022[which(dat_2021_2022$Location == "Bridgeport"), "lon"] <- -73.171353
dat_2021_2022[which(dat_2021_2022$Location == "Bridgeport"), "lat"] <- 41.205805
dat_2021_2022[which(dat_2021_2022$Location == "Bridgeport"), "Group"] <- "Deer 2-Bridgeport (Residential/Forested)"

dat_2021_2022 <- dat_2021_2022 %>% dplyr::select(lon, lat, Group)

dat_gc_mice_gu[, "lat"] <- mean(dat_gc_mice_gu$lat) # take mean coordinates of all geocoded residential mice locations
dat_gc_mice_gu[, "lon"] <- mean(dat_gc_mice_gu$lon)

dat_all <- data.frame(rbind.data.frame(dat_gc_mice_gu, dat_2022_nb, dat_2021_2022))

# Count number of samples by group
count_by_group <- ddply(dat_all, .(Group), nrow)
dat_all_unique <- dat_all %>% 
  select(lon, lat, Group) %>% 
  dplyr::distinct()
dat_all_unique_count <- left_join(dat_all_unique, count_by_group, by = "Group")
dat_all_unique_count <- dat_all_unique_count %>% dplyr::rename("Number" = "V1")

# Set labeling
dat_all_unique_count$Label <- c("Mice 1", "Mice 2", "Deer 1", "Deer 2")
dat_all_unique_count <- dat_all_unique_count %>% dplyr::rename("Sample Size" = "Number")

# Save for import into ArcGIS Online to create map
write.csv(dat_all_unique_count, file = paste0("dat_all_unique_count", ".csv"))

# Create sampling map + sVNT results (Fig 1) ------------------------------

img_context <- image_read(paste0(path_figures, "pancov_map_all.jpg"))
img_context <- image_trim(img_context)
p_context <- image_ggplot(img_context, interpolate = TRUE)
p_context <- p_context + theme(plot.margin = unit(c(10,10,10,10), "pt"))

img_nab <-  image_read(paste0(path_figures, "Fig_NAb_pos.jpg"))
img_nab <- image_trim(img_nab)
p_nab <- image_ggplot(img_nab, interpolate = TRUE)

p <- plot_grid(p_context, p_nab,
               axis = "tl",
               nrow = 2, 
               ncol = 1,
               rel_heights = c(0.45, 0.55),
               labels = c("A", "B"),
               hjust = -5)

ggsave(paste(path_figures, "Fig_1.jpg", sep=""), p)
