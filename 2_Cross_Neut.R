
# Code Description --------------------------------------------------------

# 1) Create Fig 2 (variant-specific sVNT)
# 2) Create Supp Fig 3 (VNT)
# Note: raw data available upon request

# Load libraries + set paths ---------------------------------------------------

library(readxl)
library(plyr)
library(dplyr)
library(reshape2)
library(viridis)
library(ggplot2)

file_loc <- getwd()
path_figures <- paste0(file_loc, "/Figures/")

# Variant-Specific sVNT Results (Genscript cPass) -------------------------------------------

# Load data 
sheet_date <- "053023"

# dat <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                   sheet = "NAb_Pos_Cross_Neut")

# Clean up data
dat[which((dat$`Ear Tag #` == "35") & (dat$`Recapture?` == "Yes")), "Ear Tag #"] <- "35-2"
dat[which((dat$`Ear Tag #` == "98-2") & (dat$`Recapture?` == "Yes")), "Ear Tag #"] <- "Deer Neg. Ctl."
dat <- dat %>% 
  dplyr::rename("Tag_No" = `Ear Tag #`) %>% 
  dplyr::filter(Animal == "Deer") 

# Format result data
dat_sub <- dat %>% dplyr::select(Date, 
                                 Tag_No, 
                                 WT_New_1, WT_New_2, 
                                 Alpha_1, Alpha_2, 
                                 Delta_1, Delta_2, 
                                 Omicron_BA1_1, Omicron_BA1_2)
dat_sub_melt <- reshape2::melt(dat_sub, .(Tag_No, Date), variable.name = "Variant", value.name = "Perc_Inhib")
dat_sub_melt$Duplicate <- rep(c(1, 2), each = 8, times = 4) # each sample run in duplicate on the plate
dat_sub_melt$Variant_New <- rep(c("WT_New", "Alpha", "Delta", "Omicron_BA1"), each = 16, times = 1) # add names not specific to duplicate (i.e 1 or 2)
dat_sub_melt <- dat_sub_melt %>% dplyr::select(Tag_No, Date, Perc_Inhib, Duplicate, Variant_New) # drop columns no longer need

# Format controls data
dat_controls <- dat %>% dplyr::select(Date, 
                                      Tag_No, 
                                      WT_New_PosCtl_1, WT_New_PosCtl_2, 
                                      Alpha_PosCtl_1, Alpha_PosCtl_2,
                                      Delta_PosCtl_1, Delta_PosCtl_2, 
                                      Omicron_BA1_PosCtl_1, Omicron_BA1_PosCtl_2) 
dat_controls <- dat_controls[1, ] # controls are same across all samples
dat_controls_melt <- reshape2::melt(dat_controls, .(Tag_No, Date), variable.name = "Variant", value.name = "Perc_Inhib")
dat_controls_melt <- dat_controls_melt %>% dplyr::select(Variant, Perc_Inhib) # drop columns no longer need
dat_controls_melt$Variant_New <- rep(c("WT_New", "Alpha", "Delta", "Omicron_BA1"), each = 2, times = 1) # add variant name so can match to dat_sub_melt
dat_controls_melt$Date <- rep(NA, nrow(dat_controls_melt)) # need same columsn to bind with dat_sub_melt
dat_controls_melt$Duplicate <- rep(c(1, 2), each = 1, times = 4)
dat_controls_melt$Tag_No <- rep("PosCtl", nrow(dat_controls_melt))
dat_controls_melt <- dat_controls_melt[, c("Tag_No", "Date", "Perc_Inhib", "Duplicate", "Variant_New")]

# Combine sample and controls data
dat_sub_melt <- rbind.data.frame(dat_sub_melt, dat_controls_melt)
dat_sub_melt$Perc_Inhib <- as.numeric(dat_sub_melt$Perc_Inhib)

# Update labels and determine order based on average WT percent inhibition value
dat_sub_melt$Variant_New <- factor(dat_sub_melt$Variant_New,
                                   levels = c("WT_New", "Alpha", "Delta", "Omicron_BA1"),
                                   labels = c("Wild Type", "Alpha (B.1.1.7)", "Delta (B.1.617.2)", "Omicron (BA.1)"))

dat_sub_melt_WT <- dat_sub_melt %>% dplyr::filter(Variant_New == "Wild Type")
dat_order <- ddply(dat_sub_melt_WT, .(Tag_No), summarize, Mean_Perc_Inhib = mean(Perc_Inhib), Date = Date) # calculate avg % inhibition across duplicates
sample_order <- dat_order[order(dat_order$Mean_Perc_Inhib, decreasing = TRUE), "Tag_No"]
date_order <- as.Date(dat_order[order(dat_order$Mean_Perc_Inhib, decreasing = TRUE), "Date"])
label_order <- rep(c("Kit Wild Type Pos. Ctl.", # rename b/c sample IDs don't mean anything to external audience
                     "1 (2021-06-15)", # 94
                     "2 (2021-06-15)", # 35
                     "3 (2022-07-12)", # 37
                     "3R (2022-07-28)", # 37 R
                     "4 (2022-01-12)", # 70
                     "2R (2021-07-22)", # 35R
                     "5 (2022-01-12)", # 71
                     "Deer Neg. Ctl."), each = 2)
# Select colors
viridis_names <-c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
n <- 9
par(mfrow=c(4,2), mar = c(1,1,1,1))
f <- sapply(1:8, function(x) image(matrix(1:n, n, 1), col = viridis(n=n, option = LETTERS[x]), axes =FALSE, main = viridis_names[x]))

# Set up tracker for order, color, and shape data - will need for VNT to match
dat_key <- cbind.data.frame(Tag_No = sample_order, Label = label_order)
num_deer <- nrow(dat_key) / 2 # divide by 2 due to duplicates
dat_key$Color <- rep(c("black", viridis::turbo(n=7), "gray30"), each = 2)
dat_key$Shape <- rep(seq(0, (num_deer-1), by = 1), each = 2) # subtract 1 because start at 0

# Drop Omicron BA.1 b/c did not have correct positive control due to immune escape properties
dat_sub_melt <- dat_sub_melt %>% dplyr::filter(Variant_New != "Omicron (BA.1)") 

# Plot percent inhibition relative to he negative control by variant per deer
p <- ggplot(dat_sub_melt) + 
  theme_bw() +
  geom_line(aes(x= Variant_New, y = Perc_Inhib, color = Tag_No, group = interaction(Tag_No, Duplicate), linetype = Tag_No)) +
  geom_point(aes(x= Variant_New, y = Perc_Inhib, color = Tag_No, group = interaction(Tag_No, Duplicate), shape = Tag_No)) +
  geom_hline(yintercept = 30, linetype = "dashed", color = "black", size = 0.5) +
  xlab("SARS-CoV-2 Variant") +
  ylab("Percent Signal Reduction Rel. to Negative Control") +
  annotate("text", x = 3.3, y = 32, label="≥ 30% = Positive", fontface = "bold") +
  scale_color_manual(name = "ID (Sample Date)", 
                     values = rep(c("black", viridis::turbo(n = 7), "gray30"), each = 2), 
                     breaks = sample_order, 
                     labels = label_order) +
  scale_linetype_manual(name = "ID (Sample Date)", 
                        values = c("PosCtl" = "dotted", 
                                   "94" = "solid",
                                   "35" = "solid",
                                   "37" = "solid",
                                   "37-2" = "solid",
                                   "70" = "solid",
                                   "35-2" = "solid",
                                   "71" = "solid",
                                   "Deer Neg. Ctl." = "dotted"),
                        breaks = sample_order,
                        labels = label_order) +
  scale_shape_manual(name = "ID (Sample Date)", 
                     values = rep(seq(0, (num_deer-1), by = 1), each = 2), 
                     breaks = sample_order,
                     labels = label_order) +
  labs(shape = "ID (Sample Date)", color = "ID (Sample Date)", linetype = "ID (Sample Date)") +
  theme(axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 12))

ggsave(paste(path_figures, "Fig_2.jpg", sep=""), p, height = 5, width =10)

# Variant-Specific VNT results ------------------

# Load data 
dat <- read_excel(path = paste0(file_loc, "/VNT_Results_Formatted", ".xlsx"), .name_repair = "minimal")

# Format data
dat_melt <- reshape2::melt(dat, .(variant, plate, sample), variable.name = "Deer", value.name = "Viability")
dat_melt$Percent_Viability <- dat_melt$Viability * 100
dat_melt$Deer <- as.character(dat_melt$Deer)

# Split deer name and duplicate into separate columns
dat_melt$Duplicate <- rep(NA, nrow(dat_melt))
dat_melt[which((endsWith(dat_melt$Deer, "_1")) == TRUE), "Duplicate"] <- 1
dat_melt[which((endsWith(dat_melt$Deer, "_2")) == TRUE), "Duplicate"] <- 2
dat_melt$Deer_No <- substr(dat_melt$Deer, start = 1, stop = 2) # grab first two characters

# Assign matching IDs
dat_melt$Tag_No <- rep(NA, nrow(dat_melt))
dat_melt[which(dat_melt$Deer_No == "D1"), "Tag_No"] <- "35"
dat_melt[which(dat_melt$Deer_No == "D2"), "Tag_No"] <- "94"
dat_melt[which(dat_melt$Deer_No == "D3"), "Tag_No"] <- "35-2"
dat_melt[which(dat_melt$Deer_No == "D4"), "Tag_No"] <- "70"
dat_melt[which(dat_melt$Deer_No == "D5"), "Tag_No"] <- "37-2"
dat_melt[which(dat_melt$Deer_No == "D6"), "Tag_No"] <- "Deer Neg. Ctl."
dat_melt[which(dat_melt$Deer_No == "HU"), "Tag_No"] <- "PosCtl"

# Set factors for labeling
dat_melt$sample <- factor(dat_melt$sample,
                          levels = c("1:20", "1:60", "1:180", "1:540", "1:1620", "1:4860", "1:14580", "1:43740"))
dat_melt$variant <- factor(dat_melt$variant,
                           levels = c("WA01", "Delta", "BA5"),
                           labels = c("Wild Type", "Delta (B.1.617.2)", "Omicron (BA.5)"))
dat_melt$plate <- factor(dat_melt$plate,
                         levels = c("1", "2"),
                         labels = c("Plate 1", "Plate 2"))

# Order by descending mean percent cell viability
dat_melt_WT <- dat_melt %>% dplyr::filter(variant == "Wild Type")
dat_melt_WT_mean <- ddply(dat_melt_WT, .(Tag_No), summarize, Mean_Perc_Via = mean(Percent_Viability))
dat_order <- dat_melt_WT_mean %>% dplyr::filter(Tag_No != "PosCtl") # drop controls
dat_order <- dat_order %>% dplyr::filter(Tag_No != "Deer Neg. Ctl.")
sample_order_new <- dat_order[order(dat_order$Mean_Perc_Via, decreasing = TRUE), "Tag_No"]
sample_order_new <- rep(sample_order_new, each = 2) # for duplicates
sample_order_new <- c(rep("PosCtl", each = 2), sample_order_new, rep("Deer Neg. Ctl.", each = 2)) # add controls

# Match colors and shapes from dat_key from sVNT
dat_key_new <- data.frame(Tag_No = sample_order_new)
dat_key_new <- distinct(dat_key_new) # o/w will duplicate
dat_key_new <- dat_key_new %>% dplyr::left_join(dat_key, by = "Tag_No") # join to sVNT key so can match on colors/shapes/labels
label_order_new <- dat_key_new$Label
label_order_new[1:2] <- "Human Pos. Ctl."
color_order_new <- dat_key_new$Color
shape_order_new <- dat_key_new$Shape

p <- ggplot(dat_melt) +
  theme_bw() +
  geom_line(aes(x = sample, y = Percent_Viability, color = Tag_No, group = interaction(Tag_No, Duplicate), linetype = Tag_No)) +
  geom_point(aes(x = sample, y = Percent_Viability, color = Tag_No, group = interaction(Tag_No, Duplicate), shape = Tag_No)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", size = 0.5) +
  annotate("text", x = 4.85, y = 110, label="mock", fontface = "bold", size = 3) +
  annotate("text", x = 6, y = 90, label="(no sera, no virus)", fontface = "bold", size = 3) +
  xlab("Sera Dilution") +
  ylab("Percent Cell Viability Relative to Mean Mock Control") +
  facet_wrap(~plate) +
  ylim(0, 225) +
  scale_linetype_manual(name = "ID (Sample Date)", 
                        values = c("PosCtl" = "dotted", 
                                   "94" = "solid",
                                   "35" = "solid",
                                   "70" = "solid",
                                   "37-2" = "solid",
                                   "35-2" = "solid",
                                   "Deer Neg. Ctl." = "dotted"),
                        breaks = sample_order_new,
                        labels = label_order_new) +
  scale_color_manual(name = "ID (Sample Date)",
                     values = color_order_new,
                     breaks = sample_order_new, 
                     labels = label_order_new) +
  scale_shape_manual(name = "ID (Sample Date)", 
                     values= shape_order_new, 
                     breaks = sample_order_new,
                     labels = label_order_new) +
  facet_grid(plate ~ variant) +
  labs(shape = "ID (Sample Date)", color = "ID (Sample Date)", linetype = "ID (Sample Date)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11),
        axis.title = element_text(size = 12))

ggsave(paste(path_figures, "Fig_S3.jpg", sep=""), p, height = 5, width =10)
