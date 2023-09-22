
# Code Description --------------------------------------------------------

# 1) Create Supp Fig 7 (ELISA data)
# Note: raw data available upon request

# Load libraries + set paths ---------------------------------------------------

library(readxl)
library(plyr)
library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(ggpubr)

file_loc <- getwd()
sheet_date <- "053023"
path_figures <- paste0(file_loc, "/Figures/")

## 2020-2021 data for Guilford Mice
### Note: do not need 2022 data b/c only prescreened 2020/2021 samples
# dat_2020 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                        sheet = "Mice_Guilford_2020")
# dat_2021 <- read_excel(path = paste0(file_loc, "/Metadata_Sheet/All_Sample_Metadata_Results_", sheet_date, ".xlsx"),
#                        sheet = "Mice_Guilford_2021")

# Format Data -------------------------------------------------------------

# Add dataset identifier
dat_2020$Dataset <- rep("Mice_2020", nrow(dat_2020))
dat_2021$Dataset <- rep("Mice_2021", nrow(dat_2021))

# Keep only columns of interest
keep_selected_col <- function(dat) {
  dat_new <- dat %>% dplyr::select(Date, 
                                   'Tag #', 
                                   'Recapture?',
                                   ELISA_Plate_No,
                                   ELISA_S_Protein,
                                   ELISA_N_Protein,
                                   ELISA_BSA_Protein,
                                   ELISA_S_Protein_PosCtl,
                                   ELISA_N_Protein_PosCtl,
                                   ELISA_S_Protein_NegCtl,
                                   ELISA_N_Protein_NegCtl,
                                   ELISA_BSA_Protein_NegCtl,
                                   cPass_Result_Run1,
                                   cPass_Result_Run2,
                                   Dataset)
  dat_new
}

dat_2020_new <- keep_selected_col(dat_2020)
dat_2021_new <- keep_selected_col(dat_2021)

# Combine datasets
dat_all <- rbind.data.frame(dat_2020_new, 
                            dat_2021_new)
dat_all <- as.data.frame(dat_all)

# Change protein value columns to numeric
convert_char_to_numeric <- function(dat, col_name) {
  dat[, col_name] <- as.numeric(dat[, col_name])
  dat_new <- dat
  dat_new
}

dat_all <- convert_char_to_numeric(dat_all, "ELISA_S_Protein")
dat_all <- convert_char_to_numeric(dat_all, "ELISA_N_Protein")
dat_all <- convert_char_to_numeric(dat_all, "ELISA_BSA_Protein")

dat_all <- convert_char_to_numeric(dat_all, "ELISA_S_Protein_PosCtl")
dat_all <- convert_char_to_numeric(dat_all, "ELISA_N_Protein_PosCtl")

dat_all <- convert_char_to_numeric(dat_all, "ELISA_S_Protein_NegCtl")
dat_all <- convert_char_to_numeric(dat_all, "ELISA_N_Protein_NegCtl")
dat_all <- convert_char_to_numeric(dat_all, "ELISA_BSA_Protein_NegCtl")

class(dat_all$ELISA_S_Protein) # test a few columns
class(dat_all$ELISA_BSA_Protein_NegCtl)

# Convert plate number to factor
dat_all$ELISA_Plate_No <- factor(dat_all$ELISA_Plate_No, 
                                 levels = c("1", "2", "3", "4", "5/6", "7/8", "9/10", "11/12", "13/14", "15/16"))

# Calculate relative values
dat_all$S_div_BSA <- dat_all$ELISA_S_Protein / dat_all$ELISA_BSA_Protein # for Plates 1-4 (2021 samples)
dat_all$S_div_Sctl <- dat_all$ELISA_S_Protein / dat_all$ELISA_S_Protein_PosCtl # for Plate 5/6 - 15/16 (2020 samples)
dat_all$N_div_Nctl <- dat_all$ELISA_N_Protein / dat_all$ELISA_N_Protein_PosCtl # for Plate 5/6 - 15/16 (2020 samples)


# Plot ELISA Plates 1-4 Data ----------------------------------------------
## S protein / BSA protein value w/ cutoff for cPass
## Add neg control (lack pos control)

dat_p1_to_p4 <- dat_all %>% dplyr::filter(ELISA_Plate_No %in% c("1", "2", "3", "4"))

unique_negctl_S <- unique(dat_p1_to_p4[, c("ELISA_Plate_No", "ELISA_S_Protein_NegCtl")])
unique_negctl_BSA <-  unique(dat_p1_to_p4[, c("ELISA_Plate_No", "ELISA_BSA_Protein_NegCtl")])
unique_negctl_both <-  unique_negctl_S %>% dplyr::left_join(unique_negctl_BSA, by = "ELISA_Plate_No")
unique_negctl_S_div_BSA <- ddply(unique_negctl_both, .(ELISA_Plate_No), summarize, S_div_BSA = ELISA_S_Protein_NegCtl/ELISA_BSA_Protein_NegCtl)
unique_negctl_S_div_BSA <- cbind.data.frame(unique_negctl_S_div_BSA, Type = rep("negctl", 4))
dat_p1_to_p4_sub <- dat_p1_to_p4 %>% dplyr::select(ELISA_Plate_No, S_div_BSA)
dat_p1_to_p4_sub$Type <- rep("Sample", nrow(dat_p1_to_p4_sub))
dat_p1_to_p4_ctl <- rbind.data.frame(dat_p1_to_p4_sub, unique_negctl_S_div_BSA)

p <- ggplot(dat_p1_to_p4_ctl, aes(x= ELISA_Plate_No, y = S_div_BSA)) + 
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Type), position = position_jitter(width = 0.1), alpha = 0.5) +
  geom_hline(yintercept = 1.57, linetype = "dashed", color = "black", size = 0.5) +
  labs("Plate") +
  xlab("Plate") +
  ylab("Spike Protein / BSA Protein OD Value") +
  scale_color_manual(values = c("red", "black")) +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11))
p
plate_1_4 <- p

# Plot ELISA Plates 5-16 Data ---------------------------------------------

dat_p5_to_p16 <- dat_all %>% dplyr::filter(ELISA_Plate_No %in% c("5/6", "7/8", "9/10", "11/12", "13/14", "15/16"))
unique(dat_p5_to_p16$ELISA_Plate_No)

group_labels <- c("5/6" = "Plate 5/6", 
                  "7/8" = "Plate 7/8",
                  "9/10" = "Plate 9/10",
                  "11/12" = "Plate 11/12",
                  "13/14" = "Plate 13/14",
                  "15/16" = "Plate 15/16")

unique_negctl_S <- unique(dat_p5_to_p16[, c("ELISA_Plate_No", "ELISA_S_Protein_NegCtl", "ELISA_S_Protein_PosCtl")])
unique_negctl_N <-  unique(dat_p5_to_p16[, c("ELISA_Plate_No", "ELISA_N_Protein_NegCtl", "ELISA_N_Protein_PosCtl")])
unique_negctl_both <-  unique_negctl_S %>% dplyr::left_join(unique_negctl_N, by = "ELISA_Plate_No")
unique_negctl_S_N <- ddply(unique_negctl_both, .(ELISA_Plate_No), summarize, 
                           S_div_Sctl = ELISA_S_Protein_NegCtl/ELISA_S_Protein_PosCtl,
                           N_div_Nctl = ELISA_N_Protein_NegCtl/ELISA_N_Protein_PosCtl)
unique_negctl_S_N <- cbind.data.frame(unique_negctl_S_N, Type = rep("negctl", 6))
dat_p5_to_p16_sub <- dat_p5_to_p16 %>% dplyr::select(ELISA_Plate_No, S_div_Sctl, N_div_Nctl)
dat_p5_to_p16_sub$Type <- rep("Sample", nrow(dat_p5_to_p16_sub))
dat_p5_to_p16_ctl <- rbind.data.frame(dat_p5_to_p16_sub, unique_negctl_S_N)

p <- ggplot(dat_p5_to_p16_ctl) + 
  geom_point(aes(x= N_div_Nctl, y = S_div_Sctl, color = Type), alpha = 0.5) + 
  theme_bw() +
  facet_wrap(~ELISA_Plate_No, labeller = labeller(ELISA_Plate_No = group_labels)) +
  # scale_color_discrete("Plate Number") +
  scale_color_manual(values = c("red", "black")) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  xlab("Sample / Positive Control OD Value (N Protein)") +
  ylab("Sample / Positive Control OD Value (S Protein)") +
  theme(strip.text.x = element_text(size = 12),
        legend.position = "none",  
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11)) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha = 0.2) +
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = 1, alpha = 0.2) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 2.5, alpha = 0.2) 

p
plate_5_16 <- p

# Combine All ELISA Plate Data --------------------------------------------

p <- grid.arrange(ggarrange(
  plate_1_4 + annotate("text", x=4.4, y=1.66, label="1.57x", fontface = "bold", color = "black"),
  plate_5_16,
  labels = c("A", "B")))

ggsave(paste(path_figures, "Fig_S7.jpg", sep=""), p, height = 5, width =10)

