library(tidyverse)
library(janitor)
library(lubridate)
library(here)
library(stringr)
library(psych)
library(ppcor)
library(car)
library(haven)
library(corrplot)

# =============== Rating Data Cleaning =================

setwd("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Fall_2022_rating_analyses")

rating <- read_csv("1st_batch_ratings.txt", col_names = c("trial", "sequence"))

# separate sequence number from the rest based on first space
sep_seq_num <- rating %>%
  separate(sequence, into = c("seq_num", "sequence"), sep = " ", extra = "merge") %>%
  mutate(seq_num = as.numeric(seq_num))

# check those whose seq_num weren't recorded
# those repeated consecutively don't have their seq_num recorded on the 2nd occurrance
tbd <- sep_seq_num %>% # to be deleted
  filter(seq_num < 13) %>% 
  # if seq_num > 12 we are sure that the seq_num is recorded,
  # if not we are unsure whether seq_num recorded is actual seq_num or first note of sequence,
  # so we select all seq_num < 13 and figure out manually
  arrange(seq_num) %>%
  slice(38:41)
# sequence "9 5 6 4 2 8 3 3 3 8 3" is labeled as #4, but is actually not found in the sequences list
# these 4 trials are consecutive; delete those from data & only maintain ratings on the first occurrance
del_rep <- anti_join(sep_seq_num, tbd)

# extract rating from sequence (run only once per session)
sep_seq_num$rating <- sub(".*\\s", "", sep_seq_num$sequence)

# extract rating using excel
# (https://www.extendoffice.com/excel/formulas/excel-extract-text-after-last-instance-of-character.html)

seq_rat <- sep_seq_num %>%
  mutate(rating = str_remove(rating, ";")) %>%
  mutate(rating = as.numeric(rating)) %>%
  group_by(seq_num) %>%
  summarize(mean_rat = mean(rating, na.rm = TRUE))

# seq 0-136 (137 never got rated)

# load sequence generations; delete; & add id manually in excel beforehand
sequences <- read_tsv("Summer_2022_sequences_bp.txt")

# convert text to list of numbers
joined <- left_join(sequences, seq_rat, by = "seq_num") %>%
  # convert seq chr to list of numbers for calculation purposes, keep og col for exporting
  mutate(seq_list = str_extract_all(seq, "[0-9]+"))

# check duplicated sequences
sequences %>% 
  filter(sequences %>% 
           dplyr::select(seq) %>% 
           duplicated())
# 44 duplicates

# if we were to combine duplicated sequences one person made into one seq
mean_rat <- joined %>%
  group_by(seq, seq_list, BP_Name) %>%
  summarize(mean_rating = mean(mean_rat, na.rm = TRUE))

# fill in length & # of unique pitches for each sequence
mean_rat$seq_length = 0
mean_rat$unique_pitches = 0

for (i in 1:(length(mean_rat$seq_list))) {
  mean_rat$seq_length[i] = length(mean_rat$seq_list[[i]])
  mean_rat$unique_pitches[i] = length(unique(mean_rat$seq_list[[i]]))
}

# turns out that the sequences that have 10+ unique pitches but still got rated low-creativity
# were 1-12, 12-1 or sth like that (low entropy)

# create a new var containing log2(length) to make the outlier (len=40) less influential
mean_rat$log_len <- log2(mean_rat$seq_length)

# create a new variable containing intervals in each sequence
int <- list()
for (i in 1:nrow(mean_rat)) {
  int[i] = list(diff(as.numeric(mean_rat$seq_list[[i]])))
}
mean_rat$interval <- int

# entropy - https://rstudio-pubs-static.s3.amazonaws.com/455435_30729e265f7a4d049400d03a18e218db.html
# table() - https://www.geeksforgeeks.org/create-a-tabular-representation-of-data-in-r-programming-table-function/

# compute Shannon entropy; input list of intervals, e.g. mean_rat$int[[1]]
entropy <- function(seq) {
  # calculate frequency for each interval 
  freq <- table(seq)/length(seq)
  # vectorize frequency list
  vec <- as.data.frame(freq)[[2]]
  # drop 0 to avoid NaN resulting from log2
  vec <- vec[vec > 0]
  # compute entropy
  -sum(vec * log2(vec))
}

mean_rat$entropy <- 0
for (i in 1:nrow(mean_rat)) {
  mean_rat$entropy[i] = entropy(mean_rat$interval[[i]])
}

write_csv(mean_rat, "mean_rat_1.csv")

# =============== Qualtrics Data Cleaning =================

# Qualtrics data before deleting blank & pilot: "BP_Fall_2022_Qualtrics_Original.csv"

# Qualtrics data after deleting complete blank & Psyche's pilot (deleted 5 pilot + 4 completely blank ones)
qualtrics_22 <- read_csv("BP_Fall_2022_Qualtrics.csv") %>%
  slice(-1, -2) %>%
  # change to date format
  mutate(StartDate = ymd_hms(StartDate, tz = "America/Denver"),
         # extract date and make them numeric
         date = as.numeric(gsub("-", "", date(StartDate))),
         # label semester
         batch = case_when(date < 20220501 ~ "SP22",
                           date >= 20220501 & date < 20220901 ~ "SM22",
                           date >= 20220901 ~ "FL22"),
         # create a column indicating if the participant is a rater (R) or a generator (G)
         # based on uploaded BP file name, or if they did not upload a file at all (NA)
         GeneRate = if_else(is.na(BP_Name), "NA", if_else(grepl("rating", BP_Name), "R", "G")),
         BP_Name = str_replace(BP_Name, "\\.TXT$", str_to_lower(".txt")),
         Lang_1st = str_to_title(Lang_1st), # convert 1st letter to uppercase
         Race = factor(Race, levels = c(1:8), 
                       labels = c("American Indian/Alaska Native", "Asian", "Black/African American", 
                                  "Hispanic/Latino", "Native Hawaiian/Pacific Islander", 
                                  "White", "2+ Races", "Prefer not to answer")),
         mus_train = if_else(Music == "No" | Music == "no" | Music == "No." | Music == "Ni" |
                               Music == "N/A" | is.na(Music), "No", "Yes"),
         .before = Name)

# examine those whose BP file type is wrong
qualtrics_22 %>%
  filter(BP_Type != "text/plain" | is.na(BP_Name)) %>%
  dplyr::select(StartDate, Name, BP_Name)
# 220428ALOO, 220523NA, 220613Boyuan didn't have BP file recorded, keep GeneRate as NA
# manually change EXUE and MLIU from NA to R, drop duplicated SGOT and DJAN rows & NA rows

# upon cross checking, manually edit the following classification (why did this change VPAT's GeneRate to NA?)
qualtrics_new_22 <- qualtrics_22 %>%
  mutate(GeneRate = if_else((Name == "Enxi" | Name == "Michelle Liu" | Name == "Paula"), "R", GeneRate),
         GeneRate = if_else(BP_Name == "221206_VPAT.txt", "G", GeneRate)) %>%
  rename(BMRQ_MS = SC6, BMRQ_EE = SC7, BMRQ_MR = SC8, BMRQ_SM = SC9, BMRQ_SR = SC10, BMRQ_AM = SC11, eBMRQ = SC12, 
         MSI_AE = SC13, MSI_PA = SC14, MSI_MT = SC15, MSI_SA = SC16, MSI_EM = SC17, MSI_GS = SC18)

# output each batch separately
qualtrics_1st_batch_g <- qualtrics_new_22 %>%
  filter(batch == "SM22" & GeneRate == "G")

qualtrics_1st_batch_r <- qualtrics_new_22 %>%
  filter((batch == "SM22" | batch == "FL22") 
         & GeneRate == "R" & 
           # delete 2 pilots
           Name != "Sarah sehgal" & Name != "Sophia Gitlin") %>%
  mutate(BP_Name = str_replace(BP_Name, "^(\\d{6})([^_])", "\\1_\\2"),
         # manually correct the incorrectly saved names
         BP_Name = ifelse(!grepl("_rating\\.txt$", BP_Name), paste0(BP_Name, "_rating.txt"), BP_Name),
         BP_Name = ifelse(BP_Name == "220310_PSEF.txt_rating.txt", "220310_PSEF_rating.txt", BP_Name),
         BP_Name = ifelse(BP_Name == "222909_VREM_rating.txt", "220929_VREM_rating.txt", BP_Name),
         BP_Name = ifelse(BP_Name == "222909_LOKI_rating.txt", "220929_LOKI_rating.txt", BP_Name),
         BP_Name = ifelse(BP_Name == "220310_PSEF_rating.txt", "221003_PSEF_rating.txt", BP_Name),
         BP_Name = ifelse(BP_Name == "221010_OSMI_rating.txt", "221110_OSMI_rating.txt", BP_Name)) %>%
  filter(BP_Name != "220928_XLI_rating.maxpat_rating.txt") # delete b/c we lost her BP data

qualtrics_2nd_batch_g <- qualtrics_new_22 %>%
  filter(batch == "FL22" & GeneRate == "G")

# save them into separate .csv files
write_csv(qualtrics_1st_batch_g, "qualtrics_g_1.csv")
write_csv(qualtrics_1st_batch_r, "qualtrics_r_1.csv")
write_csv(qualtrics_2nd_batch_g, "qualtrics_g_2.csv")

# ============= G ================
setwd("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_generation/Summer_2022_sequences")

# 3 people's generation data included previous participants' data, manually deleted in .txt before running

# define function to calculate entropy
entropy <- function(seq) {
  # calculate frequency for each interval 
  freq <- table(seq)/length(seq)
  # vectorize frequency list
  vec <- as.data.frame(freq)[[2]]
  # drop 0 to avoid NaN resulting from log2
  vec <- vec[vec > 0]
  # compute entropy
  -sum(vec * log2(vec))
}

qualtrics_1st_batch_g$num_seq <- 0
qualtrics_1st_batch_g$length <- 0
qualtrics_1st_batch_g$log_len <- 0
qualtrics_1st_batch_g$pitches <- 0
qualtrics_1st_batch_g$entropy <- 0

# create a loop that outputs each participant's data
for (file in qualtrics_1st_batch_g$BP_Name) {
  bp <- read_csv(file, col_names = c("trial", "seq")) %>%
    mutate(seq = str_remove(seq, ";"),
           seq = str_extract_all(seq, "[0-9]+"))
  bp$seq_length = 0
  bp$log_len = 0
  bp$unique_pitches = 0
  for (i in 1:(length(bp$seq))) {
    bp$seq_length[i] = length(bp$seq[[i]])
    bp$log_len[i] = log2(bp$seq_length[i])
    bp$unique_pitches[i] = length(unique(bp$seq[[i]]))
  }
  int <- list()
  for (i in 1:nrow(bp)) {
    int[i] = list(diff(as.numeric(bp$seq[[i]])))
  }
  bp$interval <- int
  
  bp$entropy <- 0
  for (i in 1:nrow(bp)) {
    bp$entropy[i] = entropy(bp$interval[[i]])
  }
  qualtrics_1st_batch_g$num_seq[qualtrics_1st_batch_g$BP_Name == file] <- nrow(bp)
  qualtrics_1st_batch_g$length[qualtrics_1st_batch_g$BP_Name == file] <- mean(bp$seq_length)
  qualtrics_1st_batch_g$log_len[qualtrics_1st_batch_g$BP_Name == file] <- mean(bp$log_len)
  qualtrics_1st_batch_g$pitches[qualtrics_1st_batch_g$BP_Name == file] <- mean(bp$unique_pitches)
  qualtrics_1st_batch_g$entropy[qualtrics_1st_batch_g$BP_Name == file] <- mean(bp$entropy)
}

# load Corinna's data
spss <- read_sav("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_generation/Fall_2022_analyses/BPsequencer_FULLSET.sav") %>%
  arrange(RecordedDate) %>%
  dplyr::select(`LANG#`, MAXMSP, FINGERS, EDIBLE, BRICK, CAN, CLIP, DTTtotal)
# may select other vars according to needs

# after making sure the MAXMSP names are correct for all
joined_spss_g <- left_join(qualtrics_1st_batch_g, spss, by = c("BP_Name" = "MAXMSP"), keep = FALSE)

# =============== Creativity per generator =================
creativity <- mean_rat %>%
  group_by(BP_Name) %>%
  summarize(mean_creativity_rating = mean(mean_rating, na.rm = TRUE),
            sd_creativity_rating = sd(mean_rating, na.rm = TRUE))

joined_creativity <- left_join(joined_spss_g, creativity) %>%
  dplyr::select(20:22, 24:30, 34, 212:224, 267:280) 

write_csv(joined_creativity, "~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Fall_2022_rating_analyses/joined_creativity_g_1.csv")

# ============= R ================
setwd("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Summer_Fall_2022_ratings")

qualtrics_1st_batch_r$num_rated <- 0
qualtrics_1st_batch_r$rating_avg <- 0
qualtrics_1st_batch_r$rating_sd <- 0

for (file in qualtrics_1st_batch_r$BP_Name) {
  bp <- read_csv(file, col_names = c("trial", "seq")) %>%
    mutate(seq = str_remove(seq, ";"))
  bp$rating <- sub(".*\\s", "", bp$seq)
  qualtrics_1st_batch_r$num_rated[qualtrics_1st_batch_r$BP_Name == file] <- nrow(bp) # number of sequences each rater rated
  qualtrics_1st_batch_r$rating_avg[qualtrics_1st_batch_r$BP_Name == file] <- mean(as.numeric(bp$rating), na.rm = TRUE) # mean rating given
  qualtrics_1st_batch_r$rating_sd[qualtrics_1st_batch_r$BP_Name == file] <- sd(as.numeric(bp$rating), na.rm = TRUE)
}

joined_spss_r <- left_join(qualtrics_1st_batch_r, spss, by = c("BP_Name" = "MAXMSP"), keep = FALSE) %>%
  dplyr::select(20:21, 24:30, 34, 212:224, 267:276)

write_csv(joined_spss_r, "~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Fall_2022_rating_analyses/joined_creativity_r_1.csv")