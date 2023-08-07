library(tidyverse)
library(haven)
library(psych)
library(corrplot)

setwd("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Spring_2023_rating_analyses")

rating <- read_csv("Spring_2023_all_ratings.txt", col_names = c("trial", "sequence"))

# separate sequence number from the rest based on first space
sep_seq_num <- rating %>%
  separate(sequence, into = c("seq_num", "sequence"), sep = " ", extra = "merge") %>%
  mutate(seq_num = as.numeric(seq_num)) 

sep_seq_num %>%
  group_by(seq_num) %>%
  count() %>%
  arrange(-seq_num)

# check those whose seq_num weren't recorded
# those repeated consecutively don't have their seq_num recorded on the 2nd occurrance
tbd <- sep_seq_num %>% # to be deleted
  filter(seq_num < 13) %>% 
  # if seq_num > 12 we are sure that the seq_num is recorded,
  # if not we are unsure whether seq_num recorded is actual seq_num or first note of sequence,
  # so we select all seq_num < 13 and figure out manually
  arrange(seq_num) %>%
  # delete the following ***needs to be changed once data is updated
  slice(1, 2, 31, 55, 60, 62, 73, 85, 108, 123, 124, 148, 150, 302)
# sequences 24, 36, 22, 113, 91, 72, 136, 34, 62, 99, 103 have been repeated consecutively
# delete those from data, keeping only 1st rating
del_rep <- anti_join(sep_seq_num, tbd) # by all vars

# extract rating from sequence by extracting the number before; (run only once per session)
del_rep$rating <- sub(".*\\s", "", del_rep$sequence)

# alternative: extract rating using excel
# =TRIM(RIGHT(SUBSTITUTE(C2," ",REPT(" ",LEN(C2))),LEN(C2)))
# https://www.extendoffice.com/excel/formulas/excel-extract-text-after-last-instance-of-character.html

seq_rat <- del_rep %>%
  mutate(rating = str_remove(rating, ";")) %>% # remove ;
  mutate(rating = as.numeric(rating))

# check if the frequency of selection was similar across sequences
seq_rat %>%
  ggplot(aes(seq_num)) +
  geom_bar() +
  labs(title = "Frequency of sequence selection", x = "Sequence number", y = "Count")

seq_rat %>%
  group_by(seq_num) %>%
  count() %>%
  arrange(seq_num)
# sequences 1-136 were more frequently selected, this is bc we mistakenly set maxpat to 
# select the first 136 sequences at first, and as soon as we realized that, we changed 
# maxpat so that it starts selecting the latter sequences
# seqs 660, 741 never got selected

seq_rat %>%
  filter(is.na(rating))
# no NA

# load sequence generations
sequences <- read_tsv("Fall_2022_sequences_bp.txt")

# convert text to list of numbers
joined <- left_join(seq_rat, sequences, by = "seq_num") %>%
  # convert seq chr to list of numbers
  mutate(seq_list = str_extract_all(seq, "[0-9]+")) %>%
  dplyr::select(-trial) %>%
  filter(!(seq_num %in% 310:333)) # delete the problematic sequences

# check if there are sequences with the same content but created by different people
sequences %>% 
  filter(sequences %>% 
           dplyr::select(seq) %>% 
           duplicated())
# 2 duplicates

# combine duplicated sequences into 1
mean_rat_dup <- joined %>%
  group_by(seq, seq_list, BP_Name) %>%
  summarize(mean_rating = mean(rating, na.rm = TRUE))

# fill in length & # of unique pitches for each sequence
mean_rat_dup$seq_length = 0
mean_rat_dup$unique_pitches = 0

for (i in 1:(length(mean_rat_dup$seq_list))) {
  mean_rat_dup$seq_length[i] = length(mean_rat_dup$seq_list[[i]])
  mean_rat_dup$unique_pitches[i] = length(unique(mean_rat_dup$seq_list[[i]]))
}

# create a new var containing log2(length) to make the outlier (len=80) less influential
mean_rat_dup$log_len = log2(mean_rat_dup$seq_length)

# create a new variable containing intervals in each sequence
int <- list()
for (i in 1:nrow(mean_rat_dup)) {
  int[i] = list(diff(as.numeric(mean_rat_dup$seq_list[[i]])))
}
mean_rat_dup$interval <- int

# instructions for how to calculate entropy 
# https://rstudio-pubs-static.s3.amazonaws.com/455435_30729e265f7a4d049400d03a18e218db.html
# table() - https://www.geeksforgeeks.org/create-a-tabular-representation-of-data-in-r-programming-table-function/

# compute Shannon entropy; input list of intervals, e.g. mean_rat_dup$int[[1]]
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

mean_rat_dup$entropy <- 0
for (i in 1:nrow(mean_rat_dup)) {
  mean_rat_dup$entropy[i] = entropy(mean_rat_dup$interval[[i]])
}

write_csv(mean_rat_dup, "~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Spring_2023_rating_analyses/mean_rat_2.csv")

# ========================= Qualtrics ===============================
qualtrics_23 <- read_csv("BP_Spring_2023_Qualtrics.csv") %>%
  slice(-1, -2) %>%
  # change to date format
  mutate(StartDate = ymd_hms(StartDate, tz = "America/Denver"),
         # extract date and make them numeric
         date = as.numeric(gsub("-", "", date(StartDate))),
         # label semester
         batch = "SP23",
         GeneRate = if_else(is.na(BP_Name), "NA", "R"),
         # add "rating" label if not already
         Lang_1st = str_to_title(Lang_1st), # convert 1st letter to uppercase
         Race = factor(Race, levels = c(1:8), 
                       labels = c("American Indian/Alaska Native", "Asian", "Black/African American", 
                                  "Hispanic/Latino", "Native Hawaiian/Pacific Islander", 
                                  "White", "2+ Races", "Prefer not to answer")),
         mus_train = if_else(Music == "No" | Music == "no" | Music == "none" | 
                               Music == "NO" | is.na(Music), "No", "Yes"),
         .before = Name)

# examine those whose BP file type is wrong
qualtrics_23 %>%
  filter(BP_Type != "text/plain" | is.na(BP_Name) | !grepl("_rating.txt$", BP_Name)) %>%
  dplyr::select(StartDate, Name, BP_Name, GeneRate)

# JOST: wrong survey (supposed to fill out the EEG survey)
# LHEL: nosebleed (finished up to BMRQ8); delete for BMRQ & GMSI analyses only
# JLAO: music playback for SRQ didn't work, finished up to BP rating; delete for SRQ&BMRQ&GMSI analyses
# SKUM: duplicate, delete

qualtrics_new_23 <- qualtrics_23 %>%
  mutate(BP_Name = case_when(BP_Name == "230214_XLI" ~ "230214_XLI_rating.txt",
                             BP_Name == "230221_MGON_ratings.txt" ~ "230221_MGON_rating.txt",
                             BP_Name == "230317_YHE" ~ "230317_YHE_rating.txt",
                             BP_Name == "230328JBAR" ~ "230328_JBAR_rating.txt",
                             BP_Name == "280323_ALEE.txt" ~ "280323_ALEE_rating.txt",
                             BP_Name == "230410_HHAL_rating (1).txt" ~ "230410_HHAL_rating.txt",
                             BP_Name == "230301FHE_rating.txt" ~ "230301_MRAM_rating.txt",
                             TRUE ~ BP_Name))

# ============================== DAT ================================
DAT <- qualtrics_new_23 %>%
  dplyr::select(BP_Name, starts_with("DAT"), -date) %>%
  rename(id = BP_Name, word1 = DAT_4, word2 = DAT_5, word3 = DAT_6, word4 = DAT_7, 
         word5 = DAT_8, word6 = DAT_9, word7 = DAT_10, word8 = DAT_11, word9 = DAT_12, word10 = DAT_13)

write_tsv(DAT, "DAT.tsv")
# upload to https://www.datcreativity.com/analyse and download the scored data

dat_scored <- read_tsv("dat-scored.tsv")[c(1, 12)]

qualtrics_w_dat <- left_join(qualtrics_new_23, dat_scored, by = c("BP_Name" = "id"))

# ========================================================
qualtrics_2nd_batch_r <- qualtrics_w_dat %>%
  filter(batch == "SP23", GeneRate == "R") %>%
  rename(BMRQ_MS = SC0, BMRQ_EE = SC1, BMRQ_MR = SC2, BMRQ_SM = SC3, BMRQ_SR = SC4, BMRQ_AM = SC5, eBMRQ = SC6, 
         MSI_AE = SC7, MSI_PA = SC8, MSI_MT = SC9, MSI_SA = SC10, MSI_EM = SC11, MSI_GS = SC12) %>%
  mutate(BP_Name = str_replace(BP_Name, "^(\\d{6})([^_])", "\\1_\\2"),
         BP_Name = ifelse(BP_Name == "230227_LUAN_rating.txt", "230227_LDEN_rating.txt", BP_Name),
         BP_Name = ifelse(BP_Name == "280323_ALEE_rating.txt", "230328_ALEE_rating.txt", BP_Name)) %>%
  filter(BP_Name != "022723_JOST_rating.txt") # he should fill out the EEG survey; delete

write_csv(qualtrics_2nd_batch_r, "~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Spring_2023_rating_analyses/qualtrics_r_2.csv")

# =========== Generation =============

qualtrics_2nd_batch_g <- read_csv("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Fall_2022_rating_analyses/qualtrics_g_2.csv") 

# =============== Manipulate Qualtrics data ================= 
setwd("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_generation/Fall_2022_sequences/Fall_2022_sequences")

# 9 people's data included previous participants' data, manually deleted in .txt before running

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

# manually change typo in names
qualtrics_2nd_batch_g$BP_Name[qualtrics_2nd_batch_g$BP_Name == "221013_NMON"] = "221013_NMON.txt"
qualtrics_2nd_batch_g$BP_Name[qualtrics_2nd_batch_g$BP_Name == "221017__ABUR.txt"] = "221017_ABUR.txt"
qualtrics_2nd_batch_g$BP_Name[qualtrics_2nd_batch_g$BP_Name == "221209_SKHA"] = "221209_SKHA.txt"

qualtrics_2nd_batch_g$num_seq <- 0
qualtrics_2nd_batch_g$length <- 0
qualtrics_2nd_batch_g$log_len <- 0
qualtrics_2nd_batch_g$pitches <- 0
qualtrics_2nd_batch_g$entropy <- 0

# create a loop that ouputs each participant's data
for (file in qualtrics_2nd_batch_g$BP_Name) {
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
  qualtrics_2nd_batch_g$num_seq[qualtrics_2nd_batch_g$BP_Name == file] <- nrow(bp)
  qualtrics_2nd_batch_g$length[qualtrics_2nd_batch_g$BP_Name == file] <- mean(bp$seq_length)
  qualtrics_2nd_batch_g$log_len[qualtrics_2nd_batch_g$BP_Name == file] <- mean(bp$log_len)
  qualtrics_2nd_batch_g$pitches[qualtrics_2nd_batch_g$BP_Name == file] <- mean(bp$unique_pitches)
  qualtrics_2nd_batch_g$entropy[qualtrics_2nd_batch_g$BP_Name == file] <- mean(bp$entropy)
}

# load Corinna's data (directory is /Users/eva at home but /Users/xi.wu at NEU)
spss <- read_sav("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_generation/Fall_2022_analyses/BPsequencer_FULLSET.sav") %>%
  arrange(RecordedDate) %>%
  dplyr::select(`LANG#`, MAXMSP, FINGERS, EDIBLE, BRICK, CAN, CLIP, DTTtotal)
# may select other vars according to needs

# after making sure the MAXMSP names are correct for all
joined_spss_g <- left_join(qualtrics_2nd_batch_g, spss, by = c("BP_Name" = "MAXMSP"), keep = FALSE)
# 1 less row b/c VPAT doesn't have Qualtrics survey recorded

# =============== Creativity per generator =================
creativity <- mean_rat_dup %>%
  group_by(BP_Name) %>%
  summarize(mean_creativity_rating = mean(mean_rating, na.rm = TRUE),
            sd_creativity_rating = sd(mean_rating, na.rm = TRUE))

joined_creativity_g <- left_join(joined_spss_g, creativity, by = "BP_Name") %>%
  dplyr::select(20:22, 24:30, 34, 212:224, 267:280)

write_csv(joined_creativity_g, "~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Spring_2023_rating_analyses/joined_creativity_g_2.csv")

# ============= R ================
setwd("~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Spring_2023_ratings")

qualtrics_2nd_batch_r$num_rated <- 0
qualtrics_2nd_batch_r$rating_avg <- 0
qualtrics_2nd_batch_r$rating_sd <- 0

for (file in qualtrics_2nd_batch_r$BP_Name) {
  bp <- read_csv(file, col_names = c("trial", "seq")) %>%
    mutate(seq = str_remove(seq, ";"))
  bp$rating <- sub(".*\\s", "", bp$seq)
  qualtrics_2nd_batch_r$num_rated[qualtrics_2nd_batch_r$BP_Name == file] <- nrow(bp) # number of sequences each rater rated
  qualtrics_2nd_batch_r$rating_avg[qualtrics_2nd_batch_r$BP_Name == file] <- mean(as.numeric(bp$rating), na.rm = TRUE) # mean rating given
  qualtrics_2nd_batch_r$rating_sd[qualtrics_2nd_batch_r$BP_Name == file] <- sd(as.numeric(bp$rating), na.rm = TRUE)
}

joined_creativity_r <- qualtrics_2nd_batch_r %>%
  dplyr::select(20:21, 24:30, 34, 168:184)

write_csv(joined_creativity_r, "~/Dropbox/MINDLabWes/projects/CreativityImprovisation/BP_SL_NES/BP_sequences/Sequence_rating/rating_data/Spring_2023_rating_analyses/joined_creativity_r_2.csv")