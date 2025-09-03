
cat("\014")
rm(list = ls())

pacman::p_load(dplyr, tidyr, stringr, ggplot2, ggpubr, ggbreak, scales, cowplot, multcompView)

# ────────────────────────────────────────────────────────────────────────────
# STEP 1:  Read in all six luciferase CSVs and bind into one data frame
# ────────────────────────────────────────────────────────────────────────────

luc_1_df <- read.csv("luc_0626_filtered_max_lum.csv")
luc_2_df <- read.csv("luc_0710_filtered_max_lum.csv")
luc_3_df <- read.csv("luc_0724_1_filtered_max_lum.csv")
luc_4_df <- read.csv("luc_0724_filtered_max_lum.csv")
luc_5_df <- read.csv("luc_0731_filtered_max_lum.csv")
luc_6_df <- read.csv("luc_0814_filtered_max_lum.csv")

luc_all_df <- bind_rows(
  luc_1_df, luc_2_df, luc_3_df,
  luc_4_df, luc_5_df, luc_6_df
)

# Quick sanity checks since we're working with many data points
stopifnot("luminescence" %in% colnames(luc_all_df)) 
stopifnot("promoter"      %in% colnames(luc_all_df))
stopifnot("date_tested"   %in% colnames(luc_all_df))
stopifnot("biological_replicate" %in% colnames(luc_all_df))

head(luc_all_df)

# ────────────────────────────────────────────────────────────────────────────
# STEP 2:  Keep only the four actual promoter::luciferse fusions (249, 253, 255, 256)
# ────────────────────────────────────────────────────────────────────────────
luc_four <- luc_all_df %>%
  filter(promoter %in% c("249", "253", "255", "256")) %>% 
  filter(biological_replicate != "incubation of leaf disks") # This assay was done differently from all the others; biological replicates are not clear. 

# Confirm we dropped everything else
unique(luc_four$promoter)  # This should show only 249, 253, 255, 256 and note, 256 promoter is identical to 260 promoter

normalize_255 <- luc_four %>%  # gets a mean and SD from 255 for days/biorep
  filter(promoter == "255") %>% # low control to normalize to
  group_by(date_tested, biological_replicate, promoter) %>% 
  dplyr::summarise(mean = mean(luminescence, na.rm = T), sd = sd(luminescence, na.rm = T), .groups = 'keep') 

luc_long_annotated_filtered_normalized_df <- left_join(luc_four, normalize_255, by = c("date_tested","biological_replicate")) # adds the mean and sd back to the data frame



############ Z Score Calculation to bring order to the data and see where all the other promoter fusions rank compared to PI 638858 (255)

luc_z_score <- luc_long_annotated_filtered_normalized_df %>% # makes z-scores based on 255 mean and sd for each experiment day
  group_by(date_tested, biological_replicate, promoter.x, well) %>%
  mutate(luminescence_z_score = (luminescence - mean) / sd) %>%
  ungroup() %>%
  dplyr::select(well, luminescence, promoter = promoter.x, biological_replicate, mean, sd, luminescence_z_score, date_tested) # I all of a sudden started getting an error with the select function so I had to specify dplyr

colnames(luc_z_score)

# Lets add some factors
luc_z_score$promoter <- as.factor(luc_z_score$promoter)

# now order the factors so col is first
luc_z_score$promoter = ordered(x=luc_z_score$promoter, 
                                      c('249','256','260','253','255')) # adds factor order

# use to check levels of your factors
levels(luc_z_score$promoter)

# write.csv(luc_all_df, file =  "luciferasealldata.csv")
# write.csv(luc_z_score, file =  "luciferaseallZscore.csv")

andean_plants <- luc_z_score %>% 
  filter(!promoter %in% c("35_S","EV")) %>%  # removes non promoter groups
  filter(!promoter %in% c("G19"))
andean_plants$promoter[andean_plants$promoter == "256"] <- "260" # Since 256 = 260

andean_plants <- andean_plants %>% 
  mutate(temp_promoter = promoter,promoter = case_when(temp_promoter == "249" ~ "PI_638852",
                                                       temp_promoter == "260" ~ "PI_638856",
                                                       temp_promoter == "253" ~ "PI_661803",
                                                       temp_promoter == "255" ~ "PI_638858"))

# Lets add some factors
andean_plants$promoter <- as.factor(andean_plants$promoter)

# write.csv(andean_plants, file = "andeanZscore.csv")

# now order the factors so col is first
andean_plants$promoter = ordered(x=andean_plants$promoter, 
                               c('PI_638852','PI_638856','PI_661803','PI_638858')) # adds factor order



pacman::p_load(ggpubr, rstatix, car)

# Fit model first so we test what the model assumes
m <- aov(luminescence_z_score ~ promoter, data = andean_plants)

# Normality of residuals (plus visual check)
shapiro.test(residuals(m)) # W = 0.88379, p-value = 5.732e-15 ;strong evidence against normality
ggqqplot(residuals(m)) # plot looks whack

# Homogeneity of variances (median-centered Levene = Brown–Forsythe)
leveneTest(luminescence_z_score ~ promoter, andean_plants) # p = 9.34e-14; strong evidence of unequal variances


########## Let us once again proceed to a WELCH's ANOVA

# 1) Welch's ANOVA (no equal-variances assumption)
welch <- oneway.test(luminescence_z_score ~ promoter, data = andean_plants, var.equal = FALSE)
print(welch) # p-value = 5.898e-13 Great - let's do a post-hoc test

# 2) Games–Howell post hoc (pairwise)
gh1 <- games_howell_test(luminescence_z_score ~ promoter, data = andean_plants)
print(gh1)

# 3) Build named p-value vector
pairs1 <- gh1 %>%
  transmute(
    Comparison = str_replace_all(paste(group1, group2, sep = "-"), "\\s+", ""),
    p.adj = p.adj
  ) %>%
  filter(!is.na(p.adj)) %>%
  distinct(Comparison, .keep_all = TRUE)

pvec1 <- setNames(pairs1$p.adj, pairs1$Comparison)

# 4) Compact letter display
letters1 <- multcompLetters(pvec1, threshold = 0.05)$Letters
letters1_tbl <- data.frame(Group = names(letters1), Letters = letters1, row.names = NULL)

letters1_tbl

# Now turn the letters tibble into a data frame
letters1_tbl <- data.frame(
  promoter = names(letters1),
  Letters = letters1,
  row.names = NULL
)

aov_tbl <- andean_plants %>% 
  dplyr::group_by(promoter) %>% 
  summarize(mean = mean(luminescence_z_score, na.rm = T), quant = quantile(luminescence_z_score, probs = 0.75, na.rm = T)) %>% 
  arrange(desc(mean))

# join into your summary table of means/quantiles
annot_tbl <- aov_tbl %>%
  left_join(letters1_tbl, by = "promoter")

annot_tbl

library(readr)
write_csv(andean_plants, "andean_plants.csv")


fig2c <- ggplot(data = andean_plants,
                             mapping = aes(x = promoter, y = luminescence_z_score)) +
  geom_boxplot(outlier.alpha = 0, width = 0.6) +
  geom_jitter(aes(color = paste(biological_replicate, date_tested)),
              size = 1.5, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) +
  geom_text(data = annot_tbl, aes(x = promoter, y = quant, label = Letters), 
            size = 3, hjust = -0.5, vjust = -0.5) +
  #scale_y_continuous(labels = label_comma(), limits = c(-3, 15)) +
  labs(x = "pINR by Accession",
       y = "RLU to PI 638858 mean") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none",
        # panel.grid.minor = element_blank(),
        # panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"))

print(fig2c)

ggsave("fig2c.eps", height = 3,width =2.5, units = "in")


# I want to be able to make this claim "Promoters from the two In11-responsive lines, pINR PI 638852 and pINR PI 638856, drove XX % higher luciferase activity than the pINR PI 638858 from an unresponsive line, while a second promoter pINR PI 661803 was similarly low in activity  (Fig 2C)."


# Compute fold to PI 638858 *within each day/biorep* using the `mean` column that is already joined in
pct_tbl <- luc_long_annotated_filtered_normalized_df %>%
  # Use the original promoter from the left side of the join
  transmute(
    luminescence,
    ctrl_mean = mean,                 # per-day/per-biorep PI 638858 mean
    promoter  = promoter.x
  ) %>%
  # Recode 256 -> 260, then map to accessions used in your figure
  mutate(
    promoter  = if_else(promoter == "256", "260", promoter),
    accession = case_when(
      promoter == "249" ~ "PI 638852",
      promoter == "260" ~ "PI 638856",
      promoter == "253" ~ "PI 661803",
      promoter == "255" ~ "PI 638858",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(accession)) %>%
  mutate(fold_vs_638858 = luminescence / ctrl_mean)

# Summarize the mean fold for each accession
fold_by_acc <- pct_tbl %>%
  group_by(accession) %>%
  summarise(mean_fold = mean(fold_vs_638858, na.rm = TRUE), .groups = "drop")

# Control mean fold (should be ~1, but pull it explicitly)
ctrl_fold <- fold_by_acc %>% filter(accession == "PI 638858") %>% pull(mean_fold)

# % change for the two In11‑responsive accessions COMBINED (pool all wells, then average)
pct_two_responsive <- pct_tbl %>%
  filter(accession %in% c("PI 638852","PI 638856")) %>%
  summarise(mean_fold = mean(fold_vs_638858, na.rm = TRUE)) %>%
  mutate(pct_vs_638858 = (mean_fold / ctrl_fold - 1) * 100) %>%
  pull(pct_vs_638858)

# % for PI 661803 relative to PI 638858
pct_661803 <- fold_by_acc %>%
  filter(accession == "PI 661803") %>%
  mutate(pct_vs_638858 = (mean_fold / ctrl_fold - 1) * 100) %>%
  pull(pct_vs_638858)


