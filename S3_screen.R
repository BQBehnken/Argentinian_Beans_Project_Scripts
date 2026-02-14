
# Author: TC


# Input:  qpcrdata.csv 
# Output: TonioqPCR.eps

cat("\014")
rm(list = ls())

# --- Load packages 

pacman::p_load(dplyr, ggplot2, multcompView)

# --- Import data 

df <- read.csv("qpcrdata.csv")

head(df)


# ANOVA + Tukey HSD


anova <- aov(DDCQ ~ Accession, df)
summary(anova)

tukey <- TukeyHSD(anova)

cld <- multcompLetters4(anova, tukey)

# Build summary table with means, 75th percentile, and compact letter display
table_of_anova <- df %>%
  dplyr::group_by(Accession) %>%
  summarize(
    mean  = mean(DDCQ, na.rm = TRUE),
    quant = quantile(DDCQ, probs = 0.75, na.rm = TRUE)
  ) %>%
  arrange(desc(mean))

cld <- as.data.frame.list(cld$Accession)
table_of_anova$cld <- cld$Letters


# Figure: INR expression across accessions

print(df)

# Set accession display order
df$Accession <- factor(
  df$Accession,
  levels = c("G19833", "PI 638852", "W6 17024", "PI 638856",
             "PI 661803", "PI 638858", "W6 17019"),
  ordered = TRUE
)

plot <- ggplot(df, aes(x = Accession, y = DDCQ)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(
    shape    = 16,
    alpha    = 1,
    size     = 2,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(
    data  = table_of_anova,
    aes(x = Accession, y = quant, label = cld),
    size  = 3,
    vjust = -1,
    hjust = -1
  ) +
  labs(
    x = expression(paste(italic("Phaseolus vulgaris"), " Accession")),
    y = "Relative Expression of INR"
  ) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1)
  )

print(plot)

ggsave(file = "TonioqPCR.eps", plot = plot, width = 4, height = 3, units = "in")


# Mean fold-difference: responsive vs. unresponsive accessions


responsive   <- c("PI 638852", "PI 638856", "W6 17024")
unresponsive <- c("PI 638858", "PI 661803", "W6 17019")

mean_responsive <- df %>%
  filter(Accession %in% responsive) %>%
  summarize(mean_DDCQ = mean(DDCQ, na.rm = TRUE)) %>%
  pull(mean_DDCQ)

mean_unresponsive <- df %>%
  filter(Accession %in% unresponsive) %>%
  summarize(mean_DDCQ = mean(DDCQ, na.rm = TRUE)) %>%
  pull(mean_DDCQ)

mean_for_paper <- mean_responsive / mean_unresponsive

percent_upregulation <- ((mean_responsive - mean_unresponsive) / abs(mean_unresponsive)) * 100

print(percent_upregulation)
