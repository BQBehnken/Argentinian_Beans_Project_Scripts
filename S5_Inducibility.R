
#
# Hypothesis: INR expression in PI 255 < PI 260 under undamaged conditions,
# and wounding or peptide treatment induces INR expression.
#
# Result: No statistically significant induction of INR by wounding or
# peptide treatment was detected.
#
# Input:  S5_Data.csv
# Output: plot_260.eps 
# =============================================================================

cat("\014")
rm(list = ls())

# --- Load packages

pacman::p_load(
  dplyr, ggplot2, emmeans, multcomp, multcompView,
  ggpubr, car, dunn.test
)

# --- Import data
# Set working directory: Session -> Set Working Directory -> Source File Location

df <- read.csv("S5_Data.csv")

colnames(df)
head(df)

# =============================================================================
# Statistical tests
# =============================================================================
# Performed for due diligence; no significant differences were detected.
# The plot is reported without statistical annotations.

str(df)

# Exclude the untreated PI 260 calibrator from statistical comparisons
df_stat <- df %>% dplyr::filter(!Target %in% "UnTrt260Avg")

# --- Normality: Shapiro-Wilk test 
# W = 0.9889, p = 0.9932 (data are normally distributed)
shapiro.test(df_stat$X260DDCQFC)
ggqqplot(df_stat$X260DDCQFC)

# --- Homogeneity of variance: Levene's test
# F = 4.41e+30, p = 2.2e-16 (variances are NOT homogeneous)
df_levene <- leveneTest(X260DDCQFC ~ Target, data = df_stat)
print(df_levene)

# --- Kruskal-Wallis test (nonparametric, given unequal variances)
# chi-squared = 6.0833, df = 4, p = 0.193 (not significant)
kruskal_test <- kruskal.test(X260DDCQFC ~ Target, data = df_stat)
print(kruskal_test)

# --- Dunn's post hoc test (Bonferroni correction)
dunn_260 <- dunn.test(
  x      = df_stat$X260DDCQFC,
  g      = df_stat$Target,
  method = "bonferroni"
)
print(dunn_260)

# --- Compact letter display (for reference only; not used in figure)
rank_260 <- lm(X260DDCQFC ~ Target, data = df_stat)
summary(rank_260)

emm_260 <- emmeans(rank_260, ~ Target)
cld_260 <- cld(emm_260, adjust = "bonferroni")
print(cld_260)

# =============================================================================
# Figure: INR expression (supplemental)
# =============================================================================

# Summarize mean and standard error per treatment
summary <- df %>%
  group_by(Target) %>%
  summarize(
    mean_DDCQ = mean(X260DDCQFC),
    se_DDCQ   = sd(X260DDCQFC) / sqrt(n())
  )

# Set treatment display order
summary$Target <- factor(
  summary$Target,
  levels = c("UnTrt260Avg", "UnTrt255Avg", "W6HR", "W24HR", "Pep6HR", "Pep24HR")
)

# Bar plot with error bars and individual data points
plot_260 <- ggplot() +
  geom_col(
    data = summary,
    aes(x = Target, y = mean_DDCQ, fill = Target)
  ) +
  geom_errorbar(
    data = summary,
    aes(x = Target, ymin = mean_DDCQ - se_DDCQ, ymax = mean_DDCQ + se_DDCQ),
    width = 0.2
  ) +
  geom_jitter(
    data     = df,
    aes(x = Target, y = X260DDCQFC, fill = Target),
    shape    = 16,
    alpha    = 1,
    size     = 2,
    position = position_dodge(width = 0.5)
  ) +
  theme_linedraw() +
  labs(
    title = "INR Expression in 255 Against 260",
    x     = "Treatment",
    y     = expression("Fold Change ("*Delta*Delta*"Cq FC)")
  ) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none"
  )

print(plot_260)

ggsave("plot_260.eps", width = 4, height = 3, units = "in")