cat("\014")
rm(list = ls())

pacman::p_load(dplyr, ggplot2, readr, car, emmeans)

# Testing base INR expression in Backcross

df <- readr::read_csv("Fig3a_expression2.csv", show_col_types = FALSE) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("CR4", "CR5")),
    BioRep   = as.integer(BioRep)
  )

# Decoder
# CR4 = G19833-like
# CR5 = PI 638858-like

colnames(df)

head(df)




# ==== Part A: "Normalcy" checks & omnibus test across genotypes ===============
# We compare INR expression magnitude across genotypes: response = DCQ (log2 scale)

# Testing the log data
shapiro.test(df$DCQ) # W = 0.87005, p-value = 0.1001

# 1) ANOVA residual normality for linear data
aov_fit <- aov(DCQFC ~ Genotype, data = df)
shapiro_res <- shapiro.test(residuals(aov_fit))
print(shapiro_res) # W = 0.9381, p-value = 0.5321 # No evidence against normality

# 2) Homogeneity of variance
lev <- car::leveneTest(DCQFC ~ Genotype, data = df)
lev_p <- lev[1, "Pr(>F)"]
print(lev) # 0.06422 # No evidence against group variance

head(df)

# Simple two-sided Student’s t-test since data is "normal"
t.test(DCQ ~ Genotype, data = df, var.equal = TRUE) # t = -6.9261, df = 8, p-value = 0.0001213



min_mean <- df %>%
  group_by(Genotype) %>%
  summarise(m = mean(DCQFC, na.rm = TRUE)) %>%
  summarise(min_val = min(m)) %>%
  pull(min_val)

df <- df %>%
  mutate(DCQFC_rel = DCQFC / min_mean)

y_max = 18
y_min = 0

plot2 <- ggplot(df, aes(x = Genotype, y = DCQFC_rel, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 1) +
  geom_jitter(size = 3, position = position_dodge(width = 0.5)) +
  labs(title = "",
       x = "Genotype", y = "INR Relative Expression") +
  theme_linedraw() +
  scale_y_continuous(limits = c(y_min, y_max)) +
  theme(legend.position = "none")

print(plot2)


ggsave("plot2.eps", width = 2, height = 3)

# I chose to annotate the comparison bar with statistical significance in Illustrator rather than code it in, as I am far more experience in Adobe software than I am in R. I just scaled y to leave some headroom at the top for that. 




