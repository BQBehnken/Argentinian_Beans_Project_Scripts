cat("\014")
rm(list = ls())

pacman::p_load(pacman, dplyr, ggplot2, emmeans, rio, tidyr, multcompView, BiocManager, car, readr, rstatix)

# Testing the hypothesis that in low expression P. vulgaris, marker gene like MYB will also be less expressed than in high expression plant.

df <- readr::read_csv("Fig3a_expression1.csv", show_col_types = FALSE) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("G19833", "PI638858", "C4", "C5")),
    BioRep   = as.integer(BioRep)
  )

# Decoder
# C4 = G19833-like
# C5 = PI 638858-like

colnames(df)

head(df)

# ==== Part A: "Normalcy" checks ===============
# We compare induction magnitude across genotypes: response = DCQ (Delta Cq values)
# 

# testing the log data
shapiro.test(df$DCQ) # W = 0.90312, p-value = 0.06512 # no evidence against normality

# 1) Testing ANOVA residual normality on linear data
aov_fit <- aov(DDCQFC ~ Genotype, data = df)
shapiro_res <- shapiro.test(residuals(aov_fit))
print(shapiro_res) # W = 0.92378, p-value = 0.1508 # No evidence against normality

# 2) Homogeneity of variance
lev <- car::leveneTest(DDCQFC ~ Genotype, data = df)
lev_p <- lev[1, "Pr(>F)"]
print(lev) # 0.1191 no evidence against group variations being equal



## ANOVA ##
aov <- aov(DDCQFC ~ Genotype, data = df)

summary(aov) # p =  2.9e-08; at least 1 group is different; proceed with Tukey

tukey <- TukeyHSD(aov)

# View the results
print(tukey)

# Visualizing Tukey Results in boxplot: set up model
pairwise <- lm(DDCQFC ~ Genotype, data = df)

# get (adjusted) weight means per group
emm <- emmeans(object = pairwise,
                       specs = "Genotype")

# add letters to each mean
cld <- cld(object = emm,
                       adjust = "Tukey",
                       Letters = letters,
                       alpha = 0.05)
# check output
cld

table_of_anova <- df %>% 
  dplyr::group_by(Genotype) %>% 
  summarize(mean = mean(DCQFC, na.rm = T), quant = quantile(DCQFC, probs = 0.75, na.rm = T)) %>% 
  arrange(desc(mean))

table_of_anova$cld <- cld$.group

# show output
table_of_anova




# y_max = 18
# y_min = 0

plot <- ggplot(df, aes(x = Genotype, y = DDCQFC, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 1) +
  geom_jitter(size = 2, position = position_dodge(width = 0.5)) +
  labs(title = "",
       x = "Genotype", y = "INR Relative Expression") +
  theme_linedraw() +
  # scale_y_continuous(limits = c(y_min, y_max)) +
  theme(legend.position = "none") + 
  geom_text(data = table_of_anova,
            aes(x = Genotype, y = quant, label = cld),
            vjust = -1, hjust = -1, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))

print(plot)


ggsave("plot.eps", width = 2, height = 3)


