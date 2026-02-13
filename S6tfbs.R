cat("\014")
rm(list = ls())

pacman::p_load(tidyr, dplyr, ggplot2, ggmsa, ggpubr, cowplot)

tfbs_df <- read.delim("fimo.txt") # reads in transcription factor data
snp_df <- read.csv(file = "Book1.csv") # reads in SNP data

snp_df$Family <- as.character(snp_df$Family)

typeof(tfbs_df$start)
tfbs_df <- tfbs_df %>% 
  mutate(fasta = paste0(">",pattern.name,"\n",matched.sequence,"\n"))

# Identify where the TFBS are
tfbs_longer_df <-pivot_longer(tfbs_df, cols = c("start","stop"), values_to = "region")

# Filter for only G19
tfbs_longer_df_g19 <- tfbs_longer_df %>%
  filter(sequence.name == "G19")

# Plot
tfbs_plot_g19 <- ggplot(tfbs_longer_df_g19, aes(x = region-1852, y = Family)) +
  geom_line(aes(group = interaction(pattern.name, region)), color = "grey60", size = 0.6) +
  geom_point(aes(shape = strand, color = score), size = 3, position = position_dodge(width = 0.5)) +
  # Only one facet
  facet_grid(rows = vars(sequence.name)) +
  theme_linedraw(base_size = 16) +
  labs(
    y = "TF Family",
    x = "Region (bp, relative to TSS)",
    caption = "0 = TSS"
  ) +
  scale_x_continuous(n.breaks = 18, limits = c(-1853,0)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    strip.text.y = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.caption = element_text(size = 13)
  ) +
  scale_color_gradient(low = "dodgerblue", high = "firebrick", name = "Score")

print(tfbs_plot_g19)

# Send it
ggsave("tfbs_plot_g19.eps",plot = tfbs_plot_g19, width = 15, height = 10, units = "in")