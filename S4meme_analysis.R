# P. vulgaris with the following background nucleotide frequencies:
#   A = 0.353, C = 0.147, G = 0.147, T = 0.353
#
# Logos are exported as EPS files in three styles per motif:
#   1. Base R ICM logo via universalmotif::view_motifs()
#   2. ggseqlogo bits logo (information content)
#   3. ggseqlogo probability logo
#
# A final Okabe-Ito colour-blind-friendly version (bits) is also saved.
#
# Motifs included (in order):
#   1. MYB        — Phvul.005G047400 (MP00350)
#   2. MIKC-MADS  — Phvul.006G169600 (MP00609)
#   3. Nin-like   — Phvul.007G071900 (MP00449)
#   4. bHLH       — Phvul.001G168700 (MP00074)
#   5. bZIP       — Phvul.007G229600 (MP00294)

cat("\014")
rm(list = ls())

# --- Load packages 

pacman::p_load(
  universalmotif,
  Biostrings,
  GenomicRanges,
  rtracklayer,
  ggseqlogo,
  ggplot2
)

# --- Define Okabe-Ito colour scheme (colour-blind friendly) 

oi <- c(
  A = "#0072B2",
  C = "#009E73",
  G = "#E69F00",
  T = "#D55E00"
)

oi_scheme <- make_col_scheme(chars = names(oi), cols = unname(oi))

# --- Shared background frequencies
# From file: ../../promoter_background/Pvu.bg

bg_freq <- c(A = 0.353, C = 0.147, G = 0.147, T = 0.353)

# --- Helper function 
# Builds a universalmotif object from a tab-delimited PPM string, sets
# background frequencies, and exports all four logo formats (base R ICM,
# ggseqlogo bits, ggseqlogo probability, and Okabe-Ito bits).

make_and_save_logos <- function(ppm_txt, motif_name, tf_family, nsites,
                                file_prefix) {
  
  # Parse PPM text into matrix (positions x bases)
  ppm <- as.matrix(read.table(text = ppm_txt))
  colnames(ppm) <- c("A", "C", "G", "T")
  ppm <- t(ppm)                        # transpose to bases x positions
  rownames(ppm) <- c("A", "C", "G", "T")
  
  # Create universalmotif object
  mot <- create_motif(
    ppm,
    name     = motif_name,
    nsites   = nsites,
    type     = "PPM",
    strand   = "+-",
    alphabet = "DNA"
  )
  
  # Set background frequencies
  if ("bkgr" %in% slotNames(mot)) {
    mot@bkgr <- bg_freq
  } else if ("background" %in% slotNames(mot)) {
    mot@background <- bg_freq
  }
  
  # --- 1. Base R ICM logo (EPS) 
  postscript(
    file       = paste0(file_prefix, "_logo.eps"),
    width      = 6,
    height     = 3,
    horizontal = FALSE,
    paper      = "special",
    onefile    = FALSE
  )
  view_motifs(mot, use.type = "ICM")
  dev.off()
  
  # --- 2-3. ggseqlogo bits and probability logos (EPS)
  mot_ppm <- convert_type(mot, "PPM")
  mot_icm <- convert_type(mot, "ICM")
  
  p_prob <- ggseqlogo(mot_ppm@motif, method = "prob") +
    ggtitle(paste0(tf_family, " ", motif_name, " (probabilities)"))
  
  p_bits <- ggseqlogo(mot_icm@motif, method = "bits") +
    ggtitle(paste0(tf_family, " ", motif_name, " (bits)"))
  
  print(p_bits)
  
  ggsave(
    filename = paste0(file_prefix, "_bits.eps"),
    plot     = p_bits,
    width    = 6, height = 3, dpi = 600
  )
  
  ggsave(
    filename = paste0(file_prefix, "_prob.eps"),
    plot     = p_prob,
    width    = 6, height = 3, dpi = 600
  )
  
  # --- 4. Okabe-Ito colour-blind friendly bits logo (EPS) 
  p_oi <- ggseqlogo(mot_icm@motif,
                    method     = "bits",
                    col_scheme = oi_scheme) +
    labs(
      title = paste0(tf_family, " ", motif_name, " motif (bits)"),
      x     = "Position",
      y     = "Bits"
    ) +
    theme_classic(base_size = 14)
  
  print(p_oi)
  
  ggsave(
    filename = paste0(file_prefix, "_okabe-minimal.eps"),
    plot     = p_oi,
    width    = 6, height = 3, dpi = 600
  )
}


# =============================================================================
# 1. MYB — Phvul.005G047400 (MP00350)
# =============================================================================
# PlantTFDB: http://planttfdb.cbi.pku.edu.cn/tf.php?sp=Pvu&did=Phvul.005G047400.1#bind_motif
# letter-probability matrix: alength= 4 w= 21 nsites= 599 E= 3.0e-932

myb_txt <- "0.273790	  0.150250	  0.303840	  0.272120	
  0.260434	  0.198664	  0.253756	  0.287145	
  0.322204	  0.133556	  0.225376	  0.318865	
  0.280467	  0.138564	  0.332220	  0.248748	
  0.283806	  0.160267	  0.258765	  0.297162	
  0.283806	  0.128548	  0.295492	  0.292154	
  0.287145	  0.155259	  0.307179	  0.250417	
  0.268781	  0.141903	  0.227045	  0.362270	
  0.368948	  0.101836	  0.275459	  0.253756	
  0.358932	  0.080134	  0.338898	  0.222037	
  0.191987	  0.245409	  0.116861	  0.445743	
  0.100167	  0.000000	  0.888147	  0.011686	
  0.001669	  0.000000	  0.400668	  0.597663	
  0.000000	  0.000000	  0.000000	  1.000000	
  0.634391	  0.000000	  0.041736	  0.323873	
  0.000000	  0.000000	  1.000000	  0.000000	
  0.003339	  0.000000	  0.904841	  0.091820	
  0.000000	  0.000000	  0.000000	  1.000000	
  0.444073	  0.000000	  0.549249	  0.006678	
  0.238731	  0.198664	  0.480801	  0.081803	
  0.404007	  0.118531	  0.377295	  0.100167"

make_and_save_logos(
  ppm_txt     = myb_txt,
  motif_name  = "Phvul.005G047400",
  tf_family   = "MYB",
  nsites      = 599,
  file_prefix = "Phvul_myb"
)


# =============================================================================
# 2. MIKC-MADS — Phvul.006G169600 (MP00609)
# =============================================================================
# PlantTFDB: http://planttfdb.cbi.pku.edu.cn/tf.php?sp=Pvu&did=Phvul.006G169600.1#bind_motif
# letter-probability matrix: alength= 4 w= 19 nsites= 388 E= (not recorded)

mads_txt <- "0.301546 0.146907 0.134021 0.417526
0.288660 0.172680 0.105670 0.432990
0.198454 0.314433 0.115979 0.371134
0.087629 0.085052 0.043814 0.783505
0.118557 0.033505 0.061856 0.786082
0.167526 0.061856 0.188144 0.582474
0.000000 0.951031 0.000000 0.048969
0.000000 0.682990 0.000000 0.317010
0.469072 0.108247 0.061856 0.360825
0.069588 0.257732 0.061856 0.610825
0.154639 0.000000 0.005155 0.840206
0.033505 0.000000 0.000000 0.966495
0.188144 0.105670 0.067010 0.639175
0.164948 0.090206 0.126289 0.618557
0.219072 0.002577 0.778351 0.000000
0.005155 0.046392 0.804124 0.144330
0.286082 0.252577 0.103093 0.358247
0.543814 0.134021 0.067010 0.255155
0.567010 0.074742 0.085052 0.273196"

make_and_save_logos(
  ppm_txt     = mads_txt,
  motif_name  = "Phvul.006G169600",
  tf_family   = "MIKC-MADS",
  nsites      = 388,
  file_prefix = "Phvul_MIKC-MADS"
)


# =============================================================================
# 3. Nin-like — Phvul.007G071900 (MP00449)
# =============================================================================
# PlantTFDB: http://planttfdb.cbi.pku.edu.cn/tf.php?sp=Pvu&did=Phvul.007G071900.1#bind_motif
# letter-probability matrix: alength= 4 w= 15 nsites= 595 E= 4.9e-576

nin_txt <- "0.447059	  0.109244	  0.117647	  0.326050	
  0.341176	  0.129412	  0.164706	  0.364706	
  0.299160	  0.089076	  0.183193	  0.428571	
  0.010084	  0.181513	  0.001681	  0.806723	
  0.065546	  0.000000	  0.934454	  0.000000	
  0.329412	  0.174790	  0.339496	  0.156303	
  0.000000	  1.000000	  0.000000	  0.000000	
  0.000000	  0.443697	  0.040336	  0.515966	
  0.018487	  0.389916	  0.021849	  0.569748	
  0.010084	  0.080672	  0.000000	  0.909244	
  0.038655	  0.015126	  0.000000	  0.946218	
  0.099160	  0.472269	  0.322689	  0.105882	
  0.400000	  0.112605	  0.415126	  0.072269	
  0.305882	  0.122689	  0.351261	  0.220168	
  0.346218	  0.122689	  0.280672	  0.250420"

make_and_save_logos(
  ppm_txt     = nin_txt,
  motif_name  = "Phvul.007G071900",
  tf_family   = "Nin-like",
  nsites      = 595,
  file_prefix = "Phvul_Ninlike"
)


# =============================================================================
# 4. bHLH — Phvul.001G168700 (MP00074)
# =============================================================================
# PlantTFDB: http://planttfdb.cbi.pku.edu.cn/tf.php?sp=Pvu&did=Phvul.001G168700.1#bind_motif
# letter-probability matrix: alength= 4 w= 14 nsites= 114 E= 0

bhlh_txt <- "0.324561	  0.105263	  0.324561	  0.245614	
  0.219298	  0.000000	  0.464912	  0.315789	
  0.307018	  0.245614	  0.254386	  0.192982	
  0.464912	  0.070175	  0.350877	  0.114035	
  0.236842	  0.236842	  0.280702	  0.245614	
  0.078947	  0.149123	  0.640351	  0.131579	
  0.333333	  0.122807	  0.122807	  0.421053	
  0.000000	  0.903509	  0.096491	  0.000000	
  1.000000	  0.000000	  0.000000	  0.000000	
  0.000000	  0.991228	  0.008772	  0.000000	
  0.000000	  0.000000	  1.000000	  0.000000	
  0.000000	  0.017544	  0.000000	  0.982456	
  0.000000	  0.000000	  1.000000	  0.000000	
  0.192982	  0.201754	  0.447368	  0.157895"

make_and_save_logos(
  ppm_txt     = bhlh_txt,
  motif_name  = "Phvul.001G168700",
  tf_family   = "bHLH",
  nsites      = 114,
  file_prefix = "Phvul_bHLH"
)


# =============================================================================
# 5. bZIP — Phvul.007G229600 (MP00294)
# =============================================================================
# PlantTFDB: http://planttfdb.cbi.pku.edu.cn/tf.php?sp=Pvu&did=Phvul.007G229600.1#bind_motif
# letter-probability matrix: alength= 4 w= 18 nsites= 595 E= 3.2e-1644

bzip_txt <- "0.405042	  0.119328	  0.220168	  0.255462	
  0.364706	  0.099160	  0.280672	  0.255462	
  0.376471	  0.134454	  0.265546	  0.223529	
  0.270588	  0.065546	  0.164706	  0.499160	
  0.073950	  0.031933	  0.672269	  0.221849	
  0.263866	  0.200000	  0.497479	  0.038655	
  0.332773	  0.154622	  0.000000	  0.512605	
  0.000000	  0.485714	  0.494118	  0.020168	
  1.000000	  0.000000	  0.000000	  0.000000	
  0.000000	  1.000000	  0.000000	  0.000000	
  0.000000	  0.000000	  1.000000	  0.000000	
  0.000000	  0.000000	  0.000000	  1.000000	
  0.000000	  0.000000	  1.000000	  0.000000	
  0.000000	  0.000000	  0.749580	  0.250420	
  0.045378	  0.954622	  0.000000	  0.000000	
  0.662185	  0.084034	  0.206723	  0.047059	
  0.253782	  0.208403	  0.292437	  0.245378	
  0.319328	  0.196639	  0.146218	  0.337815"

make_and_save_logos(
  ppm_txt     = bzip_txt,
  motif_name  = "Phvul.007G229600",
  tf_family   = "bZIP",
  nsites      = 595,
  file_prefix = "Phvul_bZIP"
)

