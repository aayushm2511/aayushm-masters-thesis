# ============================================================
# Response Diversity & Community Imbalance Analysis
# Western Ghats Tree Communities — CWD Gradient
# ============================================================
#
# PIPELINE://,m
#   1.  Load & inspect data
#   2.  Collapse sites sharing identical CWD values
#         — sp.occ = 1 if present in ANY of the combined sites
#         — results in 143 unique CWD values used as the analytical unit
#   3.  Fit one binomial GAM per species (sp.occ ~ s(CWD))
#   4.  Build species × CWD-group PROBABILITY matrix (183 × 143)
#   5.  Build species × CWD-group DERIVATIVE matrix  (183 × 143)
#         — central finite differences on the probability scale: dP/dCWD
#   6.  Build species × CWD-group PRESENCE-ABSENCE matrix (183 × 143)
#   7.  Mask derivative matrix to species present at each CWD group
#   8.  Calculate IMBALANCE per CWD group
#         = |mean(dP/dCWD)| across species present at that CWD group
#   9.  Visualisations
#  10.  Save outputs
# REFERENCES:
#   Polazzo et al. (2025) Ecology Letters 28: e70224
#   Ross et al. (2023) Methods in Ecology and Evolution 14: 1150-1167
# ============================================================
#AI USAGE:
# Claude Sonnet 4.6 was used to refactor and debug the code given below. User-edits were made to ensure conformation to Polazzo et al., 2025. 
rm(list = ls())
# ── 0. Packages ──────────────────────────────────────────────

library(mgcv)       
library(tidyverse)  
library(ggplot2)    
library(patchwork)  


#raw_data

dat_raw <- read.csv("/Users/aayushmanchalwar/Desktop/thesis/gams/data/geb13350-sup-0002-tree-species-occ.csv",
                    stringsAsFactors = FALSE)

cat("─── Raw data summary ────────────────────────────────────\n")
cat("Total rows       :", nrow(dat_raw), "\n")
cat("Unique sites     :", n_distinct(dat_raw$Plot_ID), "\n")
cat("Unique CWD values:", n_distinct(dat_raw$CWD), "\n")
cat("Unique species   :", n_distinct(dat_raw$Species), "\n")
cat("────────────────────────────────────────────────────────\n\n")


#cwd grouping

dat <- dat_raw %>%
  group_by(CWD, Species) %>%
  summarise(sp.occ = max(sp.occ), .groups = "drop")


cwd_groups <- dat %>%
   
  dplyr::select(CWD) %>%
  distinct() %>%
  arrange(CWD) %>%
  mutate(cwd_id = seq_len(n()))

n_cwd    <- nrow(cwd_groups)    # 143
cwd_vec  <- cwd_groups$CWD     

species_list <- sort(unique(dat$Species))
n_species    <- length(species_list)    # 183

cat("─── Collapsed data summary ──────────────────────────────\n")
cat("CWD groups (analytical units) :", n_cwd,     "\n")
cat("Species                       :", n_species, "\n")
cat("CWD range                     :",
    round(min(cwd_vec), 1), "to", round(max(cwd_vec), 1), "mm\n")

# Occurrence counts after collapsing
occ_counts <- dat %>%
  filter(sp.occ == 1) %>%
  count(Species, name = "n_presences") %>%
  arrange(n_presences)

cat("Presences per species (post-collapse):\n")
cat("  Min:", min(occ_counts$n_presences),
    "| Median:", median(occ_counts$n_presences),
    "| Max:", max(occ_counts$n_presences), "\n")
cat("────────────────────────────────────────────────────────\n\n")


#GAM fitting

h <- 0.5   # mm — finite difference step for dP/dCWD


prob_matrix  <- matrix(NA_real_, nrow = n_species, ncol = n_cwd,
                       dimnames = list(species_list, as.character(cwd_vec)))
deriv_matrix <- matrix(NA_real_, nrow = n_species, ncol = n_cwd,
                       dimnames = list(species_list, as.character(cwd_vec)))

gam_list        <- vector("list", n_species)
names(gam_list) <- species_list

cat("Fitting GAMs for", n_species, "species...\n")


for (i in seq_along(species_list)) {

  sp     <- species_list[i]
  sp_dat <- dat %>% filter(Species == sp)

  # ── Fit GAM ──────────────────────────────────────────────
  gam_fit <- gam(
    sp.occ ~ s(CWD, bs = "tp"),
    family = binomial(link = "logit"),
    data   = sp_dat,
    method = "REML"
  )
  gam_list[[sp]] <- gam_fit

  # ── P(occurrence) at each CWD group — response (probability) scale ──
  prob_matrix[i, ] <- predict(
    gam_fit,
    newdata = data.frame(CWD = cwd_vec),
    type    = "response"
  )

  # ── dP/dCWD — central finite difference on the probability scale ────
  # dP/dCWD ≈ [P(CWD + h) − P(CWD − h)] / (2h)
  P_plus  <- predict(gam_fit,
                     newdata = data.frame(CWD = cwd_vec + h),
                     type    = "response")
  P_minus <- predict(gam_fit,
                     newdata = data.frame(CWD = cwd_vec - h),
                     type    = "response")

  deriv_matrix[i, ] <- (as.numeric(P_plus) - as.numeric(P_minus)) / (2 * h)


}

cat("\nAll GAMs fitted successfully.\n\n")


# presence-absence matrix

pa_wide <- dat %>%
  dplyr::select(CWD, Species, sp.occ) %>%
  pivot_wider(names_from = CWD, values_from = sp.occ)

pa_matrix <- as.matrix(
  pa_wide[match(species_list, pa_wide$Species),
          as.character(cwd_vec)]
)
rownames(pa_matrix) <- species_list

# integrity checks
stopifnot(
  "PA matrix row count mismatch"      = nrow(pa_matrix) == n_species,
  "PA matrix column count mismatch"   = ncol(pa_matrix) == n_cwd,
  "PA matrix contains non-0/1 values" = all(pa_matrix %in% c(0L, 1L))
)
cat("Presence-absence matrix built and verified.\n")
cat("Total presences in PA matrix:", sum(pa_matrix), "\n\n")




masked_deriv <- deriv_matrix * pa_matrix


# imbalance calculation
#
# Polazzo et al. (2025) fundamental imbalance:
#   Imbalance_j = | (1/S_j) * Σ_i (dP_i/dCWD)_j |
#
# S_j = number of species present at CWD group j.
# Sum and mean are over present species only.

n_present_vec  <- integer(n_cwd)

imbalance_vec  <- numeric(n_cwd)

for (j in seq_len(n_cwd)) {

  present_idx       <- which(pa_matrix[, j] == 1)
  n_present_vec[j]  <- length(present_idx)

  if (length(present_idx) > 0) {
    mean_d             <- mean(masked_deriv[present_idx, j])
    
    imbalance_vec[j]   <- abs(mean_d)
  } else {
    
    imbalance_vec[j]   <- NA_real_
  }
}

# assemble results table
results <- cwd_groups %>%
  mutate(
    n_species_present = n_present_vec, #species richness
    
    imbalance         = imbalance_vec
  )
print(summary(results$imbalance))



#viz

# select 9 representative species for curve panels ──
# 3 rarest, 3 around the median, 3 most common (post-collapse counts)
n_sp <- nrow(occ_counts)
sampled_species <- c(
  occ_counts$Species[1:3],
  occ_counts$Species[(floor(n_sp / 2) - 1):(floor(n_sp / 2) + 1)],
  occ_counts$Species[(n_sp - 2):n_sp]
)

cwd_seq <- seq(min(cwd_vec), max(cwd_vec), length.out = 500)



# sample prob. occ. curves

prob_plot_dat <- lapply(sampled_species, function(sp) {
  pred  <- predict(gam_list[[sp]],
                   newdata = data.frame(CWD = cwd_seq),
                   type    = "response",
                   se.fit  = TRUE)
  n_occ <- occ_counts$n_presences[occ_counts$Species == sp]
  tibble(
    sp_label = paste0(sp, "  (n=", n_occ, ")"),
    CWD      = cwd_seq,
    prob     = as.numeric(pred$fit),
    se       = as.numeric(pred$se.fit)
  )
}) %>% bind_rows() %>%
  mutate(sp_label = factor(sp_label, levels = unique(sp_label)))

fig1 <- ggplot(prob_plot_dat, aes(x = CWD, y = prob)) +
  geom_ribbon(aes(ymin = pmax(0, prob - 1.96 * se),
                  ymax = pmin(1, prob + 1.96 * se)),
              fill = "steelblue", alpha = 0.20) +
  geom_line(colour = "steelblue3", linewidth = 0.9) +
  facet_wrap(~ sp_label, ncol = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "GAM-fitted probability of occurrence curves",
    subtitle = "Top row: rare | Middle: median | Bottom: common species",
    x        = "Climatic Water Deficit (mm)",
    y        = "P(occurrence)"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 7.5))


# derivative curves

deriv_plot_dat <- lapply(sampled_species, function(sp) {
  P_plus  <- as.numeric(predict(gam_list[[sp]],
                                newdata = data.frame(CWD = cwd_seq + h),
                                type    = "response"))
  P_minus <- as.numeric(predict(gam_list[[sp]],
                                newdata = data.frame(CWD = cwd_seq - h),
                                type    = "response"))
  n_occ   <- occ_counts$n_presences[occ_counts$Species == sp]
  tibble(
    sp_label = paste0(sp, "  (n=", n_occ, ")"),
    CWD      = cwd_seq,
    deriv    = (P_plus - P_minus) / (2 * h)
  )
}) %>% bind_rows() %>%
  mutate(sp_label = factor(sp_label, levels = levels(prob_plot_dat$sp_label)))

fig2 <- ggplot(deriv_plot_dat, aes(x = CWD, y = deriv)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_line(colour = "tomato3", linewidth = 0.9) +
  facet_wrap(~ sp_label, ncol = 3, scales = "free_y") +
  labs(
    title    = "First derivatives of probability curves  (dP/dCWD)",
,
    x        = "Climatic Water Deficit (mm)",
    y        = "dP/dCWD"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 7.5))


# imb vs cwd

fig3 <- ggplot(results, aes(x = CWD, y = imbalance)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method  = "gam",
              formula = y ~ s(x, bs = "tp", k = 20),
              colour  = "black", linewidth = 0.9, se = TRUE) +
  
  labs(
    title    = "Community imbalance across the CWD gradient",
    subtitle = "n = 143",
    x        = "Climatic Water Deficit (mm)",
    y        = "Imbalance  (|mean dP/dCWD|)"
  ) +
  theme_bw(base_size = 12)



# richness vs cwd

fig4<- ggplot(results, aes(x = CWD, y = n_species_present)) +
  geom_point(colour = "steelblue3", size = 2.2, alpha = 0.8) +
  geom_smooth(method  = "gam",
              formula = y ~ s(x, bs = "tp", k = 15),
              colour  = "darkorange", linewidth = 0.9, se = TRUE) +
  labs(
    title    = "Species richness across the CWD gradient",
    subtitle = "Observed richness per CWD group (n = 143)",
    x        = "Climatic Water Deficit (mm)",
    y        = "Number of species present"
  ) +
  theme_bw(base_size = 12)






print(fig1)
print(fig2)
print(fig3)
print(fig4)



# output

write.csv(results, "site_imbalance.csv", row.names = FALSE)

# imbalance objtect for thesis_analysis.R
saveRDS(
  list(
    # matrices: 183 rows (species, alphabetical) × 143 cols (CWD, low→high)
    prob_matrix   = prob_matrix,    # GAM-predicted P(occurrence)
    deriv_matrix  = deriv_matrix,   # dP/dCWD (all species)
    pa_matrix     = pa_matrix,      # binary presence-absence 
    masked_deriv  = masked_deriv,   # deriv_matrix × pa_matrix
    
    gam_list      = gam_list,       # 183 gams
    
    cwd_groups    = cwd_groups,     
    species       = species_list,   
    
    imbalance     = results         # full results table 
  ),
  file = "imbalance_objects.rds"
)

