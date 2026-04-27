# ============================================================
# THESIS ANALYSIS: Response Diversity, Richness, and Ecosystem
#                  Stability Along a CWD Gradient
#                  Western Ghats Tropical Forests
#
# Methods chapter coverage:
#   §2.1  Data preparation
#   §2.2  imbalance
#   §2.3  EVI variability
#   §2.4  Species richness 
#   §2.5  Statistical analysis:
#           - Pearson correlations + gradient scatterplots (EDA)
#           - Moran's I on OLS and GAM residuals
#           - Piecewise SEM with GLS corExp components
#             (Lefcheck 2016; Dormann et al. 2007)
# ============================================================
###Claude Sonnet 4.6 was used to refactor and debug the code given below. Unedited code for visualisation was provided by Claude.



library(tidyverse)    
library(mgcv)         
library(nlme)         
library(spdep)        
library(piecewiseSEM) 
library(vegan)        
library(ggcorrplot)   
library(patchwork)    

set.seed(42)

#joining data
imbalance_df <- read.csv("site_imbalance.csv")
evi_df       <- read.csv("evi_temporal_variability.csv")


cat("=== imbalance_df ===\n"); str(imbalance_df)
cat("\n=== evi_df ===\n");     str(evi_df)

cat("\nCWD range (imbalance_df):", paste(range(imbalance_df$CWD), collapse = " to "), "\n")
cat("CWD range (evi_df)      :", paste(range(evi_df$CWD),        collapse = " to "), "\n")

dat <- merge(imbalance_df, evi_df, by = "CWD", all = FALSE)
cat("\nRows after merge:", nrow(dat), "(expected 143)\n")
if (nrow(dat) < 140) stop("Fewer than 140 rows after merge — check CWD matching.")

dat <- dat %>%
  mutate(
    log_imbalance = log10(imbalance),
    log_richness  = log10(n_species_present),
    log_evi_var   = log10(evi_variability)
  )

cat("\n=== Key variable summaries ===\n")
dat %>%
  dplyr::select(CWD, imbalance, n_species_present, evi_variability,
                log_imbalance, log_richness, log_evi_var) %>%
  summary() %>%
  print()


#exploratory data analysis

# correlation matrix
cor_vars <- dat %>%
  dplyr::select(log_evi_var, log_imbalance, log_richness, CWD) %>%
  rename(
    "log10(EVI Var)"   = log_evi_var,
    "log10(Imbalance)" = log_imbalance,
    "log10(Richness)"  = log_richness,
    "CWD"              = CWD
  )

cor_mat <- cor(cor_vars, use = "complete.obs", method = "pearson")

cat("=== Pearson correlation matrix ===\n")
print(round(cor_mat, 3))


p_cor <- ggcorrplot(
  cor_mat,
  method   = "circle",
  type     = "lower",
  lab      = TRUE,
  lab_size = 4,
  colors   = c("#D7191C", "white", "#2C7BB6"),
  p.mat    = cor_pmat(cor_vars),
  title    = "Fig. EDA-1: Pearson correlation matrix (log-transformed variables)"
) + theme_bw(base_size = 12)
print(p_cor)


p_pairs <- ggpairs(
  cor_vars,
  lower = list(continuous = wrap("smooth", method = "loess",
                                 colour = "#2C7BB6", alpha = 0.5, se = TRUE)),
  diag  = list(continuous = wrap("densityDiag", fill = "#2C7BB6", alpha = 0.5)),
  upper = list(continuous = wrap("cor", size = 4.5))
) +
  
  theme_bw(base_size = 11)
print(p_pairs)


# imb/rich/var vs CWD

gradient_panel <- function(data, y_var, y_lab, pt_col, k = 10) {
  ggplot(data, aes_string(x = "CWD", y = y_var)) +
    geom_point(alpha = 0.5, colour = pt_col, size = 1.8) +
    geom_smooth(method  = "gam",
                formula = y ~ s(x, k = k),
                colour  = "black", se = TRUE, linewidth = 0.9) +
    labs(x        = "CWD (mm/yr)",
         y        = y_lab,
         subtitle = "GAM smoother ± 95% CI") +
    theme_bw(base_size = 12)
}

p_g1 <- gradient_panel(dat, "log_imbalance", "log10(Imbalance)",        "#D7191C")
p_g2 <- gradient_panel(dat, "log_richness",  "log10(Species Richness)", "#1A9641")
p_g3 <- gradient_panel(dat, "log_evi_var",   "log10(EVI Variability)",  "#2C7BB6")

print(
  (p_g1 | p_g2) / p_g3 +
    plot_annotation(
      title    = "Fig. 1: Gradient analyses along the CWD gradient",
      subtitle = "More negative CWD = higher drought stress"
    )
)

# richness vs imb
p_rich_imb <- ggplot(dat, aes(x = log_richness, y = log_imbalance)) +
  geom_point(aes(colour = CWD), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", colour = "black", se = TRUE, linewidth = 1) +
  scale_colour_viridis_c(name = "CWD (mm/yr)", direction = -1) +
  labs(
    title = "Fig. 2: Richness vs imbalance along the CWD gradient",
    x     = expression(log[10](Species~Richness)),
    y     = expression(log[10](Imbalance))
  ) +
  theme_bw(base_size = 12)
print(p_rich_imb)


#ols model to extract model residuals for moran's I test

M1 <- lm(
  log_evi_var ~ log_imbalance + log_richness + CWD,
  data = dat
)
print(summary(M1))




# Moran's I calculated on OLS residuals, Siginificance justifies SEM usage with GLS component models (ref. Dormann et al., 2007)

#converting coords to kms

dat <- dat %>%
  mutate(
    lon_km = Longitude_mean * 111 * cos(mean(Latitude_mean) * pi / 180),
    lat_km = Latitude_mean  * 111
  )

# spatial weights for morans test. knn used since sites have no intrinsice neighbourhood structure

coords_mat <- as.matrix(dat[, c("Longitude_mean", "Latitude_mean")])
knn8       <- knearneigh(coords_mat, k = 8)
nb8        <- knn2nb(knn8)
lw8        <- nb2listw(nb8, style = "W")


# moran's tests

print(moran.test(dat$log_evi_var, lw8))

cat("\n=== Moran's I: M1 (OLS) residuals ===\n")
print(moran.test(residuals(M1), lw8))




#piecewiseSEM
# Each component model is a GLS with an exponential spatial covariance structure (corExp) and nugget term

# component models
m_evi <- gls(
  log_evi_var ~ log_imbalance + log_richness + CWD,
  data        = dat,
  correlation = corExp(form = ~lon_km + lat_km, nugget = TRUE),
  method      = "ML"
)

m_imb <- gls(
  log_imbalance ~ CWD,
  data        = dat,
  correlation = corExp( form = ~lon_km + lat_km, nugget = TRUE),
  method      = "ML"
)

m_rich <- gls(
  log_richness ~ CWD,
  data        = dat,
  correlation = corExp( form = ~lon_km + lat_km, nugget = TRUE),
  method      = "ML"
)

# testing if SAC corrections worked

print(moran.test(residuals(m_evi, type = "normalized"), lw8))
# p > 0.05 confirms corExp successfully removed spatial autocorrelation

#fitting pSEM
psem_fit <- psem(m_evi, m_imb, m_rich, data = dat)

cat("\n=== Piecewise SEM summary ===\n")
summary(psem_fit, .progressBar = FALSE)


# d-sep test to check if richness and imbalance have any relationship

print(dSep(psem_fit, .progressBar = FALSE))

# path coefs
cat("\n=== Standardised path coefficients ===\n")
coefs(psem_fit, standardize = "scale", intercepts = FALSE) %>%
  dplyr::select(Response, Predictor, Std.Estimate, Estimate,
                Std.Error, P.Value) %>%
  mutate(
    across(where(is.numeric), ~round(., 4)),
    sig = case_when(
      P.Value < 0.001 ~ "***",
      P.Value < 0.01  ~ "**",
      P.Value < 0.05  ~ "*",
      P.Value < 0.10  ~ ".",
      TRUE            ~ "ns"
    )
  ) %>%
  print(row.names = FALSE)



#viz


p_A <- ggplot(dat, aes(x = log_imbalance, y = log_evi_var)) +
  geom_point(aes(colour = CWD), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", colour = "black",
              se = TRUE, linewidth = 1) +
  scale_colour_viridis_c(name = "CWD (mm/yr)", direction = -1) +
  labs(x = expression(log[10](Imbalance)),
       y = expression(log[10](EVI~Variability))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p_B <- ggplot(dat, aes(x = log_richness, y = log_evi_var)) +
  geom_point(aes(colour = CWD), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", colour = "black",
              se = TRUE, linewidth = 1) +
  scale_colour_viridis_c(name = "CWD (mm/yr)", direction = -1) +
  labs(x = expression(log[10](Species~Richness)),
       y = expression(log[10](EVI~Variability))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(
  p_A + p_B +
    plot_annotation(
      title      = "Fig. 4: Imbalance and richness vs EVI variability",
      subtitle   = "Points coloured by CWD. Lines are OLS fits ± 95% CI.",
      tag_levels = "A"
    )
)

