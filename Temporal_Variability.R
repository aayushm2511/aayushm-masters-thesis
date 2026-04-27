################################################################################
### Temporal Variability — White et al. (2020), grouped by CWD              ###
### Western Ghats EVI — 16-day MODIS composites, 2004–2024                  ###
################################################################################
###
### FORMULA (White et al. 2020):
###
###   v_i = sd[ EVI(i,t) − mean_{u∈m}[ EVI(i,u) ] ]²
###        = var( EVI − monthly_mean )  across all time points t
###
###
### DATA CLEANING (EVI):
###   - Year == 1899 rows: empty parsing artefacts, removed.
###   - EVI < 0: MODIS fill values (−3000) set to zero
###
###AI USAGE:
###Claude Sonnet 4.6 was used to refactor and debug the code given below. User-edits were made to ensure conformation to White et al., 2020 in variability calculations
################################################################################

rm(list = ls())

library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)



evi_raw <- read.csv("/Users/aayushmanchalwar/Desktop/thesis/stability/cleaned_evi.csv",                    stringsAsFactors = FALSE)
sp_raw  <- read.csv("/Users/aayushmanchalwar/Desktop/thesis/gams/data/geb13350-sup-0002-tree-species-occ.csv",  stringsAsFactors = FALSE)

cat("EVI raw rows:", nrow(evi_raw), "\n")
cat("Species raw rows:", nrow(sp_raw), "\n\n")




evi <- evi_raw %>%
  filter(Year != 1899,
         !is.na(ID), ID != "",
         !is.na(EVI))

evi <- evi %>% mutate(EVI = pmax(EVI, 0))


# CWD grouping


cwd_lookup <- sp_raw %>%
  group_by(Plot_ID) %>%
  summarise(CWD = first(CWD), .groups = "drop") %>%
  rename(ID = Plot_ID)

cat("CWD lookup: ", nrow(cwd_lookup), "plots,",
    n_distinct(cwd_lookup$CWD), "unique CWD values\n")

# verification
missing_cwd <- setdiff(unique(evi$ID), cwd_lookup$ID)
if (length(missing_cwd) > 0) {
  cat("WARNING: ", length(missing_cwd),
      "EVI sites have no CWD — they will be dropped:\n")
  print(missing_cwd)
} else {
  cat("All EVI sites matched to a CWD value.\n")
}
cat("\n")




evi <- evi %>% left_join(cwd_lookup, by = "ID")


group_summary <- evi %>%
  group_by(CWD) %>%
  summarise(n_sites = n_distinct(ID),
            n_obs   = n(),
            .groups = "drop")

cat("CWD groups summary:\n")
cat("  Total groups:", nrow(group_summary), "\n")
cat("  Groups with 1 site:", sum(group_summary$n_sites == 1), "\n")
cat("  Groups with 2+ sites:", sum(group_summary$n_sites > 1), "\n")
cat("  Max sites in one group:", max(group_summary$n_sites), "\n")
cat("  Obs per group — median:", median(group_summary$n_obs),
    "  range:", min(group_summary$n_obs), "–", max(group_summary$n_obs), "\n\n")

# Exclude JJAS
MONSOON_MONTHS <- c(6, 7, 8, 9)
evi <- evi %>% filter(!Month %in% MONSOON_MONTHS)


#
# Pool all EVI observations from all sites within the same CWD group.
# baseline_mean[CWD, m] = mean of all EVI values in that group × month cell,
# across all sites and all years.

baseline <- evi %>%
  group_by(CWD, Month) %>%
  summarise(
    monthly_mean = mean(EVI, na.rm = TRUE),
    n_baseline   = n(),
    .groups = "drop"
  )


# anomaly(t) = EVI(t) − monthly_mean[ CWD_group, month(t) ]

evi_anomaly <- evi %>%
  left_join(baseline, by = c("CWD", "Month")) %>%
  mutate(anomaly = EVI - monthly_mean)



# v_i  = var(anomalies)          — White et al. temporal variability
# CV_i = sd(EVI) / mean(EVI)    — coefficient of variation, raw EVI
#
# Mean lat/lon: centroid of the group's member sites

site_coords <- evi %>%
  group_by(ID, CWD) %>%
  summarise(Latitude  = first(Latitude),
            Longitude = first(Longitude),
            .groups = "drop")

group_coords <- site_coords %>%
  group_by(CWD) %>%
  summarise(Latitude_mean  = mean(Latitude),
            Longitude_mean = mean(Longitude),
            n_sites        = n(),
            member_sites   = paste(sort(ID), collapse = "; "),
            .groups = "drop")

variability <- evi_anomaly %>%
  group_by(CWD) %>%
  summarise(
    evi_variability = var(anomaly,  na.rm = TRUE),
    evi_cv          = sd(EVI,       na.rm = TRUE) / mean(EVI, na.rm = TRUE),
    n_obs           = sum(!is.na(anomaly)),
    n_months        = n_distinct(Month[!is.na(anomaly)]),
    .groups = "drop"
  ) %>%
  left_join(group_coords, by = "CWD") %>%
  mutate(
    data_flag = case_when(
      n_obs < 24      ~ "insufficient_data",
      n_months < 12   ~ "incomplete_seasonality",
      TRUE            ~ "ok"
    )
  ) %>%
  arrange(CWD)


write.csv(variability, "evi_temporal_variability.csv", row.names = FALSE)
cat("Saved: evi_temporal_variability.csv\n")
cat("Columns: CWD, evi_variability, evi_cv, n_obs, n_months, n_sites,\n")
cat("         Latitude_mean, Longitude_mean, member_sites, data_flag\n\n")



ok <- variability 



# viz

p1a <- ggplot(ok, aes(x = evi_variability)) +
  geom_histogram( bins = 30,
                 fill = viridis(1, begin = 0.4), colour = "white",
                 linewidth = 0.2) +
 
  labs(title = "Temporal Variability",
       subtitle = "var(EVI − monthly mean)  |  pooled within CWD groups",
       x = "v_i", y = "Frequency") +
  theme_classic(base_size = 13)

p1b <- ggplot(ok, aes(x = evi_cv)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = viridis(1, begin = 0.65), colour = "white",
                 linewidth = 0.2) +
  geom_density(linewidth = 0.8) +
  labs(title = "Coefficient of Variation  (CV)",
       subtitle = "sd(EVI) / mean(EVI)  |  pooled within CWD groups",
       x = "CV", y = "Density") +
  theme_classic(base_size = 13)

p1 <- plot_grid(p1a, p1b, ncol = 2, labels = c("a", "b"))
ggsave("plot1_distributions.png", p1, width = 12, height = 5, dpi = 300)
cat("Saved: plot1_distributions.png\n")


# v_i vs CV

p2a <- ggplot(ok, aes(x = CWD, y = evi_variability)) +
  geom_point(aes(size = n_sites), alpha = 0.6,
             colour = viridis(1, begin = 0.3)) +
  geom_smooth(method = "loess", span = 0.5,
              colour = "firebrick", fill = "firebrick",
              alpha = 0.2, linewidth = 0.9) +
  scale_size_continuous(name = "Sites\nin group", range = c(1.5, 6)) +
  labs(title = "Temporal Variability across the CWD gradient",
       subtitle = "Point size = number of sites pooled in that CWD group",
       x = "CWD  (more negative = more water-stressed)",
       y = "v_i  [var of EVI anomalies]") +
  theme_classic(base_size = 13)

p2b <- ggplot(ok, aes(x = CWD, y = evi_cv)) +
  geom_point(aes(size = n_sites), alpha = 0.6,
             colour = viridis(1, begin = 0.6)) +
  geom_smooth(method = "loess", span = 0.5,
              colour = "firebrick", fill = "firebrick",
              alpha = 0.2, linewidth = 0.9) +
  scale_size_continuous(name = "Sites\nin group", range = c(1.5, 6)) +
  labs(title = "CV across the CWD gradient",
       x = "CWD  (more negative = more water-stressed)",
       y = "CV  [sd/mean of raw EVI]") +
  theme_classic(base_size = 13)

p2 <- plot_grid(p2a, p2b, ncol = 1, labels = c("a", "b"))
ggsave("plot2_variability_vs_CWD.png", p2, width = 9, height = 10, dpi = 300)
cat("Saved: plot2_variability_vs_CWD.png\n")



p3 <- ggplot(ok, aes(x = evi_variability, y = evi_cv)) +
  geom_point(aes(size = n_sites), alpha = 0.5,
             colour = viridis(1, begin = 0.55)) +
  geom_smooth(method = "lm", colour = "firebrick",
              fill = "firebrick", alpha = 0.2, linewidth = 0.9) +
  scale_size_continuous(name = "Sites\nin group", range = c(1.5, 6)) +
  labs(title = "v_i vs CV across CWD groups",
       subtitle = "If correlated, both metrics capture the same signal",
       x = "Temporal Variability v_i  [var of anomalies]",
       y = "CV  [sd/mean of raw EVI]") +
  theme_classic(base_size = 13)

ggsave("plot3_vi_vs_cv.png", p3, width = 7, height = 6, dpi = 300)
cat("Saved: plot6_vi_vs_cv.png\n\n")

