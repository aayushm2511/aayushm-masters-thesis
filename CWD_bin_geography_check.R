#checking within-group spatial spread
###Claude Sonnet 4.6 was used to generate the code given below. Ouput was manually checked, no edits were made. 

library(tidyverse)
library(ggplot2)


occ <- read_csv("/Users/aayushmanchalwar/Desktop/thesis/gams/data/geb13350-sup-0002-tree-species-occ.csv", show_col_types = FALSE) %>%
  dplyr::distinct(Plot_ID, Latitude, Longitude, CWD) %>%
  dplyr::mutate(CWD_key = round(CWD, 4))

cat(sprintf("%d unique plots across %d CWD bins\n",
            n_distinct(occ$Plot_ID), n_distinct(occ$CWD_key)))


bin_summary <- occ %>%
  dplyr::group_by(CWD_key) %>%
  dplyr::summarise(
    n_plots   = n_distinct(Plot_ID),
    lat_mean  = mean(Latitude),
    lon_mean  = mean(Longitude),
    lat_sd    = sd(Latitude),
    lon_sd    = sd(Longitude),
    lat_range = diff(range(Latitude)),
    lon_range = diff(range(Longitude)),
    .groups   = "drop"
  )

cat(sprintf("\nBins with >1 plot: %d / %d\n", sum(bin_summary$n_plots > 1), nrow(bin_summary)))
cat(sprintf("Lat spread within bin (median): %.3f deg\n", median(bin_summary$lat_range, na.rm=TRUE)))
cat(sprintf("Lon spread within bin (median): %.3f deg\n", median(bin_summary$lon_range, na.rm=TRUE)))
cat(sprintf("Max lat spread within any bin:  %.3f deg\n", max(bin_summary$lat_range, na.rm=TRUE)))

# Within-bin spatial spread distribution
spread_df <- bin_summary %>%
  dplyr::filter(n_plots > 1) %>%
  dplyr::mutate(
    max_spread_deg = pmax(lat_range, lon_range),
    max_spread_km  = max_spread_deg * 111  # rough conversion
  )

fig1 <- ggplot(spread_df, aes(x = max_spread_km)) +
  geom_histogram(bins = 20, fill = "#2C7BB6", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = median(spread_df$max_spread_km),
             linetype = "dashed", colour = "#D7191C", linewidth = 0.9) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("median = %.1f km\nmax = %.1f km",
                           median(spread_df$max_spread_km),
                           max(spread_df$max_spread_km)),
           hjust = 1.1, vjust = 1.3, size = 3.5, colour = "grey20") +
  labs(title    = "Within-bin spatial spread (multi-plot bins only)",
       subtitle = sprintf("%d bins with >1 plot | max spread = max(lat_range, lon_range)",
                          nrow(spread_df)),
       x = "Maximum within-bin spread (km, approximate)",
       y = "Number of bins") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  colour = "grey40"))

ggsave("CWD_bin_map_spread.png", fig1, width = 8, height = 5, dpi = 300)
cat("Saved CWD_bin_map_spread.png\n")

