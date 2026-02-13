# Adapted from ADS
cat("\014")
rm(list = ls())

pacman::p_load(pacman, dplyr, cowplot, ggplot2, googleway, cowplot, ggrepel, ggspatial, lwgeom, sf, rnaturalearth, rnaturalearthdata, BiocManager, rmarkdown)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

beans <- data.frame(accessions=c("W6 17491","PI 640965",
"PI 638864","PI 638852",
"PI 638857","PI 661803",
"W6 17020","PI 638858",
"W6 17024","PI 638861",
"PI 661804","PI 638863",
"PI 638856","W6 17017",
"W6 17019","W6 17021",
"PI 638859","W6 17025",
"PI 638860","W6 17028",
"PI 638862","W6 13032"),
lat=c(-23.883,
-25.16138889,
-26.23305556,
-23.70388889,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65,
-23.65),
lng=c(-65.066,-65.61138889,
-65.48305556,-65.53805556,
-65.41666667,-65.41666667,
-65.41666667,-65.41666667,
-65.41666667,-65.41666667,
-65.41666667,-65.41666667,
-65.41666667,-65.41666667,
-65.41666667,-65.41666667,
-65.41666667,-65.41666667,
-65.41666667,-65.41666667,
-65.41666667,-65.41666667),
in11=c("2-NR",
"2-NR",
"1-R",
"1-R",
"2-NR",
"2-NR",
"2-NR",
"2-NR",
"1-R",
"1-R",
"1-R",
"1-R",
"1-R",
"1-R",
"2-NR",
"1-R",
"1-R",
"1-R",
"1-R",
"2-NR",
"1-R",
"1-R")
)

write.csv(beans, file =  "beans.csv")

# BUBBLE MAP WITH SCALED CIRCLES

beans_summary <- beans %>%
    group_by(lat, lng, in11) %>%
    summarize(count = n())  # Aggregate count by location and response type

bubblemap <- ggplot(data = world) +
    geom_sf(fill = "grey85") +
    geom_point(data = beans_summary, aes(x = lng, y = lat, size = count, color = in11), alpha = 1) +
    scale_color_manual(values = c("1-R" = "red", "2-NR" = "blue")) +
    scale_size(range = c(2, 8)) + # Scale sizes to make differences more visible
    coord_sf(xlim = c(-75, -45.5), ylim = c(-32, -15), expand = FALSE) +
       labs(x = "Longitude", y = "Latitude", color = "Response", size = "Count") +
       theme_bw() +
       theme(legend.position = "right")

print(bubblemap)

# ggsave("bubblemap.eps")


# install.packages(c("terra", "tidyterra",   # modern raster handling + ggplot helpers
                   #"elevatr",              # quick DEM downloads
                   #"gganimate",            # optional: turn 12 plots into a gif
                   #"viridis"))             # nice colour scales

pacman::p_load(terra, tidyterra, elevatr, gganimate, viridis)

# Keep Beans as is
# beans <- ...       # <‑‑ your existing code
world  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Get elevation once (static)
# Build a bounding box a little larger than your points
bb <- st_as_sf(st_as_sfc(st_bbox(c(xmin = -62.5, xmax = -68.5,
                                   ymin = -27, ymax = -22.5), crs = 4326))) # c(xmin = -75, xmax = -45.5, ymin = -32, ymax = -15), crs = 4326) original size. 

dem  <- elevatr::get_elev_raster(bb, z = 9, clip = "bbox")  # ~1‑km cells for 9. You'll need lots of RAM for 10 and higher. But it will get nice resolution
dem  <- terra::rast(dem)                                    # coerce to SpatRaster


# Get WorldClim v2.1 monthyl rasters for mean temperature and precipitation
# --- one‑time download, ~20 MB each ----

#Data From: Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37 (12): 4302-4315.

# 2025-05-21: I first made the larger map using 10m of resolution. But to zoom in for the supplemental, I need 30s of resolution, so I redownloaded that and changed the code accordingly. 
# Only the file names needed to be changed. (10m -> 30s)

wc_dir <- "./worldclim30s"

# --------------------------------------------------------------
# 1.  Temperature stack  (wc2.1_30s_tavg_01.tif … _12.tif)
# --------------------------------------------------------------
tavg_files <- list.files(
  wc_dir,
  pattern = "^wc2\\.1_30s_tavg_\\d{2}\\.tif$",   # <-- new pattern
  full.names = TRUE
)

tavg_files <- sort(tavg_files)                   # keeps Jan…Dec order
tavg <- terra::rast(tavg_files) |> terra::crop(bb)
names(tavg) <- month.abb                         # Jan,Feb,…,Dec


# --------------------------------------------------------------
# 2.  Precipitation stack  (wc2.1_30s_prec_01.tif … _12.tif)
# --------------------------------------------------------------
prec_files <- list.files(
  wc_dir,
  pattern = "^wc2\\.1_30s_prec_\\d{2}\\.tif$",
  full.names = TRUE
)

prec <- terra::rast(sort(prec_files)) |> terra::crop(bb)
names(prec) <- month.abb

# Load, crop, and name 12 layers

# Temperature
tavg <- terra::rast(list.files(wc_dir, pattern = "tavg_\\d{2}\\.tif$", full.names = TRUE)) |>
        terra::crop(bb)
names(tavg) <- month.abb                    # Jan–Dec

# Precipitation (mm month‑1)
prec <- terra::rast(list.files(wc_dir, pattern = "prec_\\d{2}\\.tif$", full.names = TRUE)) |>
        terra::crop(bb)
names(prec) <- month.abb

# Creating the base of the map
# Preparing DEM for plotting
dem_df <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)

names(dem_df)[3] <- "elev"            # or: dem_df <- rename(dem_df, elev = 3)

base_map <- ggplot() +
  # Elevation shading first
  geom_raster(data = dem_df,
              aes(x, y, fill = elev), interpolate = TRUE, alpha = 1) + 
  scale_fill_viridis_c(name = "Elevation (m)", option = "C", direction = -1) +

  # World country borders on top of DEM
  geom_sf(data = world, fill = NA, colour = "grey20", linewidth = 0.3) +

  # Your beans as bubbles
  geom_point(data = beans_summary,
             aes(lng, lat, size = count, colour = in11),
             stroke = 0.2) +
  scale_colour_manual(values = c("1-R" = "red", "2-NR" = "blue"),
                      name = "Response") +
  scale_size(range = c(2, 8), name = "Count") +

  coord_sf(xlim = c(-68.5, -62.5),
           ylim = c(-27, -22.5),
           expand = FALSE) +
  theme_bw() +
  theme(legend.position = "right")

base_map

# You can also save this and clip to whichever coordinates you want. Just make sure to keep the coordinates consistent throughout the script. 


pacman::p_load("patchwork", "ggnewscale", "ggplot2", "terra")

make_climate_map <- function(month_abbr,
                             clim_stack,          # tavg  *or*  prec
                             legend_title,
                             palette_option) {
  # 1. find the layer index for this month
  idx <- match(month_abbr, names(clim_stack))
  if (is.na(idx)) {
    stop("Month “", month_abbr, "” not found in raster names:\n",
         paste(names(clim_stack), collapse = ", "))
  }
  # 2. extract that single-layer SpatRaster
  layer <- clim_stack[[idx]]
  
  # 3. turn it into a simple data.frame with x,y,value
  clim_df <- terra::as.data.frame(layer, xy = TRUE, na.rm = TRUE)
  names(clim_df)[3] <- "value"
  
  # 4. build the ggplot
  ggplot() +
    # elevation (first fill scale)
    geom_raster(data = dem_df,
                aes(x, y, fill = elev),
                interpolate = TRUE, alpha = 0.40) +
    scale_fill_viridis_c(name = "Elevation (m)",
                         option = "C", direction = -1) +
    
    # reset fill for the next layer
    ggnewscale::new_scale_fill() +
    
    # climate layer (second fill scale)
    geom_raster(data = clim_df,
                aes(x, y, fill = value),
                alpha = 0.70) +
    scale_fill_viridis_c(name   = legend_title,
                         option = palette_option,
                         guide  = guide_colourbar(barwidth = 10)) +
    
    # country outlines + beans
    geom_sf(data = world, fill = NA,
            colour = "grey20", linewidth = 0.30) +
    geom_point(data = beans_summary,
               aes(lng, lat, size = count, colour = in11),
               stroke = 0.20) +
    scale_colour_manual(values = c("1-R" = "red", "2-NR" = "blue"),
                        name = "Response") +
    scale_size(range = c(2, 8), name = "Count") +
    
    coord_sf(xlim = c(-62.5, -68.5), ylim = c(-27, -22.5), expand = FALSE) +
    theme_bw() +
    theme(legend.position = "right") +
    labs(title = paste("Bean accessions –", month_abbr))
}

# programs for looping 12 months of each and displaying side by side. 
for (m in month.abb) {
  p_temp <- make_climate_map(
    month_abbr     = m,
    clim_stack     = tavg,
    legend_title   = "Mean T (°C ×10)",
    palette_option = "A"
  )
  p_prec <- make_climate_map(
    month_abbr     = m,
    clim_stack     = prec,
    legend_title   = "Rain (mm)",
    palette_option = "B"
  )
  
  combined <- p_temp + p_prec + plot_layout(guides = "collect")
  
  ggsave(
    filename = paste0("beans_", tolower(m), "_side.eps"),
    plot     = combined,
    width    = 12,
    height   = 7,
    dpi      = 300
  )
}

