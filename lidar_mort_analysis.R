# Street tree mortality estimates, models, and risk maps for New York City
# script by Dan Katz but borrows some elements from previous work by Dave Miller


#next steps from here:
# get new data from ging-yan:
# add population density
# add the land use at various scales including tree cover
# updated construction variables
#
# do model selection across taxa to decide which variants of spatial data to include
# e.g., impervious at 5 10 or 15 m
# do model selection to decide whether the pluto data should be included 
# create an SI with model selection results

# revise forest plot with coefficients
# create an SI with model diagnostics (ROC values, maybe binned residuals)



#### set up work environment and load data #####################################
#library(tidyverse)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(data.table)

library(ciTools) 
library(sf)
library(basemaps)
library(terra)
library(tidyterra)
library(ggspatial)
library(viridis)
library(scales)

# library(DHARMa)
# library(splines)
library(INLA)
library(inlabru)
library(fmesher)
library(sn) #needed by INLA::inla.posterior.sample (joint posterior draws of the fixed effects)

#library(ResourceSelection) 
library(pROC)
library(metafor)

library(patchwork)
library(knitr)
library(kableExtra)
library(gt)


data_dir <- "data"
#v4 includes the socioeconomic variables; see all_variables_v4_README.xlsx for column definitions
tree_vars_full <- readRDS(file.path(data_dir, "tree_mortality_variables/all_variables_v4.rds"))

#analysis CRS: NAD83 / New York Long Island in meters (tree coordinates arrive in EPSG:2263, US survey feet)
crs_m <- 32118

#load in nyc boundary polygon and save in the format needed by fmesher
nyc_boundary <- st_read(file.path(data_dir, "nyc_boundary_polygon/nybb.shp")) %>%
  st_union() %>% #combine the different boroughs
  st_transform(., crs = crs_m)
#the borough outline has ~80k shoreline vertices; buffer and simplify it so the mesh isn't driven by coastline detail
#(a 100 m buffer / 50 m simplification left sliver triangles that made INLA fail at max.edge 100-200 m)
nyc.bdry <- nyc_boundary %>% st_buffer(500) %>% st_simplify(dTolerance = 150) %>% fm_as_segm()

#load in basemap for plotting
nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 

### parse landuse ###################################################

  # Convert LandUse from integer to character for unique labels
  appendCharLU <- function(lu){
    if (is.na(lu)){
      lu <- paste0("LU_", lu)
    } else {
      if(nchar(as.character(lu)) < 2){
        lu <- paste0("0", lu)
      }
      lu <- paste0("LU_", lu)
    }
    return(lu)
  }
  
  # filter canopy change and relabel as canopy_endstate, relevel to LandUse_char
  tree_vars_full <- tree_vars_full %>% 
    mutate(LandUse_char = sapply(LandUse, appendCharLU))
  tree_vars_full$LandUse_char <- relevel(as.factor(tree_vars_full$LandUse_char), "LU_01") # this is Mike's recommendation for base case, may need to collapse land use classes
  
  # Collapsed land use
  tree_vars_full$LandUse_char_collapsed <- tree_vars_full$LandUse_char %>% as.character()
  x <- tree_vars_full$LandUse_char_collapsed %>%
    case_match( "LU_01" ~ "LowDensityResidential",
                c("LU_02", "LU_03") ~ "HigherDensityResidential",
                c("LU_04", "LU_05") ~ "CommercialMix",
                c("LU_06", "LU_07", "LU_10") ~ "IndustryTransport",
                "LU_08" ~ "PublicInst",
                "LU_09" ~ "OpenOutdoorRec",
                "LU_11" ~ "Vacant")
  # trees more than 10 m from any tax parcel get their own class instead of the nearest lot's land use
  pluto_dist_max_ft <- 10 / 0.3048006 # pluto_dist is in US survey feet
  x[tree_vars_full$pluto_dist > pluto_dist_max_ft] <- "NotNearParcel"
  tree_vars_full$LandUse_char_collapsed <- as.factor(x)
  tree_vars_full$LandUse_char_collapsed <- relevel(as.factor(tree_vars_full$LandUse_char_collapsed), "LowDensityResidential") # this is Mike's recommendation for base case, may need to collapse land use classes


  summary(tree_vars_full$LandUse_char_collapsed)
### create dataframe for analysis ##############################################

#record the number of trees remaining after each filter for the SI sample-flow table
sample_flow <- tibble(step = character(), n = integer())
log_n <- function(df, step) {
  force(df) #evaluate the earlier pipeline steps (and their log_n calls) before appending this row
  sample_flow <<- add_row(sample_flow, step = step, n = nrow(df))
  df
}

tree_vars_full_nona_format <- tree_vars_full %>%
  log_n("2015 census records") %>%
  filter(status == "Alive") %>%
  log_n("alive in 2015 census") %>%
  #create derived variables and filter out rows
  mutate(mort = case_when( canopy_change == 1 ~ 0,      #1 = no change; assumed survival
                           canopy_change == 2 ~ NA,     #2 = canopy gain, not to include in this analysis
                           canopy_change == 3 ~ 1)) %>% #3 = canopy loss; assumed mortality (death or removal)
  filter(!is.na(mort)) %>% #NA canopy_change = no canopy at the stem in 2017 or 2021
  log_n("canopy loss or no change at stem (2017-2021)") %>%
  mutate(dbh_cm = tree_dbh * 2.54) %>%  #convert dbh from inches to cm
  filter(dbh_cm > 7.62) %>%  # remove small trees; 7.6 cm is 3 inches listed in Bigelow
  filter(dbh_cm < 300) %>%  #removing trees with a DBH greater than the maximum ever recorded in NYC
  log_n("DBH between 7.62 and 300 cm") %>%
  filter(!is.na(LandUse_char_collapsed)) %>% #missing land use would otherwise be coded 0 in every LU dummy
  log_n("land use known (or > 10 m from a parcel)") %>%
  mutate(
         #BldgClass_fac = factor(BldgClass), 
         #BldgClass_Group_fac = factor(BldgClass_Group),
         in_sandy_zone_bool = as.logical(in_sandy_zone),
         is_B_cons_bool = as.logical(is_B_cons_10m_17_21), #construction during the 2017-2021 lidar interval
         is_DM_bool = as.logical(is_DM_10m_17_21),
         is_B_cons_1317_bool = as.logical(is_B_cons_10m_13_17), #construction 2013-2017, before the lidar interval (lagged damage)
         is_DM_1317_bool = as.logical(is_DM_10m_13_17),
         is_LU_lowdensres = case_when(LandUse_char_collapsed == "LowDensityResidential" ~ 1, .default = 0),
         is_LU_hidensres = case_when(LandUse_char_collapsed == "HigherDensityResidential" ~ 1, .default = 0),
         is_LU_comm = case_when(LandUse_char_collapsed == "CommercialMix" ~ 1, .default = 0),
         is_LU_indust = case_when(LandUse_char_collapsed == "IndustryTransport" ~ 1, .default = 0),
         is_LU_publicinst = case_when(LandUse_char_collapsed == "PublicInst" ~ 1, .default = 0),
         is_LU_outdoorrec = case_when(LandUse_char_collapsed == "OpenOutdoorRec" ~ 1, .default = 0),
         is_LU_vacant = case_when(LandUse_char_collapsed == "Vacant" ~ 1, .default = 0),
         is_LU_noparcel = case_when(LandUse_char_collapsed == "NotNearParcel" ~ 1, .default = 0),
         imp_5_2017 = BD_5m_2017 + RD_5m_2017 + OI_5m_2017 + RR_5m_2017,
         imp_10_2017 = BD_10m_2017 + RD_10m_2017 + OI_10m_2017 + RR_10m_2017, #impervious = building + road + other impervious + railroad (railroad ~0)
         imp_15_2017 = BD_15m_2017 + RD_15m_2017 + OI_15m_2017 + RR_15m_2017,
         pop_100m_log = log1p(pop_100m), #population within 100 m is right-skewed
         # impervious as % of visible (non-canopy) ground; top-down land cover hides surfaces under the crown, so imp_10_2017 tracks crown size
         vis_10_2017 = imp_10_2017 + GS_10m_2017 + SO_10m_2017 + WA_10m_2017,
         vis_15_2017 = imp_15_2017 + GS_15m_2017 + SO_15m_2017 + WA_15m_2017,
         imp_vis_10_2017 = case_when(vis_10_2017 >= 1 ~ 100 * imp_10_2017 / vis_10_2017, # use 10 m when at least 1% of the buffer is visible
                                     vis_15_2017 >= 1 ~ 100 * imp_15_2017 / vis_15_2017, # otherwise fall back to 15 m (r = 0.88 with 10 m)
                                     .default = NA), # ~180 trees fully under canopy at 15 m are dropped
         stewardship = case_when(steward == "None" ~ 0,
                                is.na(steward) ~ 0, #assuming that NA means no signs of stewardship
                                 steward == "1or2" ~ 1,
                                 steward == "3or4" ~ 1,
                                 steward == "4orMore" ~ 1)) %>%
  select(tree_id, geom,  #only include the variables that will be used in the analysis (to avoid any extra NA values)
         genus, species, 
         mort, # response: 1 = canopy loss at the stem 2017-2021, 0 = no change
         dbh_cm, 
         stewardship, # whether there were any signs of stewardship observed
         #LandUse_char_collapsed, #BldgClass_fac, BldgClass_Group_fac, LandUse_char,  # building type (detailed), building type (high level), land use
         is_LU_lowdensres, is_LU_hidensres, is_LU_comm, is_LU_indust, is_LU_publicinst, is_LU_outdoorrec, is_LU_vacant, is_LU_noparcel,
         in_sandy_zone_bool, # sandy inundation zone
         is_B_cons_bool, is_DM_bool, # building construction, building demolition (2017-2021); street construction not used (data problems)
         is_B_cons_1317_bool, is_DM_1317_bool, # same, 2013-2017
         summer_max, summer_min, # temperature: daytime max and nighttime min (r = 0.11); mean and heat-day counts omitted
        # TC_10m_diff, GS_10m_diff, SO_10m_diff, WA_10m_diff, BD_10m_diff, RD_10m_diff, OI_10m_diff, RR_10m_diff, # land cover diff, doing 10 m
        imp_5_2017, GS_5m_2017, BD_5m_2017, TC_5m_2017, 
        imp_10_2017, GS_10m_2017, BD_10m_2017,  TC_10m_2017,#ORD_10_2017, I_10_2017, RR_10m_2017, # land cover 2017, doing 10 m #WA_10m_2017, 
        imp_15_2017, GS_15m_2017, BD_15m_2017, TC_15m_2017,
        imp_vis_10_2017, # impervious share of visible ground (model term)
        near_NYCHA_10m, # descriptive only; not a model term
        near_park_10m,
        pop_100m_log, # population density
        TC_100m_2017, # neighborhood tree canopy (not the focal crown, which dominates 5-15 m buffers)
        # TC_10m_2021, GS_10m_2021, SO_10m_2021, WA_10m_2021, BD_10m_2021, RD_10m_2021, OI_10m_2021, RR_10m_2021, # land cover 2021, doing 10 m
         RPL_THEME1, RPL_THEME3) %>% # SVI socioeconomic status and racial/ethnic minority status themes; RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME3, RPL_THEME4, RPL_THEMES) %>%
            # One big thing to note is that SVI is not available for parks and certain public land uses, so will limit applicability in certain parts of the city
    separate_wider_delim(., cols = geom, delim = ",", names = c( "x", "y")) %>%  #extract coordinates from geom text string 
    mutate(y = readr::parse_number(y), #crs = 2263 
           x = readr::parse_number(x)) %>% 
    filter(!is.na(imp_vis_10_2017)) %>% 
    log_n("visible ground within 15 m (impervious share defined)") %>% 
    drop_na() %>%  #remove any rows that have an NA value (in practice, missing SVI)
    log_n("complete covariates (SVI theme 1 known)")

  #convert coordinates from US survey feet (EPSG:2263) to meters (EPSG:32118)
  xy_m <- st_as_sf(tree_vars_full_nona_format, coords = c("x", "y"), crs = 2263) %>% st_transform(crs_m) %>% st_coordinates()
  tree_vars_full_nona_format <- tree_vars_full_nona_format %>% mutate(x = xy_m[, "X"], y = xy_m[, "Y"])
  
  #spatial cross-validation folds: 2 km x 2 km blocks, randomly assigned (fixed seed) to 5 folds; the same assignment for every taxon
  cv_block_m <- 2000
  cv_n_folds <- 5
  tree_vars_full_nona_format <- tree_vars_full_nona_format %>% 
    mutate(cv_block = paste(floor(x / cv_block_m), floor(y / cv_block_m)))
  cv_blocks <- sort(unique(tree_vars_full_nona_format$cv_block))
  set.seed(20261009)
  cv_fold_of_block <- setNames(sample(rep_len(1:cv_n_folds, length(cv_blocks))), cv_blocks)
  tree_vars_full_nona_format$cv_fold <- unname(cv_fold_of_block[tree_vars_full_nona_format$cv_block])
  
  sample_flow <- sample_flow %>% mutate(n_removed = lag(n) - n)
  sample_flow
  write_csv(sample_flow, "output/SI_sample_flow.csv")
    
### including the list of species with enough individuals to analyze ##################
  cut_off_n <- 5000 #cut off number of individuals
  
  #species-level models: species (not genus-only labels such as "Prunus" or "Acer") with more than cut_off_n trees
  species_n <- 
    tree_vars_full_nona_format %>% 
    filter(grepl(" ", species)) %>% #genus-only census labels have no species epithet
    count(species) %>% 
    filter(n > cut_off_n)
  
  #genus-level models: trees not in a species model, pooled by genus, where that remainder has more than cut_off_n trees
  #(includes genus-only labels, e.g. the census "Prunus" = cherry); everything else is "other"
  genus_n <- 
    tree_vars_full_nona_format %>% 
    filter(!species %in% species_n$species) %>% 
    count(genus) %>% 
    filter(n > cut_off_n)
  
  tree_vars_full_nona_format <- tree_vars_full_nona_format %>% 
    mutate(sp_a = case_when(species %in% species_n$species ~ species,
                            genus %in% genus_n$genus & genus %in% sub(" .*", "", species_n$species) ~ paste(genus, "(other spp.)"), #genus also has species models
                            genus %in% genus_n$genus ~ paste(genus, "spp."),
                            TRUE ~ "other")) %>% 
    add_count(sp_a, name = "n") #add the n for the sp_a back to the main dataframe
  
  tree_vars_full_nona_format %>% count(sp_a, sort = TRUE) %>% print(n = Inf)
  
### create species level variables ###############################################
  
  # #create bins for tree DBH
  # tree_vars_full_nona_format <-
  #   tree_vars_full_nona_format %>% 
  #   group_by(sp_a) %>% 
  #   mutate(dbh_quintile = as.factor(ntile(dbh_cm, 5)), #calculate survival per quintile or quartile per species
  #          dbh_decile = as.factor(ntile(dbh_cm, 10))) %>%  
  #   ungroup() 
  #str(tree_vars_full_nona_format$dbh_quintile)
  
  #save scale parameters for numeric variables, pooled across all taxa so that 1 SD is the same raw change for every species
  scale_vars <- c("summer_max", "summer_min", "imp_vis_10_2017", "TC_100m_2017", "pop_100m_log", "RPL_THEME1", "RPL_THEME3")
  scale_params <-
    tree_vars_full_nona_format %>% 
    select(all_of(scale_vars)) %>% 
    summarise(across(everything(),
                     list(mean = ~mean(.x, na.rm=TRUE),
                          sd   = ~sd(.x,   na.rm=TRUE))))
  
  # Apply Z-scaling (pooled)
  tree_vars_full_nona_format <- tree_vars_full_nona_format %>% 
    mutate(across(all_of(scale_vars),
                  ~ (.x - mean(.x, na.rm=TRUE)) / sd(.x, na.rm=TRUE),
                  .names = "{.col}_z"))
  
### some last pre-analysis data exploration ########################################
  tree_vars_full_nona_format %>% sample_n(10000) %>%
    ggplot(aes(x = x, y = y, z = is_LU_vacant)) + stat_summary_hex(fun = "mean", binwidth = 1000) + theme_bw() + scale_fill_viridis_c()

  # #map of tree samples across NYC
  # tree_vars_full_nona_format %>% 
  #   filter(sp_a == "Quercus palustris") %>% #sample_n(100000) %>%
  #   ggplot(aes(x = x, y = y)) + geom_hex(binwidth = 5000) + theme_bw() + scale_fill_viridis_c()


  tree_vars_full_nona_format %>% sample_n(10000) %>%
    ggplot(aes(x = TC_5m_2017, y = TC_15m_2017)) + geom_point()

  cor(tree_vars_full_nona_format$TC_5m_2017, tree_vars_full_nona_format$TC_15m_2017)
  numeric_df <- tree_vars_full_nona_format %>% select(where(~ is.numeric(.x) | is.logical(.x)), -tree_id)
  cor_matrix <- cor(numeric_df, use = "complete.obs")
  
  
### table 1: summary of mortality by focal species  ####################################
  
  table_1 <- tree_vars_full_nona_format %>% 
    group_by(sp_a) %>% 
    summarize(n = n(),
              median_dbh = median(dbh_cm),
              mort_4yr = mean(mort),
              annual_mort = 1 - (1 - mort_4yr) ^ (1/4))  #annual mortality
  
  table_1_all_trees <- tree_vars_full_nona_format %>% 
    summarize(sp_a = "all trees", 
              n = n(),
              median_dbh = median(dbh_cm),
              mort_4yr = mean(mort),
              annual_mort = 1 - (1 - mort_4yr) ^ (1/4))
  
  table_1 <- rbind(table_1, table_1_all_trees) %>% 
    mutate(mort_4yr = mort_4yr * 100, annual_mort = annual_mort * 100)
  
  #add common names
  common_name_lookup <- read_csv(file.path(data_dir, "tree_mortality_variables/common_name_lookup.csv")) %>%
    rename(sp_a = species)
  table_1 <- left_join(table_1, common_name_lookup) %>% 
    select(sp_a, common_name, n, median_dbh, mort_4yr, annual_mort)
  
  table_1  %>% ungroup() %>% 
    gt() %>% 
    fmt_number(columns  = c(median_dbh, mort_4yr, annual_mort), decimals = 1) |>
    fmt_integer(columns = n, sep_mark = ",") %>% 
    cols_label(
      sp_a = "species",
      common_name = "common name",
      n = "n",
      median_dbh = "median DBH (cm)",
      mort_4yr = "4-year mortality (%)",
      annual_mort = "annual mortality (%)" ) |>
    tab_style(style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2),
      style = "solid"),
      locations = cells_column_labels()) 
  #%>%  gtsave( paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/table1.docx"))  
  
  
  
### Fig 2: mortality as a function of DBH #########################################
  fig2 <- tree_vars_full_nona_format %>% 
  select(sp_a, mort, dbh_cm) %>% 
  group_by(sp_a) %>% 
  mutate(quartile = ntile(dbh_cm, 4)) %>%  #calculate mortality per DBH quartile per species
  ungroup() %>% 
  group_by(sp_a, quartile) %>% 
  summarize(median_dbh = median(dbh_cm),
            mort_4yr = mean(mort),
            annual_mort = 100 * (1 - (1 - mort_4yr)^(1/4)),
            n = n()) %>% 
  ggplot(aes(x = median_dbh, y = annual_mort, color = sp_a)) + geom_point() + geom_line()+ theme_bw() + xlab("DBH (cm)") + ylab("annual mortality (%)") +
  scale_color_viridis_d(option = "turbo", name = "species") +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.text = element_text(face = "italic"))      

  # ggsave(paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/fig_2.png"),
  #        width = 7, height = 5, units = "in", dpi = 400)


  
### Fig 3: maps of observed mortality for each taxon #########################################
  #load in basemap
  nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
  nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
  
  #function to remove bins with lower than a certain n
  bin_removal <- function(x) {
    if(length(x) < 10) return(NA) 
    return((1 - (1 - mean(x)) ^ (1/4)) * 100) #annual mortality
  }
  
  #project the data to the crs of the basemap tile (3857)
  tree_vars_full_nona_format_sf_3857 <- st_as_sf(tree_vars_full_nona_format, crs = crs_m, coords = c("x", "y")) %>% 
    st_transform(., crs = 3857) %>% 
    bind_cols(st_coordinates(.) %>% as.data.frame())
  
  #create figure
  fig3_mort_map <- 
    ggplot() + ggthemes::theme_few() +   
    geom_spatraster_rgb(data = nyc_topo_spatrast) +
    stat_summary_hex(data = tree_vars_full_nona_format_sf_3857, aes(x = X, y = Y, z =  mort), fun = bin_removal, bins = 20) +
    facet_wrap(~sp_a, ncol = 6) +
    scale_fill_viridis_c(option = "turbo", name = "annual \nmortality (%)  ", 
                         limits = c(0, 6),
                         guide = guide_colorbar(barwidth = 10),
                         oob    = squish) + 
    xlab("") + ylab("") + 
    theme(strip.text = element_text(face = "italic"),
          legend.position = "bottom",
          #legend.background = element_rect(fill = "white", color = "grey80"),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # ggsave(fig3_mort_map, filename = paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/Fig3_mortalitymap.png"),
  #        width = 10, height = 6.5, units = "in", dpi = 400)
  

  
### Spatial mortality model with INLA #####################################
sp_list <- unique(tree_vars_full_nona_format$sp_a) %>% sort()
sp_df <- data.frame(sp = sp_list, model = as.character(1:length(sp_list))) 

## Run settings (environment variables; defaults = the main analysis) ------------------------------
  # CV_FOLD: 0 = fit all trees; 1-5 = hold out that spatial fold and save predictions for it (cross-validation)
  # MESH_EDGE: inner mesh max.edge in m; INLA_THREADS: INLA threads per fit (lower when running fits in parallel)
  cv_holdout <- as.integer(Sys.getenv("CV_FOLD", "0"))
  cv_mode <- cv_holdout > 0
  mesh_edge <- as.numeric(Sys.getenv("MESH_EDGE", "100"))
  inla_threads <- Sys.getenv("INLA_THREADS", "8:1")

## Mesh and SPDE (Matern with PC priors), shared by all taxa; all distances in meters ------------------------------
  mesh <- fm_mesh_2d_inla(
    boundary = nyc.bdry, #loaded at start of script
    crs      = fm_crs(crs_m), #same CRS as the tree points and the NYC boundary (needed e.g. to mask the field-map pixel grid)
    max.edge = c(mesh_edge, 2000), #100 m inside the city (~range/5 vs the 500 m range prior; fitted ranges ~150-450 m); coarse in the outer extension
    offset   = c(-0.01, 3000), #negative = relative to boundary; outer extension avoids edge effects
    cutoff   = 50
  )
  cat("Mesh has", mesh$n, "vertices\n")
  #degenerate (zero-area) triangles give non-finite FEM matrices, and INLA then aborts with an unhelpful error
  stopifnot("mesh has degenerate triangles; adjust the boundary buffer/simplification" = all(is.finite(fm_fem(mesh)$g1@x)))
  #ggplot() + gg(mesh) + theme_minimal() #plot mesh to double check it
  
  # PC priors: P(range < 500 m) = 0.05 ; P(sigma > 1) = 0.05 (see ANALYSIS_PLAN.md section 4)
  spde <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(500, 0.05),
    prior.sigma = c(1, 0.05)
  )

  lidar_interval_yr <- 4 #years between the 2017 and 2021 lidar flights; enters as an exposure offset

#each taxon's results are saved as a small summary file (coefficients, fitted values, spatial field map);
#if a taxon's summary file already exists the fit is skipped, so an interrupted run can be resumed.
#full fits on the 100 m mesh are 4-6 GB each, so they are only saved if requested
  fit_dir <- Sys.getenv("FIT_DIR", "output/fits_100m")
  dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
  fit_file <- function(sp, type) file.path(fit_dir, paste0(gsub("[^A-Za-z]+", "_", sp), "_", type, ".rds"))
  save_full_fit <- FALSE

#pixel grid over NYC for the spatial field maps: the posterior mean and SD at the mesh vertices are interpolated to each pixel
#(the mean is exact; the interpolated SD is an approximation that is adequate for mapping)
  field_grid <- fm_pixels(mesh, dims = c(300, 300), mask = nyc_boundary, format = "sf")
  field_xy <- st_coordinates(field_grid)
  field_A <- fm_evaluator(mesh, loc = field_xy)

#create empty lists to save output
  mort_model_list <- vector("list", length(sp_list))
  field_list <- vector("list", length(sp_list))
  trees_model_list <- vector("list", length(sp_list))
  roc_list <- vector("list", length(sp_list))

#run model for focal species; set the TAXA_TO_FIT environment variable (semicolon-separated) to fit only some taxa
taxa_to_fit <- if (nzchar(Sys.getenv("TAXA_TO_FIT"))) strsplit(Sys.getenv("TAXA_TO_FIT"), ";")[[1]] else sp_list
for (i in which(sp_list %in% taxa_to_fit)){
#for (i in 1:2){
 
  sp_focal <- sp_list[i] #sp_focal <- sp_list[6]
  print(paste(i, sp_focal))
  
  if (cv_mode && file.exists(fit_file(sp_focal, paste0("cv_fold", cv_holdout)))) next #this taxon and fold already done
  
  if (file.exists(fit_file(sp_focal, "summary"))) { #already fitted: reload the summary and move on
    fit_summary <- readRDS(fit_file(sp_focal, "summary"))
    mort_model_list[[i]] <- fit_summary$model
    trees_model_list[[i]] <- fit_summary$trees
    roc_list[[i]] <- fit_summary$roc
    field_list[[i]] <- fit_summary$field
    next
  }
  
  sp_sub_format <- tree_vars_full_nona_format %>% filter(sp_a == sp_focal)
    
  trees <- sp_sub_format %>% 
    mutate(log_exposure = log(lidar_interval_yr)) %>% 
    select(tree_id, mort, log_exposure, stewardship, is_B_cons_bool, is_DM_bool,
           is_B_cons_1317_bool, is_DM_1317_bool,
           is_LU_hidensres, is_LU_comm, is_LU_indust, is_LU_publicinst, is_LU_outdoorrec, is_LU_vacant, is_LU_noparcel, #LowDensityResidential is the reference
           dbh_cm, summer_max_z, summer_min_z, imp_vis_10_2017_z, TC_100m_2017_z, pop_100m_log_z, RPL_THEME1_z, RPL_THEME3_z, in_sandy_zone_bool,
           near_park_10m, #NYCHA proximity removed: too few flagged trees per taxon
           cv_fold, x, y) %>% 
    mutate(across(where(is.logical), as.numeric))
  if (cv_mode) trees <- trees %>% mutate(mort_obs = mort, mort = ifelse(cv_fold == cv_holdout, NA, mort)) #held-out trees: response unknown
  trees_sf <- st_as_sf( trees, coords = c("x", "y"), crs = crs_m, remove = FALSE) 
  
  #1D mesh for the DBH smooth (cm), spanning this taxon's range
  dbh_mesh <- fm_mesh_1d(seq(min(trees$dbh_cm), max(trees$dbh_cm), length.out = 20))
  
  ## ---- Model components and formula ---------------------------------------
  # linear terms get an N(0, 1) prior (mean.linear = 0, prec.linear = 1)
  cmp <- mort ~ Intercept(1) +
    exposure(log_exposure, model = "offset") + #intercept is then a log annual hazard
    stewardship(stewardship, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_B_cons_bool(is_B_cons_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_DM_bool(is_DM_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_B_cons_1317_bool(is_B_cons_1317_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_DM_1317_bool(is_DM_1317_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    in_sandy_zone_bool(in_sandy_zone_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    RPL_THEME1_z(RPL_THEME1_z, model = "linear", mean.linear = 0, prec.linear = 1) +
    RPL_THEME3_z(RPL_THEME3_z, model = "linear", mean.linear = 0, prec.linear = 1) +
    summer_max_z(summer_max_z, model = "linear", mean.linear = 0, prec.linear = 1) +
    summer_min_z(summer_min_z, model = "linear", mean.linear = 0, prec.linear = 1) +
    near_park_10m(near_park_10m, model = "linear", mean.linear = 0, prec.linear = 1) +

    # land use: LowDensityResidential is the reference level, so it has no term
    is_LU_hidensres(is_LU_hidensres, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_LU_comm(is_LU_comm, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_LU_indust(is_LU_indust, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_LU_publicinst(is_LU_publicinst, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_LU_outdoorrec(is_LU_outdoorrec, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_LU_vacant(is_LU_vacant, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_LU_noparcel(is_LU_noparcel, model = "linear", mean.linear = 0, prec.linear = 1) +
    
    imp_vis_10_2017_z(imp_vis_10_2017_z, model = "linear", mean.linear = 0, prec.linear = 1) +
    # grass/shrub omitted: it is effectively the complement of imp_vis_10_2017 (r = -0.78)
    TC_100m_2017_z(TC_100m_2017_z, model = "linear", mean.linear = 0, prec.linear = 1) +
    pop_100m_log_z(pop_100m_log_z, model = "linear", mean.linear = 0, prec.linear = 1) +
    
    dbh_smooth(
      dbh_cm,
      model = "rw2",
      mapper = bru_mapper(dbh_mesh),
      scale.model = TRUE,                     # makes precision comparable across resolutions
      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))) ) + # PC prior on smoothness
    spatial_field(geometry, model = spde)
  
  ## ---- Fit with bru() -------------------------------------------------------
  lik <- bru_obs(
    formula = mort ~ .,
    family  = "binomial",
    data    = trees_sf,
    Ntrials = 1,
    control.family = list(link = "cloglog")
  )
  
  fit <- bru(
    cmp,
    lik,
    options = bru_options(
      control.compute = list(dic = !cv_mode, waic = !cv_mode, cpo = !cv_mode,
                             config = TRUE), #always on: bru() fails without it, and it is needed for the posterior draws below
      control.inla = list(int.strategy = if (cv_mode) "eb" else "auto"), #CV fits: single-point integration (same estimates, ~1/3 the memory)
      verbose = FALSE, num.threads = inla_threads #default 8 threads: same estimates as 16, ~30% less memory, slightly faster
    )
  )
  
  #summary(fit)
  
  ### saving data from focal species
  fit$sp <- sp_focal
  
  if (cv_mode) { #cross-validation: save predictions for the held-out trees, then move on to the next taxon
    held <- which(is.na(trees$mort))
    eta <- fit$summary.linear.predictor[held, c("mean", "sd")]
    z <- qnorm(ppoints(50)) #average the 4-year loss probability over the posterior of the linear predictor
    p_held <- rowMeans(1 - exp(-exp(eta$mean + outer(eta$sd, z))))
    cv_pred <- tibble(sp = sp_focal, fold = cv_holdout, tree_id = trees$tree_id[held], mort = trees$mort_obs[held],
                      eta_mean = eta$mean, eta_sd = eta$sd, p = p_held,
                      p_null = mean(trees$mort, na.rm = TRUE)) #baseline: this taxon's loss rate in the training folds
    saveRDS(list(pred = cv_pred, summary.fixed = fit$summary.fixed, summary.hyperpar = fit$summary.hyperpar,
                 n_train = sum(!is.na(trees$mort)), n_held = length(held)),
            fit_file(sp_focal, paste0("cv_fold", cv_holdout)))
    rm(fit, lik, trees_sf); gc()
    next
  }
  
  if (save_full_fit) saveRDS(fit, fit_file(sp_focal, "fit")) #full fit, needed for predict()
  
  #for attribution analyses: joint posterior draws of the fixed effects (rows = terms, columns = draws)
  fixed_names <- rownames(fit$summary.fixed)
  fixed_draws <- inla.posterior.sample(n = 1000, result = fit, num.threads = inla_threads,
                                       selection = setNames(as.list(rep(1, length(fixed_names))), fixed_names))
  fixed_draws <- sapply(fixed_draws, function(d) d$latent[, 1])
  rownames(fixed_draws) <- fixed_names

  field_nodes <- fit$summary.random$spatial_field #posterior summary at each mesh vertex
  field_list[[i]] <- tibble(sp = sp_focal, x = field_xy[, "X"], y = field_xy[, "Y"],
                            mean = fm_evaluate(field_A, field = field_nodes$mean),
                            sd   = fm_evaluate(field_A, field = field_nodes$sd))

  mort_model_list[[i]] <- list(sp = sp_focal, #small version of the fit used by the figures below
                               summary.fixed = fit$summary.fixed,
                               summary.hyperpar = fit$summary.hyperpar,
                               summary.random = list(dbh_smooth = fit$summary.random$dbh_smooth),
                               waic = fit$waic[c("waic", "p.eff")],
                               dic = fit$dic[c("dic", "p.eff")],
                               #for attribution: joint posterior draws and correlation of the fixed effects; each tree's linear
                               #predictor (log 4-year hazard, including the exposure offset) is saved with the tree data below
                               fixed_draws = fixed_draws,
                               fixed_correlation = cor(t(fixed_draws)))
  
  trees$fitted_p <- fit$summary.fitted.values$mean[1:nrow(trees)] #predicted probability of canopy loss over the 4-year interval
  trees$eta_mean <- fit$summary.linear.predictor$mean[1:nrow(trees)] #linear predictor (cloglog scale), for attribution analyses
  trees$eta_sd   <- fit$summary.linear.predictor$sd[1:nrow(trees)]
  trees$resid    <- trees$mort - trees$fitted_p
  trees$sp <- sp_focal
  trees_join <- left_join(trees, sp_sub_format %>% select(tree_id, sp_a, species, genus), by = "tree_id")
  
  trees_model_list[[i]] <- trees_join
  
  roc_curve <- roc(trees$mort, trees$fitted_p)
  roc_curve$sp <- sp_focal
  roc_list[[i]] <- roc_curve
  
  saveRDS(list(model = mort_model_list[[i]], trees = trees_join, roc = roc_curve, field = field_list[[i]]), fit_file(sp_focal, "summary"))
  rm(fit, lik, trees_sf); gc() #free the full fit before the next taxon
  
} #end species loop model run
  
### SI x: ROC and binned residuals #######################################
 
  #ROC values per species
  #roc_list <- roc_list[1:5] #truncating for testing
  
  # summary(roc_list[[2]])
  # plot(roc_list[[2]])

  roc_values <- map_dbl(roc_list, auc)
  #roc_df <- data.frame(sp = sp_list[1:5], roc = roc_values)
  roc_df <- data.frame(sp = sp_list, roc = roc_values)
  
  
  # CHANGE TO DOING THIS ACROSS ALL SP
  # ## ---- Posterior predictive / binned residual check ------------------------
  # 
  # # Binned residual plot: mean residual within bins of fitted probability
  # n_bins <- 10
  # trees$bin <- cut(trees$fitted_p, breaks = quantile(trees$fitted_p,
  #                                                    probs = seq(0, 1, length.out = n_bins + 1)),
  #                  include.lowest = TRUE)
  # 
  # binned <- trees %>%
  #   group_by(bin) %>%
  #   summarise(
  #     mean_fitted = mean(fitted_p),
  #     mean_obs    = mean(mort),
  #     n           = n(),
  #     se          = sqrt(mean_fitted * (1 - mean_fitted) / n)
  #   )
  # 
  # p_binned <- ggplot(binned, aes(x = mean_fitted, y = mean_obs)) +
  #   geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  #   geom_errorbar(aes(ymin = mean_obs - 1.96 * se, ymax = mean_obs + 1.96 * se),
  #                 width = 0.01, colour = "grey60") +
  #   geom_point(size = 2, colour = "steelblue") +
  #   labs(title = "Binned residual plot",
  #        x = "Mean predicted survival probability",
  #        y = "Observed proportion surviving") +
  #   theme_minimal(base_size = 11)
  # 
  # p_binned
   
### Figs 4-6 and S1: forest plots of fixed effects (Fig 4 cross-taxa overview; Figs 5, 6, S1 by taxon) ########
  # term labels and display order, grouped by construct; continuous terms are per 1 SD (pooled across taxa)
  term_labels <- c(
    stewardship         = "Stewardship (any signs)",
    is_B_cons_1317_bool = "Building construction 2013-17",
    is_B_cons_bool      = "Building construction 2017-21",
    is_DM_1317_bool     = "Demolition 2013-17",
    is_DM_bool          = "Demolition 2017-21",
    imp_vis_10_2017_z   = "Impervious share of visible ground (+1 SD)",
    TC_100m_2017_z      = "Tree canopy within 100 m (+1 SD)",
    summer_max_z        = "Summer max temperature (+1 SD)",
    summer_min_z        = "Summer min temperature (+1 SD)",
    in_sandy_zone_bool  = "Sandy inundation zone",
    near_park_10m       = "Within 10 m of a park",
    is_LU_hidensres     = "Land use: higher-density residential",
    is_LU_comm          = "Land use: commercial / mixed",
    is_LU_indust        = "Land use: industrial / transport",
    is_LU_publicinst    = "Land use: public / institutional",
    is_LU_outdoorrec    = "Land use: open space / recreation",
    is_LU_vacant        = "Land use: vacant",
    is_LU_noparcel      = "Land use: not near a tax parcel",
    pop_100m_log_z      = "Population within 100 m (+1 SD, log)",
    RPL_THEME1_z        = "SVI socioeconomic theme (+1 SD)",
    RPL_THEME3_z        = "SVI racial & ethnic minority theme (+1 SD)")
  
  coef_df <- map_dfr(
    mort_model_list,
    ~ as_tibble(.x$summary.fixed, rownames = "term") %>% mutate(sp = .x$sp)
  ) %>%
    rename(
      estimate = mean,
      lower    = `0.025quant`,
      upper    = `0.975quant`
    ) %>%
    select(sp, term, estimate, lower, upper, sd) %>% 
    filter(term != "Intercept")
  
  # random-effects meta-analysis across taxa for each term (see ANALYSIS_PLAN.md section 5):
  # each taxon's posterior mean and SD (log hazard ratio) are treated as an estimate and its standard error;
  # REML for the between-taxon variance, Knapp-Hartung adjusted CI (only ~19 taxa); "other" is included
  meta_fit <- function(d) {
    m <- rma(yi = d$estimate, sei = d$sd, method = "REML", test = "knha")
    pr <- predict(m)
    tibble(estimate = as.numeric(m$b), lower = m$ci.lb, upper = m$ci.ub,
           tau = sqrt(m$tau2), I2 = m$I2, pi_lower = pr$pi.lb, pi_upper = pr$pi.ub, k = m$k)
  }
  meta_df <- coef_df %>% group_by(term) %>% group_modify(~ meta_fit(.x)) %>% ungroup() %>% mutate(sp = "All taxa (meta-analysis)")
  meta_df
  
  # SI sensitivity check: related taxa may respond alike, so add a genus random effect above taxon
  meta_genus_fit <- function(d) {
    m <- rma.mv(yi = estimate, V = sd^2, random = ~ 1 | genus/sp, data = d, method = "REML", test = "t")
    tibble(estimate_genus = as.numeric(m$b), lower_genus = m$ci.lb, upper_genus = m$ci.ub,
           sd_genus = sqrt(m$sigma2[1]), sd_taxon_within_genus = sqrt(m$sigma2[2]))
  }
  meta_genus_df <- coef_df %>% mutate(genus = sub(" .*", "", sp)) %>% #"other" is its own genus group
    group_by(term) %>% group_modify(~ meta_genus_fit(.x)) %>% ungroup()
  meta_sensitivity <- meta_df %>% select(term, estimate, lower, upper, tau) %>% left_join(meta_genus_df, by = "term") %>% 
    mutate(across(c(estimate, lower, upper, estimate_genus, lower_genus, upper_genus), exp, .names = "HR_{.col}"))
  meta_sensitivity
  write_csv(meta_sensitivity, "output/SI_meta_genus_sensitivity.csv")
  
  # y-axis order (same in every figure): meta-analytic mean at the top, taxa alphabetically (so genera group together), "other" last
  meta_label <- "All taxa (meta-analysis)"
  taxa_alpha <- sort(setdiff(unique(coef_df$sp), "other"))
  sp_levels <- c("other", rev(taxa_alpha), meta_label) #ggplot draws the first level at the bottom
  forest_df <- bind_rows(coef_df %>% mutate(type = "Taxon"), meta_df %>% mutate(type = "Mean across taxa")) %>% 
    mutate(sp = factor(sp, levels = sp_levels),
           type = factor(type, levels = c("Taxon", "Mean across taxa")))
  
  # italicize genus and species names only, e.g. "Acer (other spp.)" -> italic Acer + plain "(other spp.)"
  taxon_labels <- function(x) parse(text = case_when(
    x %in% c("other", meta_label)     ~ paste0("'", x, "'"),
    grepl(" \\(other spp\\.\\)$", x)  ~ paste0("italic('", sub(" .*", "", x), "')~'(other spp.)'"),
    grepl(" spp\\.$", x)              ~ paste0("italic('", sub(" .*", "", x), "')~'spp.'"),
    TRUE                              ~ paste0("italic('", x, "')")))
  
  #readable hazard-ratio ticks within each panel's range: doubling steps on wide log axes,
  #intermediate steps on medium axes, evenly spaced on narrow ones
  hr_breaks <- function(lims) {
    ratio <- lims[2] / lims[1]
    b <- if (ratio > 6) c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16)
         else if (ratio > 1.8) c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6)
         else pretty(lims, n = 4)
    b[b >= lims[1] & b <= lims[2] & b > 0]
  }
  hr_labels <- function(x) format(x, drop0trailing = TRUE, trim = TRUE)
  
  forest_caption <- "Hazard ratio (log scale) with 95% credible interval from one spatial Bernoulli GLMM (cloglog link) per taxon. Mean across taxa: random-effects\nmeta-analysis (REML, Knapp-Hartung 95% CI). Continuous drivers: per 1 SD, pooled across taxa."
  
  # one forest plot for a set of drivers; facets follow the order of `terms`
  plot_forest <- function(terms, title, ncol = 4, shared_x = FALSE, note = NULL) {
    forest_df %>% 
      filter(term %in% terms) %>% 
      mutate(term = factor(term_labels[term], levels = term_labels[terms])) %>% 
      ggplot(aes(x = exp(estimate), y = sp, colour = type, shape = type)) +
      geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
      geom_linerange(aes(xmin = exp(lower), xmax = exp(upper)), linewidth = 0.5) +
      geom_point(size = 1.6) +
      scale_x_log10(breaks = hr_breaks, labels = hr_labels) +
      scale_y_discrete(labels = taxon_labels) +
      scale_colour_manual(values = c("Taxon" = "#2a78d6", "Mean across taxa" = "#eb6834"), name = NULL) +
      scale_shape_manual(values = c("Taxon" = 16, "Mean across taxa" = 18), name = NULL) +
      facet_wrap(~term, ncol = ncol, scales = if (shared_x) "fixed" else "free_x", labeller = label_wrap_gen(width = 30)) +
      labs(title = title, caption = paste(c(forest_caption, note), collapse = "\n"), x = "Hazard ratio", y = NULL) +
      theme_minimal(base_size = 9) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour = "#e6e5e0", linewidth = 0.3),
            strip.text = element_text(face = "bold", hjust = 0),
            panel.spacing.x = unit(1.5, "lines"), #keeps tick labels of neighboring panels from running together
            axis.text.y = element_text(size = 7),
            legend.position = "top", legend.justification = "left",
            plot.caption = element_text(hjust = 0, colour = "#52514e"),
            plot.background = element_rect(fill = "white", colour = NA))
  }
  
  # Fig 5, construction and demolition: rows = period (lagged 2013-17, concurrent 2017-21), columns = building construction, demolition;
  # shared x-axis (all binary). Street construction is not analyzed (data problems that could not be resolved)
  p_forest_construction <- plot_forest(
    c("is_B_cons_1317_bool", "is_DM_1317_bool", "is_B_cons_bool", "is_DM_bool"),
    title = "Figure 5. Effects of nearby building construction and demolition on street tree mortality", ncol = 2, shared_x = TRUE)
  p_forest_construction
  ggsave(p_forest_construction, filename = "output/Fig5_forest_construction.png", width = 7.5, height = 8.5, dpi = 300)
  
  # Fig 6, site, environment, stewardship and social context
  p_forest_site <- plot_forest(
    c("stewardship", "imp_vis_10_2017_z", "TC_100m_2017_z", "in_sandy_zone_bool", "summer_max_z", "summer_min_z", "pop_100m_log_z", "RPL_THEME1_z", "RPL_THEME3_z"),
    title = "Figure 6. Effects of stewardship, site, environment and social context on street tree mortality", ncol = 3)
  p_forest_site
  ggsave(p_forest_site, filename = "output/Fig6_forest_site_env_social.png", width = 9.5, height = 12, dpi = 300)
  
  # Fig S1 (SI): land use and park proximity
  p_forest_landuse <- plot_forest(
    c("is_LU_hidensres", "is_LU_comm", "is_LU_indust", "is_LU_publicinst", "is_LU_outdoorrec", "is_LU_vacant", "is_LU_noparcel", "near_park_10m"),
    title = "Figure S1. Effects of adjacent land use and park proximity on street tree mortality", ncol = 4,
    note = "Land use is relative to low-density (one/two-family) residential.")
  p_forest_landuse
  ggsave(p_forest_landuse, filename = "output/FigS1_forest_landuse_park.png", width = 11, height = 8.5, dpi = 300)
  
  # Fig 4, overview (precedes the by-taxon figures): mean across taxa only, one row per driver
  p_forest_mean <- meta_df %>% 
    mutate(term = factor(term_labels[term], levels = rev(term_labels))) %>% 
    ggplot(aes(x = exp(estimate), y = term)) +
    geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
    geom_linerange(aes(xmin = exp(lower), xmax = exp(upper)), colour = "#eb6834", linewidth = 0.6) +
    geom_point(colour = "#eb6834", shape = 18, size = 2.5) +
    scale_x_log10(breaks = c(0.8, 0.9, 1, 1.25, 1.5, 2, 2.5)) +
    labs(title = "Figure 4. Mean effect of each driver across taxa",
         subtitle = sprintf("Hazard ratio with 95%% CI; random-effects meta-analysis (REML, Knapp-Hartung) of %d taxa", length(unique(coef_df$sp))),
         x = "Hazard ratio (log scale)", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
  
  p_forest_mean
  ggsave(p_forest_mean, filename = "output/Fig4_cross_taxa_means.png", width = 7, height = 5.5, dpi = 300)
  
  write_csv(bind_rows(coef_df, meta_df), "output/forest_coefficients.csv")
  
### SI X: Spatial random field map ---------------------------------------------
  
  #posterior mean of each taxon's spatial field on the pixel grid saved in the loop; log hazard ratio relative to the city-wide
  #baseline after accounting for the drivers (> 0: higher mortality than the drivers explain)
  field_df <- bind_rows(field_list) %>% mutate(sp = factor(sp, levels = c(sort(setdiff(unique(sp), "other")), "other")))

  p_field <- ggplot(field_df) +
    geom_tile(aes(x = x, y = y, fill = mean)) + #tile rather than raster: pixel centers are not exactly evenly spaced
    geom_sf(data = nyc_boundary, fill = NA, colour = "grey40", linewidth = 0.15) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0,
                         limits = c(-2, 2), oob = squish, name = "spatial effect\n(log hazard ratio)") +
    facet_wrap(~sp, ncol = 6) +
    coord_sf(crs = crs_m, datum = NA) +
    labs(title = "Estimated spatial random effect (SPDE) by taxon", x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(face = "italic"), legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = NA))

  p_field
  ggsave(p_field, filename = "output/SI_spatial_field_by_taxon.png", width = 11, height = 9, dpi = 300)
  
### SI: spatial block cross-validation (2 km blocks, 5 folds; see ANALYSIS_PLAN.md section 6) ###############
  # fits are run separately with CV_FOLD = 1-5 (output/run_cv.sh); each saved the predictions for its held-out trees
  cv_dir <- "output/cv"
  cv_pred <- map_dfr(list.files(cv_dir, pattern = "_cv_fold[0-9]+\\.rds$", full.names = TRUE), ~ readRDS(.x)$pred)
  
  log_score <- function(y, p) y * log(p) + (1 - y) * log(1 - p)
  #calibration-in-the-large (intercept with the predictions as an offset; 0 = no overall bias) and calibration slope (1 = ideal),
  #both on the logit scale
  cal_intercept <- function(y, p) unname(coef(glm(y ~ offset(qlogis(p)), family = binomial))[1])
  cal_slope     <- function(y, p) unname(coef(glm(y ~ qlogis(p), family = binomial))[2])
  cv_summarise <- function(d) d %>% summarise(
    n = n(), folds = n_distinct(fold), observed = mean(mort), predicted = mean(p),
    log_score = mean(log_score(mort, p)), log_score_null = mean(log_score(mort, p_null)),
    brier = mean((mort - p)^2), brier_null = mean((mort - p_null)^2),
    auc_cv = as.numeric(pROC::auc(mort, p, quiet = TRUE)),
    cal_intercept = cal_intercept(mort, p), cal_slope = cal_slope(mort, p), .groups = "drop") %>% 
    mutate(log_score_skill = 1 - log_score / log_score_null, brier_skill = 1 - brier / brier_null)
  
  auc_insample <- tibble(sp = map_chr(roc_list, "sp"), auc_insample = map_dbl(roc_list, ~ as.numeric(pROC::auc(.x))))
  cv_metrics <- bind_rows(cv_pred %>% group_by(sp) %>% cv_summarise(),
                          cv_pred %>% cv_summarise() %>% mutate(sp = meta_label)) %>% 
    left_join(auc_insample, by = "sp") %>% 
    mutate(sp = factor(sp, levels = sp_levels)) %>% arrange(desc(sp))
  cv_metrics
  cv_pooled_label <- "All taxa (pooled)" #the CV summary pools held-out trees; it is not a meta-analysis
  write_csv(cv_metrics %>% mutate(sp = ifelse(sp == meta_label, cv_pooled_label, as.character(sp))), "output/SI_cv_metrics.csv")
  
  # calibration: observed vs predicted 4-year loss in deciles of predicted risk, by taxon and pooled
  cv_bins <- bind_rows(cv_pred, cv_pred %>% mutate(sp = meta_label)) %>% 
    group_by(sp) %>% mutate(bin = ntile(p, 10)) %>% 
    group_by(sp, bin) %>% 
    summarise(predicted = mean(p), observed = mean(mort), n = n(), .groups = "drop") %>% 
    mutate(se = sqrt(observed * (1 - observed) / n),
           sp = factor(sp, levels = rev(sp_levels))) #meta label first, then taxa alphabetically, "other" last
  
  p_cv_calibration <- ggplot(cv_bins, aes(x = predicted, y = observed)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_linerange(aes(ymin = pmax(observed - 1.96 * se, 0), ymax = observed + 1.96 * se), colour = "#2a78d6", linewidth = 0.4) +
    geom_point(colour = "#2a78d6", size = 1.2) +
    facet_wrap(~sp, ncol = 6, labeller = as_labeller(function(x) ifelse(x == meta_label, cv_pooled_label, x))) +
    coord_equal(xlim = c(0, 0.4), ylim = c(0, 0.4)) +
    scale_x_continuous(labels = label_percent(accuracy = 1), breaks = c(0, 0.2, 0.4)) +
    scale_y_continuous(labels = label_percent(accuracy = 1), breaks = c(0, 0.2, 0.4)) +
    labs(title = "Calibration of held-out predictions (spatial block cross-validation)",
         subtitle = "Observed vs predicted 4-year canopy loss in deciles of predicted risk; 2 km blocks, 5 folds; line = perfect calibration",
         x = "Predicted 4-year loss", y = "Observed 4-year loss") +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(face = "italic"), panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
  p_cv_calibration
  ggsave(p_cv_calibration, filename = "output/SI_cv_calibration.png", width = 10, height = 8.5, dpi = 300)
  
  # discrimination: in-sample vs cross-validated AUC by taxon
  p_cv_auc <- cv_metrics %>% 
    select(sp, `In-sample` = auc_insample, `Cross-validated` = auc_cv) %>% 
    tidyr::pivot_longer(-sp, names_to = "type", values_to = "auc") %>% 
    filter(!is.na(auc)) %>% 
    mutate(type = factor(type, levels = c("In-sample", "Cross-validated"))) %>% 
    ggplot(aes(x = auc, y = sp)) +
    geom_vline(xintercept = 0.5, colour = "grey60", linewidth = 0.3) +
    geom_line(aes(group = sp), colour = "#d4d3cd", linewidth = 0.6) +
    geom_point(aes(colour = type, shape = type), size = 2) +
    scale_colour_manual(values = c("In-sample" = "#a3a29c", "Cross-validated" = "#2a78d6"), name = NULL) +
    scale_shape_manual(values = c("In-sample" = 16, "Cross-validated" = 16), name = NULL) +
    scale_y_discrete(labels = function(x) { lab <- taxon_labels(x); lab[x == meta_label] <- cv_pooled_label; lab }) +
    scale_x_continuous(limits = c(0.5, 1)) +
    labs(title = "Discrimination: in-sample vs spatially cross-validated AUC",
         subtitle = "Held-out trees are in 2 km blocks not used for fitting; the spatial field cannot help there",
         x = "AUC", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          legend.position = "top", legend.justification = "left",
          plot.background = element_rect(fill = "white", colour = NA))
  p_cv_auc
  ggsave(p_cv_auc, filename = "output/SI_cv_auc.png", width = 6.5, height = 6, dpi = 300)
  
### Fig 8: losses attributable to modifiable drivers (Shapley-averaged attributable fractions) ##########
  # Associational, not causal: the share of observed 2017-21 canopy losses associated with exposure beyond a reference level,
  # under each taxon's fitted model. Exposure is only reduced, never increased: trees already at or better than the reference
  # keep their observed values. Land use, canopy within 100 m, heat, SVI, population and park proximity are held at observed
  # values (context; neighborhood canopy was left out of the attribution because its interpretation is complicated).
  # Groups interact multiplicatively (cloglog), so each group's share of the joint attributable losses is its Shapley value:
  # its average marginal contribution over all orders of removing groups (average attributable fraction; Eide & Gefeller 1995).
  # Uncertainty: posterior draws of the fixed effects; each tree's spatial field and other terms are held at the posterior mean.
  
  attr_ref <- c(imp = quantile(tree_vars_full_nona_format$imp_vis_10_2017_z, 0.10, names = FALSE)) #10th percentile, pooled
  attr_groups <- list(
    "Construction & demolition" = c("is_B_cons_bool", "is_DM_bool", "is_B_cons_1317_bool", "is_DM_1317_bool"),
    "Lack of stewardship"       = "stewardship",
    "Impervious surface"        = "imp_vis_10_2017_z",
    "Sandy inundation zone"     = "in_sandy_zone_bool")
  #counterfactual value of each term (only moves exposure toward the reference)
  attr_counterfactual <- function(trees) {
    tibble(is_B_cons_bool = 0, is_DM_bool = 0, is_B_cons_1317_bool = 0, is_DM_1317_bool = 0,
           stewardship = 1,
           imp_vis_10_2017_z = pmin(trees$imp_vis_10_2017_z, attr_ref["imp"]),
           in_sandy_zone_bool = 0, .rows = nrow(trees))
  }
  attr_n_draws <- 500
  k <- length(attr_groups)
  subsets <- as.matrix(expand.grid(rep(list(0:1), k))); colnames(subsets) <- names(attr_groups) #all 2^k combinations
  shapley_w <- function(s) factorial(s) * factorial(k - s - 1) / factorial(k)
  #for each group: the subsets without it, the matching subsets with it added, and the Shapley weights (computed once)
  subset_key <- apply(subsets, 1, paste, collapse = "")
  shapley_pairs <- lapply(seq_len(k), function(g) {
    without <- which(subsets[, g] == 0)
    with <- match(apply(subsets[without, , drop = FALSE], 1, function(x) paste(replace(x, g, 1), collapse = "")), subset_key)
    list(without = without, with = with, w = shapley_w(rowSums(subsets[without, , drop = FALSE])))
  })
  
  attr_one_taxon <- function(m, trees) {
    cf <- attr_counterfactual(trees)
    terms <- unlist(attr_groups)
    dx <- as.matrix(cf[, terms]) - as.matrix(trees[, terms]) #change in each term under the counterfactual
    draws <- m$fixed_draws[terms, seq_len(attr_n_draws), drop = FALSE]
    res <- lapply(seq_len(attr_n_draws), function(s) {
      #change in each tree's linear predictor from switching off each group, for this draw
      d_group <- sapply(attr_groups, function(g) dx[, g, drop = FALSE] %*% draws[g, s])
      losses_S <- colSums(1 - exp(-exp(trees$eta_mean + d_group %*% t(subsets)))) #expected losses with each subset of groups removed
      v <- losses_S[1] - losses_S #reduction in expected losses for each subset (row 1 = nothing removed)
      phi <- sapply(shapley_pairs, function(sp_g) sum(sp_g$w * (v[sp_g$with] - v[sp_g$without]))) #Shapley value of each group
      c(setNames(phi / losses_S[1], names(attr_groups)), joint = v[nrow(subsets)] / losses_S[1]) #attributable fractions
    })
    af <- do.call(rbind, res)
    as_tibble(af) %>% mutate(draw = row_number(), sp = m$sp, observed_losses = sum(trees$mort))
  }
  
  attr_draws <- parallel::mclapply(seq_along(mort_model_list), function(i) {
    attr_one_taxon(mort_model_list[[i]], trees_model_list[[i]])
  }, mc.cores = 7) %>% bind_rows()
  
  #citywide: attributable losses summed over taxa (per draw), as a share of all observed losses
  attr_long <- attr_draws %>% 
    tidyr::pivot_longer(c(all_of(names(attr_groups)), joint), names_to = "group", values_to = "af") %>% 
    mutate(losses = af * observed_losses)
  attr_city <- attr_long %>% group_by(group, draw) %>% 
    summarise(losses = sum(losses), observed_losses = sum(observed_losses), .groups = "drop") %>% 
    mutate(af = losses / observed_losses, sp = meta_label)
  attr_summary <- bind_rows(attr_long, attr_city) %>% 
    group_by(sp, group) %>% 
    summarise(observed_losses = first(observed_losses),
              losses = median(losses), losses_lower = quantile(af * observed_losses, 0.025), losses_upper = quantile(af * observed_losses, 0.975),
              pct = 100 * median(af), pct_lower = 100 * quantile(af, 0.025), pct_upper = 100 * quantile(af, 0.975), .groups = "drop")
  attr_summary %>% filter(sp == meta_label)
  write_csv(attr_summary %>% mutate(sp = ifelse(sp == meta_label, "All taxa (summed)", sp)), "output/attribution_summary.csv")
  
  attr_plot_df <- attr_summary %>% 
    mutate(group = factor(group, levels = c(names(attr_groups), "joint"), labels = c(names(attr_groups), "All four groups (joint)")),
           sp = factor(sp, levels = sp_levels), 
           type = factor(ifelse(sp == meta_label, "All taxa", "Taxon"), levels = c("Taxon", "All taxa")))
  p_attr <- ggplot(attr_plot_df, aes(x = pct, y = sp, colour = type, shape = type)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_linerange(aes(xmin = pct_lower, xmax = pct_upper), linewidth = 0.5) +
    geom_point(size = 1.6) +
    scale_colour_manual(values = c("Taxon" = "#2a78d6", "All taxa" = "#eb6834"), name = NULL) +
    scale_shape_manual(values = c("Taxon" = 16, "All taxa" = 18), name = NULL) +
    scale_y_discrete(labels = function(x) { lab <- taxon_labels(x); lab[x == meta_label] <- "All taxa (summed)"; lab }) +
    facet_wrap(~group, ncol = 3, scales = "free_x") +
    labs(title = "Figure 8. Share of 2017-21 canopy losses associated with modifiable drivers",
         subtitle = "Shapley-averaged attributable fraction (%) with 95% interval; exposure reduced to: no construction/demolition,\nstewardship present, 10th percentile impervious share, outside the Sandy zone",
         caption = "Associational, not causal. Land use, canopy within 100 m, heat, SVI, population and park proximity held at observed values.\nIntervals reflect uncertainty in the driver effects.",
         x = "% of observed canopy losses", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0), panel.spacing.x = unit(1.5, "lines"),
          axis.text.y = element_text(size = 7), legend.position = "top", legend.justification = "left",
          plot.caption = element_text(hjust = 0, colour = "#52514e"),
          plot.background = element_rect(fill = "white", colour = NA))
  p_attr
  ggsave(p_attr, filename = "output/Fig8_attribution.png", width = 10, height = 9, dpi = 300)
  
### Fig 7: map of predicted mortality by individual tree #######################################################
 tree_preds <- bind_rows(trees_model_list)
  

   #project the predictions to the crs of the basemap tile (3857)
  fitted_preds_sf <- st_as_sf(tree_preds, crs = crs_m, coords = c("x", "y")) %>% 
    st_transform(., crs = 3857) %>% 
    bind_cols(st_coordinates(.) %>% as.data.frame())


  
  #zoom in on individual trees
  ggplot() + ggthemes::theme_few() +   
    geom_spatraster_rgb(data = nyc_topo_spatrast) +
    geom_point(data = fitted_preds_sf, aes(x = X, y = Y, color = fitted_p)) +
    coord_sf(xlim = c(-8233000, (-8233000 + 3000)), ylim = c(4955000, 4955000 + 3000)) +
    #facet_wrap(~sp_a, ncol = 6) +
    scale_color_viridis_c(option = "turbo", name = "predicted \n4-year \nmortality", limits = c(0, 0.3), oob = squish) + xlab("") + ylab("") + 
    theme(strip.text = element_text(face = "italic"),
          legend.position = c(0.14, 0.7), #legend.position = c(0.92, 0.14),
          legend.background = element_rect(fill = "white", color = "grey80"),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) 
  
  
  # #create figure for predicted survival across the city
  # ggplot() + ggthemes::theme_few() +   
  #   geom_spatraster_rgb(data = nyc_topo_spatrast) +
  #   stat_summary_hex(data = fitted_preds_sf, aes(x = X, y = Y, z = fitted_p), fun = median, bins = 40) +
  #   facet_wrap(~sp_a, ncol = 6) +
  #   scale_fill_viridis_c(option = "turbo", direction = -1, name = "predicted \nmedian \nsurvival (%)") + xlab("") + ylab("") + 
  #   theme(strip.text = element_text(face = "italic"),
  #         legend.position = c(0.14, 0.7), #legend.position = c(0.92, 0.14),
  #         legend.background = element_rect(fill = "white", color = "grey80"),
  #         panel.grid.major = element_blank(),  
  #         panel.grid.minor = element_blank(),
  #         axis.text = element_blank(),
  #         axis.ticks = element_blank())
