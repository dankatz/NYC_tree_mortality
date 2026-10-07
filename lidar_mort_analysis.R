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
         # TEMPORARY: is_S_cons_10m_17_21 appears inverted in v4 (TRUE for 88% of trees vs 6.8% for 2013-17); flipped pending confirmation from Ging-Yan
         is_S_cons_bool = !as.logical(is_S_cons_10m_17_21),
         is_DM_bool = as.logical(is_DM_10m_17_21),
         is_B_cons_1317_bool = as.logical(is_B_cons_10m_13_17), #construction 2013-2017, before the lidar interval (lagged damage)
         is_S_cons_1317_bool = as.logical(is_S_cons_10m_13_17),
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
         is_B_cons_bool, is_S_cons_bool, is_DM_bool, # building construction, street construction, building demolition (2017-2021)
         is_B_cons_1317_bool, is_S_cons_1317_bool, is_DM_1317_bool, # same, 2013-2017
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
         RPL_THEME1) %>% # RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME3, RPL_THEME4, RPL_THEMES) %>%
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
  scale_vars <- c("summer_max", "summer_min", "imp_vis_10_2017", "TC_100m_2017", "pop_100m_log", "RPL_THEME1")
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

## Mesh and SPDE (Matern with PC priors), shared by all taxa; all distances in meters ------------------------------
  mesh <- fm_mesh_2d_inla(
    boundary = nyc.bdry, #loaded at start of script
    max.edge = c(200, 2000), #200 m inside the city (~range/2.5 vs the 500 m range prior); 100 m (217k vertices) failed in INLA; coarse in the outer extension
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

#full fits are large (latent field + config), so each is saved to disk and only a small summary is kept in memory;
#if a taxon's summary file already exists the fit is skipped, so an interrupted run can be resumed
  fit_dir <- "output/fits"
  dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
  fit_file <- function(sp, type) file.path(fit_dir, paste0(gsub("[^A-Za-z]+", "_", sp), "_", type, ".rds"))

#create empty lists to save output
  mort_model_list <- vector("list", length(sp_list))
  trees_model_list <- vector("list", length(sp_list))
  roc_list <- vector("list", length(sp_list))

#run model for focal species; set the TAXA_TO_FIT environment variable (semicolon-separated) to fit only some taxa
taxa_to_fit <- if (nzchar(Sys.getenv("TAXA_TO_FIT"))) strsplit(Sys.getenv("TAXA_TO_FIT"), ";")[[1]] else sp_list
for (i in which(sp_list %in% taxa_to_fit)){
#for (i in 1:2){
 
  sp_focal <- sp_list[i] #sp_focal <- sp_list[6]
  print(paste(i, sp_focal))
  
  if (file.exists(fit_file(sp_focal, "summary"))) { #already fitted: reload the summary and move on
    fit_summary <- readRDS(fit_file(sp_focal, "summary"))
    mort_model_list[[i]] <- fit_summary$model
    trees_model_list[[i]] <- fit_summary$trees
    roc_list[[i]] <- fit_summary$roc
    next
  }
  
  sp_sub_format <- tree_vars_full_nona_format %>% filter(sp_a == sp_focal)
    
  trees <- sp_sub_format %>% 
    mutate(log_exposure = log(lidar_interval_yr)) %>% 
    select(tree_id, mort, log_exposure, stewardship, is_B_cons_bool, is_S_cons_bool, is_DM_bool,
           is_B_cons_1317_bool, is_S_cons_1317_bool, is_DM_1317_bool,
           is_LU_hidensres, is_LU_comm, is_LU_indust, is_LU_publicinst, is_LU_outdoorrec, is_LU_vacant, is_LU_noparcel, #LowDensityResidential is the reference
           dbh_cm, summer_max_z, summer_min_z, imp_vis_10_2017_z, TC_100m_2017_z, pop_100m_log_z, RPL_THEME1_z, in_sandy_zone_bool,
           near_park_10m, #NYCHA proximity removed: too few flagged trees per taxon
           x, y) %>% 
    mutate(across(where(is.logical), as.numeric))
  trees_sf <- st_as_sf( trees, coords = c("x", "y"), crs = crs_m, remove = FALSE) 
  
  #1D mesh for the DBH smooth (cm), spanning this taxon's range
  dbh_mesh <- fm_mesh_1d(seq(min(trees$dbh_cm), max(trees$dbh_cm), length.out = 20))
  
  ## ---- Model components and formula ---------------------------------------
  # linear terms get an N(0, 1) prior (mean.linear = 0, prec.linear = 1)
  cmp <- mort ~ Intercept(1) +
    exposure(log_exposure, model = "offset") + #intercept is then a log annual hazard
    stewardship(stewardship, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_B_cons_bool(is_B_cons_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_S_cons_bool(is_S_cons_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_DM_bool(is_DM_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_B_cons_1317_bool(is_B_cons_1317_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_S_cons_1317_bool(is_S_cons_1317_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    is_DM_1317_bool(is_DM_1317_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    in_sandy_zone_bool(in_sandy_zone_bool, model = "linear", mean.linear = 0, prec.linear = 1) +
    RPL_THEME1_z(RPL_THEME1_z, model = "linear", mean.linear = 0, prec.linear = 1) +
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
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
      verbose = FALSE, num.threads = "8:1" #8 threads: same estimates as 16, ~30% less memory, slightly faster
    )
  )
  
  #summary(fit)
  
  ### saving data from focal species
  fit$sp <- sp_focal
  saveRDS(fit, fit_file(sp_focal, "fit")) #full fit, needed for predict() (e.g., spatial field maps)
  mort_model_list[[i]] <- list(sp = sp_focal, #small version of the fit used by the figures below
                               summary.fixed = fit$summary.fixed,
                               summary.hyperpar = fit$summary.hyperpar,
                               summary.random = list(dbh_smooth = fit$summary.random$dbh_smooth),
                               waic = fit$waic[c("waic", "p.eff")],
                               dic = fit$dic[c("dic", "p.eff")])
  
  trees$fitted_p <- fit$summary.fitted.values$mean[1:nrow(trees)] #predicted probability of canopy loss over the 4-year interval
  trees$resid    <- trees$mort - trees$fitted_p
  trees$sp <- sp_focal
  trees_join <- left_join(trees, sp_sub_format %>% select(tree_id, sp_a, species, genus), by = "tree_id")
  
  trees_model_list[[i]] <- trees_join
  
  roc_curve <- roc(trees$mort, trees$fitted_p)
  roc_curve$sp <- sp_focal
  roc_list[[i]] <- roc_curve
  
  saveRDS(list(model = mort_model_list[[i]], trees = trees_join, roc = roc_curve), fit_file(sp_focal, "summary"))
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
   
### Fig 4: Forest plot of fixed effects ########################################
  # term labels and display order, grouped by construct; continuous terms are per 1 SD (pooled across taxa)
  term_labels <- c(
    stewardship         = "Stewardship (any signs)",
    is_B_cons_1317_bool = "Building construction 2013-17",
    is_B_cons_bool      = "Building construction 2017-21",
    is_DM_1317_bool     = "Demolition 2013-17",
    is_DM_bool          = "Demolition 2017-21",
    is_S_cons_1317_bool = "Street construction 2013-17",
    is_S_cons_bool      = "Street construction 2017-21 *",
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
    RPL_THEME1_z        = "SVI socioeconomic theme (+1 SD)")
  
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
  meta_df <- coef_df %>% group_by(term) %>% group_modify(~ meta_fit(.x)) %>% ungroup() %>% mutate(sp = "All taxa (RE mean)")
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
  
  # taxa ordered by sample size (largest at top), with the meta-analytic mean above them
  sp_n <- tree_vars_full_nona_format %>% count(sp_a) %>% arrange(n)
  sp_levels <- c(sp_n$sp_a, "All taxa (RE mean)")
  forest_df <- bind_rows(coef_df %>% mutate(type = "Taxon"), meta_df %>% mutate(type = "Mean across taxa")) %>% 
    mutate(sp = factor(sp, levels = sp_levels),
           term = factor(term_labels[term], levels = term_labels),
           type = factor(type, levels = c("Taxon", "Mean across taxa")))
  
  hr_breaks <- function(lims) { b <- pretty(lims, n = 4); b[b > 0] } #readable hazard-ratio ticks within each panel's range
  
  p_forest <- ggplot(forest_df, aes(x = exp(estimate), y = sp, colour = type, shape = type)) +
    geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
    geom_linerange(aes(xmin = exp(lower), xmax = exp(upper)), linewidth = 0.5) +
    geom_point(size = 1.6) +
    scale_x_log10(breaks = hr_breaks) +
    scale_y_discrete(labels = function(x) parse(text = ifelse(x %in% c("other", "All taxa (RE mean)"), #italicize species names only
                                                              paste0("'", x, "'"), paste0("italic('", x, "')")))) +
    scale_colour_manual(values = c("Taxon" = "#2a78d6", "Mean across taxa" = "#eb6834"), name = NULL) +
    scale_shape_manual(values = c("Taxon" = 16, "Mean across taxa" = 18), name = NULL) +
    facet_wrap(~term, ncol = 4, scales = "free_x", labeller = label_wrap_gen(width = 32)) +
    labs(
      title = "Figure 4. Effects of drivers on street tree mortality",
      subtitle = "Hazard ratio (log scale) with 95% credible interval; one spatial Bernoulli GLMM (cloglog link) per taxon",
      caption = "Mean across taxa: random-effects meta-analysis (REML, Knapp-Hartung 95% CI). Continuous drivers: per 1 SD, pooled across taxa.\nLand use is relative to low-density residential. * 2017-21 street construction indicator flipped pending confirmation (see ANALYSIS_PLAN.md).",
      x = "Hazard ratio", y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#e6e5e0", linewidth = 0.3),
          strip.text = element_text(face = "bold", hjust = 0),
          axis.text.y = element_text(size = 6.5),
          legend.position = "top", legend.justification = "left",
          plot.caption = element_text(hjust = 0, colour = "#52514e"),
          plot.background = element_rect(fill = "white", colour = NA))
  
  p_forest
  ggsave(p_forest, filename = "output/Fig4_forest_by_taxon.png", width = 11, height = 15, dpi = 300)
  
  # compact version: mean across taxa only, one row per driver
  p_forest_mean <- meta_df %>% 
    mutate(term = factor(term_labels[term], levels = rev(term_labels))) %>% 
    ggplot(aes(x = exp(estimate), y = term)) +
    geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
    geom_linerange(aes(xmin = exp(lower), xmax = exp(upper)), colour = "#eb6834", linewidth = 0.6) +
    geom_point(colour = "#eb6834", shape = 18, size = 2.5) +
    scale_x_log10(breaks = c(0.8, 0.9, 1, 1.25, 1.5, 2, 2.5)) +
    labs(title = "Mean effect across taxa (random-effects meta-analysis)",
         subtitle = sprintf("Hazard ratio with 95%% CI (REML, Knapp-Hartung); %d taxa", length(unique(coef_df$sp))),
         x = "Hazard ratio (log scale)", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
  
  p_forest_mean
  ggsave(p_forest_mean, filename = "output/Fig4_forest_mean.png", width = 7, height = 5.5, dpi = 300)
  
  write_csv(bind_rows(coef_df, meta_df), "output/Fig4_coefficients.csv")
  
### SI X: Spatial random field map ---------------------------------------------
  
  #all taxa share one mesh; i_field picks the taxon to map
  i_field <- 1
  field_fit <- readRDS(fit_file(sp_list[i_field], "fit"))
  field_pred <- predict(
    field_fit,
    fm_pixels(mesh, dims = c(300, 300), mask = nyc_boundary),
    ~ spatial_field
  ) %>% 
    cbind(., st_coordinates(.))
  
  
  p_field <- ggplot() +
    gg(field_pred, aes(color = mean)) +
    
    scale_color_viridis_c() +
    #coord_equal() +
    labs(title = paste("Estimated spatial random effect (SPDE):", sp_list[i_field]),
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal(base_size = 11)
  
  p_field + geom_sf(data = nyc_boundary, fill = NA) # + geom_point(data = trees, aes(x = x, y = y), size = 0.4, colour = "white", alpha = 0.5) +
  
### Fig 5: map of predicted mortality by individual tree #######################################################
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
