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

library(patchwork)
library(knitr)
library(kableExtra)
library(gt)


your_path_for_box <- "C:/Users/dsk273/Box/Katz lab/NYC/"
tree_vars <- fread(paste0(your_path_for_box, '/tree_mortality_variables/all_variables_v2.csv'))
socioeco_vars <- fread(paste0(your_path_for_box,'/tree_mortality_variables/model_df_socioeco.csv'))


# merge tree and socioeconomic variables
socioeco_vars_select <- socioeco_vars %>% select(tree_id, #medincomeE has too many NAs to reasonably include
                                                 estimate_c_perc_unoccupied, estimate_c_perc_non_white, estimate_c_perc_unemployed, estimate_c_perc_poverty,
                                                 RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME4, RPL_THEMES)
tree_vars_full <- left_join(tree_vars, socioeco_vars_select) 

#load in nyc boundary polygon and save in the format needed by fmesher
nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>%
  st_union() %>% #combine the different boroughs
  st_transform(., crs = 2263)
nyc.bdry <-  as(nyc_boundary, "Spatial") %>% fm_as_segm()

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
  tree_vars_full$LandUse_char_collapsed <- as.factor(x)
  tree_vars_full$LandUse_char_collapsed <- relevel(as.factor(tree_vars_full$LandUse_char_collapsed), "LowDensityResidential") # this is Mike's recommendation for base case, may need to collapse land use classes


  summary(tree_vars_full$LandUse_char_collapsed)
### create dataframe for analysis ##############################################

tree_vars_full_nona_format <- tree_vars_full %>%
  #create derived variables and filter out rows
  mutate(canopy_endstate = case_when( canopy_change == 1 ~ 1,      #1 = no change; assumed survival
                                      canopy_change == 2 ~ NA,     #2 = canopy gain, not to include in this analysis
                                      canopy_change == 3 ~ 0)) %>% #3 = canopy loss; assumed mortality
  filter(!is.na(canopy_endstate)) %>%
  mutate(dbh_cm = tree_dbh * 2.54) %>%  #convert dbh from inches to cm
  filter(dbh_cm > 7.62) %>%  # remove small trees; 7.6 cm is 3 inches listed in Bigelow
  filter(dbh_cm < 300) %>%  #removing trees with a DBH greater than the maximum ever recorded in NYC 
  filter(pluto_dist < 20) %>% # tree distance 
        #   # reasonable distance threshold for land use proximity? This is in survey foot (effectively feet)
        #   quantile(tree_vars_full_nona$pluto_dist, seq(0, 1, 0.01))
        # # picking 20, this is approx 95% of the remaining trees. Will apply this in the mutate part
  mutate(
         #BldgClass_fac = factor(BldgClass), 
         #BldgClass_Group_fac = factor(BldgClass_Group),
         in_sandy_zone_bool = as.logical(in_sandy_zone),
         is_B_cons_bool = as.logical(is_B_cons),
         is_S_cons_bool = as.logical(is_S_cons),
         is_DM_bool = as.logical(is_DM),
         is_LU_lowdensres = case_when(LandUse_char_collapsed == "LowDensityResidential" ~ 1, .default = 0),
         is_LU_hidensres = case_when(LandUse_char_collapsed == "HigherDensityResidential" ~ 1, .default = 0),
         is_LU_comm = case_when(LandUse_char_collapsed == "CommercialMix" ~ 1, .default = 0),
         is_LU_indust = case_when(LandUse_char_collapsed == "IndustryTransport" ~ 1, .default = 0),
         is_LU_publicinst = case_when(LandUse_char_collapsed == "PublicInst" ~ 1, .default = 0),
         is_LU_outdoorrec = case_when(LandUse_char_collapsed == "OpenOutdoorRec" ~ 1, .default = 0),
         is_LU_vacant = case_when(LandUse_char_collapsed == "Vacant" ~ 1, .default = 0),
         imp_5_2017 = BD_5_2017 + RD_5_2017 + OI_5_2017 + RR_5_2017,
         imp_10_2017 = BD_10_2017 + RD_10_2017 + OI_10_2017 + RR_10_2017,
         imp_15_2017 = BD_15_2017 + RD_15_2017 + OI_15_2017 + RR_15_2017,
         stewardship = case_when(steward == "None" ~ 0,
                                is.na(steward) ~ 0, #assuming that NA means no signs of stewardship
                                 steward == "1or2" ~ 1,
                                 steward == "3or4" ~ 1,
                                 steward == "4orMore" ~ 1)) %>%
  select(tree_id, geom,  #only include the variables that will be used in the analysis (to avoid any extra NA values)
         genus, species, 
         canopy_endstate, # dependent variable 
         dbh_cm, 
         stewardship, # whether there were any signs of stewardship observed
         #LandUse_char_collapsed, #BldgClass_fac, BldgClass_Group_fac, LandUse_char,  # building type (detailed), building type (high level), land use
         is_LU_lowdensres, is_LU_hidensres, is_LU_comm, is_LU_indust, is_LU_publicinst, is_LU_outdoorrec, is_LU_vacant,
         in_sandy_zone_bool, # sandy inundation zone
         is_B_cons_bool, is_S_cons_bool, is_DM_bool, # building construction, street construction, building demolition
         summer_mean, #summer_max, summer_min, days_max_27, days_max_32, days_max_35, # temperature
        # TC_10_diff, GS_10_diff, SO_10_diff, WA_10_diff, BD_10_diff, RD_10_diff, OI_10_diff, RR_10_diff, # land cover diff, doing 10 m
        imp_5_2017, GS_5_2017, BD_5_2017, TC_5_2017, 
        imp_10_2017, GS_10_2017, BD_10_2017,  TC_10_2017,#ORD_10_2017, I_10_2017, RR_10_2017, # land cover 2017, doing 10 m #WA_10_2017, 
        imp_15_2017, GS_15_2017, BD_15_2017, TC_15_2017,
        near_NYCHA_15m,
        near_park_15m,
        # TC_10_2021, GS_10_2021, SO_10_2021, WA_10_2021, BD_10_2021, RD_10_2021, OI_10_2021, RR_10_2021, # land cover 2021, doing 10 m
         RPL_THEME1, RPL_THEME3) %>% # RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME3, RPL_THEME4, RPL_THEMES) %>%
            # One big thing to note is that SVI is not available for parks and certain public land uses, so will limit applicability in certain parts of the city
    separate_wider_delim(., cols = geom, delim = ",", names = c( "x", "y")) %>%  #extract coordinates from geom text string 
    mutate(y = readr::parse_number(y), #crs = 2263 
           x = readr::parse_number(x)) %>% 
    drop_na()  #remove any rows that have an NA value
    
### including the list of species with enough individuals to analyze ##################
  cut_off_n <- 5000 #cut off number of individuals
  
  #species with a n above that number
  species_n <- 
    tree_vars_full_nona_format %>% 
    group_by(species) %>% 
    summarize(n = n()) %>% 
    filter(n > cut_off_n) %>% 
    filter(species != "Acer") #remove the unidentified Acer individuals
  
  #label non-focal species as "other"
  tree_vars_full_nona_format <- tree_vars_full_nona_format %>% 
    mutate(sp_a = case_when(species %in% species_n$species ~ species,
                            TRUE ~ "other")) %>% 
    left_join(., species_n) #add the n for the sp_a back to the main dataframe
  
### create species level variables ###############################################
  
  # #create bins for tree DBH
  # tree_vars_full_nona_format <-
  #   tree_vars_full_nona_format %>% 
  #   group_by(sp_a) %>% 
  #   mutate(dbh_quintile = as.factor(ntile(dbh_cm, 5)), #calculate survival per quintile or quartile per species
  #          dbh_decile = as.factor(ntile(dbh_cm, 10))) %>%  
  #   ungroup() 
  #str(tree_vars_full_nona_format$dbh_quintile)
  
  #save scale numeric variables within species
  scale_params <-
    tree_vars_full_nona_format %>% 
    group_by(sp_a) %>% 
    select(x, y, dbh_cm, summer_mean, imp_10_2017, GS_10_2017, BD_10_2017, RPL_THEME1, RPL_THEME3) %>% 
    summarise(across(everything(),
                     list(mean = ~mean(.x, na.rm=TRUE),
                          sd   = ~sd(.x,   na.rm=TRUE))))
  
  # Apply Z-scaling
  tree_vars_full_nona_format <- tree_vars_full_nona_format %>% 
    group_by(sp_a) %>% 
    mutate(across(c(x, y, dbh_cm, summer_mean, imp_10_2017, GS_10_2017, BD_10_2017, RPL_THEME1, RPL_THEME3),
                  ~ (.x - mean(.x, na.rm=TRUE)) / sd(.x, na.rm=TRUE),
                  .names = "{.col}_z")) %>% 
    ungroup()
  
### some last pre-analysis data exploration ########################################
  tree_vars_full_nona_format %>% sample_n(10000) %>%
    ggplot(aes(x = x, y = y, z = is_LU_vacant)) + stat_summary_hex(fun = "mean", binwidth = 1000) + theme_bw() + scale_fill_viridis_c()

  # #map of tree samples across NYC
  # tree_vars_full_nona_format %>% 
  #   filter(sp_a == "Quercus palustris") %>% #sample_n(100000) %>%
  #   ggplot(aes(x = x, y = y)) + geom_hex(binwidth = 5000) + theme_bw() + scale_fill_viridis_c()


  tree_vars_full_nona_format %>% sample_n(10000) %>%
    ggplot(aes(x = TC_5_2017, y = TC_15_2017)) + geom_point()

  cor(tree_vars_full_nona_format$TC_5_2017, tree_vars_full_nona_format$TC_15_2017)
  cor(tree_vars_full_nona_format)
  # 
  # numeric_df <- tree_vars_full_nona_format[, sapply(tree_vars_full_nona_format, is.numeric)]
  # cor_matrix <- cor(numeric_df, use = "complete.obs")
  
  
### table 1: summary of survival by focal species  ####################################
  
  table_1 <- tree_vars_full_nona_format %>% 
    group_by(sp_a) %>% 
    summarize(n = n(),
              median_dbh = median(dbh_cm),
              sum_alive = sum(canopy_endstate),
              mean_alive = sum_alive/n,
              annual_surv = mean_alive ^ (1/4))  #annual survival
  
  table_1_all_trees <- tree_vars_full_nona_format %>% 
    summarize(sp_a = "all trees", 
              n = n(),
              median_dbh = median(dbh_cm),
              sum_alive = sum(canopy_endstate),
              mean_alive = sum_alive/n,
              annual_surv = mean_alive ^ (1/4))
  
  table_1 <- rbind(table_1, table_1_all_trees) %>% 
    select(-sum_alive, -mean_alive) %>% 
    mutate(annual_surv = annual_surv * 100)
  
  #add common names
  common_name_lookup <- read_csv(paste0(your_path_for_box,'/tree_mortality_variables/common_name_lookup.csv')) %>% 
    rename(sp_a = species)
  table_1 <- left_join(table_1, common_name_lookup) %>% 
    select(sp_a, common_name, n, median_dbh, annual_surv)
  
  table_1  %>% ungroup() %>% 
    gt() %>% 
    fmt_number(columns  = c(median_dbh, annual_surv), decimals = 1) |>
    fmt_integer(columns = n, sep_mark = ",") %>% 
    cols_label(
      sp_a = "species",
      common_name = "common name",
      n = "n",
      median_dbh = "median DBH (cm)",
      annual_surv = "annual survival (%)" ) |>
    tab_style(style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2),
      style = "solid"),
      locations = cells_column_labels()) 
  #%>%  gtsave( paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/table1.docx"))  
  
  
  
### Fig 2: survival as a function of DBH #########################################
  fig2 <- tree_vars_full_nona_format %>% 
  select(sp_a, canopy_endstate, dbh_cm) %>% 
  group_by(sp_a) %>% 
  mutate(quintile = ntile(dbh_cm, 4)) %>%  #calculate survival per quintile or quartile per species
  ungroup() %>% 
  group_by(sp_a, quintile) %>% 
  summarize(median_dbh = median(dbh_cm),
            mean_mort = mean(canopy_endstate),
            mean_annual_mort = 100 * (mean_mort^(1/4)),
            n = n()) %>% 
  ggplot(aes(x = median_dbh, y = mean_annual_mort, color = sp_a)) + geom_point() + geom_line()+ theme_bw() + xlab("DBH (cm)") + ylab("annual survival (%)") +
  scale_color_viridis_d(option = "turbo", name = "species") +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.text = element_text(face = "italic"))      

  # ggsave(paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/fig_2.png"),
  #        width = 7, height = 5, units = "in", dpi = 400)


  
### Fig 3: maps of observed survival for each taxon #########################################
  #load in basemap
  nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
  nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
  
  #function to remove bins with lower than a certain n
  bin_removal <- function(x) {
    if(length(x) < 10) return(NA) 
    return((mean(x) ^ (1/4)) * 100)
  }
  
  #project the data to the crs of the basemap tile (3857)
  tree_vars_full_nona_format_sf_3857 <- st_as_sf(tree_vars_full_nona_format, crs = 2263, coords = c("x", "y")) %>% 
    st_transform(., crs = 3857) %>% 
    bind_cols(st_coordinates(.) %>% as.data.frame())
  
  #create figure
  fig3_mort_map <- 
    ggplot() + ggthemes::theme_few() +   
    geom_spatraster_rgb(data = nyc_topo_spatrast) +
    stat_summary_hex(data = tree_vars_full_nona_format_sf_3857, aes(x = X, y = Y, z =  canopy_endstate), fun = bin_removal, bins = 20) +
    facet_wrap(~sp_a, ncol = 6) +
    scale_fill_viridis_c(option = "turbo", direction = -1, name = "mean \nsurvival (%)  ", 
                         limits = c(90, 100),
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

#create empty lists to save output
  mort_model_list <- vector("list", length(sp_list))
  trees_model_list <- vector("list", length(sp_list))
  roc_list <- vector("list", length(sp_list))

#run model for focal species
for (i in 1:length(sp_list)){      
#for (i in 1:2){
 
  sp_focal <- sp_list[i] #sp_focal <- sp_list[6]
  print(paste(i, sp_focal))
  
  sp_sub_format <- tree_vars_full_nona_format %>% filter(sp_a == sp_focal) %>%  #
     group_by(x, y) %>% filter(n() == 1) %>%  ungroup() %>%  #remove some rows that have duplicated location information
    filter(!is.na(y))
    
  trees <- sp_sub_format %>% 
    select(canopy_endstate, stewardship, is_B_cons_bool, is_S_cons_bool, is_DM_bool,
           is_LU_lowdensres, is_LU_hidensres, is_LU_comm, is_LU_indust, is_LU_publicinst, is_LU_outdoorrec, is_LU_vacant,
           dbh_cm_z, summer_mean_z, imp_10_2017_z, GS_10_2017_z, BD_10_2017_z, RPL_THEME1_z, in_sandy_zone_bool,
           near_NYCHA_15m, near_park_15m,
           x, y)
  trees_sf <- st_as_sf( trees, coords = c("x", "y"), remove = FALSE) 
  
  #create covariates
  dbh_grid <-  seq(min(sp_sub_format$dbh_cm_z), max(sp_sub_format$dbh_cm_z), length.out = 10)
  
  ## Mesh and SPDE (Matern with PC priors) ------------------------------
  max.edge <- diff(range(st_coordinates(trees_sf)[,1]))/(3*5) #https://rpubs.com/jafet089/886687
  bound.outer <- diff(range(st_coordinates(trees_sf)[,1]))/10

  mesh <- fm_mesh_2d_inla(
    loc      = as.matrix(trees[, c("x", "y")]),
    max.edge = c(1000, 5000), #c(1,2)*max.edge,
    offset   = c(max.edge, bound.outer),
    cutoff = max.edge/10,
    boundary = nyc.bdry #loaded at start of script
  )
  
  #plot mesh to double check it
  ggplot() + gg(mesh) + geom_sf(data= trees_sf,col='purple',size=1.7,alpha=0.5) + theme_minimal()
  cat("Mesh has", mesh$n, "vertices\n")
  
  # PC priors: P(range < 15) = 0.05 ; P(sigma > 1.5) = 0.05
  spde <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(15, 0.05),
    prior.sigma = c(1.5, 0.05)
  )
  
  ## ---- Model components and formula ---------------------------------------
  cmp <- canopy_endstate ~ Intercept(1) +
    stewardship(stewardship, model = "linear") +
    is_B_cons_bool(is_B_cons_bool, model = "linear") +
    is_S_cons_bool(is_S_cons_bool, model = "linear") +
    is_DM_bool(is_DM_bool, model = "linear") +
    in_sandy_zone_bool(in_sandy_zone_bool, model = "linear") +
    RPL_THEME1_z(RPL_THEME1_z, model = "linear") +
    summer_mean_z(summer_mean_z, model = "linear") +
    near_NYCHA_15m(near_NYCHA_15m, model = "linear") +
    near_park_15m(near_park_15m, model = "linear") +

    is_LU_lowdensres(is_LU_lowdensres, model = "linear") +
    is_LU_hidensres(is_LU_hidensres, model = "linear") +
    is_LU_comm(is_LU_comm, model = "linear") +
    is_LU_indust(is_LU_indust, model = "linear") +
    is_LU_publicinst(is_LU_publicinst, model = "linear") +
    is_LU_outdoorrec(is_LU_outdoorrec, model = "linear") +
    is_LU_vacant(is_LU_vacant, model = "linear") +
  
    # imp_5_2017_z(imp_5_2017_z, model = "linear") +
    # GS_5_2017_z(GS_5_2017_z, model = "linear") +
    # BD_5_2017_z(BD_5_2017_z, model = "linear") +
    
    imp_10_2017_z(imp_10_2017_z, model = "linear") +
    GS_10_2017_z(GS_10_2017_z, model = "linear") +
    BD_10_2017_z(BD_10_2017_z, model = "linear") +
    
    # imp_15_2017_z(imp_15_2017_z, model = "linear") +
    # GS_15_2017_z(GS_15_2017_z, model = "linear") +
    # BD_15_2017_z(BD_15_2017_z, model = "linear") +
    
    
    
    dbh_smooth(
      dbh_cm_z,
      model = "rw2",
      values = dbh_grid,
      scale.model = TRUE,                     # makes precision comparable across resolutions
      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))) ) + # PC prior on smoothness
    spatial_field(geometry, model = spde)
  
  ## ---- Fit with bru() -------------------------------------------------------
  lik <- bru_obs(
    formula = canopy_endstate ~ .,
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
      verbose = FALSE
    )
  )
  
  #summary(fit)
  
  ### saving data from focal species
  fit$sp <- sp_focal
  mort_model_list[[i]] <- fit
  
  trees$fitted_p <- fit$summary.fitted.values$mean[1:nrow(trees)]
  trees$resid    <- trees$canopy_endstate - trees$fitted_p
  trees$sp <- sp_focal
  trees_join <- left_join(trees, sp_sub_format)
  
  trees_model_list[[i]] <- trees_join
  
  roc_curve <- roc(trees$canopy_endstate, trees$fitted_p)
  roc_curve$sp <- sp_focal
  roc_list[[i]] <- roc_curve
  
  
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
  #     mean_obs    = mean(canopy_endstate),
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
  #fixed_summary <- fit$summary.fixed
  #mort_model_list <- mort_model_list[1:5]
  #mort_model_list[[1]]$summary.fixed
  
  coef_df <- map_dfr(
    mort_model_list,
    ~ as_tibble(.x$summary.fixed, rownames = "term"),
    .id = "model"
  )%>%
    rename(
      estimate = mean,
      lower    = `0.025quant`,
      upper    = `0.975quant`
    ) %>%
    select(model, term, estimate, lower, upper, sd) %>% 
    left_join(., sp_df)
  
  coef_df
  
  
  p_forest <- ggplot(coef_df, aes(x = estimate, y = sp)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   height = 0.15, colour = "grey30") +
    geom_point(size = 3) +
    facet_wrap(~term, scales = "free")+
    labs(
      title = "Figure 1. Effects of covariates on tree survival",
      subtitle = "Binomial GLM with cloglog link (posterior mean \u00b1 95% CrI)",
      x = "Coefficient (cloglog scale)", y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(face = "italic"))
  
  p_forest
  
  
  
  
### SI X: Spatial random field map ---------------------------------------------
  
  field_pred <- predict(
    mort_model_list[[1]],
    fm_pixels(mesh, mask = NULL),
    ~ spatial_field
  ) %>% 
    cbind(., st_coordinates(.)) %>% 
    st_set_crs(2263)
  
  
  p_field <- ggplot() +
    gg(field_pred, aes(color = mean)) +
    
    scale_color_viridis_c() +
    #coord_equal() +
    labs(title = "Figure 2. Estimated spatial random effect (SPDE)",
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal(base_size = 11)
  
  p_field + geom_sf(data = nyc_boundary, fill = NA) # + geom_point(data = trees, aes(x = x, y = y), size = 0.4, colour = "white", alpha = 0.5) +
  
### Fig 5: map of predicted survival by individual tree #######################################################
 tree_preds <- bind_rows(trees_model_list)
  

   #project the predictions to the crs of the basemap tile (3857)
  fitted_preds_sf <- st_as_sf(tree_preds, crs = 2263, coords = c("x", "y")) %>% 
    st_transform(., crs = 3857) %>% 
    bind_cols(st_coordinates(.) %>% as.data.frame())


  
  #zoom in on individual trees
  ggplot() + ggthemes::theme_few() +   
    geom_spatraster_rgb(data = nyc_topo_spatrast) +
    geom_point(data = fitted_preds_sf, aes(x = X, y = Y, color = fitted_p)) +
    coord_sf(xlim = c(-8233000, (-8233000 + 3000)), ylim = c(4955000, 4955000 + 3000)) +
    #facet_wrap(~sp_a, ncol = 6) +
    scale_color_viridis_c(option = "turbo", direction = -1, name = "predicted \nsurvival (%)", limits = c(0.7, 1)) + xlab("") + ylab("") + 
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




# ### Fig 4: map of predicted survival by species #######################################################
# 
# #project the predictions to the crs of the basemap tile (3857)
# fitted_preds_sf <- st_as_sf(trees, crs = 2263, coords = c("x", "y")) %>% 
#             st_transform(., crs = 3857) %>% 
#             bind_cols(st_coordinates(.) %>% as.data.frame())
#   
# # #load in nyc boundary polygon
# # nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
# #   st_union() %>% #combine the different boroughs
# #   st_transform(., crs = 3857)
# 
# #load in basemap
# nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
# nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 
# 
# #create figure
# fig_pred_fit <- ggplot() + ggthemes::theme_few() +   
#   geom_spatraster_rgb(data = nyc_topo_spatrast) +
#   stat_summary_hex(data = fitted_preds_sf, aes(x = X, y = Y, z = fitted_p), fun = median, bins = 40) +
#   #facet_wrap(~sp_a, ncol = 6) +
#   scale_fill_viridis_c(option = "turbo", direction = -1, name = "predicted \nmedian \nsurvival (%)") + xlab("") + ylab("") + 
#   theme(strip.text = element_text(face = "italic"),
#         legend.position = c(0.14, 0.7), #legend.position = c(0.92, 0.14),
#         legend.background = element_rect(fill = "white", color = "grey80"),
#         panel.grid.major = element_blank(),  
#         panel.grid.minor = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#   # annotation_scale(location = "br",  # "bl" for bottom-left, other options exist
#   #                  bar_cols = c("black", "white"), # Colors of the scale bar segments
#   #                  style = "ticks",
#   #                  text_cex = 0.8) +  # Text size for the scale bar label
#   # annotation_north_arrow(location = "br", height = unit( 0.8, "cm"), style = north_arrow_minimal,
#   #                        pad_x = unit(1, "cm"), pad_y = unit(1, "cm")) +
#   # theme(  legend.position = c(0.1, 0.9),  # Places the legend at the top-left corner
#   #         legend.justification = c(0.1, 0.9)) # Aligns the legend box to its top-left corner)+
# 
# ggsave(paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/fig_4.png"),
#        width = 10, height = 6.5, units = "in", dpi = 400)


### map of model residuals #####################################################
#project the predictions to the crs of the basemap tile (3857)
fitted_preds_sf <- st_as_sf(fitted_preds, crs = 2263, coords = c("lon", "lat")) %>% 
  st_transform(., crs = 3857) %>% 
  bind_cols(st_coordinates(.) %>% as.data.frame())

# #load in nyc boundary polygon
# nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
#   st_union() %>% #combine the different boroughs
#   st_transform(., crs = 3857)

#load in basemap
nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 

#function to remove bins with lower than a certain n
bin_removal <- function(x) {
  if(length(x) < 10) return(NA) 
  return(mean(x))
}

#create figure
fig_map_resid <- 
  ggplot() + ggthemes::theme_few() +   
  geom_spatraster_rgb(data = nyc_topo_spatrast) +
  stat_summary_hex(data = fitted_preds_sf, aes(x = X, y = Y, z =  canopy_endstate - pred_fit), fun = bin_removal, bins = 20) +
  facet_wrap(~sp_a, ncol = 6) +
  scale_fill_viridis_c(option = "turbo", direction = -1, name = "residuals (%)") + xlab("") + ylab("") + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = c(0.92, 0.14),
        legend.background = element_rect(fill = "white", color = "grey80"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
# annotation_scale(location = "br",  # "bl" for bottom-left, other options exist
#                  bar_cols = c("black", "white"), # Colors of the scale bar segments
#                  style = "ticks",
#                  text_cex = 0.8) +  # Text size for the scale bar label
# annotation_north_arrow(location = "br", height = unit( 0.8, "cm"), style = north_arrow_minimal,
#                        pad_x = unit(1, "cm"), pad_y = unit(1, "cm")) +
# theme(  legend.position = c(0.1, 0.9),  # Places the legend at the top-left corner
#         legend.justification = c(0.1, 0.9)) # Aligns the legend box to its top-left corner)+

ggsave(fig_map_resid, filename = paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/SI_resid_mapb.png"),
       width = 10, height = 6.5, units = "in", dpi = 400)


fitted_preds_sf %>% 
  mutate(resid = canopy_endstate - pred_fit) %>% 
  group_by(sp_a) %>% summarize(mean_resid = mean(resid))
  ggplot(aes(x= resid)) + geom_histogram() + facet_wrap(~sp_a)
future_preds %>% 
  group_by(sp_a) %>% 
  summarize(pred_surv_mean = mean(pred_surv, na.rm = TRUE),
            pred_surv_median = median(pred_surv, na.rm = TRUE))
  

ggplot(future_preds, aes(x = x, y = lat)) + geom_hex() + facet_wrap(~sp_a) + theme_bw() + scale_fill_viridis_c()
ggplot(future_preds, aes(x = x, y = lat, z = pred_surv)) + stat_summary_hex(fun = median, bins = 30) + facet_wrap(~sp_a) + theme_bw() + scale_fill_viridis_c()

# 
# 
# 
# #predicted future mortality risk - should be pretty similar to previous figure
# ggplot(future_preds, aes(x = x, y = lat, z = pred_surv)) + 
#   stat_summary_hex(fun = median, bins = 150) + #facet_wrap(~sp_a) + 
#   scale_fill_viridis_c(option = "turbo", direction = -1, name = "median survival (%)") + xlab("longitude") + ylab("latitude") + 
#   theme_bw() + 
#   theme(strip.text = element_text(face = "italic"),
#         legend.position = c(0.8, 0.1),
#         legend.background = element_rect(fill = "white", color = "grey80"),
#         panel.grid.major = element_blank(),  
#         panel.grid.minor = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())

### SI X: map of tree sample size by species #########################################################
ggplot(future_preds, aes(x = x, y = lat, z = pred_surv)) + 
  geom_hex() + facet_wrap(~sp_a) + theme_bw() + 
  scale_fill_viridis_c(option = "turbo", name = "n") + xlab("longitude") + ylab("latitude") + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = c(0.8, 0.1),
        legend.background = element_rect(fill = "white", color = "grey80"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())



### OLD VERSION FROM DAVE coefficients for select logistic regression variables ######################################
# # Setup plot with facets by variable, with species down the y-axis
# 
# for (i in 1:length(sp_list)){  #for (i in 12:12){
#   print(i)
#   coeff_names <- names(mort_model_list[[i]]$coefficients)
#   beta <- mort_model_list[[i]]$coefficients
#   model_summary <- summary(mort_model_list[[i]])
#   p_value <- model_summary$coefficients[,4]
#   
#   #ci <- confint(mort_model_list[[i]]) # maybe include this inside the model? 
#   # this fails for many models
#   #Error in approx(sp$y, sp$x, xout = cutoff) : 
#   #  need at least two non-NA values to interpolate
#   # leaving out for now, not required for the moment. Will need to revise
#   
#   # merge in case p-value df is too small
#   p_value_df <- data.frame(names(p_value), p_value)
#   coeff_df <- data.frame(coeff_names)
#   p_value_df_merge <- merge(p_value_df, coeff_df, by.x = "names.p_value.", by.y = "coeff_names", all = TRUE, sort = FALSE) # do not sort
#   p_value <- p_value_df_merge$p_value # this has NA values built in?
#   
#   if (i == 1){
#     df_mort_models <- cbind.data.frame(rep(sp_list[i], length(coeff_names)), coeff_names, beta, p_value) #, ci)
#   } else {
#     df_mort_model_sub <- cbind.data.frame(rep(sp_list[i], length(coeff_names)), coeff_names, beta, p_value) #, ci)
#     df_mort_models <- rbind.data.frame(df_mort_models, df_mort_model_sub)
#   }
# }
# 
# colnames(df_mort_models)[1] <- "species"
# df_mort_models$OR <- exp(df_mort_models$beta)
# df_mort_models$sig <- FALSE
# df_mort_models$sig[which(df_mort_models$p_value < 0.05)] <- TRUE
# df_mort_models$direction <- "Insignificant"
# df_mort_models$direction[which(df_mort_models$sig & df_mort_models$OR > 1)] <- "Positive"
# df_mort_models$direction[which(df_mort_models$sig & df_mort_models$OR < 1)] <- "Negative"
# 
# df_mort_models %>% 
#   filter(coeff_names %in% c("imp_10_2021", "steward_level", "in_sandy_zone_boolTRUE", "RPL_THEME1_bin(0.25,0.5]", "RPL_THEME1_bin(0.5,0.75]", "RPL_THEME1_bin(0.75,1]")) %>%
#   ggplot(aes(color = direction)) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   geom_point(aes(x = OR, y = species), size = 2) +
#   #geom_linerange(aes(xmin = lower_95, xmax = upper_95, y = vars_2)) +
#   scale_y_discrete(limits = rev) +
#   scale_color_manual(values = c("gray80", "red", "blue")) +
#   facet_wrap(~coeff_names, nrow = 1) +
#   labs(x = "Odds Ratio", y = "", color = "",
#        title = "Odds ratios by tree species\nNote: genus-only labels are not inclusive of the related species labels") +
#   theme_bw() +
#   theme(strip.text = element_text(face = "italic"))
# # note that genus-only labels here are those trees that did not have a species label, and are not inclusive of the species labels


