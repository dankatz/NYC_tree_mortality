# Street tree mortality estimates, models, and risk maps for New York City
# this script is a logistic regression analysis of tree mortality
# original script by David Miller with revision by Dan Katz


#### set up work environment and load data #####################################
library(tidyverse)
library(data.table)
library(ResourceSelection) # hoslem.test
library(pROC)
library(ciTools) #install.packages("here")
library(sf)
library(basemaps)
library(terra)
library(here)
library(gt)

your_path_for_box <- "C:/Users/dsk273/Box/Katz lab/NYC/"
tree_vars <- fread(paste0(your_path_for_box, '/tree_mortality_variables/all_variables.csv'))
socioeco_vars <- fread(paste0(your_path_for_box,'/tree_mortality_variables/model_df_socioeco.csv'))


# merge tree and socioeconomic variables
socioeco_vars_select <- socioeco_vars %>% select(tree_id, bg_area_sqft, medincomeE, estimate_c_perc_unoccupied, estimate_c_perc_non_white, estimate_c_perc_unemployed, estimate_c_perc_poverty,
                                                 RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME4, RPL_THEMES)
tree_vars_full <- left_join(tree_vars, socioeco_vars_select)



### create derived variables ###################################################

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
  tree_vars_full <- tree_vars_full %>% filter(canopy_change %in% c(1, 3)) %>%
    mutate(canopy_endstate = (canopy_change - 3)*(-1/2),
           LandUse_char = sapply(LandUse, appendCharLU))
  tree_vars_full$LandUse_char <- relevel(as.factor(tree_vars_full$LandUse_char), "LU_01") # this is Mike's recommendation for base case, may need to collapse land use classes
  
  # Checking SVI variables
  # need to include SVI variables in here, work with the first three
  # RPL_THEME1: socioeconomic status
  # RPL_THEME2: household characteristics
  # RPL_THEME3: racial and ethnic minority status
  # RPL_THEME4: housing type and transportation (don't include this one)
  # RPL_THEMES: summary of all the above
  
  # sp <- "Quercus palustris" #"Acer platanoides"
  # tree_vars_full %>% filter(species == sp) %>%
  #   ggplot(aes(x = RPL_THEME1, y = RPL_THEMES)) +
  #   #ggplot(aes(x = RPL_THEME2, y = RPL_THEME3)) +
  #   geom_point() +
  #   geom_smooth(method = "lm")
  
  # RPL_THEME1 is very correlated with RPL_THEMES
  cor(tree_vars_full$RPL_THEME1, tree_vars_full$RPL_THEMES, use = "complete.obs") # R = 0.90
  # others are not as correlated
  cor(tree_vars_full$RPL_THEME2, tree_vars_full$RPL_THEMES, use = "complete.obs") # R = 0.64
  cor(tree_vars_full$RPL_THEME3, tree_vars_full$RPL_THEMES, use = "complete.obs") # R = 0.66
  cor(tree_vars_full$RPL_THEME4, tree_vars_full$RPL_THEMES, use = "complete.obs") # R = 0.64
  # check relationship of 1, 2, and 3
  cor(tree_vars_full$RPL_THEME1, tree_vars_full$RPL_THEME2, use = "complete.obs") # R = 0.49
  cor(tree_vars_full$RPL_THEME1, tree_vars_full$RPL_THEME3, use = "complete.obs") # R = 0.63
  cor(tree_vars_full$RPL_THEME2, tree_vars_full$RPL_THEME3, use = "complete.obs") # R = 0.37
  # 1 and 3 are somewhat correlated and may need to remove race & ethnicity but that's tbd
  
  tree_vars_full_nona <- na.omit(tree_vars_full) # This is a *big* hit, from 684k to 496k. Will have to see if RPL_THEMEx is worth it for predictive power.
  # One big thing to note is that SVI is not available for parks and certain public land uses, so will limit applicability in certain parts of the city
  
  # reasonable distance threshold for land use proximity? This is in survey foot (effectively feet)
  quantile(tree_vars_full_nona$pluto_dist, seq(0, 1, 0.01))
  # picking 20, this is approx 95% of the remaining trees. Will apply this in the mutate part
  
  # Collapsed land use
  tree_vars_full_nona$LandUse_char_collapsed <- tree_vars_full_nona$LandUse_char %>% as.character()
  x <- tree_vars_full_nona$LandUse_char_collapsed %>%
    case_match( "LU_01" ~ "LowDensityResidential",
                c("LU_02", "LU_03") ~ "HigherDensityResidential",
                c("LU_04", "LU_05") ~ "CommercialMix",
                c("LU_06", "LU_07", "LU_10") ~ "IndustryTransport",
                "LU_08" ~ "PublicInst",
                "LU_09" ~ "OpenOutdoorRec",
                "LU_11" ~ "Vacant")
  tree_vars_full_nona$LandUse_char_collapsed <- as.factor(x)
  tree_vars_full_nona$LandUse_char_collapsed <- relevel(as.factor(tree_vars_full_nona$LandUse_char_collapsed), "LowDensityResidential") # this is Mike's recommendation for base case, may need to collapse land use classes


### create dataframe for analysis ##############################################
# Format input df as needed
tree_vars_full_nona_format <- tree_vars_full_nona %>%
  mutate(dbh_cm = tree_dbh * 2.54) %>%  #convert dbh from inches to cm
  filter(dbh_cm > 7.62) %>%  # remove small trees; 7.6 cm is 3 inches listed in Bigelow
  filter(tree_dbh < 118) %>%  #removing trees with a DBH greater than the maximum ever recorded in NYC 
  filter(pluto_dist < 20) %>% # tree distance 
  mutate(BldgClass_fac = factor(BldgClass), 
         BldgClass_Group_fac = factor(BldgClass_Group),
         in_sandy_zone_bool = as.logical(in_sandy_zone),
         is_B_cons_bool = as.logical(is_B_cons),
         is_S_cons_bool = as.logical(is_S_cons),
         is_DM_bool = as.logical(is_DM)) %>%
  select(tree_id, geom, 
         genus, species, 
         canopy_endstate, # dependent variable 
         tree_dbh, dbh_cm, steward_level, # tree diameter at breast height, number of signs of stewardship
         BldgClass_fac, BldgClass_Group_fac, LandUse_char, LandUse_char_collapsed, # building type (detailed), building type (high level), land use
         in_sandy_zone_bool, # sandy inundation zone
         is_B_cons_bool, is_S_cons_bool, is_DM_bool, # building construction, street construction, building demolition
         summer_max, summer_min, summer_mean, days_max_27, days_max_32, days_max_35, # temperature
         TC_10_diff, GS_10_diff, SO_10_diff, WA_10_diff, BD_10_diff, RD_10_diff, OI_10_diff, RR_10_diff, # land cover diff, doing 10 m
         TC_10_2017, GS_10_2017, SO_10_2017, WA_10_2017, BD_10_2017, RD_10_2017, OI_10_2017, RR_10_2017, # land cover 2017, doing 10 m
         TC_10_2021, GS_10_2021, SO_10_2021, WA_10_2021, BD_10_2021, RD_10_2021, OI_10_2021, RR_10_2021, # land cover 2021, doing 10 m
         RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME3, RPL_THEME4, RPL_THEMES) %>%
  mutate(imp_10_2021 = BD_10_2017 + RD_10_2017 + OI_10_2017 + RR_10_2017) %>% 
    separate_wider_delim(., cols = geom, delim = ",", names = c( "lon", "lat")) %>%  #extract coordinates from geom text string 
    mutate(lat = readr::parse_number(lat),
           lon = readr::parse_number(lon))
   
    
# do bins of tree dbh
#quantile(tree_vars_full_nona_format$tree_dbh, seq(0,1,0.1))
# Remove small trees
tree_vars_full_nona_format$tree_dbh_bin <- cut(tree_vars_full_nona_format$tree_dbh, c(3, 6, 9, 12, 18, 24, 500))
tree_vars_full_nona_format$tree_dbh_bin %>% table()
# 3 inches is 7.6 cm listed in Bigelow

# Social vulnerability can also be binned
tree_vars_full_nona_format$RPL_THEME1_bin <- cut(tree_vars_full_nona_format$RPL_THEME1, c(0, 0.25, 0.5, 0.75, 1), include.lowest = TRUE) # just doing these bins for now
tree_vars_full_nona_format$RPL_THEME2_bin <- cut(tree_vars_full_nona_format$RPL_THEME2, c(0, 0.25, 0.5, 0.75, 1), include.lowest = TRUE) # just doing these bins for now
tree_vars_full_nona_format$RPL_THEME3_bin <- cut(tree_vars_full_nona_format$RPL_THEME3, c(0, 0.25, 0.5, 0.75, 1), include.lowest = TRUE) # just doing these bins for now
tree_vars_full_nona_format$RPL_THEME4_bin <- cut(tree_vars_full_nona_format$RPL_THEME4, c(0, 0.25, 0.5, 0.75, 1), include.lowest = TRUE) # *won't use this but for completeness
tree_vars_full_nona_format$RPL_THEMES_bin <- cut(tree_vars_full_nona_format$RPL_THEMES, c(0, 0.25, 0.5, 0.75, 1), include.lowest = TRUE) # just doing these bins for now

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
                          TRUE ~ "other"))



### analyzing per species mortality with a glm #####################################
sp_list <- unique(tree_vars_full_nona_format$sp_a) %>% sort()

#where to save output
mort_model_list <- vector("list", length(sp_list))
fitted_preds_list <- vector("list", length(sp_list))
future_preds_list <- vector("list", length(sp_list))

for (i in 1:length(sp_list)){
 # for (i in 12:12){
  print(i)
  sp <- sp_list[i] #sp <- sp_list[12]
  sp_sub_format <- tree_vars_full_nona_format %>% filter(species == sp_a)
  mort_model <- glm( canopy_endstate ~ imp_10_2021 + tree_dbh_bin + is_DM_bool + steward_level + in_sandy_zone_bool + LandUse_char_collapsed + 
                       RPL_THEME1_bin + RPL_THEME2_bin + RPL_THEME3_bin +
                       imp_10_2021 + steward_level + is_DM_bool, #tree_dbh_bin:is_DM_bool,
                     family = binomial(link = "cloglog"),
                     na.action = na.exclude,
                     data = sp_sub_format)
  mort_model$species <- sp
  
  # not every tree gets modeled, may be an issue with the input predictor variables? Using default g value of 10
  # p.value, want it to be > 0.05 (NOT significant), example: mort_model_list[[2]]$hoslem$p.value
    mort_model$hoslem <- hoslem.test(mort_model$model$canopy_endstate, mort_model$fitted.values)

  # area under the curve, example: mort_model_list[[2]]$roc_curve$auc
    mort_model$roc_curve <- roc(mort_model$model$canopy_endstate, mort_model$fitted.values)

  #save mortality model
  mort_model_list[[i]] <- mort_model
  
  # can save this as an RDS file (if needed, this runs quickly)
  # saveRDS()
  
  ## extract fitted values from model -------------------------------
  summary(mort_model)
  mort_model$fitted.values
  
  tree_vars_fitted_focal <- 
  tree_vars_full_nona_format %>% 
    filter( sp_a == sp) %>% 
    select(-canopy_endstate) %>% 
    mutate(pred_fit = predict.glm(object = mort_model, newdata = ., type = "response")) %>% 
    select(tree_id, lat, lon, species, sp_a, tree_dbh, pred_fit)
  
      #save predictions for the focal species
      fitted_preds_list[[i]] <- tree_vars_fitted_focal
  
  
  ## extract predictions of future survival from model -------------
  tree_vars_predict_focal <- tree_vars_full_nona_format %>% 
    filter(canopy_endstate == 1) %>% 
    filter( sp_a == sp) %>% 
    select(-canopy_endstate) %>% 
    mutate(pred_surv = predict.glm(object = mort_model, newdata = ., type = "response")) %>% 
    select(tree_id, lat, lon, species, sp_a, tree_dbh, pred_surv)
    #pred_surv = add_ci(data = ., model =  mort_model, type = "response", alpha = 0.05)
  
    #save predictions for the focal species
    future_preds_list[[i]] <- tree_vars_predict_focal
  
  # ## extract residuals from model
  #   tree_vars_resid_focal <-
  #     tree_vars_full_nona_format %>% 
  #     mutate(resid_response = residuals(mort_model, type = "response"))
  #   
  #   #save predictions for the focal species
  #   resids_preds_list[[i]] <- tree_vars_resid_focal
    
}

#save survival predictions
fitted_preds <- bind_rows(fitted_preds_list)
future_preds <- bind_rows(future_preds_list) 



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
    locations = cells_column_labels()) %>% 
  gtsave( paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/table1.docx"))  




### Fig 2: survival by species #######################################################
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
  ggplot(aes(x = median_dbh, y = mean_annual_mort, color = sp_a)) + geom_point() + geom_line()+ theme_bw() + xlab("DBH (cm)") + ylab("annual mortality (%)") +
  scale_color_viridis_d(option = "turbo", name = "species") +
  theme(panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),
      legend.text = element_text(face = "italic"))      

ggsave(paste0(your_path_for_box, "tree_mortality/NYC_st_tree_results/fig_2.png"),
       width = 7, height = 5, units = "in", dpi = 400)



### Fig 4: #######################################################

#load in nyc boundary polygon
nyc_boundary <- st_read( "C:/Users/dsk273/Box/Katz lab/NYC/nyc_boundary_polygon/nybb.shp") %>% 
  st_union() %>% #combine the different boroughs
  st_transform(., crs = 2263)
#nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)), crs= 32618)
# nyc_boundary_invert <- st_difference(nyc_boundary_box, nyc_boundary)

nyc_topo_rast <- basemap_raster(nyc_boundary, map_service = "carto", map_type = "light_no_labels") #basemap_raster(nyc_boundary, map_service = "esri", map_type = "world_hillshade")
nyc_topo_spatrast <- rast(nyc_topo_rast) #convert to spatrast for plotting 

ggplot() + ggthemes::theme_few() +   
  geom_spatraster_rgb(data = nyc_topo_spatrast) +
  #geom_spatraster(data = prod_400m_focal_sum/1000, alpha = 0.6) +
  scale_fill_viridis_c(na.value = "transparent", 
                       #option = "magma",
                       name = "pollen production \n(trillions of grains within 400 m)",
                       labels = scales::label_comma()) +
  annotation_scale(location = "br",  # "bl" for bottom-left, other options exist
                   bar_cols = c("black", "white"), # Colors of the scale bar segments
                   style = "ticks",
                   text_cex = 0.8) +  # Text size for the scale bar label
  annotation_north_arrow(location = "br", height = unit( 0.8, "cm"), style = north_arrow_minimal,
                         pad_x = unit(1, "cm"), pad_y = unit(1, "cm")) +
  theme(  legend.position = c(0.1, 0.9),  # Places the legend at the top-left corner
          legend.justification = c(0.1, 0.9)) # Aligns the legend box to its top-left corner)+



future_preds %>% 
  group_by(sp_a) %>% 
  summarize(pred_surv_mean = mean(pred_surv, na.rm = TRUE),
            pred_surv_median = median(pred_surv, na.rm = TRUE))
  
  
ggplot(future_preds, aes(x = tree_dbh, y = pred_surv)) + geom_hex() + facet_wrap(~sp_a) + theme_bw()


ggplot(future_preds, aes(x = lon, y = lat)) + geom_hex() + facet_wrap(~sp_a) + theme_bw() + scale_fill_viridis_c()

ggplot(future_preds, aes(x = lon, y = lat, z = pred_surv)) + stat_summary_hex(fun = median, bins = 30) + facet_wrap(~sp_a) + theme_bw() + scale_fill_viridis_c()

#modeled mortality risk for trees
ggplot(fitted_preds, aes(x = lon, y = lat, z = pred_fit)) + 
  stat_summary_hex(fun = median, bins = 40) + facet_wrap(~sp_a) + theme_bw() + 
  scale_fill_viridis_c(option = "turbo", direction = -1, name = "median survival (%)") + xlab("longitude") + ylab("latitude") + 
  theme(strip.text = element_text(face = "italic"),
    legend.position = c(0.8, 0.1),
    legend.background = element_rect(fill = "white", color = "grey80"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank())


#predicted future mortality risk - should be pretty similar to previous figure
ggplot(future_preds, aes(x = lon, y = lat, z = pred_surv)) + 
  stat_summary_hex(fun = median, bins = 150) + #facet_wrap(~sp_a) + 
  scale_fill_viridis_c(option = "turbo", direction = -1, name = "median survival (%)") + xlab("longitude") + ylab("latitude") + 
  theme_bw() + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = c(0.8, 0.1),
        legend.background = element_rect(fill = "white", color = "grey80"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

#SI X: map of tree sample size by species
ggplot(future_preds, aes(x = lon, y = lat, z = pred_surv)) + 
  geom_hex() + facet_wrap(~sp_a) + theme_bw() + 
  scale_fill_viridis_c(option = "turbo", name = "n") + xlab("longitude") + ylab("latitude") + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = c(0.8, 0.1),
        legend.background = element_rect(fill = "white", color = "grey80"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())



### Fig 3: coefficients for select logistic regression variables ######################################
# Setup plot with facets by variable, with species down the y-axis

for (i in 1:length(sp_list)){  #for (i in 12:12){
  print(i)
  coeff_names <- names(mort_model_list[[i]]$coefficients)
  beta <- mort_model_list[[i]]$coefficients
  model_summary <- summary(mort_model_list[[i]])
  p_value <- model_summary$coefficients[,4]
  
  #ci <- confint(mort_model_list[[i]]) # maybe include this inside the model? 
  # this fails for many models
  #Error in approx(sp$y, sp$x, xout = cutoff) : 
  #  need at least two non-NA values to interpolate
  # leaving out for now, not required for the moment. Will need to revise
  
  # merge in case p-value df is too small
  p_value_df <- data.frame(names(p_value), p_value)
  coeff_df <- data.frame(coeff_names)
  p_value_df_merge <- merge(p_value_df, coeff_df, by.x = "names.p_value.", by.y = "coeff_names", all = TRUE, sort = FALSE) # do not sort
  p_value <- p_value_df_merge$p_value # this has NA values built in?
  
  if (i == 1){
    df_mort_models <- cbind.data.frame(rep(sp_list[i], length(coeff_names)), coeff_names, beta, p_value) #, ci)
  } else {
    df_mort_model_sub <- cbind.data.frame(rep(sp_list[i], length(coeff_names)), coeff_names, beta, p_value) #, ci)
    df_mort_models <- rbind.data.frame(df_mort_models, df_mort_model_sub)
  }
}

colnames(df_mort_models)[1] <- "species"
df_mort_models$OR <- exp(df_mort_models$beta)
df_mort_models$sig <- FALSE
df_mort_models$sig[which(df_mort_models$p_value < 0.05)] <- TRUE
df_mort_models$direction <- "Insignificant"
df_mort_models$direction[which(df_mort_models$sig & df_mort_models$OR > 1)] <- "Positive"
df_mort_models$direction[which(df_mort_models$sig & df_mort_models$OR < 1)] <- "Negative"

df_mort_models %>% 
  filter(coeff_names %in% c("imp_10_2021", "steward_level", "in_sandy_zone_boolTRUE", "RPL_THEME1_bin(0.25,0.5]", "RPL_THEME1_bin(0.5,0.75]", "RPL_THEME1_bin(0.75,1]")) %>%
  ggplot(aes(color = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(aes(x = OR, y = species), size = 2) +
  #geom_linerange(aes(xmin = lower_95, xmax = upper_95, y = vars_2)) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("gray80", "red", "blue")) +
  facet_wrap(~coeff_names, nrow = 1) +
  labs(x = "Odds Ratio", y = "", color = "",
       title = "Odds ratios by tree species\nNote: genus-only labels are not inclusive of the related species labels") +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
# note that genus-only labels here are those trees that did not have a species label, and are not inclusive of the species labels


#####

#cut(tree_vars_full$tree_dbh, 10) %>% table()

#####
# Cross validation steps?


# 
# #####
# 
# 
# # Acer platanoides has significant negative relationship with RPL_THEME3 (race)
# # Quercus palustris has significant positive relationship with RPL_THEME1 (socioeconomic), which is backwards? Sig negative with RPL_THEME2 household characteristics
# 
# # can add in possible interactions too
# 
# tree_vars_full_nona_format %>% filter(species %in% c("Acer platanoides", "Quercus palustris")) %>%
# ggplot() +
#   geom_point(aes(x = RPL_THEME1, y = RPL_THEME3, color = canopy_endstate), shape = 1) +
#   facet_wrap(~species)
# 
# 
# 
# # OLD CODE #
# #####
# for(i in 1:length(sp_list)){
#   
#   sp <- sp_list[i]
#   print(sp)
#   
#   sp_sub <- tree_vars %>% 
#     filter(species == sp & canopy_change %in% c(1, 3)) %>%
#     mutate(canopy_endstate = (canopy_change - 3)*(-1/2), # convert 1 and 3 to 1 and 0, where 1 is alive and 0 is dead. This is because the models in the literature are focused on likelihood of survival
#            LandUse_char = sapply(LandUse, appendCharLU))
#   
#   sp_sub <- sp_sub %>% 
#     filter(LandUse_char != "LU_NA")
#   
#   #sp_sub$LandUse_char <- relevel(as.factor(sp_sub$LandUse_char), "LU_NA") # make NA the default level, is this the best approach?
#   sp_sub$LandUse_char <- relevel(as.factor(sp_sub$LandUse_char), "LU_09") # "LU09_Open_Space_and_Outdoor_Recreation" is base case
#   
#   #####
#   # Can prefilter for pluto_dist and xis_dist if needed
#   
#   # Format input df as needed
#   sp_sub_format <- sp_sub %>% 
#     mutate(BldgClass_fac = factor(BldgClass), 
#            BldgClass_Group_fac = factor(BldgClass_Group),
#            in_sandy_zone_bool = as.logical(in_sandy_zone),
#            is_B_cons_bool = as.logical(is_B_cons),
#            is_S_cons_bool = as.logical(is_S_cons),
#            is_DM_bool = as.logical(is_DM)) %>%
#     select(canopy_endstate, # dependent variable - does this need to be factor??
#            tree_dbh, steward_level, # tree diameter at breast height, number of signs of stewardship
#            BldgClass_fac, BldgClass_Group_fac, LandUse_char, # building type (detailed), building type (high level), land use
#            in_sandy_zone_bool, # sandy inundation zone
#            is_B_cons_bool, is_S_cons_bool, is_DM_bool, # building construction, street construction, building demolition
#            summer_max, summer_min, summer_mean, days_max_27, days_max_32, days_max_35, # temperature
#            TC_10_diff, GS_10_diff, SO_10_diff, WA_10_diff, BD_10_diff, RD_10_diff, OI_10_diff, RR_10_diff, # land cover diff, doing 10 m
#            TC_10_2017, GS_10_2017, SO_10_2017, WA_10_2017, BD_10_2017, RD_10_2017, OI_10_2017, RR_10_2017, # land cover 2017, doing 10 m
#            TC_10_2021, GS_10_2021, SO_10_2021, WA_10_2021, BD_10_2021, RD_10_2021, OI_10_2021, RR_10_2021) # land cover 2021, doing 10 m
#   
#   # Subset and organize similar to the variables we would want based on 
#   sp_sub_format <- sp_sub_format %>% mutate(imp_10_2021 = BD_10_2017 + RD_10_2017 + OI_10_2017 + RR_10_2017) # building + road + other impervious + railroad
#   
#   # null_model <- glm( canopy_endstate ~ 1,
#   #                    family = binomial(link = "logit"),
#   #                    data = sp_sub_format)
#   # 
#   # mort_model_acpl <- glm( canopy_endstate ~ imp_10_2021 + tree_dbh + is_B_cons_bool,
#   #                    family = binomial(link = "logit"),
#   #                    data = sp_sub_format)
#   # 
#   # glm( canopy_endstate ~ imp_10_2021 + tree_dbh + is_DM_bool, #is_B_cons_bool,
#   #      family = binomial(link = "logit"),
#   #      data = sp_sub_format) %>% summary()
#   # # is_DM_bool has slightly smaller p value (0.11) than is_B_cons_bool (0.17)
#   # 
#   # glm( canopy_endstate ~ imp_10_2021 + tree_dbh + is_DM_bool + steward_level,
#   #      family = binomial(link = "logit"),
#   #      data = sp_sub_format) %>% summary()
#   # # steward level is highly significant too
#   # 
#   # glm( canopy_endstate ~ imp_10_2021 + tree_dbh + is_DM_bool + steward_level + LandUse_char,
#   #      family = binomial(link = "logit"),
#   #      data = sp_sub_format) %>% summary()
#   # # sig land use classes, 04 mixed residential & commercial, 05 commercial and office buildings
#   # # almost sig (<0.10): 01 one and two family buildings, 06 industrial and manufacturing, 08 public facilities & institutions 
#   # 
#   # sp_sub_format$canopy_endstate %>% table()
#   
#   mort_model <- glm( canopy_endstate ~ imp_10_2021 + tree_dbh + is_DM_bool + steward_level + in_sandy_zone_bool + LandUse_char,
#                      family = binomial(link = "logit"),
#                      #family = binomial(link = "cloglog"),
#                      data = sp_sub_format)
#   
#   mort_model_summary <- summary(mort_model)
#   
#   coef <- data.frame(mort_model_summary$coefficients)
#   coef$OR <- exp(coef$Estimate)
#   ci <- confint(mort_model)
#   ci_or <- exp(ci)
#   
#   coef <- cbind.data.frame(coef, ci_or)
#   colnames(coef) <- c("beta", "se", "Z", "p_value", "OR", "lower_95", "upper_95")
#   coef$vars <- rownames(coef)
#   # coef$vars_2 <- c("Intercept", "Impervious_10m_Pct", "DBH", "Demolition", "Stewardship_Level", "In_Sandy_Zone",
#   #                  "LU01_One_and_Two_Family_Buildings", "LU02_MultiFamily_WalkUp_Buildings", "LU03_MultiFamily_Elevator_Buildings", "LU04_Mixed_Residential_Commercial_Buildings",
#   #                  "LU05_Commercial_and_Office_Buildings", "LU06_Industrial_and_Manufacturing", "LU07_Transportation_and_Utility", "LU08_Public_Facilities_and_Institutions",
#   #                  "LU09_Open_Space_and_Outdoor_Recreation", "LU10_Parking_Facilities", "LU11_Vacant_Land")
#   coef$vars_2 <- c("Intercept", "Impervious_10m_Pct", "DBH", "Demolition", "Stewardship_Level", "In_Sandy_Zone",
#                    "LU01_One_and_Two_Family_Buildings", "LU02_MultiFamily_WalkUp_Buildings", "LU03_MultiFamily_Elevator_Buildings", "LU04_Mixed_Residential_Commercial_Buildings",
#                    "LU05_Commercial_and_Office_Buildings", "LU06_Industrial_and_Manufacturing", "LU07_Transportation_and_Utility", "LU08_Public_Facilities_and_Institutions",
#                     "LU10_Parking_Facilities", "LU11_Vacant_Land") # "LU09_Open_Space_and_Outdoor_Recreation" is base case
#   coef$sig <- FALSE
#   coef$sig[which(coef$p_value < 0.05)] <- TRUE
#   coef$direction <- "Insignificant"
#   coef$direction[which(coef$sig & coef$OR > 1)] <- "Positive"
#   coef$direction[which(coef$sig & coef$OR < 1)] <- "Negative"
#   coef$species <- sp
#   if (i == 1) {
#     coef_all <- coef
#   } else {
#     coef_all <- rbind.data.frame(coef_all, coef)
#   }
# }
# 
# coef_all_plotting <- coef_all %>% filter(vars_2 != "Intercept")
# 
# coef_all_plotting %>%
#   ggplot(aes(color = direction)) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   geom_point(aes(x = OR, y = vars_2), size = 2) +
#   geom_linerange(aes(xmin = lower_95, xmax = upper_95, y = vars_2)) +
#   scale_y_discrete(limits = rev) +
#   scale_color_manual(values = c("gray80", "red", "blue")) +
#   facet_wrap(~species) +
#   labs(x = "Odds Ratio", y = "", color = "") +
#   theme_bw() +
#   theme(strip.text = element_text(face = "italic"))
# 
# ggsave("/Users/dsk273/dsk273_files/nyc_tree_mortality/figures/logistic_regression_odds_ratio_by_species_example_logit.jpg",
#       width = 12, height = 7, units = "in")
# 
# # ggsave("/Users/dsk273/dsk273_files/nyc_tree_mortality/figures/logistic_regression_odds_ratio_by_species_example_cloglog.jpg",
# #        width = 12, height = 7, units = "in")
# 
# # Need to get model fit parameters for these as well
# # Record: AIC, AUC, Hosmer-Lemeshow
# 
# #####
# # Breakdown of the dataset by species and genus, with mortality rates
# sp_change_counts <- tree_vars %>% 
#   group_by(species, canopy_change) %>% 
#   dplyr::summarize(sp_cc_counts = length(canopy_change)) %>% 
#   pivot_wider(id_cols = species, names_from = canopy_change, values_from = sp_cc_counts)
# sp_change_counts <- sp_change_counts %>% 
#   mutate(total_alive_dead = `1` + `3`,
#          pct_alive = round(`1`/(`1` + `3`)*100,4))
# 
# sp_change_counts <- sp_change_counts[order(sp_change_counts$total_alive_dead, decreasing = TRUE),] %>% drop_na()
# 
# library(ggrepel)
# 
# ggplot() +
#   geom_point(data = sp_change_counts, aes(x = total_alive_dead, y = pct_alive), color = "gray80") +
#   geom_point(data = sp_change_counts[1:12,], aes(x = total_alive_dead, y = pct_alive, color = species)) +
#   geom_text_repel(data = sp_change_counts[1:12,], aes(x = total_alive_dead, y = pct_alive, color = species, label = species)) +
#   labs(x = "Total Alive + Dead", y = "Survival %") +
#   theme_bw() +
#   theme(legend.position = "none")
# ggsave("/Users/dsk273/dsk273_files/nyc_tree_mortality/figures/top_species_survival_counts.jpg",
#        width = 6, height = 5, units = "in")
