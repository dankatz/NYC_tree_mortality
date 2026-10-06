library(INLA)
library(inlabru)
library(fmesher)
library(sf)
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
library(knitr)
library(kableExtra)

sp_focal <- sp_list[6]

sp_sub_format <- tree_vars_full_nona_format %>% filter(sp_a == sp_focal) %>%  #
  group_by(x, y) %>% filter(n() == 1) %>%  ungroup() %>%  #remove some rows that have duplicated location information
  filter(!is.na(y)) %>% 
  sample_n(10000) 
  
  #create covariates
dbh_grid <-  seq(min(sp_sub_format$dbh_cm_z), max(sp_sub_format$dbh_cm_z), length.out = 10)


trees <- sp_sub_format %>% 
  select(canopy_endstate, dbh_cm_z, RPL_THEME1_z, in_sandy_zone_bool, x, y)
trees_sf <- st_as_sf( sp_sub_format, coords = c("x", "y"), remove = FALSE) %>% 
  select(canopy_endstate, dbh_cm_z, RPL_THEME1_z, in_sandy_zone_bool, x, y)




## ---- 2.2 Mesh and SPDE (Matern with PC priors) ------------------------------
#load in nyc boundary polygon
nyc_boundary <- st_read("data/nyc_boundary_polygon/nybb.shp") %>%
  st_union() %>% #combine the different boroughs
  st_transform(., crs = 2263)
nyc.bdry <-  as(nyc_boundary, "Spatial") %>% fm_as_segm()

max.edge <- diff(range(st_coordinates(trees_sf)[,1]))/(3*5) #https://rpubs.com/jafet089/886687
bound.outer <- diff(range(st_coordinates(trees_sf)[,1]))/10


# mesh <- inla.mesh.2d(
#   loc      = as.matrix(trees[, c("x", "y")]),
#   max.edge = c(2000, 10000), #c(1,2)*max.edge,
#   offset   = c(max.edge, bound.outer),
#   cutoff = max.edge/5
# )

mesh <- fm_mesh_2d_inla(
  loc      = as.matrix(trees[, c("x", "y")]),
  max.edge = c(1000, 5000), #c(1,2)*max.edge,
  offset   = c(max.edge, bound.outer),
  cutoff = max.edge/10,
  boundary = nyc.bdry
)


ggplot() + gg(mesh) + geom_sf(data= trees_sf,col='purple',size=1.7,alpha=0.5) + theme_minimal()
cat("Mesh has", mesh$n, "vertices\n")

# PC priors: P(range < 15) = 0.05 ; P(sigma > 1.5) = 0.05
spde <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(15, 0.05),
  prior.sigma = c(1.5, 0.05)
)

## ---- 2.3 Model components and formula ---------------------------------------

cmp <- canopy_endstate ~ Intercept(1) +
  RPL_THEME1_z(RPL_THEME1_z, model = "linear") +
  #dbh_cm_z(dbh_cm_z, model = "rw2") +
  dbh_smooth(
    dbh_cm_z,
    model = "rw2",
    values = dbh_grid,
    scale.model = TRUE,                     # makes precision comparable across resolutions
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))  # PC prior on smoothness
  ) +
  in_sandy_zone_bool(in_sandy_zone_bool, model = "linear") + 
  #habitat(habitat, model = "factor_contrast") +
  spatial_field(geometry, model = spde)

## ---- 2.4 Fit with bru() -------------------------------------------------------
## NOTE: recent versions of inlabru expect the likelihood to be specified via
## like(), and the `options` argument to be a bru_options() object (not a
## plain list). Passing family/data/Ntrials/control.family directly to bru()
## (older syntax) can trigger an internal coercion error of the form:
##   "Error in missing(x) || is.null(x) || is.na(x) :
##    'length = 2' in coercion to 'logical(1)'"
## The like()/bru_options() syntax below avoids this.

lik <- like(
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

summary(fit)



## ============================================================================
## 3. MODEL DIAGNOSTICS
## ============================================================================

## ---- 3.1 Fixed-effect posterior summaries -----------------------------------

fixed_summary <- fit$summary.fixed
print(fixed_summary)

## ---- 3.2 Fixed-effect posterior density plots --------------------------------

fixed_names <- rownames(fit$summary.fixed)

fixed_plots <- lapply(fixed_names, function(nm) {
  df <- as.data.frame(fit$marginals.fixed[[nm]])
  ggplot(df, aes(x = x, y = y)) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(title = nm, x = expression(beta), y = "Density") +
    theme_minimal(base_size = 10)
})

wrap_plots(fixed_plots, ncol = 2)

## Predict the smooth term over the covariate grid
dbh_pred <- predict(
  fit,
  data.frame(dbh_cm_z = dbh_grid),
  ~ dbh_smooth
)

p_smooth <- ggplot(dbh_pred, aes(x = dbh_cm_z, y = mean)) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "steelblue", alpha = 0.25) +
  geom_line(colour = "steelblue", linewidth = 1) +
  labs(title = "Estimated smooth effect of soil moisture",
       x = "DBH (z-scored)", y = "Effect (cloglog scale)") +
  theme_minimal(base_size = 11)

p_smooth

## ---- 3.3 Hyperparameter posteriors (spatial range & SD) -----------------------

range_post <- inla.tmarginal(function(x) x, fit$marginals.hyperpar$`Range for spatial_field`)
sigma_post <- inla.tmarginal(function(x) x, fit$marginals.hyperpar$`Stdev for spatial_field`)

p_range <- ggplot(as.data.frame(range_post), aes(x, y)) +
  geom_line(colour = "darkorange", linewidth = 0.8) +
  #geom_vline(xintercept = true_range, linetype = "dashed") +
  labs(title = "Posterior: spatial range", x = "Range (m)", y = "Density") +
  theme_minimal(base_size = 10)

p_sigma <- ggplot(as.data.frame(sigma_post), aes(x, y)) +
  geom_line(colour = "darkorange", linewidth = 0.8) +
  #geom_vline(xintercept = true_sigma, linetype = "dashed") +
  labs(title = "Posterior: spatial field SD", x = "Sigma", y = "Density") +
  theme_minimal(base_size = 10)

p_range + p_sigma

## ---- 3.4 Model fit statistics: WAIC, DIC, CPO --------------------------------

cat("WAIC:", fit$waic$waic, "\n")
cat("DIC :", fit$dic$dic,  "\n")

# CPO / PIT-based checks: failure values and PIT histogram
cpo_fail <- sum(fit$cpo$failure > 0, na.rm = TRUE)
cat("Number of CPO failure flags:", cpo_fail, "\n")

pit_df <- data.frame(pit = fit$cpo$pit)
p_pit <- ggplot(pit_df, aes(x = pit)) +
  geom_histogram(bins = 20, fill = "grey60", colour = "white") +
  labs(title = "PIT histogram (should be approx. uniform)",
       x = "Probability Integral Transform", y = "Count") +
  theme_minimal(base_size = 10)

p_pit




## ---- 3.5 Posterior predictive / binned residual check ------------------------

trees$fitted_p <- fit$summary.fitted.values$mean[1:10000]
trees$resid    <- trees$canopy_endstate - trees$fitted_p

#area under the curve, example: mort_model_list[[2]]$roc_curve$auc
roc_curve <- roc(trees$canopy_endstate, trees$fitted_p)
roc_curve 



# Binned residual plot: mean residual within bins of fitted probability
n_bins <- 10
trees$bin <- cut(trees$fitted_p, breaks = quantile(trees$fitted_p,
                                                   probs = seq(0, 1, length.out = n_bins + 1)),
                 include.lowest = TRUE)

binned <- trees %>%
  group_by(bin) %>%
  summarise(
    mean_fitted = mean(fitted_p),
    mean_obs    = mean(canopy_endstate),
    n           = n(),
    se          = sqrt(mean_fitted * (1 - mean_fitted) / n)
  )

p_binned <- ggplot(binned, aes(x = mean_fitted, y = mean_obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = mean_obs - 1.96 * se, ymax = mean_obs + 1.96 * se),
                width = 0.01, colour = "grey60") +
  geom_point(size = 2, colour = "steelblue") +
  labs(title = "Binned residual plot",
       x = "Mean predicted survival probability",
       y = "Observed proportion surviving") +
  theme_minimal(base_size = 11)

p_binned




## ---- 3.6 Mesh diagnostic plot -------------------------------------------------

p_mesh <- ggplot() +
  gg(mesh) +
  geom_point(data = trees, aes(x = x, y = y), size = 0.6, colour = "red", alpha = 0.6) +
  coord_equal() +
  labs(title = "SPDE mesh and tree locations") +
  theme_minimal(base_size = 10)

p_mesh

## ============================================================================
## 4. PUBLICATION-QUALITY RESULTS
## ============================================================================

## ---- 4.1 Coefficient table (log-cloglog scale, with 95% CrI) ------------------

coef_table <- fixed_summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term") %>%
  transmute(
    Term,
    Estimate = round(mean, 3),
    `95% CrI lower` = round(`0.025quant`, 3),
    `95% CrI upper` = round(`0.975quant`, 3)
  ) %>%
  mutate(Term = recode(Term,
                       Intercept = "Intercept",
                       soil_moisture = "Soil moisture (z-scored)",
                       dbh_std = "DBH (z-scored)",
                       `habitat1` = "Habitat: Mixed vs Conifer",
                       `habitat2` = "Habitat: Broadleaf vs Conifer"
  ))

# Render as a nicely formatted table for a manuscript
kable(coef_table, format = "html", caption =
        "Table 1. Posterior summaries (mean and 95% credible intervals) for fixed effects from the binomial cloglog spatial model of tree survival.") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))

## ---- 4.2 Forest plot of fixed effects ----------------------------------------

p_forest <- ggplot(coef_table[coef_table$Term != "Intercept", ],
                   aes(x = Estimate, y = reorder(Term, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = `95% CrI lower`, xmax = `95% CrI upper`),
                 height = 0.15, colour = "grey30") +
  geom_point(size = 3, colour = "steelblue") +
  labs(
    title = "Figure 1. Effects of covariates on tree survival",
    subtitle = "Binomial GLM with cloglog link (posterior mean \u00b1 95% CrI)",
    x = "Coefficient (cloglog scale)", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p_forest

## ---- 4.3 Spatial random field map ---------------------------------------------

field_pred <- predict(
  fit,
  fm_pixels(mesh, mask = NULL),
  ~ spatial_field
) %>% 
  cbind(., st_coordinates(.))


p_field <- ggplot() +
  gg(field_pred, aes(color = mean)) +
  scale_color_viridis_c(name = "Spatial\neffect") + #scale_color_viridis_c() +
  geom_point(data = trees, aes(x = x, y = y), size = 0.4, colour = "white", alpha = 0.5) +
  #coord_equal() +
  labs(title = "Figure 2. Estimated spatial random effect (SPDE)",
       x = "Easting (m)", y = "Northing (m)") +
  theme_minimal(base_size = 11)

p_field

# ggplot() +
#   gg(field_pred, aes(color = mean)) +
#   scale_color_viridis_c(name = "Spatial\neffect") + #scale_color_viridis_c() +
#   #geom_point(data = trees, aes(x = x, y = y), size = 0.4, colour = "white", alpha = 0.5) +
#   coord_sf(xlim = c(1000000, 1090000), ylim = c(200000, 290000)) + 
#   labs(title = "Figure 2. Estimated spatial random effect (SPDE)",
#        x = "Easting (m)", y = "Northing (m)") +
#   theme_minimal(base_size = 11)

## ---- 4.4 Predicted survival probability map (for a reference habitat/DBH) ----

# Predict the full linear predictor including fixed effects, holding
# soil_moisture and dbh_std at observed surface values is not possible on a
# pixel grid (those are point covariates), so here we map the combined
# spatial contribution to the survival probability holding other covariates
# at their mean (0) and habitat at the reference level (Conifer).

pred_p <- predict(
  fit,
  fm_pixels(mesh, mask = NULL),
  ~ 1 - exp(-exp(Intercept + spatial_field))
)

p_prob <- ggplot() +
  gg(pred_p, aes(color = mean)) +
  scale_color_viridis_c(name = "Predicted\nP(survival)") + #limits = c(0, 1)
  #coord_equal() +
  labs(title = "Figure 3. Predicted survival probability surface",
       subtitle = "Conifer habitat, mean soil moisture and DBH",
       x = "Easting (m)", y = "Northing (m)") +
  theme_minimal(base_size = 11)

p_prob

## ---- 4.5 Map of observed survival outcomes -------------------------------------

p_obs <- ggplot(trees, aes(x = x, y = y, colour = factor(canopy_endstate))) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = c("0" = "firebrick", "1" = "forestgreen"),
                      labels = c("Dead", "Alive"), name = "Status") +
  #coord_equal() +
  labs(title = "Observed tree survival status", x = "Easting (m)", y = "Northing (m)") +
  theme_minimal(base_size = 11)

p_obs

## ============================================================================
## 5. SAVE OUTPUTS (optional)
## ============================================================================

# ggsave("forest_plot.png", p_forest, width = 7, height = 4, dpi = 300)
# ggsave("spatial_field.png", p_field, width = 6, height = 5, dpi = 300)
# ggsave("predicted_survival.png", p_prob, width = 6, height = 5, dpi = 300)
# ggsave("binned_residuals.png", p_binned, width = 5, height = 5, dpi = 300)
# saveRDS(fit, "tree_survival_inla_fit.rds")

## ============================================================================
## END OF SCRIPT
## ============================================================================