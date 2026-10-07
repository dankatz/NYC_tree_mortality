# Analysis plan: NYC street tree mortality

Status: **draft, updated 2026-10-06**. Items marked **[ASK G-Y]** are questions about the data for Ging-Yan. The model specification in section 3 is implemented in `lidar_mort_analysis.R`.

This plan sets out how the analysis in `lidar_mort_analysis.R` will be reworked. It is also the source for rewriting the Methods. Data: `data/tree_mortality_variables/all_variables_v4.rds`, with columns documented in `all_variables_v4_README.xlsx`.

## 1. Questions

1. How do size, stewardship, disturbance, site, heat, land use and social context affect street tree mortality, and how do those effects differ among species?
2. Where in the city are individual trees at the highest risk of dying?

**Decided:** Q1 (driver effects) is the main result. The Q2 risk maps are an application of the fitted models.

## 2. Response and sample

**Response.** The response is `mort = 1` if the tree had canopy loss at its stem point between the 2017 and 2021 lidar flights (`canopy_change == 3`), and `mort = 0` if there was no change (`canopy_change == 1`). Trees with canopy gain (`== 2`) are excluded. Mortality is the modeled event because it is the rare outcome (roughly 5–10% over 4 years). Modeling it directly lets the cloglog link give log hazard ratios (section 4).

**What "loss" measures.** Canopy loss at the stem can mean any of the following:

- death/removal of a living tree, for example for construction
- heavy pruning
- loss of a neighboring crown that overhung the point

We will not try to untangle natural causes of death from intentional removal of trees, because that is not possible with this data.

**Sample.** Filters are applied in this order. The number of trees removed at each step will be reported in an SI table, along with how the losses are distributed across boroughs and land-use types.

1. All 2015 census records: 683,788
2. Tree recorded as `Alive` in the 2015 census. Dead trees and stumps are currently removed only indirectly, because their species is missing; this should be an explicit filter.
3. `canopy_change` is 1 or 3. NA (79,970 records; 59,572 of them trees alive in 2015) means there was no canopy at the stem point in either 2017 or 2021 (confirmed by Ging-Yan). These trees are excluded. The sample is therefore trees with canopy at the stem in 2017. Trees that died or were removed between the 2015 census and the 2017 flight, and trees whose crowns don't register at the stem point, are outside it.
4. DBH greater than 3 in (7.62 cm) and less than 300 cm.
5. Land use from the nearest tax lot. There is no longer a distance filter. Trees more than 10 m (32.8 US survey ft) from any lot get their own land-use class, "not near a tax parcel" (about 16,000 trees, 3.3%). Trees within 10 m of a lot whose land use is missing are excluded: 2,027 trees (0.4%), mostly on miscellaneous (Z9), easement (Z7) or condo-billing (R0) lots. Their loss rate is the same as other trees (9.2% vs 9.1%), but they are spatially clustered, mostly in lower Manhattan, probably Battery Park City. **Provisional; may revisit.**
6. Complete covariates. Missing SVI theme 1 removes 1.4% of trees, and 55% of those are within 10 m of a park (vs 6% overall), so this filter does not remove a random subset. Note it in the SI.
7. Impervious share of visible ground is undefined when canopy covers the whole buffer. 184 trees with less than 1% visible ground even at 15 m are excluded (see note B).
8. Drop the duplicate-location filter (`group_by(x, y) %>% filter(n() == 1)`). It removes every tree at a shared coordinate, and the SPDE model handles points at the same location without problems. Report how many trees share coordinates.

**Taxa (decided 2026-10-07).** Taxa are modeled at species level where possible, with a genus level where the data allow and a single cutoff of more than 5,000 trees, with no exceptions:
1. **Species models:** species (binomial names) with more than 5,000 trees. That gives **17** species.
2. **Genus models:** trees not already in a species model, pooled by genus, where that remainder has more than 5,000 trees. Genus-only census labels go here. That gives **3**:
   - *Prunus* spp. (23,428 trees): the census genus-only "cherry" (16,432), *P. cerasifera*, *P. virginiana* and *P. serotina*.
   - *Acer* (other spp.) (14,643): unidentified *Acer*, sycamore, sugar, hedge, Amur maple and others.
   - *Quercus* (other spp.) (11,891): swamp white, willow, sawtooth, Shumard, white oak and others.
3. **"Other":** everything else (36,959 trees, 57 genera, including the non-American elms and the remaining ash).

That makes **21 models**; the fewest canopy losses in any model is 366 (*Liquidambar*). A 2,500 cutoff (26 models, with a Chinese elm model and a split *Prunus*) was considered and rejected to keep this broad paper simple. The genus/taxon multilevel meta-analysis (section 5) now has several genera with multiple taxa: *Acer* 4, *Quercus* 3, *Tilia* 3.

## 3. Drivers and how each is measured

The drivers are fixed **a priori** from the literature (Hilbert et al. 2019, plus heat and Sandy for NYC). Every taxon gets the same set of drivers, and a driver stays in the model even when its effect is near zero for a particular taxon. There is no variable selection per species.

The "direction" column gives the hypothesized effect on mortality. All candidate measures are in v4.

| Construct | Direction | Candidate measures | **Model term (decided)** |
|---|---|---|---|
| Size | U-shaped | `tree_dbh` | RW2 smooth on DBH (cm) |
| Stewardship | − | `steward` (recorded 2015) | **Binary:** any signs vs. none. Only 3% of trees have 3 or more signs (`output/SI_stewardship_levels.png`) |
| Building construction | + | `is_B_cons_{10,15,20,30}m_13_17`; `is_B_cons_10m_17_21` | **Both windows, 10 m,** as separate terms. See note A |
| Demolition | + | `is_DM_{10,15,20,30}m_13_17`; `is_DM_10m_17_21` | Both windows, 10 m. **Kept as a variable of interest** despite being rare (1.1% and 0.5% of trees). Every demolition tree is also flagged for building construction in the same window, so the coefficient is the extra effect of demolition beyond construction |
| Street construction | + | `is_S_cons_{5,10,15}m_13_17`; `is_S_cons_10m_17_21` | Both windows, 10 m. **The 2017–21 indicator is flipped (TEMPORARY):** as delivered it flags 88% of trees, which looks inverted. Flipped, it flags 11.7%. **[ASK G-Y]** Dan has emailed to confirm; remove the flip if a corrected file arrives |
| Surface cover around the stem (2017) | + | `{RD,OI,RR,BD,GS,SO,WA}_{10,15}m_2017` | **Impervious share of visible ground,** 10 m, with a 15 m fallback. See note B. No separate grass term |
| Neighborhood canopy (2017) | ? | `TC_{25,50,75,100}m_2017` | **100 m.** Not 5–15 m, which mostly measures the focal tree's own crown |
| Heat | + | `summer_mean/max/min`, `days_max_{27,32,35}` | **`summer_max` + `summer_min`.** See note C |
| Flooding legacy (Sandy 2012) | + | `in_sandy_zone` | As is |
| Land use (PLUTO) | varies | `LandUse`, 11 codes collapsed to 7 classes, plus "not near a tax parcel" (> 10 m) | **8 classes;** low-density residential is the reference, so **7 terms** |
| Public housing | + | `near_NYCHA_{5,10,15}m` | **Removed.** Only 0.64% of trees are flagged at 10 m, and 12 of 19 taxa have fewer than 10 canopy losses among flagged trees |
| Park adjacency | − | `near_park_{5,10,15}m` | **10 m.** Trees 10–15 m away are often across the street. A possible measurement-choice comparison across all taxa (5/10/15 m) |
| Human activity | + | `pop_{50,100,150,200}m` | **100 m,** log(1 + x) |
| Socioeconomic vulnerability | + | `RPL_THEME1–4`, `RPL_THEMES`, ACS % variables | **`RPL_THEME1`** (socioeconomic status). Median income is a possible SI sensitivity check only (missing for 9.8% of trees) |

This gives **23 fixed-effect terms** plus the DBH smooth and the spatial field. The smallest taxon (*Ulmus americana*, about 480 canopy losses) has about 20 losses per term.

**Excluded as predictors, because they are partly the outcome:** all `*_2021` land cover, all `*_diff` change variables, and `TC_{5,10,15}m_2017`. Including them would let the outcome leak into the predictors.

**Note A, construction timing.** The 2017–21 construction indicators are exactly the v2 variables the script uses now. They overlap the window in which mortality is measured, so they capture both construction damage and removal of trees to make way for construction. The 2013–17 indicators capture lagged damage only. **Decided:** both windows go in the main model. The 2013–17 coefficient is interpreted as lagged construction damage. The 2017–21 coefficient is interpreted as damage or intentional removal, since the two can't be separated.

**Note B, surface cover.** The eight land-cover shares sum to 100% for each tree. The land cover is classified from above, so surfaces under a crown are classified as canopy. As a result, impervious cover as a share of the whole buffer partly measures crown size (r = −0.64 with DBH). **Decided:** a single term, `imp_vis_10_2017` = 100 × impervious ÷ (impervious + grass/shrub + soil + water), where impervious = building + road + other impervious + railroad. This is the paved share of the visible ground (r = −0.08 with DBH).
- It uses the 10 m buffer when at least 1% of it is visible, and falls back to 15 m otherwise (r = 0.88 between the 10 m and 15 m shares). The 184 trees still under 1% visible at 15 m are excluded.
- These trees are large and have low loss (2.4%), so dropping all trees without visible ground at 10 m would have biased the sample.
- Grass/shrub is not a separate term. Soil and water are about 0, so grass is essentially the complement of the impervious share (r = −0.78 with the current grass term, −1.00 for grass's share of visible ground).
- The variable is bunched near the top: the median is 95%, and at least a quarter of trees have 100%.

**Note C, heat.** Heat varies along two separate dimensions:
- Daytime: `summer_max` and the heat-day counts are nearly interchangeable (r = 0.90–0.98).
- Nighttime: `summer_min` is almost unrelated to the daytime measures (r = 0.11 with `summer_max`).

`summer_mean` blends the two. **Decided:** `summer_max` and `summer_min` as separate terms; the mean and the heat-day counts are not used. Temperature varies little across the city (SD under 0.5°C) and comes on a grid of roughly 100 m.

### Redundancy, resolved from predictor correlations only (never from the outcome)

Rule: drop or combine a term if pooled |r| > 0.7 or VIF > 5, choosing which on scientific grounds. The check was run on 2026-10-06 (`output/SI_predictor_correlations.csv`, `output/SI_predictor_vif.csv`), before the impervious redefinition and the removal of grass and NYCHA.
- **VIFs:** the highest pooled VIF was 2.4. The highest in any single taxon was 4.96 (park adjacency in *Ulmus americana*).
- **Correlated pairs:** none exceeded 0.7 pooled.
  - DBH vs. impervious: −0.66. Resolved by the visible-ground redefinition, now −0.08.
  - Near park vs. open space/recreation land use: 0.60.
  - Summer min vs. population: 0.57.
- **Pairs over 0.7 in a single taxon:** near park vs. not-near-a-parcel (0.75), and near park vs. open space/recreation (0.71). Both are kept.
- Grass vs. the new impervious share: −0.78. Grass was dropped.
- **Spatial confounding:** heat, SVI and the Sandy zone vary smoothly across the city, so they compete with the spatial field. Expect wider intervals for these three. The Discussion should say this explicitly rather than reading a weak effect as no effect.
- **Rerun the check** on the final term set before the final fits.

### Measurement choices (buffer sizes, which heat or SVI metric)

These are made **once for all taxa**. The defaults in the table are fixed a priori. The SI then shows that the conclusions hold under the alternatives: 5 m and 15 m cover, the other construction buffers, the 5 m and 15 m park buffers, the other neighborhood-canopy and population scales, and `RPL_THEMES` or median income. If a data-driven choice is wanted instead, use WAIC or the log score summed over all taxa, never chosen per taxon.

## 4. Model, fit separately for each taxon with the same structure

- **Likelihood:** Bernoulli, with a cloglog link on `mort`, and offset log(Δt) where Δt is the years between the lidar flights (about 4). With this setup, exp(β) is a hazard ratio and the intercept is an annual hazard.
- **Fixed effects:** the drivers in section 3.
  - **Decided:** continuous predictors are z-scored using the **pooled mean and SD across all taxa**, so a "1 SD" effect is the same raw change for every species. This replaces the current scaling within each species.
  - Binary predictors are left as 0/1.
  - Prior: N(0, 1) on each coefficient, in place of INLA's nearly flat default.
- **DBH:** RW2 on an explicit `fm_mesh_1d` over DBH, with a PC prior on its precision.
- **Spatial field:** Matérn SPDE with PC priors.
  - Work in **meters**: project to EPSG:32118 (NAD83 / New York Long Island, m) instead of EPSG:2263 (feet).
  - **Decided:** range prior P(range < 500 m) = 0.05. The covariates capture most effects at the scale of a block face, so the field should mainly capture neighborhood-scale processes. The SI shows sensitivity at 250 m and 1 km.
  - Proposed SD prior: P(σ > 1) = 0.05 on the cloglog scale, which allows roughly a 7-fold difference in hazard between neighborhoods. The current prior is P(σ > 1.5) = 0.05.
  - The mesh inner `max.edge` is about 100 m (about range/5), with `cutoff` and the outer extension also in meters. `trees_sf` gets an explicit CRS. **Correction:** a 100 m mesh over NYC's 783 km² has about 217,000 vertices, not the ~90k estimated earlier. A test fit on that mesh failed inside INLA (not out of memory; cause not yet diagnosed). **Current setting: `max.edge` = 200 m** (about range/2.5), to be revisited later, for example a finer mesh or INLA settings that make 100 m work.
- **Software:** inlabru/INLA, with WAIC, CPO and `config` computed.

## 5. Synthesis across species

- **Primary:** random-effects meta-analysis of each driver's coefficient across taxa, using `metafor::rma(yi, sei, method = "REML", test = "knha")`.
  - Inputs: each taxon's posterior mean log hazard ratio (`yi`) and its posterior SD (`sei`), treating the posteriors as approximately normal.
  - REML estimates the between-taxon variance. The Knapp–Hartung adjustment gives honest CIs with only ~19 taxa.
  - "Other" is **included** for now (may revisit).
  - Reported for each driver: the mean (hazard ratio and CI), τ and I², and a 95% prediction interval (in `output/Fig4_coefficients.csv`). The framing emphasizes the individual species, which the audience knows well; the mean summarizes across them, and generalizing to other species is secondary. Prediction intervals are therefore not drawn on the figures.
- **SI sensitivity check:** a multilevel model with genus and taxon random effects (`rma.mv(..., random = ~ 1 | genus/sp)`), because related taxa may respond alike (three *Acer*, three *Tilia*, two *Quercus*). Output: `output/SI_meta_genus_sensitivity.csv`.
- **Possible extension:** meta-regression on species traits (e.g., salt or drought tolerance) to explain differences between species.
- **Possible upgrade:** one hierarchical model with species-specific slopes partially pooled toward a common mean. This depends on computing cost and on whether taxa can share a spatial field. Not the starting point.

## 6. Validation

- **Spatial block cross-validation:** about 2 km blocks and 5–10 folds. Report pooled log score, AUC and calibration for each taxon. This is the out-of-sample evidence behind the risk maps.
- **Diagnostics in the fitted models:** LOO-PIT from CPO, binned-residual plots and calibration plots.
- In-sample AUC is reported only alongside cross-validated AUC, never on its own.
- Hosmer–Lemeshow is dropped; with n ≈ 40k per taxon it rejects almost any model.

## 7. Checking that canopy loss means mortality

- Draw a random sample of loss trees, about 200 stratified by taxon and by construction status.
- Check each against Parks Forestry work orders and removal records, and against Street View imagery from before and after.
- Report the share confirmed dead or removed, and how much of the construction signal is removal.

## 8. Outputs

| Item | Content |
|---|---|
| Table 1 | Taxa, n, observed 4-year and annual mortality |
| Fig 2 | Mortality versus DBH, by taxon |
| Fig 3 | Map of observed mortality |
| Fig 4 | Forest plot by driver and taxon, with the meta-analytic mean |
| Fig 5 | Maps of risk for individual trees, from cross-validated or full-model predictions |
| SI | Sample flow; predictor correlations; sensitivity to measurement choices; spatial field maps; cross-validation and calibration; the loss-check results |

## 9. Code changes needed before refitting

- [x] Read v4 `.rds` and drop the socioeconomic join; update renamed columns. Change to logical types and text `LandUse` codes.
- [x] Model terms as in section 3: land use (7 terms plus not-near-parcel), heat pair, impervious share of visible ground, 100 m canopy and population, both construction windows, park 10 m; NYCHA and grass removed.
- [x] Unidentified *Acer* go into "other" (the script already did this; only the comment said "remove").
- [x] Drop `RPL_THEME3` from the kept columns (it isn't a model term).
- [ ] Remove the TEMPORARY street-construction flip if Ging-Yan sends a corrected file.
- [x] Add the explicit `status == "Alive"` filter, and build the sample-flow table (`output/SI_sample_flow.csv`).
- [x] Recode the response to `mort` and add the exposure offset.
- [x] Convert coordinates to meters with an explicit CRS. Re-specify the mesh and priors in meters. One mesh is shared by all taxa, built on a buffered, simplified borough boundary.
- [x] Remove the duplicate-location filter.
- [x] Use pooled z-scaling, N(0,1) fixed-effect priors, and `fm_mesh_1d` for DBH (20 knots over each taxon's DBH range, in cm).
- [x] Join on `tree_id`, not on all shared columns.
- [x] Field map: shared mesh, masked to NYC, axes in meters.
- [x] Fix `cor()` on non-numeric columns. The missing `%>%` was in the leftover code, which is now removed.
- [x] Remove the leftover GLM code after Fig 5: `fitted_preds`, `future_preds` and the residual maps.
- [x] Table 1 and Figs 2–5 report mortality rather than survival; the Fig 2 variables are renamed.
- [ ] Choose the mesh resolution (see section 4) and run all taxa.
- [ ] Spatial block CV, meta-analysis, sensitivity runs.

## Decision log

| Date | Decision |
|---|---|
| 2026-10-06 | v4 `.rds` is the canonical input data |
| 2026-10-06 | No per-species variable selection; the same a priori structure for every taxon |
| 2026-10-06 | Main result is driver effects (Q1); risk maps are an application |
| 2026-10-06 | Construction: both 2017–21 (damage or removal) and 2013–17 (lagged damage) |
| 2026-10-06 | ~~Land cover: road, building, other impervious and grass/shrub as separate terms~~ (superseded below) |
| 2026-10-06 | Stewardship: binary, any signs vs. none |
| 2026-10-06 | Predictors standardized across all taxa, not within each taxon |
| 2026-10-06 | Spatial range prior P(range < 500 m) = 0.05; sensitivity at 250 m and 1 km |
| 2026-10-06 | `canopy_change` NA = no canopy at the stem in either year; excluded (confirmed by Ging-Yan) |
| 2026-10-06 | No PLUTO distance filter; trees > 10 m from a lot form a "not near a tax parcel" land-use class |
| 2026-10-06 | Land use: low-density residential is the reference (7 terms) |
| 2026-10-06 | Trees within 10 m of a lot with missing land use excluded (2,027); provisional |
| 2026-10-06 | NYCHA proximity removed (too few flagged trees per taxon) |
| 2026-10-06 | Heat: `summer_max` + `summer_min`; mean and heat-day counts not used |
| 2026-10-06 | Park adjacency at 10 m |
| 2026-10-06 | Population and neighborhood canopy at 100 m (population as log(1 + x)) |
| 2026-10-06 | Construction and demolition, both windows, at 10 m; demolition kept as a variable of interest |
| 2026-10-06 | Land cover: single term, impervious share of visible ground (10 m, 15 m fallback); no separate grass term |
| 2026-10-06 | 2017–21 street construction flipped, pending Ging-Yan's confirmation (TEMPORARY) |
| 2026-10-06 | Mesh `max.edge` 200 m for now (100 m mesh failed in INLA); revisit later |
| 2026-10-07 | Mesh failures traced to sliver triangles from the coastline; boundary now buffered 500 m and simplified 150 m, `cutoff` 50 m |
| 2026-10-07 | INLA run with 8 threads (same estimates as 16, less memory) |
| 2026-10-07 | Meta-analysis: `metafor` REML + Knapp–Hartung, "other" included; genus/taxon multilevel model as SI; prediction intervals reported in the table but not on the figures |
| 2026-10-07 | Taxa: species models (> 5,000 trees) plus genus models for the remaining trees of a genus (> 5,000); 17 species + 3 genus models + "other" = 21; no exceptions (Chinese elm not separated) |
