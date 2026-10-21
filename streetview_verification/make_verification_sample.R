# Street View verification of canopy loss: draw the sample and build the review spreadsheet
# Draws 200 trees recorded as lost (canopy loss at the stem, 2017-2021) from the analysis sample, stratified
# proportionally by model taxon, plus 50 randomly ordered backups (used only to replace trees that cannot be verified).
# Run from the project root after the 100 m model fits (output/fits_100m/*_summary.rds, which hold the analysis sample).
# See ANALYSIS_PLAN.md section 7 and streetview_verification/README.md.

library(dplyr)
library(purrr)
library(readr)
library(sf)
library(openxlsx)

n_primary <- 200
n_backup <- 50
set.seed(20261020)
out_dir <- "streetview_verification"

#analysis sample (one row per tree) from the saved model summaries
trees <- map(list.files("output/fits_100m", pattern = "_summary.rds$", full.names = TRUE), ~ readRDS(.x)$trees) %>%
  bind_rows() %>%
  select(tree_id, taxon = sp_a, census_species = species, genus, dbh_cm, mort, x, y,
         construction_2013_17 = is_B_cons_1317_bool, construction_2017_21 = is_B_cons_bool)
losses <- trees %>% filter(mort == 1)

#proportional allocation of the 200 primary trees across taxa (largest-remainder rounding so the total is exactly 200)
alloc <- losses %>% count(taxon) %>%
  mutate(exact = n_primary * n / sum(n), n_sample = floor(exact), remainder = exact - n_sample)
alloc$n_sample[order(-alloc$remainder)[seq_len(n_primary - sum(alloc$n_sample))]] <-
  alloc$n_sample[order(-alloc$remainder)[seq_len(n_primary - sum(alloc$n_sample))]] + 1
primary <- losses %>%
  inner_join(alloc %>% select(taxon, n_sample), by = "taxon") %>%
  group_by(taxon) %>% group_modify(~ slice_sample(.x, n = .x$n_sample[1])) %>% ungroup() %>%
  select(-n_sample) %>% slice_sample(prop = 1) %>% #shuffle so reviewers do not work taxon by taxon
  mutate(role = "primary")
backup <- losses %>% filter(!tree_id %in% primary$tree_id) %>% slice_sample(n = n_backup) %>% mutate(role = "backup")
sample_df <- bind_rows(primary, backup)

#coordinates (analysis CRS is EPSG:32118, meters) to latitude/longitude, and borough
pts <- st_as_sf(sample_df, coords = c("x", "y"), crs = 32118, remove = FALSE) %>% st_transform(4326)
boro <- st_read("data/nyc_boundary_polygon/nybb.shp", quiet = TRUE) %>% st_transform(4326) %>% select(borough = BoroName)
pts <- st_join(pts, boro, join = st_nearest_feature)
ll <- st_coordinates(pts)

common <- read_csv("data/tree_mortality_variables/common_name_lookup.csv", show_col_types = FALSE)

#neighbors on either side, so the reviewer can tell which tree is the focal one: all 2015 census records (alive, dead, stumps).
#first = nearest census tree within 25 m; second = nearest one in roughly the opposite direction (more than 90 degrees away)
census <- readRDS("data/tree_mortality_variables/all_variables_v4.rds") %>%
  select(tree_id, status, species, tree_dbh, geom) %>%
  tidyr::separate_wider_delim(geom, delim = ",", names = c("xf", "yf")) %>%
  mutate(x = parse_number(xf) * 0.3048006, y = parse_number(yf) * 0.3048006) #US survey feet (EPSG:2263) to meters (same scale as EPSG:32118)
compass <- function(dx, dy) c("N", "NE", "E", "SE", "S", "SW", "W", "NW")[(round((atan2(dx, dy) * 180 / pi) %% 360 / 45) %% 8) + 1]
describe <- function(nb) {
  what <- case_when(nb$status == "Stump" ~ "stump",
                    nb$status == "Dead" ~ paste0("dead tree, ", nb$tree_dbh, " in"),
                    TRUE ~ paste0(coalesce(common$common_name[match(nb$species, common$species)], nb$species), ", ", nb$tree_dbh, " in"))
  sprintf("%s %.0f m: %s", compass(nb$dx, nb$dy), nb$d, what)
}
neighbors_text <- function(fx, fy, focal_id) {
  near <- census %>% filter(abs(x - fx) < 25, abs(y - fy) < 25, tree_id != focal_id) %>%
    mutate(dx = x - fx, dy = y - fy, d = sqrt(dx^2 + dy^2)) %>% filter(d <= 25) %>% arrange(d)
  if (nrow(near) == 0) return("no census trees within 25 m")
  first <- near[1, ]
  angle_diff <- abs(((atan2(near$dy, near$dx) - atan2(first$dy, first$dx) + pi) %% (2 * pi)) - pi)
  other <- near[angle_diff > pi / 2, ]
  paste0(describe(first), "; ", if (nrow(other) > 0) describe(other[1, ]) else "no census tree on the other side within 25 m")
}

#focal tree coordinates in the census CRS converted to meters, to match the census neighbors above
pts$x <- NULL; pts$y <- NULL
focal_xy <- census %>% select(tree_id, x, y)
pts <- pts %>% left_join(focal_xy, by = "tree_id")

sheet <- pts %>% st_drop_geometry() %>%
  mutate(latitude = round(ll[, "Y"], 6), longitude = round(ll[, "X"], 6),
         sample_id = sprintf("%s%03d", ifelse(role == "primary", "P", "B"), ave(seq_along(role), role, FUN = seq_along)),
         dbh_inches_2015 = round(dbh_cm / 2.54),
         nearby_building_construction = case_when(construction_2013_17 == 1 & construction_2017_21 == 1 ~ "2013-17 and 2017-21",
                                                  construction_2013_17 == 1 ~ "2013-17", construction_2017_21 == 1 ~ "2017-21",
                                                  TRUE ~ "none"),
         #Google Maps link that drops a red pin on the tree's recorded location, in satellite view from ~80 m up
         #(!3m1!1e3 = satellite imagery; tested 2026-10-21)
         map_pin = sprintf("https://www.google.com/maps/place/%.6f,%.6f/@%.6f,%.6f,80m/data=!3m1!1e3", latitude, longitude, latitude, longitude)) %>%
  left_join(common, by = c("census_species" = "species")) %>%
  mutate(neighbors_2015 = mapply(neighbors_text, x, y, tree_id)) %>%
  transmute(sample_id, role, tree_id, census_species, common_name, model_taxon = taxon, dbh_inches_2015, borough,
            map_pin, neighbors_2015, latitude, longitude, nearby_building_construction,
            #to be filled in by the reviewer
            before_image_date = "", tree_present_before = "", after_image_date = "", outcome = "", confidence = "",
            notes = "", reviewer = "")

#workbook: Trees (to fill in), Codes (allowed values and definitions), About
outcome_codes <- tibble::tribble(
  ~code,                     ~meaning,
  "GONE",                    "Tree present before; no tree after (empty pit, stump, paved over)",
  "GONE_REPLANTED",          "Tree present before; replaced by a clearly different, young tree after",
  "DEAD_STANDING",           "Tree still standing after but dead (no leaves in a leaf-on image, or clearly dead)",
  "ALIVE_MAJOR_CROWN_LOSS",  "Same tree alive after, but a large part of the crown is gone (heavy pruning, storm damage)",
  "ALIVE_NO_CLEAR_CHANGE",   "Same tree alive after with a similar crown",
  "NO_TREE_BEFORE",          "No tree at this location in the before imagery",
  "NOT_VERIFIABLE",          "No usable imagery in either window, tree hidden, or cannot tell which tree it is")

wb <- createWorkbook()
addWorksheet(wb, "Trees"); addWorksheet(wb, "Codes"); addWorksheet(wb, "About")
writeData(wb, "Trees", sheet)
#clickable map links
links <- sheet$map_pin; class(links) <- "hyperlink"; names(links) <- rep("Map with pin", nrow(sheet))
writeData(wb, "Trees", x = links, startCol = which(names(sheet) == "map_pin"), startRow = 2, colNames = FALSE)
writeData(wb, "Codes", outcome_codes, startCol = 1)
writeData(wb, "Codes", data.frame(yes_no = c("Yes", "No", "Unclear")), startCol = 4)
writeData(wb, "Codes", data.frame(confidence = c("High", "Low")), startCol = 5)
writeData(wb, "About", data.frame(About = c(
  "Street View verification of recorded canopy loss, NYC street tree mortality study.",
  "See README.md in this folder for full instructions.",
  sprintf("Sample: %d primary trees (P###; proportional to taxon) and %d backups (B###; use in order only to replace NOT_VERIFIABLE primaries).", n_primary, n_backup),
  sprintf("Drawn %s with seed 20261020 from %d trees recorded as lost (canopy loss at the stem point between the 2017 and 2021 lidar flights).", Sys.Date(), nrow(losses)))))
n <- nrow(sheet) + 1
col_of <- function(name) which(names(sheet) == name)
dataValidation(wb, "Trees", cols = col_of("outcome"), rows = 2:n, type = "list", value = sprintf("'Codes'!$A$2:$A$%d", nrow(outcome_codes) + 1))
dataValidation(wb, "Trees", cols = col_of("tree_present_before"), rows = 2:n, type = "list", value = "'Codes'!$D$2:$D$4")
dataValidation(wb, "Trees", cols = col_of("confidence"), rows = 2:n, type = "list", value = "'Codes'!$E$2:$E$3")
#make the columns to fill in stand out, freeze the header and ID columns
fill_cols <- col_of("before_image_date"):col_of("reviewer")
addStyle(wb, "Trees", createStyle(fgFill = "#FFF8DC"), rows = 1:n, cols = fill_cols, gridExpand = TRUE)
addStyle(wb, "Trees", createStyle(textDecoration = "bold"), rows = 1, cols = seq_along(sheet), gridExpand = TRUE, stack = TRUE)
freezePane(wb, "Trees", firstRow = TRUE, firstCol = TRUE)
setColWidths(wb, "Trees", cols = seq_along(sheet), widths = "auto")
setColWidths(wb, "Codes", cols = 1:2, widths = c(26, 95))
setColWidths(wb, "About", cols = 1, widths = 140)
saveWorkbook(wb, file.path(out_dir, "verification_sheet.xlsx"), overwrite = TRUE)

cat("primary:", sum(sheet$role == "primary"), " backup:", sum(sheet$role == "backup"), " taxa in primary:", n_distinct(primary$taxon), "\n")
print(alloc %>% select(taxon, losses = n, n_sample) %>% arrange(desc(losses)), n = Inf)
