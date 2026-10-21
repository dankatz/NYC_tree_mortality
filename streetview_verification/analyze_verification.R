# Street View verification of canopy loss: summarize the completed review spreadsheet
# Reads streetview_verification/verification_sheet.xlsx after review and reports how often a recorded canopy loss
# corresponds to a tree that is actually gone or dead (precision of the loss classification), with exact 95% CIs.
# Backups (B###) are used in order, one for each primary tree (P###) marked NOT_VERIFIABLE.
# See streetview_verification/README.md.

library(dplyr)
library(readxl)
library(readr)

sheet <- read_excel("streetview_verification/verification_sheet.xlsx", sheet = "Trees")
stopifnot("some reviewed rows have no outcome" = all(!is.na(sheet$outcome[sheet$role == "primary"])))

primary <- sheet %>% filter(role == "primary")
n_replace <- sum(primary$outcome == "NOT_VERIFIABLE")
backup_used <- sheet %>% filter(role == "backup", !is.na(outcome), outcome != "NOT_VERIFIABLE") %>%
  arrange(sample_id) %>% slice_head(n = n_replace)
if (nrow(backup_used) < n_replace) warning("not enough verified backups to replace all NOT_VERIFIABLE primary trees")
verified <- bind_rows(primary %>% filter(outcome != "NOT_VERIFIABLE"), backup_used)

#loss confirmed = the tree is gone, replaced, or dead; ALIVE_MAJOR_CROWN_LOSS is reported separately (canopy really was lost)
verified <- verified %>% mutate(
  confirmed = outcome %in% c("GONE", "GONE_REPLANTED", "DEAD_STANDING"),
  canopy_lost = confirmed | outcome == "ALIVE_MAJOR_CROWN_LOSS")
share_ci <- function(x) { b <- binom.test(sum(x), length(x)); c(pct = 100 * mean(x), lower = 100 * b$conf.int[1], upper = 100 * b$conf.int[2]) }

cat("Trees verified:", nrow(verified), " (primary:", sum(verified$role == "primary"), ", backups used:", nrow(backup_used),
    ", primary NOT_VERIFIABLE:", n_replace, ")\n\n")
cat("Outcomes:\n"); print(verified %>% count(outcome) %>% mutate(pct = round(100 * n / sum(n), 1)))
cat("\nRecorded loss confirmed (tree gone, replaced, or dead), % with exact 95% CI:\n"); print(round(share_ci(verified$confirmed), 1))
cat("Canopy actually lost (confirmed, or alive with major crown loss), %:\n"); print(round(share_ci(verified$canopy_lost), 1))
cat("\nBy nearby building construction (does accuracy differ with exposure?):\n")
print(verified %>% group_by(nearby_building_construction) %>% summarise(n = n(), pct_confirmed = round(100 * mean(confirmed), 1)))
cat("\nBy borough:\n")
print(verified %>% group_by(borough) %>% summarise(n = n(), pct_confirmed = round(100 * mean(confirmed), 1)))
cat("\nLow-confidence calls:", sum(verified$confidence == "Low", na.rm = TRUE), "\n")

write_csv(verified, "streetview_verification/verification_results.csv")
