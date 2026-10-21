# Street View verification of recorded tree losses

## Purpose

The NYC street tree mortality study detects tree loss from lidar. A tree counts as "lost" when the tree canopy at its stem location in 2017 is gone by 2021. That signal can sometimes be wrong:
- the tree may only have been heavily pruned
- a neighboring tree's crown that used to overhang the spot may be the thing that disappeared

This check estimates **how often a recorded loss really is a tree that is gone or dead**. We do **not** try to work out *why* a tree was lost (natural death vs. removal); only whether it is gone.

## What is in this folder

| File | What it is |
|---|---|
| `verification_sheet.xlsx` | **The file to fill in.** Sheet *Trees*: one row per tree. Sheet *Codes*: the allowed outcome codes and their meanings. |
| `README.md` | These instructions |
| `make_verification_sample.R` | Code that drew the sample (for the record; no need to run it) |
| `analyze_verification.R` | Code that summarizes the completed sheet (run by the analysis team) |

## The sample

- **200 primary trees** (`P001`–`P200`): a random sample of trees recorded as lost, in proportion to each species' share of all losses.
- **50 backup trees** (`B001`–`B050`): use these **only** to replace a primary tree you mark `NOT_VERIFIABLE`, and use them **in order**. For example, if 6 primaries are not verifiable, review `B001`–`B006`, plus more if some backups are not verifiable either.

The rows are already in random order. Please work through them **in the order given** and **do not sort, delete or rename columns**.

## How to review one tree (about 3–5 minutes)

1. **Open the location.** Click **Map with pin** in the row. Google Maps opens in **satellite view**, looking straight down at the block, with a red pin on the tree's recorded location. Use the satellite view to see which tree crown or pit is at the pin and how it sits relative to its neighbors. Then drag the yellow Street View person onto the street right beside the pin. The pinned location also shows inside Street View, which tells you which tree pit to look at.
2. **Find the tree.** The recorded location can be off by a few meters. Use these to pick the right street tree (usually in a sidewalk tree pit):
   - the species and trunk diameter (`dbh_inches_2015`, measured in 2015)
   - **`neighbors_2015`**: the nearest 2015 census tree on each side of the focal tree, with compass direction, distance, species (or "stump" / "dead tree") and diameter. For example, "NE 9 m: Callery pear, 11 in; SW 17 m: Callery pear, 7 in" means the focal tree sits between those two.

   If two or more trees could still be the one, use your best judgment and set **confidence = Low**, or mark `NOT_VERIFIABLE` if you really cannot tell.
3. **Look at the "before" imagery.** On a desktop browser, open the imagery dates (click the date or clock icon in the top-left Street View box, or "See more dates"). Choose an image from **2015–2017**, ideally in leaf-on season (May–October) and as close to 2017 as possible. Record its month and year in `before_image_date` (e.g., `2016-08`). Record whether the tree is there in `tree_present_before` (Yes / No / Unclear).
4. **Look at the "after" imagery.** Choose an image from **2021 or later**, as close to 2021 as possible and ideally leaf-on. Record its month and year in `after_image_date`.
5. **Choose the outcome** from the drop-down list:

| Code | Use when |
|---|---|
| `GONE` | The tree was there before and is not there after: an empty pit, a stump, or the pit paved over |
| `GONE_REPLANTED` | The tree was there before, and after there is a clearly different, young tree in its place |
| `DEAD_STANDING` | The tree is still standing after, but dead (bare in a summer image, or clearly dead) |
| `ALIVE_MAJOR_CROWN_LOSS` | The same tree is alive after, but a large part of its crown is gone (heavy pruning or storm damage) |
| `ALIVE_NO_CLEAR_CHANGE` | The same tree is alive after, with a similar crown |
| `NO_TREE_BEFORE` | There is no tree at this spot even in the before imagery |
| `NOT_VERIFIABLE` | No usable imagery in one of the windows, the tree is hidden (trucks, scaffolding), or you cannot tell which tree it is |

6. **Set `confidence`** to High or Low, and add anything useful in `notes` (e.g., "construction fence in 2021", "neighboring tree removed instead"). Put your initials in `reviewer`.

**Tips**
- Winter images can make a living deciduous tree look dead. Prefer leaf-on images, and use `DEAD_STANDING` only when you are fairly sure.
- If the only after-image is from much later (e.g., 2024), it is still useful. Note the date; the analysis team can judge.
- A newly planted young tree in the same pit means the original tree is gone (`GONE_REPLANTED`).

## When you are done

Save the workbook (keep the name `verification_sheet.xlsx` and the *.xlsx* format) and send it back. If more than one person reviews, please split the rows rather than reviewing the same trees. A second, independent review of about 20 trees is helpful for checking agreement, but optional.

## For the analysis team

Put the completed `verification_sheet.xlsx` back in this folder and run, from the project root:

```
Rscript streetview_verification/analyze_verification.R
```

It reports:
- the share of recorded losses confirmed as gone, replaced or dead, with an exact 95% CI
- the share with real canopy loss, which also counts `ALIVE_MAJOR_CROWN_LOSS`
- breakdowns by nearby construction and by borough

It also writes `verification_results.csv`. See `ANALYSIS_PLAN.md`, section 7.
