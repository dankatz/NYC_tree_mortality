# Package repositories: Posit Package Manager (CRAN snapshot) plus the INLA repository.
# Inside the dev container (Ubuntu 24.04) use P3M's prebuilt Linux binaries; elsewhere
# (e.g., RStudio on Windows or macOS) the standard P3M URL serves the right binaries.
in_container <- file.exists("/usr/local/lib/R/site-library") && Sys.info()[["sysname"]] == "Linux"
options(repos = c(
  CRAN = if (in_container) "https://p3m.dev/cran/__linux__/noble/latest" else "https://p3m.dev/cran/latest",
  INLA = "https://inla.r-inla-download.org/R/stable"
))

# In the container, keep editor tooling (languageserver, httpgd, vscDebugger) from the
# site library visible; project packages in renv/library still take precedence.
if (in_container) Sys.setenv(RENV_CONFIG_EXTERNAL_LIBRARIES = "/usr/local/lib/R/site-library")
rm(in_container)

source("renv/activate.R")

# R reads only one .Rprofile, so also run the user profile (VS Code httpgd hook).
if (file.exists("~/.Rprofile")) source("~/.Rprofile")
