# NECESSARY FUNCTIONS AND DATA SETS OF FULLER ET AL. 2026
## This script is called at the start of the R Markdown files.
## It loads all items necessary.
## Finalized by Sam Schurkamp on July 1, 2026

# ==============================================================================================
# Step 0.1: Load (or load + install) all required packages 
# ==============================================================================================
# Vector of all packages required
packages <- c(
  ## organizational functions
  "tidyr", "dplyr", "stringr", "lubridate", "kableExtra", "data.table", "plyr","ggh4x","reshape",
  # used during the second lollipop graph prep
  "forcats",
  # improve the import functions for data within repositories
  "readr", "here",
  # graphing functions
  "ggplot2", "patchwork",
  # packages used in summary data sets
  "flextable", "officer", "purrr",
  # running statistical tests
  "car",
  # for the Bayesian summaries
  "AER", "vcd", "skellam",
  # for the hurdle model
  "glmmTMB"
  )

# Function to load packages
## (will also install and load any uninstalled packages)
install_and_load <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Run the function (load packages and install if necessary)
install_and_load(packages)

# ==============================================================================================
# Step 0.2: Load and modify richness data set
# ==============================================================================================

# read the appropriate CSV from GitHub into R
all_pc_aru_IMPORT <- readr::read_csv(here::here("data", "Fuller2026_cleaned_data.csv"))

# calculate richness for ARU and PC separately
pc_aru_richness <- all_pc_aru_IMPORT %>%
  filter(!CommonName %in% c("----", "Trumpeter Swan", "Sandhill Crane", "Snowy Egret", "Forster's Tern", "Sedge Wren", "Common Tern", "Wilson's Snipe")) %>% # had to add the sedge wren and common tern
  dplyr::filter(IndivCount >=1) %>%
  dplyr::group_by(Protocol, Route, Year, Point) %>%
  dplyr::summarize(Richness = length(unique(BirdCd)))

# calculate richness when ARU + PC are added together
pc_aru_richness_added_together <- all_pc_aru_IMPORT %>%
  filter(!CommonName %in% c("----", "Trumpeter Swan", "Sandhill Crane", "Snowy Egret", "Forster's Tern", "Sedge Wren", "Common Tern", "Wilson's Snipe")) %>%
  dplyr::filter(IndivCount >=1) %>%
  dplyr::group_by(Route, Point, Year) %>% # exclude Protocol to get combined richness
  dplyr::summarize(Richness = length(unique(BirdCd))) %>%
  dplyr::mutate(Protocol = "ARU + PC")

# combine the resulting DFs into one DF that has 3 protocol (ARU alone, PC alone, ARU + PC combined)
all_richness_no0 <- rbind(pc_aru_richness, pc_aru_richness_added_together)

# zero-fill richness summary  
all_data_plot <- all_richness_no0 %>%
  dplyr::select(Point, Protocol, Richness, Route, Year) %>%
  pivot_wider(names_from = Protocol, values_from = Richness) %>%
  pivot_longer(cols = c(ARU, PC, `ARU + PC`), names_to = "Protocol", values_to = "Richness") %>%
  dplyr::mutate(Richness = ifelse(is.na(Richness), 0, Richness)) 

# create a wide-form version
all_data_long <- cast(all_data_plot, Point+Route+Year~Protocol, value = "Richness")

# ==============================================================================================
# Step 0.3: De-clutter R environment
# ==============================================================================================
remove(
  packages,
  install_and_load,
  pc_aru_richness_added_together,
  all_richness_no0,
  pc_aru_richness
  )
















