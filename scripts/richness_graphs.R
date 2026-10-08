# Code for richness graphs displayed in Fuller et al. 2026
# Finalized by Sam Schurkamp on July 1, 2026

# ==============================================================================================
# load preamble (required functions & data sets) 
# ==============================================================================================
source(here::here("scripts", "preamble.R"), echo = TRUE)

# ==============================================================================================
# lollipop graph: total detections by method per site
# ==============================================================================================
# create a pointYear column, which fixes the lack of a facet wrap grouping the right points
all_data_plot <- all_data_plot %>%
  mutate(pointYear = paste0(Point, Year))

# manually mutate the data to allow for a better jitter
all_data_plot <- all_data_plot %>%
  dplyr::mutate(
    Protocol = factor(Protocol, levels = c("PC", "ARU", "ARU + PC")),
    protocol_order = as.integer(Protocol)
  ) %>%
  dplyr::group_by(Point, pointYear, Richness) %>%   # replace pointID if needed
  dplyr::mutate(
    n_tied = dplyr::n(),
    tie_rank = rank(protocol_order, ties.method = "first"),
    jitter_offset = dplyr::if_else(
      n_tied > 1,
      (tie_rank - (n_tied + 1) / 2) * 0.1,
      0
    ),
    richness_adjusted = Richness + jitter_offset
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n_tied, -tie_rank, -jitter_offset, -protocol_order)

# final data tweaks before graphing

# create dataset for segment lines
segment_lines <- all_data_plot %>%
  dplyr::arrange(Route, Year, Point, richness_adjusted) %>%
  dplyr::group_by(Route, Year, Point) %>%
  dplyr::mutate(next_richness = lead(richness_adjusted)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(next_richness))

# order the pointYear column so the graph is easier to interpret
pointYear_order <- all_data_plot %>%
  dplyr::group_by(.data$pointYear) %>%
  dplyr::summarise(
    sort_value = min(.data$richness_adjusted, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$sort_value)) %>%
  dplyr::pull("pointYear")

all_data_plot <- all_data_plot %>%
  dplyr::mutate(pointYear = factor(pointYear, levels = pointYear_order))

segment_lines <- segment_lines %>%
  dplyr::mutate(pointYear = factor(pointYear, levels = pointYear_order))

# define graph aesthetics ==========
Protocol_colors <- c(
  "ARU" = "#F8766D",
  "PC"  = "#00BFC4",
  "ARU + PC" = "#6F3959"
)

Protocol_shapes <- c(
  "ARU" = 18,
  "PC"  = 18,
  "ARU + PC" = 16 # diamond
)


# create the plot ====
long_richness_plot <- ggplot(
  all_data_plot,
  aes(x = richness_adjusted, y = pointYear, # x is adjusted for jitter, y is adjusted for discrete point-years
      color = Protocol, shape = Protocol)
) +
  # add the lines
  geom_segment( 
    data = segment_lines,
    aes(x = richness_adjusted, xend = next_richness,
        y = pointYear, yend = pointYear),
    inherit.aes = FALSE,
    color = "gray60",
    linewidth = 0.7
  ) +
  # add the points
  geom_point(size = 5) +
  scale_color_manual(
    values = Protocol_colors,
    name = "Detection Method",
    breaks = c("PC", "ARU", "ARU + PC"),
    labels = c("Point Count only", "ARU only", "ARU + Point Count Total")
  ) +
  scale_shape_manual(
    values = Protocol_shapes,
    name = "Detection Method",
    breaks = c("PC", "ARU", "ARU + PC"),
    labels = c("Point Count only", "ARU only", "ARU + Point Count Total")
  ) +
  scale_x_continuous(breaks = scales::breaks_width(2)) + # each vertical line is a whole number (without, the lines are every 1.5)
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_blank(), # y axis is blank since point-years are not important to identify
    legend.position = "bottom" # way too much blank space when it's to the side
  ) +
  labs(y = "Unique Point-Year", x = "Species Detected")

# generate the plot ====
long_richness_plot

# ==============================================================================================
# lollipop graph: species detections
# ==============================================================================================

# data preparation ====

# delete the misc species irrelevant to our study
pc_aru_species <- all_pc_aru_IMPORT %>%
  filter(!CommonName %in% c("----", "Trumpeter Swan", "Sandhill Crane", "Snowy Egret", "Forster's Tern", "Sedge Wren", "Common Tern", "Wilson's Snipe"))

# create a pointYear column (I like to have these be one item)
pc_aru_species <- pc_aru_species %>%
  mutate(pointYear = paste0(Point, Year))

# remove all instances of a non-observation
pc_aru_species_detections <- subset(pc_aru_species, IndivCount > 0)

# drop rows where the same species was observed in the same point-year under the same protocol
# (e.g., dropping all but one Sora per protocol per site per year)
speciesUnique <- pc_aru_species_detections %>%
  distinct(pointYear, Protocol, CommonName)

# pivot to wide:
# for each species at each site, specify which protocol found it
speciesUnique_wide <- speciesUnique %>%
  dplyr::mutate(detected = 1) %>%
  pivot_wider(names_from = Protocol, values_from = detected, values_fill = 0) %>%
  dplyr::mutate(Detection_Category = case_when(
    ARU == 1 & PC == 0 ~ "Only ARU",
    ARU == 0 & PC == 1 ~ "Only PC",
    ARU == 1 & PC == 1 ~ "ARU + PC",
    TRUE ~ "Neither"
  ))

# create list of species wholly missed by a detection category
# and the number of times a given detection category did find them
missed_species <- speciesUnique_wide %>%
  dplyr::group_by(CommonName, Detection_Category) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  filter(Detection_Category != "Neither")

# Define the order of categories for consistent positioning
category_levels <- c("Only ARU", "Only PC", "ARU + PC")
missed_species <- missed_species %>%
  dplyr::mutate(Detection_Category = factor(Detection_Category, levels = category_levels))

# Reorder species by total detections
species_order <- missed_species %>%
  dplyr::group_by(CommonName) %>%
  dplyr::summarise(total = sum(count)) %>%
  dplyr::arrange(total) %>%
  pull(CommonName)

# make sure missed_species is understood as factors 
missed_species <- missed_species %>%
  dplyr::mutate(CommonName = factor(CommonName, levels = species_order))

# graph components ==============================================================================================

# change x variable: percent of points with the detection
# (there are 60 point-years)
missed_species$count_per <- (missed_species$count /60)*100


# step 0: build the dataset for the line segments of the graph
segmentsSpecies <- missed_species %>%
  group_by(CommonName) %>%
  arrange(Detection_Category) %>%
  dplyr::reframe(
    x_start = count_per[-length(count_per)],
    x_end   = count_per[-1],
    y       = CommonName[1]
  )

# Step 1: Calculate max n per CommonName and get corresponding Detection_Category
max_n_per_common <- missed_species %>%
  dplyr::group_by(CommonName) %>%
  dplyr::mutate(max_n = max(count_per)) %>%
  filter(count_per == max_n) %>%
  slice(1) %>%  # in case of ties
  ungroup() %>%
  select(CommonName, Detection_Category_max = Detection_Category, max_n)

# Step 2: Identify CommonNames that appear only with "Only ARU"
only_aru_common <- missed_species %>%
  dplyr::group_by(CommonName) %>%
  dplyr::summarize(only_aru = all(Detection_Category == "Only ARU")) %>%
  ungroup()

# Step 3: Merge and assign custom order
ordering_info <- max_n_per_common %>%
  dplyr::left_join(only_aru_common, by = "CommonName") %>%
  dplyr::mutate(order_group = case_when(
    only_aru ~ 1,
    Detection_Category_max == "Only ARU" ~ 2,
    Detection_Category_max == "ARU + PC" ~ 3,
    TRUE ~ 4  # optional, in case of unexpected categories
  ))

# Step 4: Join back to original dataframe and arrange
species_miss_counts_ordered <- missed_species %>%
  left_join(ordering_info, by = "CommonName") %>%
  arrange(order_group, max_n, CommonName)

# re-order
row_id <- seq(1, nrow(species_miss_counts_ordered))
species_miss_counts_ordered$row_id = row_id
species_miss_counts_ordered <- species_miss_counts_ordered %>%
  arrange(desc(row_id))
species_order <- unique(species_miss_counts_ordered$CommonName)

# build the plot! ==============================================================================================

# Plot
reverse_detplot <- ggplot(missed_species,
                          aes(x = count_per, y = CommonName,
                              color = Detection_Category,
                              shape = Detection_Category)) +
  geom_segment(data = segmentsSpecies,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               inherit.aes = FALSE,
               color = "grey50", linewidth = 1) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    "Only ARU" = "coral1",
    "Only PC" = "dodgerblue",
    "ARU + PC" = "goldenrod1"),
    name = "Detected by:",
    breaks = c("Only ARU", "Only PC", "ARU + PC"),
    labels = c("Only ARU", "Only PC", "Both Methods")
  ) +
  scale_shape_manual(values = c(
    "Only ARU" = 16,
    "Only PC" = 16,
    "ARU + PC" = 16),
    name = "Detected by:",
    breaks = c("Only ARU", "Only PC", "ARU + PC"),
    labels = c("Only ARU", "Only PC", "Both Methods")
  ) +
  labs(
    #title = "Species Detection by Survey Protocol",
    x = "Percent of Point-Years With Detection",
    y = "Species",
    color = "Detected by:", 
    shape = "Detected by:"
  ) +
  scale_y_discrete(limits = species_order) +
  theme_minimal(base_size = 14)

# run the plot
reverse_detplot







