#!/usr/bin/env Rscript

# ==============================================================================
# MSA_plotter_with_N-sequon_markers.R
#
# Description:
#   This script aligns sequences from a FASTA file and generates a readable,
#   text-like Multiple Sequence Alignment (MSA) plot using 'ggplot2'.
#   The plot is wrapped, highlights N-sequons, and adds labeled markers for them.
#
# Usage:
#   Rscript create_msa_plot.R -i your_sequences.fasta -o output_plot.pdf
#
# Arguments:
#   -i, --input        (Required) Path to the input FASTA file (.fa, .fasta).
#   -o, --output       (Optional) Path for the output PDF file.
#                      Defaults to 'MSA_Plot.pdf'.
#   --no-align         (Optional) Flag to indicate that the input file is
#                      already aligned.
#   --wrap-width       (Optional) Number of characters per line before wrapping.
#                      Defaults to 100.
# ==============================================================================


# --- 1. Package Management ---
# Check for and install required packages if they are not already present.

install_if_missing <- function(pkg_name, repo = "CRAN") {
  if (!require(pkg_name, character.only = TRUE)) {
    message(paste("Installing package:", pkg_name))
    if (repo == "BiocManager") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg_name, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg_name)
    }
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
  }
}

# Install necessary packages
install_if_missing("argparse")
install_if_missing("msa", repo = "BiocManager")
install_if_missing("ggmsa", repo = "BiocManager") # Still needed for tidy_msa
install_if_missing("ggplot2")
install_if_missing("dplyr")       # For data manipulation
install_if_missing("patchwork")   # For combining plots
install_if_missing("tidyr")       # For unnest()
install_if_missing("purrr")       # For map()


# --- 2. Command-Line Argument Parsing ---

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = "Align sequences and generate a readable, wrapped MSA plot from a FASTA file.")

parser$add_argument("-i", "--input",
                    type = "character",
                    required = TRUE,
                    help = "Path to the input FASTA file (.fa, .fasta)")

parser$add_argument("-o", "--output",
                    type = "character",
                    default = "MSA_Plot.pdf",
                    help = "Path for the output PDF file [default: %(default)s]")

parser$add_argument("--no-align",
                    action = "store_true",
                    default = FALSE,
                    help = "Skip the alignment step (assumes input file is already aligned)")

parser$add_argument("--wrap-width",
                    type = "integer",
                    default = 100,
                    help = "Number of characters per line before wrapping [default: %(default)s]")

args <- parser$parse_args()


# --- 3. Main Script Logic ---

cat("Input file:", args$input, "\n")
cat("Output file:", args$output, "\n")

if (!file.exists(args$input)) {
  stop("Error: Input file not found at path: ", args$input)
}

msa_data <- NULL

if (args$no_align) {
  cat("Skipping alignment as per --no-align flag. Reading aligned file...\n")
  msa_data <- Biostrings::readAAMultipleAlignment(args$input)
} else {
  cat("Reading sequences from file...\n")
  mySequences <- tryCatch({
    Biostrings::readAAStringSet(args$input)
  }, error = function(e) {
    stop("Failed to read the FASTA file. Please ensure it is a valid amino acid FASTA file.\nError message: ", e$message)
  })

  cat("Performing multiple sequence alignment (this may take a moment)...\n")
  alignment_result <- msa::msa(mySequences)

  cat("Converting alignment object for plotting...\n")
  msa_data <- as(alignment_result, "AAMultipleAlignment")
}


# --- 4. Tidy data, find motifs, determine conservation, and sort ---

cat("Tidying MSA data for plotting...\n")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("ggmsa"))

# Convert the MSA object into a "tidy" data frame
tidy_data <- tidy_msa(msa_data)

# Set factor levels on the tidied data to ensure ggplot respects alphabetical order
sorted_names <- sort(unique(tidy_data$name))
tidy_data$name <- factor(tidy_data$name, levels = sorted_names)

# Find N-glycosylation sequons (N-X-S/T where X is not P)
sequon_pattern <- "N[^P][ST]"
sequon_positions <- tidy_data %>%
  group_by(name) %>%
  summarise(sequence_str = paste(character, collapse = ""), .groups = "drop") %>%
  # Create a list-column containing a data.frame with a predictable column name 'starts'
  mutate(matches = purrr::map(sequence_str, function(s) {
    match_starts <- gregexpr(sequon_pattern, s)[[1]]
    return(data.frame(starts = match_starts))
  })) %>%
  unnest(matches) %>% # This now safely creates a 'starts' column
  filter(starts > 0) %>%
  group_by(name, starts) %>%
  summarise(position = list(c(starts, starts + 1, starts + 2)), .groups = "drop") %>%
  unnest(position) %>%
  mutate(is_sequon = TRUE)

# **NEW**: Create a data frame for the unique marker positions
unique_marker_positions <- sequon_positions %>%
  distinct(starts) %>%
  rename(position = starts) %>%
  mutate(facet_group = floor((position - 1) / args$wrap_width))


# Join sequon info and determine conservation
tidy_data <- tidy_data %>%
  left_join(sequon_positions, by = c("name", "position")) %>%
  group_by(position) %>%
  mutate(is_conserved = n_distinct(character[character != "-"]) == 1) %>%
  ungroup() %>%
  mutate(
    font_face = ifelse(!is.na(is_sequon) & is_sequon, "bold", "plain"),
    facet_group = floor((position - 1) / args$wrap_width)
  )


# --- 5. Generate and Save the Readable MSA Plot ---
cat("Generating readable, wrapped MSA plot using ggplot2 & patchwork...\n")
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("patchwork"))

conservation_colors <- c("conserved" = "#a6cee3", "variable"  = "#fbb4ae", "gap" = "white")

num_seqs <- length(levels(tidy_data$name))

plot_list <- lapply(unique(tidy_data$facet_group), function(group) {
  
  subset_data <- tidy_data %>% filter(facet_group == group)
  start_pos <- group * args$wrap_width
  end_pos <- start_pos + args$wrap_width
  
  p <- ggplot(subset_data, aes(x = position, y = name)) +
    geom_tile(aes(fill = ifelse(character == "-", "gap", ifelse(is_conserved, "conserved", "variable"))),
              color = "gray85", height = 0.7, width = 1) +
    geom_text(aes(label = character, fontface = font_face), color = "black", size = 2.5) +
    scale_fill_manual(values = conservation_colors) +
    # **FIX**: Add 0.5 padding and turn clipping off to prevent edges being cut
    coord_cartesian(xlim = c(start_pos + 0.5, end_pos + 0.5), clip = "off") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.title = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      # Add margin to the top of the plot area to make space for markers
      plot.margin = margin(t = 25, r = 5, b = 5, l = 5, unit = "pt")
    )
    
  # Add the marker layer only for the top row of plots
  marker_data_for_facet <- unique_marker_positions %>% filter(facet_group == group)
  if(nrow(marker_data_for_facet) > 0) {
      p <- p +
        # Add the upside-down triangle marker
        geom_point(
          data = marker_data_for_facet,
          aes(x = position, y = num_seqs + 0.9), # Position above the top sequence
          shape = 6, # upside-down triangle
          size = 3,
          color = "black"
        ) +
        # **NEW**: Add the position number above the triangle
        geom_text(
          data = marker_data_for_facet,
          aes(x = position, y = num_seqs + 1.6, label = position), # Position above the triangle
          size = 2.5,
          color = "black",
          vjust = 0.5
        )
  }
  
  return(p)
})

# Combine all plot segments vertically
combined_plot <- wrap_plots(plot_list, ncol = 1)


# --- 6. Dynamically calculate PDF dimensions and save ---
num_pos <- max(tidy_data$position)
num_facets <- ceiling(num_pos / args$wrap_width)

plot_height <- (num_seqs * num_facets * 0.25) + 2
if (plot_height > 50) {
  plot_height <- 50
  cat("Warning: Plot height capped at 50 inches to prevent excessively large PDF.\n")
}

ggsave(
  args$output,
  plot = combined_plot,
  width = 11,
  height = plot_height,
  units = "in",
  dpi = 300,
  device = "pdf",
  limitsize = FALSE
)

cat("Success! Plot saved to:", args$output, "\n")
