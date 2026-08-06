#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else getwd()
need <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "purrr")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), ". Install them in a non-base Conda environment.")

library(readr); library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(purrr)

species_order <- c(
  "Bacteroides thetaiotaomicron",
  "Cutibacterium acnes",
  "Fusobacterium nucleatum",
  "Prevotella melaninogenica",
  "Veillonella parvula"
)
species_keys <- setNames(str_replace_all(species_order, " ", "_"), species_order)

path <- function(...) file.path(project_root, ...)
dir.create(path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("results", "combined"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("reports"), recursive = TRUE, showWarnings = FALSE)

read_tsv_if <- function(p) if (file.exists(p) && file.info(p)$size > 0) readr::read_tsv(p, show_col_types = FALSE) else tibble()

custom_all <- read_tsv_if(path("results", "combined", "custom_emu_expected_species_abundance_all_usable_reads.tsv"))
custom_target <- read_tsv_if(path("results", "combined", "custom_emu_expected_species_abundance_renormalized.tsv"))
broad_species <- read_tsv_if(path("results", "combined", "broad_emu_species_abundance.tsv"))
broad_unexpected <- read_tsv_if(path("results", "combined", "broad_emu_unexpected_species.tsv"))
broad_genus <- read_tsv_if(path("results", "combined", "broad_emu_genus_abundance.tsv"))
custom_accounting <- read_tsv_if(path("results", "combined", "custom_emu_read_accounting.tsv"))
broad_accounting <- read_tsv_if(path("results", "combined", "broad_emu_read_accounting.tsv"))
minimap <- read_tsv_if(path("results", "combined", "minimap2_assignment_summary.tsv"))
fastq_manifest <- read_tsv_if(path("metadata", "fastq_manifest.tsv"))
qc_summary <- read_tsv_if(path("qc", "filtered", "filter_summary.tsv"))
sample_metadata <- read_tsv_if(path("metadata", "sample_metadata.tsv"))
expected <- read_tsv_if(path("metadata", "expected_ratios_template.tsv"))

custom_long <- custom_target %>%
  mutate(species = factor(species, levels = species_order))
if (nrow(sample_metadata) > 0 && "sample_id" %in% names(sample_metadata)) {
  custom_long <- custom_long %>% left_join(sample_metadata, by = "sample_id")
}
readr::write_tsv(custom_long, path("results", "combined", "custom_emu_tidy_long.tsv"))
readr::write_csv(custom_long, path("results", "combined", "custom_emu_tidy_long.csv"))

if (nrow(custom_long) > 0) {
  p <- ggplot(custom_long, aes(sample_id, emu_abundance_pct_renormalized_among_expected, fill = species)) +
    geom_col(width = 0.85) +
    scale_fill_brewer(palette = "Set2", drop = FALSE) +
    labs(x = NULL, y = "Renormalized abundance among five expected species (%)", fill = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(path("results", "figures", "custom_emu_stacked_abundance.png"), p, width = 9, height = 5, dpi = 300)
  ggsave(path("results", "figures", "custom_emu_stacked_abundance.pdf"), p, width = 9, height = 5)
}

expected_long <- tibble()
if (nrow(expected) > 0 && "sample_id" %in% names(expected)) {
  allowed_ratio_bases <- c(
    "purified_amplicon_molar", "purified_amplicon_mass", "gDNA_mass",
    "estimated_16S_copy", "cell_count", "CFU", "culture_volume", "unknown"
  )
  expected <- expected %>% filter(!is.na(sample_id), sample_id != "")
  if (anyDuplicated(expected$sample_id)) {
    stop("metadata/expected_ratios_template.tsv contains duplicate sample_id values")
  }
  invalid_basis <- setdiff(na.omit(unique(expected$ratio_basis)), allowed_ratio_bases)
  if (length(invalid_basis)) {
    stop("Invalid ratio_basis value(s): ", paste(invalid_basis, collapse = ", "))
  }
  ratio_columns <- names(expected)[str_detect(names(expected), "_expected_pct$")]
  ratio_checks <- expected %>%
    mutate(expected_sum = rowSums(across(all_of(ratio_columns)), na.rm = TRUE),
           supplied_values = rowSums(!is.na(across(all_of(ratio_columns)))))
  missing_basis <- ratio_checks %>%
    filter(supplied_values > 0, is.na(ratio_basis) | ratio_basis == "")
  if (nrow(missing_basis) > 0) {
    stop(
      "ratio_basis is required for sample(s) with expected percentages: ",
      paste(missing_basis$sample_id, collapse = ", ")
    )
  }
  ratio_sums <- ratio_checks %>%
    filter(supplied_values > 0, abs(expected_sum - 100) > 0.5)
  if (nrow(ratio_sums) > 0) {
    stop(
      "Expected percentages must sum to 100 (+/- 0.5) for sample(s): ",
      paste(ratio_sums$sample_id, collapse = ", ")
    )
  }
  expected_long <- expected %>%
    select(sample_id, ratio_basis, ends_with("_expected_pct"), known_absent_species) %>%
    pivot_longer(ends_with("_expected_pct"), names_to = "species_key", values_to = "expected_pct") %>%
    mutate(species = str_remove(species_key, "_expected_pct$"),
           species = str_replace_all(species, "_", " "),
           species = factor(species, levels = species_order))
  readr::write_tsv(expected_long, path("results", "combined", "expected_ratios_tidy_long.tsv"))
  readr::write_csv(expected_long, path("results", "combined", "expected_ratios_tidy_long.csv"))
}

if (nrow(expected_long) > 0 && nrow(custom_long) > 0) {
  comp <- custom_long %>%
    select(sample_id, species, observed_pct = emu_abundance_pct_renormalized_among_expected) %>%
    left_join(expected_long %>% select(sample_id, species, ratio_basis, expected_pct), by = c("sample_id", "species")) %>%
    filter(!is.na(expected_pct)) %>%
    mutate(error_pct_points = observed_pct - expected_pct,
           abs_error_pct_points = abs(error_pct_points),
           recovery_ratio = if_else(expected_pct > 0, observed_pct / expected_pct, NA_real_))
  metrics <- comp %>%
    summarise(
      mean_absolute_error_pct_points = mean(abs_error_pct_points, na.rm = TRUE),
      root_mean_squared_error = sqrt(mean(error_pct_points^2, na.rm = TRUE)),
      pearson_correlation = suppressWarnings(cor(expected_pct, observed_pct, method = "pearson", use = "complete.obs")),
      spearman_correlation = suppressWarnings(cor(expected_pct, observed_pct, method = "spearman", use = "complete.obs"))
    )
  species_bias <- comp %>%
    group_by(species) %>%
    summarise(systematic_bias_pct_points = mean(error_pct_points, na.rm = TRUE),
              mean_recovery_ratio = mean(recovery_ratio, na.rm = TRUE), .groups = "drop")
  readr::write_tsv(comp, path("results", "combined", "expected_vs_observed_custom_emu.tsv"))
  readr::write_tsv(metrics, path("results", "combined", "expected_vs_observed_metrics.tsv"))
  readr::write_tsv(species_bias, path("results", "combined", "organism_specific_bias.tsv"))
  readr::write_csv(comp, path("results", "combined", "expected_vs_observed_custom_emu.csv"))
  readr::write_csv(metrics, path("results", "combined", "expected_vs_observed_metrics.csv"))
  readr::write_csv(species_bias, path("results", "combined", "organism_specific_bias.csv"))
  p <- ggplot(comp, aes(expected_pct, observed_pct, color = species)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
    geom_point(size = 2) +
    facet_wrap(~ species) +
    scale_color_brewer(palette = "Set2", drop = FALSE) +
    labs(x = "Expected (%)", y = "Observed custom Emu (%)", color = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
  ggsave(path("results", "figures", "expected_vs_observed_custom_emu.png"), p, width = 8, height = 6, dpi = 300)
  ggsave(path("results", "figures", "expected_vs_observed_custom_emu.pdf"), p, width = 8, height = 6)
}

if (nrow(broad_unexpected) > 0) {
  readr::write_csv(broad_unexpected, path("results", "combined", "broad_emu_unexpected_species.csv"))
  unexpected_plot <- broad_unexpected %>%
    filter(abundance_pct_all_usable_reads > 0) %>%
    mutate(
      taxon_label = if_else(
        is.na(species) | species == "",
        paste0("Unclassified ", coalesce(genus, "taxon")),
        species
      )
    ) %>%
    group_by(sample_id) %>%
    slice_max(abundance_pct_all_usable_reads, n = 10, with_ties = FALSE) %>%
    ungroup()
  if (nrow(unexpected_plot) > 0) {
    p <- ggplot(
      unexpected_plot,
      aes(taxon_label, abundance_pct_all_usable_reads, fill = taxon_category)
    ) +
      geom_col() +
      coord_flip() +
      facet_wrap(~ sample_id, scales = "free_y", ncol = 3) +
      labs(
        x = NULL,
        y = "Abundance among all QC-passing reads (%)",
        fill = NULL
      ) +
      theme_bw(base_size = 9)
    ggsave(path("results", "figures", "broad_emu_unexpected_species.png"), p, width = 11, height = 13, dpi = 300)
    ggsave(path("results", "figures", "broad_emu_unexpected_species.pdf"), p, width = 11, height = 13)
  }
}

if (nrow(broad_genus) > 0 && nrow(custom_target) > 0) {
  custom_genus <- custom_target %>%
    mutate(genus = str_extract(species, "^[^ ]+")) %>%
    group_by(sample_id, genus) %>%
    summarise(
      custom_emu_pct = sum(emu_abundance_pct_renormalized_among_expected),
      .groups = "drop"
    )
  broad_expected_genus <- broad_genus %>%
    filter(genus %in% str_extract(species_order, "^[^ ]+")) %>%
    group_by(sample_id) %>%
    mutate(
      broad_expected_total = sum(estimated_counts),
      broad_emu_pct = if_else(
        broad_expected_total > 0,
        100 * estimated_counts / broad_expected_total,
        0
      )
    ) %>%
    ungroup() %>%
    select(sample_id, genus, broad_emu_pct)
  genus_comparison <- full_join(
    custom_genus,
    broad_expected_genus,
    by = c("sample_id", "genus")
  )
  readr::write_tsv(genus_comparison, path("results", "combined", "custom_emu_vs_broad_emu_genus.tsv"))
  readr::write_csv(genus_comparison, path("results", "combined", "custom_emu_vs_broad_emu_genus.csv"))
  p <- ggplot(genus_comparison, aes(custom_emu_pct, broad_emu_pct, color = genus)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
    geom_point(size = 2, na.rm = TRUE) +
    facet_wrap(~ genus) +
    scale_color_brewer(palette = "Set2", drop = FALSE) +
    labs(
      x = "Custom Emu: expected-species composition (%)",
      y = "Broad Emu: expected-genus composition (%)",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
  ggsave(path("results", "figures", "custom_emu_vs_broad_emu_genus.png"), p, width = 8, height = 6, dpi = 300)
  ggsave(path("results", "figures", "custom_emu_vs_broad_emu_genus.pdf"), p, width = 8, height = 6)
}

if (nrow(expected) > 0 && nrow(custom_all) > 0 && "known_absent_species" %in% names(expected)) {
  known_absent <- expected %>%
    select(sample_id, known_absent_species) %>%
    filter(!is.na(known_absent_species), known_absent_species != "") %>%
    separate_rows(known_absent_species, sep = "\\s*[;,]\\s*") %>%
    transmute(sample_id, species = known_absent_species)
  false_positive <- known_absent %>%
    left_join(
      custom_all %>%
        select(sample_id, species, observed_pct = emu_abundance_pct_all_usable_reads),
      by = c("sample_id", "species")
    ) %>%
    mutate(observed_pct = coalesce(observed_pct, 0))
  if (nrow(false_positive) > 0) {
    readr::write_tsv(false_positive, path("results", "combined", "known_absence_false_positive_background.tsv"))
    readr::write_csv(false_positive, path("results", "combined", "known_absence_false_positive_background.csv"))
    p <- ggplot(false_positive, aes(sample_id, observed_pct, fill = species)) +
      geom_col(position = "dodge") +
      labs(
        x = NULL,
        y = "Observed abundance for known-absent species (%)",
        fill = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(path("results", "figures", "known_absence_false_positive.png"), p, width = 9, height = 5, dpi = 300)
    ggsave(path("results", "figures", "known_absence_false_positive.pdf"), p, width = 9, height = 5)
  }
}

if (nrow(custom_accounting) > 0 && "total_unassigned_or_unclassified_pct" %in% names(custom_accounting)) {
  p <- ggplot(custom_accounting, aes(sample_id, total_unassigned_or_unclassified_pct)) +
    geom_col(fill = "#6B6B6B") +
    labs(x = NULL, y = "Custom Emu unassigned or unclassified (%)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(path("results", "figures", "custom_emu_unassigned_fraction.png"), p, width = 9, height = 4.5, dpi = 300)
  ggsave(path("results", "figures", "custom_emu_unassigned_fraction.pdf"), p, width = 9, height = 4.5)
}

if (nrow(minimap) > 0 && nrow(custom_long) > 0) {
  minimap_comp <- minimap %>%
    filter(assignment_status == "confident") %>%
    mutate(species = str_replace(reference, "_ref$", ""),
           species = str_replace_all(species, "_", " ")) %>%
    group_by(sample_id, species) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    group_by(sample_id) %>%
    mutate(minimap_pct = 100 * reads / sum(reads)) %>%
    ungroup() %>%
    full_join(custom_long %>% select(sample_id, species, custom_emu_pct = emu_abundance_pct_renormalized_among_expected), by = c("sample_id", "species"))
  readr::write_tsv(minimap_comp, path("results", "combined", "custom_emu_vs_minimap2.tsv"))
  p <- ggplot(minimap_comp, aes(custom_emu_pct, minimap_pct, color = species)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
    geom_point(size = 2, na.rm = TRUE) +
    facet_wrap(~ species) +
    scale_color_brewer(palette = "Set2", drop = FALSE) +
    labs(x = "Custom Emu (%)", y = "Minimap2 confident assignments (%)", color = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
  ggsave(path("results", "figures", "custom_emu_vs_minimap2.png"), p, width = 8, height = 6, dpi = 300)
}

if (nrow(qc_summary) > 0) {
  readr::write_tsv(qc_summary, path("results", "combined", "qc_filter_summary.tsv"))
  readr::write_csv(qc_summary, path("results", "combined", "qc_filter_summary.csv"))
  qc_included <- qc_summary
  if ("analysis_status" %in% names(qc_included)) {
    qc_included <- qc_included %>% filter(analysis_status == "included")
  }
  p <- ggplot(qc_included, aes(sample_id, percent_retained)) +
    geom_col(fill = "#3C78A8") +
    labs(x = NULL, y = "Reads retained after filtering (%)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(path("results", "figures", "qc_percent_retained.png"), p, width = 9, height = 4.5, dpi = 300)
  ggsave(path("results", "figures", "qc_percent_retained.pdf"), p, width = 9, height = 4.5)
  p_count <- ggplot(qc_included, aes(sample_id, reads_after_quality_filter)) +
    geom_col(fill = "#3C8D75") +
    labs(x = NULL, y = "QC-passing reads") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(path("results", "figures", "qc_retained_read_count.png"), p_count, width = 9, height = 4.5, dpi = 300)
  ggsave(path("results", "figures", "qc_retained_read_count.pdf"), p_count, width = 9, height = 4.5)
}

session <- capture.output(sessionInfo())
writeLines(session, path("reports", "R_sessionInfo.txt"))
message("Analysis tables and figures written under results/combined and results/figures")
