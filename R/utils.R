from_here <- function(...) {
  # Resolve paths relative to repository root (where config/ exists).
  # Works even if scripts are run from another working directory.
  script_dir <- normalizePath(dirname(sys.frames()[[1]]$ofile %||% "."), winslash = "/", mustWork = FALSE)
  # Go up one level from R/ to repo root
  root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
  normalizePath(file.path(root, ...), winslash = "/", mustWork = FALSE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

log_info <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] ", ts), sprintf(...), "\n", sep = "")
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stopf("Missing R package '%s'. Please install it (see scripts_r/00_install_packages.R).", pkg)
  }
}

read_yaml <- function(path) {
  require_pkg("yaml")
  if (!file.exists(path)) stopf("Config not found: %s", path)
  yaml::read_yaml(path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_gene_list <- function(path) {
  if (!file.exists(path)) return(NULL)
  genes <- readLines(path, warn = FALSE)
  genes <- unique(trimws(genes))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) return(NULL)
  genes
}

safe_write_csv <- function(df, path) {
  ensure_dir(dirname(path))
  utils::write.csv(df, path, row.names = FALSE, quote = TRUE)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(config = "config/config.yaml")
  if (length(args) == 0) return(out)
  for (i in seq_along(args)) {
    if (args[i] == "--config" && i < length(args)) out$config <- args[i + 1]
  }
  out
}

load_cfg <- function(config_path) {
  cfg <- read_yaml(from_here("..", config_path))
  # materialize key paths
  cfg$paths$data_dir <- from_here("..", cfg$paths$data_dir)
  cfg$paths$results_dir <- from_here("..", cfg$paths$results_dir)
  cfg$paths$cache_dir <- from_here("..", cfg$paths$cache_dir)

  # output dirs
  cfg$outputs$tcga_dir <- from_here("..", cfg$outputs$tcga_dir)
  cfg$outputs$deg_dir <- from_here("..", cfg$outputs$deg_dir)
  cfg$outputs$enrichment_dir <- from_here("..", cfg$outputs$enrichment_dir)
  cfg$outputs$clustering_dir <- from_here("..", cfg$outputs$clustering_dir)
  cfg$outputs$immune_dir <- from_here("..", cfg$outputs$immune_dir)

  cfg$inputs$dna_metabolism_genes <- from_here("..", cfg$inputs$dna_metabolism_genes)
  cfg$inputs$immune_gmt <- from_here("..", cfg$inputs$immune_gmt)

  cfg
}
