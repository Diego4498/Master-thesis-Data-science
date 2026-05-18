# scripts/00_setup_paths.R
suppressPackageStartupMessages({
  if (requireNamespace("here", quietly = TRUE)) library(here)
})

detect_project_root <- function() {
  env_root <- Sys.getenv("PROJECT_ROOT", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }
  
  if (exists("here") && is.function(here::here)) {
    root <- tryCatch(here::here(), error = function(e) NA_character_)
    if (!is.na(root) && dir.exists(root)) {
      return(normalizePath(root, winslash = "/", mustWork = TRUE))
    }
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

PROJECT_ROOT <- detect_project_root()
path_root <- function(...) file.path(PROJECT_ROOT, ...)

DIR_DATA    <- path_root("data")
DIR_SCRIPTS <- path_root("scripts")
DIR_RESULTS <- path_root("results")
DIR_FIGURES <- path_root("figures")
DIR_LOGS    <- path_root("logs")

required_dirs <- c(DIR_DATA, DIR_SCRIPTS)
missing_required <- required_dirs[!vapply(required_dirs, dir.exists, logical(1))]

if (length(missing_required) > 0) {
  stop(
    "Missing required folders:\n- ",
    paste(missing_required, collapse = "\n- "),
    "\nPROJECT_ROOT = ", PROJECT_ROOT
  )
}

dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_LOGS, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  write(msg, file = file.path(DIR_LOGS, "pipeline.log"), append = TRUE)
}

set.seed(1)
log_msg("OK: paths initialized | PROJECT_ROOT =", PROJECT_ROOT)