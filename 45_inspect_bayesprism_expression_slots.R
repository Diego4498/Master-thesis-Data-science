# scripts/45_inspect_bayesprism_expression_slots.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(BayesPrism)
})

# Cambia el fold si quieres mirar otro
fold_name <- "mPCa_M1"

rds_path <- file.path(
  DIR_RESULTS, "validation", "pseudobulk_loso_level2", "bayesprism", "folds",
  paste0("bayesprism_level2_result_", fold_name, ".rds")
)

stopifnot(file.exists(rds_path))

bp_res <- readRDS(rds_path)

cat("\n===== bp_res =====\n")
print(class(bp_res))
print(slotNames(bp_res))

for (sl in slotNames(bp_res)) {
  cat("\n===== SLOT:", sl, "=====\n")
  obj <- slot(bp_res, sl)
  print(class(obj))
  
  # list-like
  if (is.list(obj)) {
    cat("Names:\n")
    print(names(obj))
  }
  
  # S4
  if (isS4(obj)) {
    cat("S4 slotNames:\n")
    print(slotNames(obj))
    for (ssl in slotNames(obj)) {
      cat("\n--- subslot:", ssl, "---\n")
      subobj <- slot(obj, ssl)
      print(class(subobj))
      
      if (is.matrix(subobj) || is.data.frame(subobj)) {
        cat("dim:\n")
        print(dim(subobj))
        cat("colnames head:\n")
        print(head(colnames(subobj)))
        cat("rownames head:\n")
        print(head(rownames(subobj)))
      }
      
      if (is.list(subobj)) {
        cat("names:\n")
        print(names(subobj))
      }
    }
  }
  
  # matrix/data.frame directly
  if (is.matrix(obj) || is.data.frame(obj)) {
    cat("dim:\n")
    print(dim(obj))
    cat("colnames head:\n")
    print(head(colnames(obj)))
    cat("rownames head:\n")
    print(head(rownames(obj)))
  }
}