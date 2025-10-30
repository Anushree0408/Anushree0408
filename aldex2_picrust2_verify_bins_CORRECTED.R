#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ALDEx2)
})

# ---- Absolute paths (yours) ----
unstrat_fp <- "/media/hdd/14Tb_1/anushreejoshi/Filtered_7001_qc/picrust_7001/path_abun_unstrat.tsv.gz"
mapping_fp <- "/media/hdd/14Tb_1/anushreejoshi/Filtered_7001_qc/picrust_7001/Final_mapping_pathways.csv"
meta_fp    <- "/media/hdd/14Tb_1/anushreejoshi/Filtered_7001_qc/metadata.filtered.txt"
out_dir    <- "/media/hdd/14Tb_1/anushreejoshi/Filtered_7001_qc/picrust_7001/aldex2_out"

denom      <- "all"
mc_samples <- 128L

# ---------- helpers ----------
force_numeric_matrix <- function(df_mat) {
  # df_mat: data.frame with rownames=features, columns=samples
  # return: numeric matrix with same dimnames
  m <- as.matrix(df_mat)
  # if any column is not numeric, coerce safely
  if (!all(apply(m, 2, is.numeric))) {
    m <- suppressWarnings(matrix(
      as.numeric(m),
      nrow = nrow(m),
      ncol = ncol(m),
      dimnames = list(rownames(df_mat), colnames(df_mat))
    ))
  }
  storage.mode(m) <- "numeric"
  m[is.na(m)] <- 0
  return(m)
}

# ---- Check files ----
if (!all(file.exists(unstrat_fp, mapping_fp, meta_fp))) {
  stop(paste("Missing file(s). Checked:\n ", unstrat_fp, "\n ", mapping_fp, "\n ", meta_fp))
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("[INFO] Inputs:\n")
cat("  unstrat :", unstrat_fp, "\n")
cat("  mapping :", mapping_fp, "\n")
cat("  meta    :", meta_fp, "\n")
cat("  out_dir :", out_dir, "\n\n")

# ---- Mapping (CSV: path_id,pathway) ----
cat("[INFO] Reading mapping...\n")
map_dt <- fread(mapping_fp)
setnames(map_dt, tolower(names(map_dt)))
if (!"path_id" %in% names(map_dt)) setnames(map_dt, names(map_dt)[1], "path_id")
if (!"pathway" %in% names(map_dt)) {
  if (ncol(map_dt) > 1) setnames(map_dt, names(map_dt)[2], "pathway") else map_dt[, pathway := path_id]
}
map_dt[, path_id := as.character(path_id)]
map_lookup <- map_dt[, .(path_id, pathway)]
cat(sprintf("[OK] mapping rows: %d  (unique path_id: %d)\n\n",
            nrow(map_dt), uniqueN(map_dt$path_id)))

# ---- PICRUSt2 unstrat (TSV.GZ, first col = path_id) ----
cat("[INFO] Reading unstrat table...\n")
un_dt <- fread(unstrat_fp)                # data.table
setnames(un_dt, names(un_dt)[1], "path_id")
un_dt[, path_id := trimws(as.character(path_id))]
samp_cols <- setdiff(names(un_dt), "path_id")

# Coerce *strictly* to numeric using type.convert then as.numeric
for (cc in samp_cols) {
  # convert like "1e+03" etc.; keep as character->numeric to avoid factor
  if (!is.numeric(un_dt[[cc]])) {
    un_dt[[cc]] <- suppressWarnings(as.numeric(as.character(un_dt[[cc]])))
  }
}
un_dt[is.na(un_dt)] <- 0
cat(sprintf("[OK] unstrat: %d pathways x %d samples\n\n", nrow(un_dt), length(samp_cols)))

# Build counts data.frame with rownames
counts_df <- as.data.frame(un_dt)
rownames(counts_df) <- counts_df$path_id
counts_df$path_id <- NULL

# ---- Metadata (TAB: SampleID, Group, AGE) ----
cat("[INFO] Reading metadata...\n")
meta_dt <- fread(meta_fp, sep = "\t")
setnames(meta_dt, tolower(names(meta_dt)))
req <- c("sampleid","group","age")
if (!all(req %in% names(meta_dt))) stop("Metadata must contain: SampleID, Group, AGE")
setnames(meta_dt, c("sampleid","group","age"), c("SampleID","Group","AGE"))
meta_dt[, SampleID := as.character(SampleID)]
meta_dt[, Group    := toupper(trimws(as.character(Group)))]
meta_dt <- meta_dt[Group %in% c("CONTROL","IBD")]
meta_dt[, AGE := suppressWarnings(as.numeric(AGE))]
meta_dt <- meta_dt[!is.na(AGE)]
meta_dt[AGE < 0, AGE := 0]; meta_dt[AGE > 12, AGE := 12]
meta_dt[, AgeBin := ifelse(AGE < 5, "0-4", "5-12")]
cat(sprintf("[OK] metadata rows (CONTROL/IBD only): %d  | bins: %s\n\n",
            nrow(meta_dt),
            paste(sprintf("%s=%d", names(table(meta_dt$AgeBin)), as.integer(table(meta_dt$AgeBin))), collapse="; ")))

# ---- Align samples & force numeric matrix ----
common <- intersect(colnames(counts_df), meta_dt$SampleID)
if (!length(common)) stop("No overlapping SampleIDs.")
counts_df <- counts_df[, common, drop = FALSE]
meta_dt   <- meta_dt[SampleID %in% common]
meta_dt   <- meta_dt[match(colnames(counts_df), meta_dt$SampleID)]  # order by columns

# Make sure final matrix is numeric
counts_mat <- force_numeric_matrix(counts_df)

cat(sprintf("[OK] overlap: %d samples | CONTROL=%d, IBD=%d\n\n",
            length(common),
            sum(meta_dt$Group=="CONTROL"), sum(meta_dt$Group=="IBD")))

# ---- Run ALDEx2 per bin ----
run_bin <- function(bin_name){
  cat(sprintf("[INFO] ALDEx2 AgeBin=%s (CONTROL vs IBD)\n", bin_name))

  mb <- meta_dt[AgeBin == bin_name]
  if (!all(c("CONTROL","IBD") %in% unique(mb$Group))) {
    cat(sprintf("[WARN] Skipping %s: need both CONTROL and IBD.\n\n", bin_name)); return(invisible(NULL))
  }

  samples <- mb$SampleID
  samples <- samples[samples %in% colnames(counts_mat)]
  if (length(samples) < 4) {
    cat(sprintf("[WARN] Skipping %s: too few samples (n=%d)\n\n", bin_name, length(samples))); return(invisible(NULL))
  }

  # Subset matrix and drop all-zero pathways
  X <- counts_mat[, samples, drop = FALSE]
  keep <- rowSums(X) > 0
  X <- X[keep, , drop = FALSE]

  # Reorder metadata to match columns exactly
  mb <- mb[match(colnames(X), mb$SampleID), ]

  # FIXED: Convert to integers (ALDEx2 requirement)
  # PICRUSt2 outputs floating-point abundances, but ALDEx2 needs integers
  # Scale by 1000 and round to preserve precision while converting to integers
  X <- round(X * 1000)
  storage.mode(X) <- "integer"
  
  # FIXED: Condition as CHARACTER VECTOR (not factor)
  # ALDEx2 expects character vector, not factor
  cond <- as.character(mb$Group)
  
  # Ensure the order is CONTROL first, then IBD for proper comparison
  # ALDEx2 will use the first unique value as reference
  cond <- factor(cond, levels = c("CONTROL", "IBD"))
  cond <- as.character(cond)

  # Debug: Print first few samples and their conditions
  cat(sprintf("  [DEBUG] Samples: %d | Conditions: CONTROL=%d, IBD=%d\n", 
              length(cond), sum(cond=="CONTROL"), sum(cond=="IBD")))

  # Run ALDEx2
  # Step 1: Generate CLR values
  clr <- ALDEx2::aldex.clr(X, cond, mc.samples = mc_samples, denom = denom, verbose = FALSE)
  
  # Step 2: Calculate effect sizes (this gives us diff.btw, diff.win, effect, overlap)
  eff <- ALDEx2::aldex.effect(clr, CI = TRUE, verbose = FALSE)
  
  # Step 3: Calculate t-test statistics (this gives us we.ep, we.eBH, wi.ep, wi.eBH)
  tt  <- ALDEx2::aldex.ttest(clr, paired = FALSE, verbose = FALSE)

  # Merge effect sizes and t-test results
  res <- merge(as.data.frame(eff), as.data.frame(tt), by = "row.names", all = TRUE)
  res$path_id <- res$Row.names
  res$Row.names <- NULL
  rownames(res) <- NULL

  # Debug: Print column names to see what ALDEx2 actually returns
  cat(sprintf("  [DEBUG] ALDEx2 output columns: %s\n", paste(names(res), collapse=", ")))

  # FIXED: Check which effect size column exists and use it
  # ALDEx2 typically returns 'effect' from aldex.effect()
  effect_col <- NULL
  if ("effect" %in% names(res)) {
    effect_col <- "effect"
  } else if ("diff.btw" %in% names(res)) {
    effect_col <- "diff.btw"
  } else {
    stop(sprintf("Neither 'diff.btw' nor 'effect' found in ALDEx2 output. Columns are: %s", 
                 paste(names(res), collapse=", ")))
  }

  # Attach pathway names + direction
  res <- res %>%
    dplyr::left_join(map_lookup, by = "path_id") %>%
    dplyr::mutate(higher_in = ifelse(.data[[effect_col]] > 0, "IBD", "CONTROL"))
  
  # Sort by significance and effect size
  if ("we.eBH" %in% names(res)) {
    res <- dplyr::arrange(res, we.eBH, dplyr::desc(abs(.data[[effect_col]])))
  }

  tag <- gsub("-", "", bin_name)
  out_csv <- file.path(out_dir, sprintf("aldex2_%s_CONTROL_vs_IBD.csv", tag))
  data.table::fwrite(res, out_csv)
  cat("  [OK] Wrote:", out_csv, "\n")

  if ("we.eBH" %in% names(res)) {
    sig <- res[res$we.eBH < 0.05, ]
    sig_csv <- file.path(out_dir, sprintf("aldex2_%s_significant_weBH_lt_0.05.csv", tag))
    data.table::fwrite(sig, sig_csv)
    cat(sprintf("  [OK] Wrote: %s (n=%d significant pathways)\n\n", sig_csv, nrow(sig)))
  } else {
    cat("  [WARN] No we.eBH column found.\n\n")
  }
}

run_bin("0-4")
run_bin("5-12")
cat("[DONE]\n")
