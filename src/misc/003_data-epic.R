#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
usage <- function() {
  cat(
    "Usage: Rscript src/003_data-epic.R [options]\n",
    "\n",
    "Build a tab-delimited EPIC study phenotype file from:\n",
    "  genetics_caco.sas7bdat, genetics_id.sas7bdat, genetics.sas7bdat,\n",
    "  and optionally Subj_Id_2015.txt for supplementary IDs.\n",
    "\n",
    "Country-based exclusion (Country in 6, 8, B) is applied ONCE, after\n",
    "aggregation and the supplementary-ID merge. A sample is excluded if it is\n",
    "flagged for an excluded country in ANY source, and a sanity check aborts\n",
    "the run if any flagged sample survives into the output.\n",
    "\n",
    "Options:\n",
    "  --input-dir DIR            Directory containing the SAS files.\n",
    "  --output FILE              Output txt path.\n",
    "  --id-column NAME           ID column to use; default: Idepic_Bio.\n",
    "  --missing-phenotype VALUE  PLINK missing phenotype code; default: -9.\n",
    "  --subj-id-file FILE        Path to Subj_Id_2015.txt for supplementary IDs\n",
    "                             not in genetics_id.sas7bdat. Defaults to\n",
    "                             $SUBJ_ID_FILE or <input-dir>/Subj_Id_2015.txt.\n",
    "  --help                     Show this help.\n",
    "\n",
    sep = ""
  )
}
script_path <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    path <- sub("^--file=", "", file_arg[[1]])
    if (!startsWith(path, "/")) {
      submit_dir <- Sys.getenv("SLURM_SUBMIT_DIR", "")
      if (nchar(submit_dir) > 0) path <- file.path(submit_dir, path)
    }
    return(normalizePath(path, mustWork = FALSE))
  }
  submit_dir <- Sys.getenv("SLURM_SUBMIT_DIR", "")
  base <- if (nchar(submit_dir) > 0) submit_dir else "."
  normalizePath(file.path(base, "src/003_data-epic.R"), mustWork = FALSE)
}
parse_args <- function(args, defaults) {
  opts <- defaults
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help" || arg == "-h") {
      usage()
      quit(status = 0)
    }
    if (grepl("^--[^=]+=", arg)) {
      key <- sub("=.*$", "", arg)
      value <- sub("^[^=]+=", "", arg)
    } else if (grepl("^--", arg)) {
      key <- arg
      i <- i + 1
      if (i > length(args)) {
        stop("Missing value for ", key, call. = FALSE)
      }
      value <- args[[i]]
    } else {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }
    if (key == "--input-dir") {
      opts$input_dir <- value
    } else if (key == "--output") {
      opts$output <- value
    } else if (key == "--id-column") {
      opts$id_column <- value
    } else if (key == "--missing-phenotype") {
      opts$missing_phenotype <- as.integer(value)
      if (is.na(opts$missing_phenotype)) {
        stop("--missing-phenotype must be an integer.", call. = FALSE)
      }
    } else if (key == "--subj-id-file") {
      opts$subj_id_file <- value
    } else {
      stop("Unknown option: ", key, call. = FALSE)
    }
    i <- i + 1
  }
  opts
}
require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required R package is not installed: ", pkg,
      "\nInstall it before running this script.",
      call. = FALSE
    )
  }
}
read_sas_df <- function(path) {
  as.data.frame(haven::zap_labels(haven::read_sas(path)), stringsAsFactors = FALSE)
}
as_clean_character <- function(x) {
  x <- as.character(x)
  x[x == ""] <- NA_character_
  x
}
as_clean_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}
# ── Country-based exclusion helpers ──────────────────────────────────────────
# Exclusion is applied exactly once (in main), after aggregation and the
# supplementary-ID merge. The raw sources are intentionally NOT row-filtered on
# read, so that (a) per-ID country conflicts are still surfaced by
# warn_collapsed_conflicts(), and (b) the single exclusion below can see a
# sample's country in *every* source and drop it if flagged anywhere.
country_excluded_id_set <- function(sources, excluded_countries) {
  # sources: list of data.frames, each with $ID and $Country (NULL entries and
  #          sources lacking either column are skipped).
  # Returns the set of IDs carrying an excluded Country in ANY source. A sample
  # flagged in one source is excluded even if another source reports a
  # non-excluded or missing country for the same ID.
  excluded_ids <- character(0)
  for (src in sources) {
    if (is.null(src) || is.null(src$ID) || is.null(src$Country) || nrow(src) == 0) {
      next
    }
    ids <- as_clean_character(src$ID)
    country <- trimws(as_clean_character(src$Country))
    hit <- !is.na(ids) & !is.na(country) & country %in% excluded_countries
    excluded_ids <- union(excluded_ids, unique(ids[hit]))
  }
  excluded_ids
}
assert_no_excluded_countries <- function(out, excluded_id_set, excluded_countries,
                                         label = "final output") {
  # Sanity check: abort the run if any sample flagged for country-based
  # exclusion is still present, either by its (collapsed) country code or by
  # membership in the flagged-ID set. A passing run prints a confirmation.
  final_country <- trimws(as_clean_character(out$country))
  bad_country <- !is.na(final_country) & final_country %in% excluded_countries
  bad_id <- out$ID %in% excluded_id_set
  bad <- bad_country | bad_id
  if (any(bad)) {
    offending <- unique(out$ID[bad])
    stop(
      "SANITY CHECK FAILED: ", sum(bad),
      " sample(s) flagged for country-based exclusion remain in ", label,
      ". Example ID(s): ", paste(utils::head(offending, 5L), collapse = ", "),
      call. = FALSE
    )
  }
  message(
    "Sanity check passed: 0 country-excluded samples in ", label,
    " (", nrow(out), " rows)."
  )
}
collapse_character <- function(values) {
  values <- unique(stats::na.omit(as_clean_character(values)))
  if (length(values) == 0) {
    return(NA_character_)
  }
  values[[1]]
}
collapse_binary_any <- function(values) {
  values <- unique(stats::na.omit(as_clean_numeric(values)))
  if (length(values) == 0) {
    return(NA_integer_)
  }
  as.integer(max(values))
}
collapse_single_numeric <- function(values) {
  values <- unique(stats::na.omit(as_clean_numeric(values)))
  if (length(values) == 0) {
    return(NA_integer_)
  }
  if (length(values) == 1) {
    return(as.integer(values[[1]]))
  }
  NA_integer_
}
aggregate_by_id <- function(df, id_col, char_cols = character(), binary_cols = character(),
                            single_numeric_cols = character()) {
  df[[id_col]] <- as_clean_character(df[[id_col]])
  df <- df[!is.na(df[[id_col]]), , drop = FALSE]
  split_df <- split(df, df[[id_col]], drop = TRUE)
  out <- data.frame(ID = names(split_df), stringsAsFactors = FALSE)
  for (col in char_cols) {
    out[[col]] <- vapply(split_df, function(part) collapse_character(part[[col]]), character(1))
  }
  for (col in binary_cols) {
    out[[col]] <- vapply(split_df, function(part) collapse_binary_any(part[[col]]), integer(1))
  }
  for (col in single_numeric_cols) {
    out[[col]] <- vapply(split_df, function(part) collapse_single_numeric(part[[col]]), integer(1))
  }
  out
}
assert_columns <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}
warn_unexpected_values <- function(df, cols, allowed, label) {
  for (col in cols) {
    values <- unique(stats::na.omit(as_clean_numeric(df[[col]])))
    unexpected <- setdiff(values, allowed)
    if (length(unexpected) > 0) {
      warning(
        label, " column ", col, " has unexpected value(s): ",
        paste(unexpected, collapse = ", "),
        call. = FALSE
      )
    }
  }
}
warn_collapsed_conflicts <- function(df, id_col, cols, label, numeric = FALSE) {
  df[[id_col]] <- as_clean_character(df[[id_col]])
  df <- df[!is.na(df[[id_col]]), , drop = FALSE]
  split_df <- split(df, df[[id_col]], drop = TRUE)
  for (col in cols) {
    conflict <- vapply(
      split_df,
      function(part) {
        values <- if (numeric) as_clean_numeric(part[[col]]) else as_clean_character(part[[col]])
        length(unique(stats::na.omit(values))) > 1
      },
      logical(1)
    )
    n_conflict <- sum(conflict)
    if (n_conflict > 0) {
      warning(
        label, " column ", col, " has conflicting non-missing values for ",
        n_conflict, " ID(s); those collapsed values are set to missing.",
        call. = FALSE
      )
    }
  }
}
collapse_numeric_by_id <- function(df, id_col, value_col) {
  df[[id_col]] <- as_clean_character(df[[id_col]])
  df <- df[!is.na(df[[id_col]]), , drop = FALSE]
  split_df <- split(df, df[[id_col]], drop = TRUE)
  data.frame(
    ID = names(split_df),
    value = vapply(split_df, function(part) collapse_single_numeric(part[[value_col]]), integer(1)),
    stringsAsFactors = FALSE
  )
}
plink_case_status <- function(in_study, case_control, missing_phenotype) {
  in_study <- as_clean_numeric(in_study)
  case_control <- as_clean_numeric(case_control)
  out <- rep.int(missing_phenotype, length(in_study))
  # A sample is included if flagged in genetics_id (GW_* = 1) OR if it has a
  # caco record — the latter handles supplementary IDs from Subj_Id_2015 that
  # are not in genetics_id but do appear in genetics_caco for this study.
  included <- (!is.na(in_study) & in_study == 1) | !is.na(case_control)
  out[included & !is.na(case_control) & case_control == 0] <- 1L
  out[included & !is.na(case_control) & case_control == 1] <- 2L
  as.integer(out)
}
read_subj_id_file <- function(path, id_column) {
  # Subj_Id_2015.txt: no header, comma-separated.
  # Columns (1-indexed): Country(1), Center(2), ID_short(3), Sex(4), DOB(5),
  #   ID_long(6), [spare](7-8), Idepic(9), Idepic_Bio(10)
  # Only Idepic_Bio is supported as id_column here.
  if (id_column != "Idepic_Bio") {
    warning("read_subj_id_file: id_column '", id_column,
            "' is not Idepic_Bio; supplementary IDs will not be loaded.",
            call. = FALSE)
    return(NULL)
  }
  lines <- readLines(path, encoding = "latin1", warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]
  result <- lapply(lines, function(line) {
    parts <- strsplit(line, ",", fixed = TRUE)[[1]]
    if (length(parts) < 10) return(NULL)
    bio     <- trimws(parts[10])
    country <- trimws(parts[1])
    centre  <- trimws(parts[2])
    sex_raw <- trimws(parts[4])
    if (is.na(bio) || bio == "") return(NULL)
    sex <- if (sex_raw %in% c("1", "2")) as.integer(sex_raw) else NA_integer_
    list(ID = bio, Country = country, Center = centre, Sex = sex)
  })
  result <- result[!vapply(result, is.null, logical(1))]
  if (length(result) == 0) return(NULL)
  df <- data.frame(
    ID      = vapply(result, `[[`, character(1), "ID"),
    Country = vapply(result, `[[`, character(1), "Country"),
    Center  = vapply(result, `[[`, character(1), "Center"),
    Sex     = vapply(result, `[[`, integer(1),   "Sex"),
    stringsAsFactors = FALSE
  )
  # Keep first occurrence of each Idepic_Bio
  df[!duplicated(df$ID), , drop = FALSE]
}
study_case_control <- function(caco_raw, ids, included, study_id_col, caco_col, id_column) {
  caco_raw[[id_column]] <- as_clean_character(caco_raw[[id_column]])
  caco_raw$Proj_Acronym <- as_clean_character(caco_raw$Proj_Acronym)
  included <- as_clean_numeric(included)
  study_rows <- caco_raw[
    !is.na(caco_raw$Proj_Acronym) &
      caco_raw$Proj_Acronym == study_id_col &
      !is.na(caco_raw[[id_column]]),
    ,
    drop = FALSE
  ]
  out <- rep(NA_integer_, length(ids))
  if (nrow(study_rows) > 0) {
    warn_collapsed_conflicts(
      study_rows,
      id_column,
      caco_col,
      paste0("genetics_caco.sas7bdat Proj_Acronym=", study_id_col),
      numeric = TRUE
    )
    study_status <- collapse_numeric_by_id(study_rows, id_column, caco_col)
    out <- study_status$value[match(ids, study_status$ID)]
  } else {
    if (any(!is.na(included) & included == 1)) {
      warning(
        "genetics_caco.sas7bdat has no rows with Proj_Acronym=", study_id_col,
        "; included samples for this study will have missing phenotype.",
        call. = FALSE
      )
    }
  }
  out
}
study_map <- base::data.frame(
  study = c(
    "Brea_01_Erneg",
    "Brea_02_Onco",
    "Clrt_01_Gecco",
    "Ecvd_01",
    "Ecvd_02",
    "Ecvd_03",
    "Glbd_01",
    "Inte_01",
    "Inte_02",
    "Inte_03",
    "Kidn_01",
    "Kidn_02",
    "Lung_01",
    "Lymp_01",
    "Neuro_01",
    "Panc_01_PS1",
    "Panc_02_PS3",
    "Pros_01_Bpc3",
    "Pros_03_Onco",  
    "Pros_04_P160555",
    "Ovar_01",
    "Stom_01",
    "Uadt_01"
  ),
  id_col = c(
    "GW_BREA_01",
    "GW_BREA_02",
    "GW_CLRT_01",
    "GW_ECVD_01",
    "GW_ECVD_02",
    "GW_ECVD_03",
    "GW_GLBD_01",
    "GW_INTE_01",
    "GW_INTE_02",
    "GW_INTE_03",
    "GW_KIDN_01",
    "GW_KIDN_02",
    "GW_LUNG_01",
    "GW_LYMP_01",
    "GW_NEUR_01",
    "GW_PANC_01",
    "GW_PANC_02",
    "GW_PROS_01",
    "GW_PROS_03",  
    "GW_PROS_04",
    "GW_OVAR_01",
    "GW_STOM_01",
    "GW_UADT_01"
  ),
  caco_col = c(
    "Cncr_Caco_Brea_Orig",
    "Cncr_Caco_Brea_Orig",
    "Cncr_Caco_Clrt_Orig",
    "Ecvd_Case_Orig",
    "Ecvd_Case_Orig",
    "Ecvd_Case_Orig",
    "Cncr_Caco_Live_Orig",
    "InterAct_T2D_Status_Orig",
    "InterAct_T2D_Status_Orig",
    "InterAct_T2D_Status_Orig",
    "Cncr_Caco_Kidn_Orig",
    "Cncr_Caco_Kidn_Orig",
    "Cncr_Caco_Lung_Orig",
    "Cncr_Caco_Lymp_Orig",
    "Cncr_Caco_Neur_Orig",
    "Cncr_Caco_Panc_Orig",
    "Cncr_Caco_Panc_Orig",
    "Cncr_Caco_Pros_Orig",
    "Cncr_Caco_Pros_Orig", 
    "Cncr_Caco_Pros_Orig",
    "Cncr_Caco_Ovar_Orig",
    "Cncr_Caco_Stom_Orig",
    "Cncr_Caco_Uadt_Orig"
  ),
  stringsAsFactors = FALSE
)
main <- function() {
  require_namespace("haven")
  root <- normalizePath(file.path(dirname(script_path()), "..", ".."), mustWork = FALSE)
  default_input_dir <- file.path(root, "data", "reference", "Epic")
  defaults <- list(
    input_dir = Sys.getenv("EPIC_REF_DIR", default_input_dir),
    output = Sys.getenv(
      "EPIC_CASE_STATUS_FILE",
      file.path(default_input_dir, "EPIC_study_case_status.txt")
    ),
    id_column = Sys.getenv("EPIC_CASE_STATUS_ID_COLUMN", "Idepic_Bio"),
    missing_phenotype = as.integer(Sys.getenv("EPIC_MISSING_PHENOTYPE", "-9")),
    subj_id_file = Sys.getenv(
      "SUBJ_ID_FILE",
      file.path(default_input_dir, "Subj_Id_2015.txt")
    )
  )
  opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
  if (is.na(opts$missing_phenotype)) {
    stop("Missing phenotype code must be an integer.", call. = FALSE)
  }
  input_dir <- normalizePath(opts$input_dir, mustWork = FALSE)
  output <- normalizePath(opts$output, mustWork = FALSE)
  caco_file <- file.path(input_dir, "genetics_caco.sas7bdat")
  id_file <- file.path(input_dir, "genetics_id.sas7bdat")
  sex_file <- file.path(input_dir, "genetics.sas7bdat")
  for (path in c(caco_file, id_file, sex_file)) {
    if (!file.exists(path)) {
      stop("Input file not found: ", path, call. = FALSE)
    }
  }
  subj_id_file <- normalizePath(opts$subj_id_file, mustWork = FALSE)
  message("==========================================")
  message(" Building EPIC Study Case Status File")
  message(" Input dir:     ", input_dir)
  message(" Output:        ", output)
  message(" ID column:     ", opts$id_column)
  message(" Subj_Id file:  ", subj_id_file)
  message("==========================================")
  caco_raw <- read_sas_df(caco_file)
  id_raw <- read_sas_df(id_file)
  sex_raw <- read_sas_df(sex_file)
  excluded_countries <- c("6", "8", "B")
  # NOTE: Country-based exclusion is intentionally NOT applied per-source here.
  # It is applied ONCE, after aggregation and the supplementary-ID merge (see
  # the "Country-based exclusion" block below), so that a sample flagged in any
  # source cannot slip back in via the supplementary path or a per-ID country
  # conflict.
  assert_columns(id_raw, c("Country", "Center", "Idepic", opts$id_column, study_map$id_col), "genetics_id.sas7bdat")
  assert_columns(sex_raw, c(opts$id_column, "Sex"), "genetics.sas7bdat")

  # ── Idepic <-> Idepic_Bio map (authoritative EPIC IDs) ──────────────────────
  # Stage 4 SAMPLE_MANIFEST joins the finalised genotype sample IIDs to this map
  # to populate the manifest's Idepic and Idepic_Bio columns. The genotype IID is
  # a composite <Idepic>_<Idepic_Bio> for most samples and a single <Idepic_Bio>
  # for others, so string-parsing cannot always recover both IDs; this map does.
  idepic_map <- unique(data.frame(
    Idepic     = trimws(as_clean_character(id_raw$Idepic)),
    Idepic_Bio = trimws(as_clean_character(id_raw[[opts$id_column]])),
    stringsAsFactors = FALSE
  ))
  idepic_map <- idepic_map[idepic_map$Idepic != "" & idepic_map$Idepic_Bio != "", , drop = FALSE]
  idepic_map_path <- file.path(dirname(output), "EPIC_Idepic_map.tsv")
  write.table(idepic_map, idepic_map_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote Idepic<->Idepic_Bio map (", nrow(idepic_map), " rows) to: ", idepic_map_path)
  assert_columns(caco_raw, c(opts$id_column, "Proj_Acronym", unique(study_map$caco_col)), "genetics_caco.sas7bdat")
  # ── Sanity diagnostic: how Country codes compare as strings ────────────────
  # Confirms that the excluded codes ("6", "8", "B") match the actual values in
  # the data once coerced to character (catches e.g. a numeric column that
  # reads back as "6.0" instead of "6", which would silently miss exclusions).
  message("Country code distribution in genetics_id.sas7bdat (string form):")
  print(table(trimws(as_clean_character(id_raw$Country)), useNA = "ifany"))
  message("Excluded country codes: ", paste(excluded_countries, collapse = ", "))
  warn_unexpected_values(id_raw, study_map$id_col, c(0, 1), "genetics_id.sas7bdat")
  warn_unexpected_values(caco_raw, unique(study_map$caco_col), c(0, 1), "genetics_caco.sas7bdat")
  warn_unexpected_values(sex_raw, "Sex", c(1, 2), "genetics.sas7bdat")
  warn_collapsed_conflicts(id_raw, opts$id_column, c("Country", "Center", "Idepic"), "genetics_id.sas7bdat")
  warn_collapsed_conflicts(sex_raw, opts$id_column, "Sex", "genetics.sas7bdat", numeric = TRUE)
  id_agg <- aggregate_by_id(
    id_raw,
    opts$id_column,
    char_cols = c("Country", "Center", "Idepic"),
    binary_cols = study_map$id_col
  )
  sex_agg <- aggregate_by_id(
    sex_raw,
    opts$id_column,
    single_numeric_cols = "Sex"
  )
  # ── Supplementary IDs from Subj_Id_2015.txt ────────────────────────────────
  subj_df <- NULL
  if (file.exists(subj_id_file)) {
    message("Reading supplementary IDs from: ", subj_id_file)
    subj_df <- read_subj_id_file(subj_id_file, opts$id_column)
    if (!is.null(subj_df)) {
      # Not filtered here: country exclusion is applied once, post-merge. The
      # subj_df Country column still feeds that single exclusion as evidence.
      new_ids <- subj_df$ID[!subj_df$ID %in% id_agg$ID]
      message("  IDs in Subj_Id_2015 (pre country-exclusion): ", nrow(subj_df))
      message("  Already in genetics_id:        ", nrow(subj_df) - length(new_ids))
      message("  Supplementary IDs to add:      ", length(new_ids))
      if (length(new_ids) > 0) {
        subj_new <- subj_df[subj_df$ID %in% new_ids, , drop = FALSE]
        # Build supplement rows aligned with id_agg columns
        supp_agg <- data.frame(
          ID      = subj_new$ID,
          Country = subj_new$Country,
          Center  = subj_new$Center,
          Idepic  = NA_character_,
          stringsAsFactors = FALSE
        )
        for (gw_col in study_map$id_col) {
          supp_agg[[gw_col]] <- NA_integer_
        }
        id_agg <- rbind(id_agg, supp_agg)
      }
    }
  } else {
    message("Subj_Id_2015 file not found (", subj_id_file, "); skipping supplementary IDs.")
  }
  # ── Sex: genetics.sas7bdat primary, Subj_Id_2015 fallback ─────────────────
  sex_match <- match(id_agg$ID, sex_agg$ID)
  sex <- sex_agg$Sex[sex_match]
  if (!is.null(subj_df)) {
    subj_sex_lookup <- stats::setNames(subj_df$Sex, subj_df$ID)
    fallback_mask <- (is.na(sex) | !(sex %in% c(1L, 2L))) & id_agg$ID %in% names(subj_sex_lookup)
    sex[fallback_mask] <- subj_sex_lookup[id_agg$ID[fallback_mask]]
  }
  sex[is.na(sex) | !(sex %in% c(1L, 2L))] <- 0L
  out <- data.frame(
    ID = id_agg$ID,
    country = id_agg$Country,
    centre = id_agg$Center,
    sex = as.integer(sex),
    stringsAsFactors = FALSE
  )
  # ── Country-based exclusion (single point, post-merge) ─────────────────────
  # Collect every ID that carries an excluded Country in ANY source, then drop
  # it from BOTH id_agg and out (kept row-aligned) so that the per-study
  # counts and the written file reflect the exclusion. This is the only place
  # country exclusion happens, and it closes the leaks identified earlier:
  #   * supplementary re-admission (excluded in genetics_id, re-added via
  #     Subj_Id_2015 with a different/blank country),
  #   * within-ID country conflict resolved toward inclusion by the collapse,
  #   * a flag present only in caco/genetics but not genetics_id.
  excluded_id_set <- country_excluded_id_set(
    list(
      data.frame(ID = caco_raw[[opts$id_column]], Country = caco_raw[["Country"]],
                 stringsAsFactors = FALSE),
      data.frame(ID = id_raw[[opts$id_column]],   Country = id_raw[["Country"]],
                 stringsAsFactors = FALSE),
      data.frame(ID = sex_raw[[opts$id_column]],  Country = sex_raw[["Country"]],
                 stringsAsFactors = FALSE),
      if (!is.null(subj_df)) {
        data.frame(ID = subj_df$ID, Country = subj_df$Country,
                   stringsAsFactors = FALSE)
      } else {
        NULL
      }
    ),
    excluded_countries
  )
  # Write the comprehensive country-exclusion ID list (flagged in ANY source)
  # alongside the case-status file. Stage 1 reads this as EPIC_country_excluded.txt
  # and removes any genotype sample whose EPIC ID appears here, so the exclusion
  # carries through to the pipeline output.
  excluded_list_path <- file.path(dirname(output), "EPIC_country_excluded.txt")
  writeLines(excluded_id_set, excluded_list_path)
  message(
    "Wrote ", length(excluded_id_set), " country-excluded ID(s) to: ",
    excluded_list_path
  )
  final_country <- trimws(as_clean_character(out$country))
  exclude_mask <- (out$ID %in% excluded_id_set) |
    (!is.na(final_country) & final_country %in% excluded_countries)
  n_excluded <- sum(exclude_mask)
  if (n_excluded > 0) {
    message(
      "Country-based exclusion: removing ", n_excluded,
      " sample(s) flagged in any source (Country in ",
      paste(excluded_countries, collapse = ", "), ")."
    )
  } else {
    message("Country-based exclusion: no samples flagged.")
  }
  id_agg <- id_agg[!exclude_mask, , drop = FALSE]
  out <- out[!exclude_mask, , drop = FALSE]
  # Sanity check: abort if any flagged sample survived into the output.
  assert_no_excluded_countries(out, excluded_id_set, excluded_countries)
  summary_rows <- vector("list", nrow(study_map))
  for (i in seq_len(nrow(study_map))) {
    study <- study_map$study[[i]]
    id_col <- study_map$id_col[[i]]
    caco_col <- study_map$caco_col[[i]]
    included <- id_agg[[id_col]]
    case_control <- study_case_control(caco_raw, id_agg$ID, included, id_col, caco_col, opts$id_column)
    out[[study]] <- plink_case_status(included, case_control, opts$missing_phenotype)
    summary_rows[[i]] <- data.frame(
      study = study,
      included = sum(!is.na(included) & included == 1),
      controls = sum(out[[study]] == 1),
      cases = sum(out[[study]] == 2),
      missing_included = sum(!is.na(included) & included == 1 & out[[study]] == opts$missing_phenotype),
      stringsAsFactors = FALSE
    )
  }
  order_idx <- order(out$country, out$centre, out$ID, na.last = TRUE)
  out <- out[order_idx, , drop = FALSE]
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  write.table(out, output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  summary_df <- do.call(rbind, summary_rows)
  message("Rows written: ", nrow(out))
  message("Study columns: ", nrow(study_map))
  message("Case/control summary:")
  print(summary_df, row.names = FALSE)
  message("Done.")
}
if (sys.nframe() == 0) {
  main()
}