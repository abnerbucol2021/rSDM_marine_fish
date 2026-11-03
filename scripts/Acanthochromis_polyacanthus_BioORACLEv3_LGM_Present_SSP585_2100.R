###############################################################################
# Acanthochromis polyacanthus: three-time-slice marine species distribution model
# A: Last Glacial Maximum (~21 ka)
# B: Present (Bio-ORACLE v3, 2010-2020 decade)
# C: Future (Bio-ORACLE v3 SSP5-8.5, 2090-2100 decade)
#
# Predictors (THREE): mean sea-surface temperature (SST, deg C), mean
# sea-surface salinity (SSS, practical salinity), signed bathymetry (m).
#
# DATA STRATEGY: Bio-ORACLE v3 does not provide LGM data. For the LGM, compute
# paleo-MARSPEC anomalies (LGM minus MARSPEC present) for SST and SSS, then add
# these anomalies to the Bio-ORACLE baseline. Use paleo-MARSPEC LGM bathymetry.
# This is an exploratory delta-method cross-product hindcast; assess biases.
#
# FUTURE BATHYMETRY: Terrain remains fixed. As a simplified sensitivity
# assumption, sea-level rise increases water depth by SEA_LEVEL_RISE_M (default
# 0.77 m: IPCC AR6 SSP5-8.5 global-median rise by 2100 relative to 1995-2014).
# NO new flooded terrestrial cells are generated without a coastal DEM.
#
# Run in RStudio: source("D:/rSDM_marine_fish/<THIS_SCRIPT_NAME>.R")
# NOTE: Internet access is required on first run to retrieve environmental data.
# v4 FIXES (October 21, 2026):
#   * Windows permission-denied CSV when Excel has older output open:
#     create a unique output RUN directory on each invocation.
#   * Write CSV via a temporary file, then rename it; check access early.
#   * Do not overwrite prior figures, CSVs, or analysis outputs.
# v3 FIXES (October 21, 2026):
#   * RStudio 2026.07 subset-expression IDE bug: script runs in one function;
#     avoids df[, columns] code forms in interactive sections.
#   * terra::extract(SpatRaster, matrix) does NOT accept ID=FALSE.
#   * Defensive input/geometry checks and stage messages.
#   * Fully defined A/B/C predictions and plotting before the figure calls.
#   * Area totals corrected using cellSize(..., mask=TRUE) on predicted
#     presence only, with shallow-water area sanity checks.
#   * Two geographic figure sets (overall + Philippines), two maps each.
#   * Coastlines, LGM shelf reconstruction, bathymetric contours,
#     outward ticks, visible frames, larger typography and area changes.
# Author: Revised from user-supplied MARSPEC present/LGM SDM, October 2026.
###############################################################################

# -----------------------------------------------------------------------------
# 0. CONFIGURATION. EDIT THESE FIRST IF REQUIRED.
# -----------------------------------------------------------------------------
# IMPORTANT: DO NOT select all lines and run with Ctrl+Enter one-by-one.
# Save this file in D:/rSDM_marine_fish and call source(...) ONCE.
# Putting all calculations in a single function also avoids a known RStudio
# console parser regression affecting expressions like dataframe[, cols].
run_acanthochromis_sdm <- function() {
message("[SDM v4] Starting A=LGM, B=present, C=SSP5-8.5 (2090s)")
ROOT <- "D:/rSDM_marine_fish"
SPECIES <- "Acanthochromis polyacanthus"
OCCURRENCE_FILE <- file.path(ROOT, "data", "Acanthochromis_polyacanthus.csv")
# IMPORTANT (Windows): if a previous results CSV is open in Excel, overwriting
# that file will fail with "Permission denied". Every run now writes into a
# NEW dated folder under the SAME project output directory, with no data loss.
# Original cached downloads continue to be reused from DATA_DIR.
OUTPUT_BASE_DIR <- file.path(ROOT, "output", "Acanthochromis_polyacanthus_BioORACLEv3_LGM_Present_SSP585")
DATA_DIR <- file.path(ROOT, "data")
BO_CACHE <- file.path(DATA_DIR, "BioORACLE_v3_cached")
MARSPEC_CACHE <- DATA_DIR  # reuse prior sdmpredictors/MARSPEC downloads
TEMP_DIR <- file.path(ROOT, "temp")

MODEL_EXTENT <- terra::ext(60, 179.95, -40, 45)  # main Indo-Pacific calibration window
PH_EXTENT <- terra::ext(116, 127.5, 4, 21.5)    # map window (NOT Philippine EEZ)
MIN_DEPTH_M <- 0
MAX_DEPTH_M <- 20
SEA_LEVEL_RISE_M <- 0.77        # edit to 0.63, 1.01, or 0 to test sensitivity
BACKGROUND_N <- 2000
RANDOM_SEED <- 20210707

# False keeps the sea floor fixed AND omits SLR-induced additional water depth.
# TRUE (default) applies -SEA_LEVEL_RISE_M to existing underwater bathymetry.
ADJUST_FUTURE_WATER_DEPTH <- TRUE

# -----------------------------------------------------------------------------
# 1. PACKAGES AND DIRECTORIES
# -----------------------------------------------------------------------------
# Install once, in RStudio's Console (uncomment or copy these commands):
# install.packages(c("terra", "sdmpredictors", "pROC", "sf",
#                    "rnaturalearth", "rnaturalearthdata", "remotes"))
# remotes::install_github("bio-oracle/biooracler")
needed <- c("terra", "biooracler", "sdmpredictors", "pROC", "sf",
            "rnaturalearth", "rnaturalearthdata")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if ("biooracler" %in% missing) {
  message("'biooracler' is available from GitHub. Install once with:")
  message("  install.packages('remotes')")
  message("  remotes::install_github('bio-oracle/biooracler')")
}
if (length(missing)) stop("Install required packages before running: ",
                          paste(missing, collapse=", "))
suppressPackageStartupMessages({
  library(terra)
  library(sf)
})
for (d in c(DATA_DIR, BO_CACHE, MARSPEC_CACHE, OUTPUT_BASE_DIR, TEMP_DIR)) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

# Reserve a new session folder without touching CSVs opened by Excel/Windows.
# Dated folder avoids both CSV locks and overwriting earlier model runs.
new_output_dir <- function(parent) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  for (attempt in seq_len(999L)) {
    suffix <- if (attempt == 1L) "" else paste0("_", attempt)
    candidate <- file.path(parent, paste0("run_", stamp, suffix))
    if (dir.exists(candidate)) next
    if (dir.create(candidate, showWarnings=FALSE, recursive=TRUE)) {
      return(candidate)
    }
    stop("Cannot create a writable output directory: ", candidate,
         "\nCheck Windows permissions, OneDrive sync, or folder security.")
  }
  stop("Unable to allocate a unique output directory under: ", parent)
}
OUTPUT_DIR <- new_output_dir(OUTPUT_BASE_DIR)

# Check write permissions immediately rather than after lengthy predictions.
write_probe <- tempfile(pattern="write_permission_test_", tmpdir=OUTPUT_DIR,
                        fileext=".txt")
tryCatch({
  writeLines("Write test", write_probe, useBytes=TRUE)
  unlink(write_probe)
}, error=function(e) stop("Output folder is not writable: ", OUTPUT_DIR,
                          "\n", conditionMessage(e)))
message("This run will write NEW outputs to: ", OUTPUT_DIR)

# Save tables safely: first write to a temporary CSV within the same folder.
# Rename to final filename once complete. This also avoids truncated files.
write_csv_safely <- function(data, final_path) {
  folder <- dirname(final_path)
  staged_path <- tempfile(pattern=".csv_stage_", tmpdir=folder,
                          fileext=".csv")
  on.exit(if (file.exists(staged_path)) unlink(staged_path), add=TRUE)
  tryCatch(utils::write.csv(data, file=staged_path, row.names=FALSE),
    error=function(e) stop("Cannot write CSV in: ", folder, "\n",
                            conditionMessage(e)))
  if (file.exists(final_path)) {
    stop("Refusing to replace an existing output file: ", final_path,
         "\nUse a fresh run folder or close the file in Excel.")
  }
  if (!file.rename(staged_path, final_path)) {
    stop("Failed to finalize CSV: ", final_path,
         "\nCheck file locks and Windows permissions.")
  }
  invisible(final_path)
}

options(sdmpredictors_datadir=MARSPEC_CACHE)
terra::terraOptions(tempdir=TEMP_DIR, memfrac=0.65, progress=1)
set.seed(RANDOM_SEED)
message("[1/10] Packages and working directories checked")
log_path <- file.path(OUTPUT_DIR, "run_configuration.txt")

# Some modelling extents cross longitude 180. The present study intentionally
# stops at 179.95 E to avoid the antimeridian wrap in this initial model.

# -----------------------------------------------------------------------------
# 2. FETCH AND CACHE INDIVIDUAL Bio-ORACLE v3 VARIABLES FROM ERDDAP
# -----------------------------------------------------------------------------
# Verified IDs and variables:
#   thetao_baseline_2000_2019_depthsurf   / thetao_mean (deg C)
#   so_baseline_2000_2019_depthsurf      / so_mean (practical salinity)
#   thetao_ssp585_2020_2100_depthsurf    / thetao_mean (deg C)
#   so_ssp585_2020_2100_depthsurf       / so_mean (practical salinity)
#   terrain_characteristics             / bathymetry_mean (m, signed)
# NOTE: Dataset IDs have '2000_2019', but the layer's 2010 timestamp denotes
# the 2010-2020 decade. Future 2090 timestamp denotes the 2090-2100 decade.
# To save bandwidth, request only one decade and crop to the study region.

bo_get <- function(dataset_id, variable, time_string, out_name) {
  out_file <- file.path(BO_CACHE, paste0(out_name, ".tif"))
  if (file.exists(out_file)) {
    message("Reading cached Bio-ORACLE: ", out_name)
    x <- terra::rast(out_file)
    names(x) <- out_name
    return(x)
  }
  message("Downloading Bio-ORACLE dataset ", dataset_id,
          ", variable ", variable, ", time ", time_string)
  limits <- list(
    time=c(time_string, time_string),
    latitude=c(terra::ymin(MODEL_EXTENT), terra::ymax(MODEL_EXTENT)),
    longitude=c(terra::xmin(MODEL_EXTENT), terra::xmax(MODEL_EXTENT))
  )
  raw <- tryCatch(
    biooracler::download_layers(
      dataset_id=dataset_id,
      variables=variable,
      constraints=limits,
      fmt="raster",
      directory=BO_CACHE
    ),
    error=function(e) stop("Bio-ORACLE download failed for ", out_name, ": ",
                            conditionMessage(e),
                            "\nCheck Internet/ERDDAP access and dataset metadata, or retry.")
  )
  x <- if (inherits(raw, "SpatRaster")) raw else terra::rast(raw)
  # Single timestamp and single variable expected; do not average or select
  # accidentally across decades (an unintentional average is a serious error).
  if (terra::nlyr(x) != 1L) {
    stop("Expected exactly ONE Bio-ORACLE layer for ", out_name,
         " but downloaded ", terra::nlyr(x),
         ". Check time selection and metadata.")
  }
  x <- terra::crop(x, MODEL_EXTENT, snap="near")
  if (!terra::is.lonlat(x)) stop("Bio-ORACLE output is not lon/lat: ", out_name)
  names(x) <- out_name
  terra::writeRaster(x, out_file, overwrite=TRUE, datatype="FLT4S")
  x
}

message("[2/10] Bio-ORACLE present/future layers (cached files reused)")
present_sst <- bo_get("thetao_baseline_2000_2019_depthsurf", "thetao_mean",
                      "2010-01-01T00:00:00Z", "present_sst")
present_sss <- bo_get("so_baseline_2000_2019_depthsurf", "so_mean",
                      "2010-01-01T00:00:00Z", "present_salinity")
present_bathy <- bo_get("terrain_characteristics", "bathymetry_mean",
                        "1970-01-01T00:00:00Z", "present_bathymetry")
future_sst <- bo_get("thetao_ssp585_2020_2100_depthsurf", "thetao_mean",
                     "2090-01-01T00:00:00Z", "future_sst_SSP585_2090_2100")
future_sss <- bo_get("so_ssp585_2020_2100_depthsurf", "so_mean",
                     "2090-01-01T00:00:00Z", "future_salinity_SSP585_2090_2100")

# Align all grids to present SST exactly. Bilinear for continuous predictors.
# Resampling is only performed when necessary, e.g. slight ERDDAP extent shifts.
match_grid <- function(x, template, method="bilinear") {
  if (terra::compareGeom(x, template, stopOnError=FALSE)) return(x)
  if (!terra::same.crs(x, template)) {
    terra::project(x, template, method=method)
  } else terra::resample(x, template, method=method)
}
present_sss <- match_grid(present_sss, present_sst)
present_bathy <- match_grid(present_bathy, present_sst)
future_sst <- match_grid(future_sst, present_sst)
future_sss <- match_grid(future_sss, present_sst)

# Ocean depths should be NEGATIVE (signed elevation). Stop rather than silently
# flipping any unverified bathymetry dataset's sign convention.
mean_bathy <- as.numeric(terra::global(present_bathy, "mean", na.rm=TRUE)[1,1])
if (!is.finite(mean_bathy) || mean_bathy >= 0) {
  stop("Check the Bio-ORACLE bathymetry sign convention; ocean-depth mean is ",
       mean_bathy, ". Expected a negative signed seabed elevation.")
}

# -----------------------------------------------------------------------------
# 3. PAST (LGM) MARSPEC DATA: FIVE LAYERS ONLY
# -----------------------------------------------------------------------------
# Scaling: paleo-MARSPEC SST and SSS source TIFFs may store integer values x100.
# Scale correction is applied to both modern and LGM MARSPEC before differencing.
ms_codes <- c(
  "MS_bathy_5m", "MS_bathy_21kya",
  "MS_biogeo08_sss_mean_5m", "MS_biogeo08_sss_mean_21kya_adjCCSM",
  "MS_biogeo13_sst_mean_5m", "MS_biogeo13_sst_mean_21kya_adjCCSM"
)
message("[3/10] MARSPEC present + LGM layers (cached files reused)")
ms <- lapply(ms_codes, function(code) {
  z <- sdmpredictors::load_layers(code, rasterstack=FALSE,
                                   datadir=MARSPEC_CACHE)
  # load_layers with rasterstack=FALSE returns a list.
  if (is.list(z) && !inherits(z, "SpatRaster")) z <- z[[1]]
  z <- terra::rast(z)
  terra::crop(z, MODEL_EXTENT)
})
names(ms) <- ms_codes

ms_fix_units <- function(x, variable) {
  v <- as.numeric(terra::global(x, "mean", na.rm=TRUE)[1,1])
  if (!is.finite(v)) stop("MARSPEC layer is completely NA: ", variable)
  if (abs(v) > 100) {
    message("Applying /100 MARSPEC scaling for ", variable,
            "; raw mean was ", round(v,2))
    x <- x / 100
  }
  v2 <- as.numeric(terra::global(x, "mean", na.rm=TRUE)[1,1])
  if (variable == "sst" && (v2 < -5 || v2 > 40)) {
    stop("Unreasonable MARSPEC SST mean after unit check: ", v2)
  }
  if (variable == "salinity" && (v2 < 0 || v2 > 48)) {
    stop("Unreasonable MARSPEC SSS mean after unit check: ", v2)
  }
  x
}
ms_sst_present <- ms_fix_units(ms[["MS_biogeo13_sst_mean_5m"]], "sst")
ms_sst_LGM <- ms_fix_units(ms[["MS_biogeo13_sst_mean_21kya_adjCCSM"]], "sst")
ms_sss_present <- ms_fix_units(ms[["MS_biogeo08_sss_mean_5m"]], "salinity")
ms_sss_LGM <- ms_fix_units(ms[["MS_biogeo08_sss_mean_21kya_adjCCSM"]], "salinity")
ms_bathy_present <- ms[["MS_bathy_5m"]]
ms_bathy_LGM <- ms[["MS_bathy_21kya"]]

# The LGM - MARSPEC-present difference on the common Bio-ORACLE raster grid.
# Do NOT simply replace Bio-ORACLE baseline values with MARSPEC absolute values.
lgm_sst_delta <- match_grid(ms_sst_LGM, present_sst) -
                 match_grid(ms_sst_present, present_sst)
lgm_sss_delta <- match_grid(ms_sss_LGM, present_sst) -
                 match_grid(ms_sss_present, present_sst)
lgm_bathy <- match_grid(ms_bathy_LGM, present_sst)

# Paleo-MARSPEC LGM ocean bathymetry already incorporates ~120 m lower sea level;
# exposed shallow shelf cells can be NA. Use a NEAREST resampled paleo sea mask.
lgm_native_ocean <- terra::ifel(ms_bathy_LGM < 0, 1, NA)
lgm_ocean_mask <- match_grid(lgm_native_ocean, present_sst, "near")
lgm_ocean_mask <- terra::ifel(lgm_ocean_mask == 1, 1, NA)

# -----------------------------------------------------------------------------
# 4. ASSEMBLE EXACTLY THREE PREDICTORS AT EACH TIME SLICE
# -----------------------------------------------------------------------------
# Signed bathymetry (negative elevation) acts as the depth predictor.
# For future, rising sea level causes seabed elevation relative to the NEW
# sea surface to become more negative: -10.0 -> -10.77 m, for +0.77 m SLR.
future_bathy <- if (ADJUST_FUTURE_WATER_DEPTH) {
  present_bathy - SEA_LEVEL_RISE_M
} else present_bathy
names(future_bathy) <- "bathymetry"

mk_stack <- function(b, s, t) {
  v <- c(b, s, t)
  names(v) <- c("bathymetry", "salinity", "sst")
  v
}
present_env <- mk_stack(present_bathy, present_sss, present_sst)
future_env <- mk_stack(future_bathy, future_sss, future_sst)
lgm_env <- mk_stack(lgm_bathy,
                    present_sss + lgm_sss_delta,
                    present_sst + lgm_sst_delta)

present_mask <- terra::ifel(present_env$bathymetry < 0, 1, NA)
future_mask <- terra::ifel(present_bathy < 0, 1, NA)
# IMPORTANT: no newly flooded present-land grid cells can be predicted in 2100.
lgm_mask <- terra::ifel(lgm_env$bathymetry < 0 & !is.na(lgm_ocean_mask), 1, NA)
present_env <- terra::mask(present_env, present_mask)
future_env <- terra::mask(future_env, future_mask)
lgm_env <- terra::mask(lgm_env, lgm_mask)
message("[4/10] Built matching bathymetry, salinity, SST rasters")
scenario_rasters <- list(present=present_env, LGM=lgm_env, future=future_env)
for (scenario_name in names(scenario_rasters)) {
  scenario_raster <- scenario_rasters[[scenario_name]]
  # Check that there are marine cells with three non-NA predictor values.
  valid_three <- terra::ifel(!is.na(scenario_raster$bathymetry) &
                            !is.na(scenario_raster$salinity) &
                            !is.na(scenario_raster$sst), 1, NA)
  valid_count <- terra::global(valid_three, "sum", na.rm=TRUE)[1, 1]
  if (!is.finite(valid_count) || valid_count < 100) {
    stop("Insufficient valid predictor cells for ", scenario_name,
         " (", valid_count, "). Check layer coverage.")
  }
}

if (!identical(names(present_env), names(future_env)) ||
    !identical(names(present_env), names(lgm_env))) {
  stop("Predictor names/order differ between scenarios.")
}

# Optional audit and output rasters.
terra::writeRaster(lgm_sst_delta, file.path(OUTPUT_DIR, "LGM_minus_modern_MARSPEC_SST_delta_C.tif"), overwrite=TRUE)
terra::writeRaster(lgm_sss_delta, file.path(OUTPUT_DIR, "LGM_minus_modern_MARSPEC_salinity_delta.tif"), overwrite=TRUE)
terra::writeRaster(future_env$bathymetry, file.path(OUTPUT_DIR, "future_2100_signed_bathymetry_SLRSensitivity.tif"), overwrite=TRUE)

# -----------------------------------------------------------------------------
# 5. OCCURRENCE CLEANING AND BACKGROUND POINTS
# -----------------------------------------------------------------------------
if (!file.exists(OCCURRENCE_FILE)) {
  stop("Occurrence file missing: ", OCCURRENCE_FILE,
       "\nExpected 'longitude' and 'latitude' column headings.")
}
message("[5/10] Cleaning occurrences and sampling background")
obs <- utils::read.csv(OCCURRENCE_FILE)
if (!all(c("longitude","latitude") %in% names(obs))) {
  stop("Occurrence CSV must contain 'longitude' and 'latitude' columns.")
}
obs <- obs[c("longitude", "latitude")]
obs$longitude <- suppressWarnings(as.numeric(as.character(obs$longitude)))
obs$latitude <- suppressWarnings(as.numeric(as.character(obs$latitude)))
obs <- obs[stats::complete.cases(obs), , drop=FALSE]
inside_model <- obs$longitude >= terra::xmin(MODEL_EXTENT) &
                obs$longitude <= terra::xmax(MODEL_EXTENT) &
                obs$latitude >= terra::ymin(MODEL_EXTENT) &
                obs$latitude <= terra::ymax(MODEL_EXTENT)
obs <- obs[inside_model, , drop=FALSE]
if (nrow(obs) == 0L) stop("No occurrences fall inside MODEL_EXTENT.")
obs <- unique(obs)
# Thin duplicate presences within the same 0.05-degree modelling raster cell.
obs$cell <- terra::cellFromXY(present_env[[1]], as.matrix(obs[c("longitude", "latitude")]))
obs <- obs[!is.na(obs$cell) & !duplicated(obs$cell), , drop=FALSE]
# Matrix (x,y) extraction never takes ID=FALSE. That argument is specific
# to extraction by a SpatVector and triggered the error in the provided log.
obs_vals <- terra::extract(present_env, as.matrix(obs[c("longitude", "latitude")]))
obs <- obs[stats::complete.cases(obs_vals), , drop=FALSE]
if (nrow(obs) < 10) stop("Fewer than 10 usable marine occurrences (",nrow(obs),
                         "). Check the CSV, extent, and environmental raster.")
message("Marine occurrence cells for modelling: ", nrow(obs))

# Sample candidate pseudo-absences from present marine grid with all predictors.
complete_present <- terra::ifel(!is.na(present_env$bathymetry) &
                                !is.na(present_env$salinity) &
                                !is.na(present_env$sst), 1, NA)
set.seed(RANDOM_SEED)
bg_xy <- terra::spatSample(complete_present,
                           size=max(4000, BACKGROUND_N*3),
                           method="random", values=FALSE,
                           na.rm=TRUE, xy=TRUE, as.df=TRUE)
bg_xy <- as.data.frame(bg_xy)
if (!all(c("x", "y") %in% names(bg_xy))) {
  stop("terra::spatSample returned coordinates with unexpected column names: ",
       paste(names(bg_xy),collapse=", "))
}
bg_xy$cell <- terra::cellFromXY(present_env[[1]], as.matrix(bg_xy[c("x", "y")]))
bg_xy <- bg_xy[!duplicated(bg_xy$cell) & !(bg_xy$cell %in% obs$cell), , drop=FALSE]
if (nrow(bg_xy) < 100) stop("Too few non-overlapping background cells.")
bg_xy <- head(bg_xy, BACKGROUND_N)

presence_xy <- obs[c("longitude", "latitude")]
background_xy <- data.frame(longitude=bg_xy$x,latitude=bg_xy$y)
model_xy <- rbind(presence_xy, background_xy)
pa <- c(rep(1L,nrow(presence_xy)), rep(0L,nrow(background_xy)))
# y is a two-column longitude/latitude matrix; no ID= argument here.
vals <- as.data.frame(terra::extract(present_env, as.matrix(model_xy)))
if (!identical(names(vals), c("bathymetry", "salinity", "sst"))) {
  stop("Unexpected predictor columns from terra::extract: ",
       paste(names(vals), collapse=", "))
}
if (nrow(vals) != nrow(model_xy)) {
  stop("terra::extract returned a different number of rows than the points.")
}
model_df <- data.frame(pa=pa, vals)
model_df <- model_df[stats::complete.cases(model_df), , drop=FALSE]
if (length(unique(model_df$pa)) < 2) stop("Only one response class after extraction.")

# -----------------------------------------------------------------------------
# 6. FIT PRESENT-DAY BINOMIAL GLM AND EVALUATE ON A HOLDOUT
# -----------------------------------------------------------------------------
message("[6/10] GLM fitting and holdout AUC")
set.seed(RANDOM_SEED)
# Stratified 80/20 random split. Non-spatial split may overestimate transferability.
train_ids <- unlist(lapply(split(seq_len(nrow(model_df)),model_df$pa),
                           function(ii) sample(ii, size=max(1, floor(0.8*length(ii))))))
training <- model_df[train_ids, , drop=FALSE]
testing <- model_df[-train_ids, , drop=FALSE]
if (length(unique(training$pa)) < 2 || length(unique(testing$pa)) < 2) {
  stop("Need BOTH classes in train and test; check occurrence/background counts.")
}
model <- stats::glm(pa ~ bathymetry + salinity + sst,
                    data=training, family=stats::binomial())
message("GLM coefficients:")
print(stats::coef(summary(model)))
if (any(!is.finite(stats::coef(model)))) {
  warning("GLM has non-finite coefficients. Risk of quasi/complete separation; ",
          "predictions may be unreliable. Check presences/background sampling.")
}
# A high random-split AUC can be optimistic when background sampling separates
# environmental niches or observations are spatially autocorrelated.
# Do not treat this validation as independent spatial cross-validation.

# Select threshold ONLY using training data (avoid tuning on holdout testing).
train_p <- stats::predict(model, newdata=training, type="response")
train_roc <- pROC::roc(training$pa, train_p, levels=c(0,1),
                       direction="<", quiet=TRUE)
threshold <- as.numeric(pROC::coords(train_roc, x="best", best.method="youden",
                                     ret="threshold", transpose=FALSE))[1]
if (!is.finite(threshold)) stop("Could not calculate a finite GLM threshold.")
test_p <- stats::predict(model, newdata=testing, type="response")
test_roc <- pROC::roc(testing$pa,test_p, levels=c(0,1),
                      direction="<",quiet=TRUE)
auc_test <- as.numeric(pROC::auc(test_roc))
message(sprintf("Random holdout AUC = %.3f; train-derived threshold = %.3f",
                auc_test, threshold))
write_csv_safely(
  data.frame(variable=names(coef(model)), coefficient=as.numeric(coef(model))),
  file.path(OUTPUT_DIR,"model_coefficients.csv"))
saveRDS(model,file.path(OUTPUT_DIR,"GLM_BioORACLE_v3_3predictors.rds"))

# Predictor correlation diagnostic on non-NA present marine cells.
cor_samples <- terra::spatSample(present_env,size=10000,method="random",
                                 values=TRUE,na.rm=TRUE,as.df=TRUE)
cor_mat <- stats::cor(as.data.frame(cor_samples),use="pairwise.complete.obs")
write_csv_safely(as.data.frame(cor_mat),
  file.path(OUTPUT_DIR,"present_predictor_correlations.csv"))

# -----------------------------------------------------------------------------
# 7. PROJECT SAME TRAINED MODEL TO A/B/C (UNMASKED DEPTH RANGE FIRST)
# -----------------------------------------------------------------------------
project_one <- function(x, tag) {
  destination <- file.path(OUTPUT_DIR,paste0(tag,"_full_extent_suitability.tif"))
  terra::predict(x, model, type="response", filename=destination,
                 overwrite=TRUE, wopt=list(datatype="FLT4S",gdal="COMPRESS=LZW"))
}
message("[7/10] Projecting GLM to all three periods")
p_lgm <- project_one(lgm_env,"A_LGM_21ka")
p_present <- project_one(present_env,"B_Present_2010_2020")
p_future <- project_one(future_env,"C_Future_SSP585_2090_2100")

# -----------------------------------------------------------------------------
# 8. BINARY PROJECTIONS, DEPTH MASKS AND CORRECTED AREA STATISTICS
# -----------------------------------------------------------------------------
# ALL overall analyses cover MODEL_EXTENT, and are not limited to 0-20 m.
# ONLY the Philippines comparison applies the post-model 0-20 m depth mask.
# Signed bathymetry convention: -20 <= depth < 0 metres is 0-20 m water.
message("[8/10] Binary habitat predictions and area corrections")

binary_raster <- function(suitability, cutpoint) {
  terra::ifel(is.na(suitability), NA,
              terra::ifel(suitability >= cutpoint, 1, 0))
}
shallow_mask <- function(bathymetry) {
  terra::ifel(!is.na(bathymetry) &
              bathymetry < -MIN_DEPTH_M &
              bathymetry >= -MAX_DEPTH_M, 1, NA)
}

# COUNT non-NA cells, including zero-valued absence cells.
count_valid_cells <- function(x) {
  value <- as.numeric(terra::global(
    terra::ifel(!is.na(x), 1, 0), "sum", na.rm=TRUE)[1, 1])
  if (!is.finite(value)) 0 else value
}

# Convert *only* selected cells to a raster mask; others MUST be NA.
# cellSize(...,mask=TRUE) is essential: mask=FALSE would count the
# entire raster rectangle, which generated impossible Philippines totals.
area_of_non_na_km2 <- function(x) {
  if (count_valid_cells(x) == 0) return(0)
  areas <- terra::cellSize(x, mask=TRUE, unit="km")
  total <- as.numeric(terra::global(areas, "sum", na.rm=TRUE)[1, 1])
  if (!is.finite(total)) 0 else total
}

suitable_area_km2 <- function(binary) {
  only_presence <- terra::ifel(!is.na(binary) & binary == 1, 1, NA)
  area_of_non_na_km2(only_presence)
}

# Full model extent suitability and binary maps.
periods <- list(
  A=list(tag="A_LGM_21ka", label="A. LGM (~21 ka)",
         suitability=p_lgm,
         binary=binary_raster(p_lgm, threshold),
         bathy=lgm_env$bathymetry),
  B=list(tag="B_Present_2010_2020", label="B. Present (2010-2020)",
         suitability=p_present,
         binary=binary_raster(p_present, threshold),
         bathy=present_env$bathymetry),
  C=list(tag="C_Future_SSP585_2090_2100", label="C. SSP5-8.5 (2090-2100)",
         suitability=p_future,
         binary=binary_raster(p_future, threshold),
         bathy=future_env$bathymetry)
)

stat_rows <- list()
for (id in names(periods)) {
  z <- periods[[id]]
  z$shallow <- shallow_mask(z$bathy)
  z$ph_suitability <- terra::mask(
    terra::crop(z$suitability, PH_EXTENT),
    terra::crop(z$shallow, PH_EXTENT))
  z$ph_binary <- terra::mask(
    terra::crop(z$binary, PH_EXTENT),
    terra::crop(z$shallow, PH_EXTENT))
  z$ph_bathy <- terra::crop(z$bathy, PH_EXTENT)
  z$ph_only_presence <- terra::ifel(
    !is.na(z$ph_binary) & z$ph_binary == 1, 1, NA)
  z$overall_only_presence <- terra::ifel(
    !is.na(z$binary) & z$binary == 1, 1, NA)
  z$ph_valid_mask <- terra::ifel(!is.na(z$ph_suitability), 1, NA)
  z$overall_valid_mask <- terra::ifel(!is.na(z$suitability), 1, NA)

  overall_area <- suitable_area_km2(z$binary)
  overall_valid_area <- area_of_non_na_km2(z$overall_valid_mask)
  ph_area <- suitable_area_km2(z$ph_binary)
  ph_valid_area <- area_of_non_na_km2(z$ph_valid_mask)
  if (overall_area > overall_valid_area + 0.01 ||
      ph_area > ph_valid_area + 0.01) {
    stop("Area sanity check FAILED for ", z$tag,
         ": predicted area exceeds all valid modelling cells.\n",
         "Overall: ", overall_area, " > ", overall_valid_area,
         "; Philippines: ", ph_area, " > ", ph_valid_area,
         ". Check masked raster and terra versions.")
  }
  stat_rows[[id]] <- data.frame(
    period=z$label,
    full_extent_predicted_area_km2=overall_area,
    full_extent_valid_marine_area_km2=overall_valid_area,
    full_extent_predicted_presence_cells=count_valid_cells(z$overall_only_presence),
    Philippines_window_0_20m_predicted_area_km2=ph_area,
    Philippines_window_0_20m_valid_area_km2=ph_valid_area,
    Philippines_window_0_20m_cells=count_valid_cells(z$ph_suitability),
    Philippines_window_0_20m_predicted_presence_cells=count_valid_cells(z$ph_only_presence),
    stringsAsFactors=FALSE
  )
  terra::writeRaster(z$binary,
    file.path(OUTPUT_DIR,paste0(z$tag,"_overall_range_binary.tif")),
    overwrite=TRUE, wopt=list(gdal="COMPRESS=LZW"))
  terra::writeRaster(z$ph_suitability,
    file.path(OUTPUT_DIR,paste0(z$tag,"_Philippines_0_20m_suitability.tif")),
    overwrite=TRUE, wopt=list(gdal="COMPRESS=LZW"))
  terra::writeRaster(z$ph_binary,
    file.path(OUTPUT_DIR,paste0(z$tag,"_Philippines_0_20m_binary.tif")),
    overwrite=TRUE, wopt=list(gdal="COMPRESS=LZW"))
  periods[[id]] <- z
}

stats_out <- do.call(rbind, stat_rows)
row.names(stats_out) <- NULL
base_full <- stats_out$full_extent_predicted_area_km2[2]
base_ph <- stats_out$Philippines_window_0_20m_predicted_area_km2[2]
if (base_full <= 0 || base_ph <= 0) {
  stop("Present reference suitable area is zero for one of the two regions. ",
       "Percentage differences cannot be computed.")
}
# Positive means more suitable area than PRESENT; negative means less.
stats_out$change_vs_present_full_percent <- 100 * (
  stats_out$full_extent_predicted_area_km2 - base_full) / base_full
stats_out$change_vs_present_PH_0_20m_percent <- 100 * (
  stats_out$Philippines_window_0_20m_predicted_area_km2 - base_ph) / base_ph

write_csv_safely(stats_out,
  file.path(OUTPUT_DIR,"A_B_C_suitable_area_summary.csv"))

# Independent, easier-to-read comparisons with present as the denominator.
# Past:     100 * (LGM - present) / present
# Future:   100 * (future - present) / present
comparison_out <- do.call(rbind, lapply(
  list(
    list(name="Overall model extent (all modelled marine depths)",
         area_col="full_extent_predicted_area_km2", pct_col="change_vs_present_full_percent"),
    list(name="Philippines map window (0-20 m water depth)",
         area_col="Philippines_window_0_20m_predicted_area_km2",pct_col="change_vs_present_PH_0_20m_percent")
  ), function(region) {
    data.frame(
      region=region$name,
      comparison=c("LGM relative to present", "Future relative to present"),
      present_reference_area_km2=c(stats_out[[region$area_col]][2],stats_out[[region$area_col]][2]),
      comparison_area_km2=c(stats_out[[region$area_col]][1],stats_out[[region$area_col]][3]),
      difference_km2=c(stats_out[[region$area_col]][1]-stats_out[[region$area_col]][2],
                       stats_out[[region$area_col]][3]-stats_out[[region$area_col]][2]),
      percent_difference_vs_present=c(stats_out[[region$pct_col]][1],stats_out[[region$pct_col]][3]),
      stringsAsFactors=FALSE
    )
  }
))
write_csv_safely(comparison_out,
  file.path(OUTPUT_DIR,"Present_vs_LGM_and_Future_percent_comparisons.csv"))
message("Corrected suitable area estimates (km2) and percent change from PRESENT:")
print(stats_out)
print(comparison_out)
if (any(stats_out$Philippines_window_0_20m_cells < 20)) {
  warning("Some Philippines panels have fewer than 20 valid 0-20 m cells. ",
          "Bio-ORACLE's coarse mean bathymetry does not resolve reef-scale habitat.")
}

# -----------------------------------------------------------------------------
# 9. PUBLICATION-STYLE A/B/C FIGURES: TWO REGIONS, TWO MAP TYPES
# -----------------------------------------------------------------------------
message("[9/10] Building overall + Philippines A/B/C figure sets")
world <- rnaturalearth::ne_countries(scale="medium",returnclass="sf")
world_overall <- suppressWarnings(sf::st_crop(world,
  xmin=terra::xmin(MODEL_EXTENT), xmax=terra::xmax(MODEL_EXTENT),
  ymin=terra::ymin(MODEL_EXTENT), ymax=terra::ymax(MODEL_EXTENT)))
world_ph <- suppressWarnings(sf::st_crop(world,
  xmin=terra::xmin(PH_EXTENT), xmax=terra::xmax(PH_EXTENT),
  ymin=terra::ymin(PH_EXTENT), ymax=terra::ymax(PH_EXTENT)))

# LGM land is MODERN land plus shelf exposed by lower glacial sea level.
# The exposed shelf is approximated by modern MARSPEC bathymetry between -120/0 m.
# This is not a fully reconstructed paleocoastline; the bathymetry used for
# LGM suitability is the period-specific paleo-MARSPEC bathymetry above.
lgm_exposed_shelf <- terra::ifel(
  !is.na(ms_bathy_present) & ms_bathy_present >= -120 &
  ms_bathy_present < 0, 1, NA)

bathy_palette <- grDevices::gray.colors(140,start=0.22,end=0.96)
suitability_palette <- grDevices::colorRampPalette(c(
  "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272",
  "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d"))(120)
red_overlay <- grDevices::adjustcolor(suitability_palette, alpha.f=0.84)
binary_red <- grDevices::adjustcolor("#d7191c", alpha.f=0.87)
land_fill <- "#ebebeb"
coast_border <- "#686868"

pct_text <- function(id, region) {
  index <- match(id, names(periods))
  if (id == "B") return("Reference: 0.00%")
  col <- if (region == "overall")
    "change_vs_present_full_percent" else
    "change_vs_present_PH_0_20m_percent"
  sprintf("%+.2f%% vs present", stats_out[[col]][index])
}

# Choose geographic ticks that work well in the two plot windows.
lon_ticks <- function(region) if (region == "overall") {
  seq(60, 180, by = 20)
} else {
  c(116, 118, 120, 122, 124, 126)
}
lat_ticks <- function(region) if (region == "overall") {
  seq(-40, 40, by = 20)
} else {
  c(5, 10, 15, 20)
}

lon_tick_labels <- function(values) {
  parse(text = paste0(values, "*degree*E"))
}
lat_tick_labels <- function(values) {
  parse(text = paste0(abs(values), "*degree*", ifelse(values < 0, "S", "N")))
}

# Region-specific layout targets used to avoid stretched maps while also
# keeping the raster image well matched to the frame.
plot_specs <- list(
  overall = list(
    triptych_width = 11200,
    triptych_height = 3200,
    single_width = 3900,
    single_height = 2600,
    layout_widths = c(1, 1, 1, 0.20),
    mar_triptych = c(4.5, 4.9, 4.4, 0.65),
    mar_single = c(4.6, 5.0, 4.5, 0.9),
    legend_mar = c(3.4, 0.6, 4.1, 0.8),
    contour_levels = c(-6000, -4000, -2000, -1000, -500, -200, -20),
    contour_lwd = 0.58,
    point_cex = 0.34,
    coast_lwd = 0.45
  ),
  philippines = list(
    triptych_width = 9000,
    triptych_height = 4800,
    single_width = 3000,
    single_height = 4200,
    layout_widths = c(1, 1, 1, 0.18),
    mar_triptych = c(4.6, 5.0, 4.5, 0.65),
    mar_single = c(4.7, 5.0, 4.6, 0.9),
    legend_mar = c(3.4, 0.6, 4.1, 0.8),
    contour_levels = c(-4000, -2000, -1000, -500, -200, -100, -50, -20),
    contour_lwd = 0.72,
    point_cex = 0.65,
    coast_lwd = 0.65
  )
)

panel_map <- function(id, region = c("overall", "philippines"),
                      mode = c("continuous", "binary")) {
  region <- match.arg(region)
  mode <- match.arg(mode)
  spec <- plot_specs[[region]]
  z <- periods[[id]]
  is_ph <- identical(region, "philippines")
  bath <- if (is_ph) z$ph_bathy else z$bathy
  suitability <- if (is_ph) z$ph_suitability else z$suitability
  presence_only <- if (is_ph) z$ph_only_presence else z$overall_only_presence
  ex <- if (is_ph) PH_EXTENT else MODEL_EXTENT
  land <- if (is_ph) world_ph else world_overall

  terra::plot(
    bath,
    col = bathy_palette,
    zlim = c(-6000, 0),
    legend = FALSE,
    axes = FALSE,
    box = FALSE,
    main = "",
    maxcell = 350000,
    xlim = c(terra::xmin(ex), terra::xmax(ex)),
    ylim = c(terra::ymin(ex), terra::ymax(ex))
  )

  # Overlay only marine predictions (and only 0-20 m in Philippines).
  if (mode == "continuous") {
    if (count_valid_cells(suitability) > 0) {
      terra::plot(
        suitability,
        add = TRUE,
        col = red_overlay,
        zlim = c(0, 1),
        legend = FALSE,
        axes = FALSE,
        maxcell = 350000
      )
    }
  } else if (count_valid_cells(presence_only) > 0) {
    terra::plot(
      presence_only,
      add = TRUE,
      col = binary_red,
      legend = FALSE,
      axes = FALSE,
      maxcell = 350000
    )
  }

  # Bathymetric contours are added AFTER red overlays for visibility.
  suppressWarnings(try(
    terra::contour(
      bath,
      levels = spec$contour_levels,
      add = TRUE,
      maxcells = 110000,
      drawlabels = FALSE,
      col = grDevices::adjustcolor("#525252", alpha.f = 0.68),
      lwd = spec$contour_lwd
    ),
    silent = TRUE
  ))

  # Modern land polygons in all three time periods.
  # Additional approximate shelf exposed during LGM only.
  if (id == "A") {
    shelf <- terra::crop(lgm_exposed_shelf, ex)
    terra::plot(
      shelf,
      add = TRUE,
      col = land_fill,
      legend = FALSE,
      axes = FALSE,
      maxcell = 250000
    )
  }
  plot(sf::st_geometry(land), add = TRUE, col = land_fill,
       border = coast_border, lwd = spec$coast_lwd)

  # Occurrences only on the present panel; black cross markers.
  if (id == "B") {
    inside <- obs$longitude >= terra::xmin(ex) &
      obs$longitude <= terra::xmax(ex) &
      obs$latitude >= terra::ymin(ex) &
      obs$latitude <= terra::ymax(ex)
    if (any(inside)) {
      graphics::points(
        obs$longitude[inside], obs$latitude[inside],
        pch = 3,
        cex = spec$point_cex,
        col = grDevices::adjustcolor("black", alpha.f = 0.65),
        lwd = 0.75
      )
    }
  }

  lonv <- lon_ticks(region)
  latv <- lat_ticks(region)
  graphics::axis(
    1,
    at = lonv,
    labels = lon_tick_labels(lonv),
    cex.axis = 1.03,
    lwd = 0,
    lwd.ticks = 1.1
  )
  graphics::axis(
    2,
    at = latv,
    labels = lat_tick_labels(latv),
    las = 1,
    cex.axis = 1.03,
    lwd = 0,
    lwd.ticks = 1.1
  )
  graphics::title(xlab = "Longitude", ylab = "Latitude",
                  cex.lab = 1.15, line = 2.55)
  graphics::title(main = z$label, cex.main = 1.24, font.main = 2, line = 2.2)
  graphics::mtext(
    pct_text(id, region),
    side = 3,
    line = 0.55,
    cex = 1.03,
    font = 2,
    col = if (id == "B") "#444444" else "#991b1b"
  )
  invisible(TRUE)
}

shared_legend <- function(mode) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  if (mode == "continuous") {
    graphics::text(0.02, 0.95, "Suitability", adj = 0, font = 2, cex = 1.08)
    graphics::text(0.02, 0.91, "(0-1)", adj = 0, cex = 0.95)
    y0 <- 0.13
    y1 <- 0.84
    for (i in seq_along(suitability_palette)) {
      graphics::rect(
        0.04, y0 + (i - 1) * (y1 - y0) / length(suitability_palette),
        0.42, y0 + i * (y1 - y0) / length(suitability_palette),
        col = suitability_palette[i], border = NA
      )
    }
    graphics::rect(0.04, y0, 0.42, y1, border = "grey45", lwd = 0.6)
    graphics::text(
      0.53, seq(y0, y1, length.out = 6),
      labels = sprintf("%.1f", seq(0, 1, length.out = 6)),
      adj = 0, cex = 0.92
    )
  } else {
    graphics::text(0.02, 0.91, "Binary", adj = 0, font = 2, cex = 1.08)
    graphics::text(0.02, 0.86, "prediction", adj = 0, font = 2, cex = 1.00)
    graphics::rect(0.05, 0.71, 0.39, 0.76, col = "#d7191c", border = NA)
    graphics::text(0.48, 0.735, "Presence", adj = 0, cex = 0.90)
    graphics::text(0.03, 0.58, "Red = suitable", adj = 0, cex = 0.88)
    graphics::text(0.03, 0.52, "Gray = land", adj = 0, cex = 0.88)
  }
  graphics::text(0.02, 0.29, "Gray contours", adj = 0, cex = 0.86)
  graphics::text(0.02, 0.25, "= depth (m)", adj = 0, cex = 0.86)
  invisible(TRUE)
}

save_triptych <- function(destination, region, mode) {
  spec <- plot_specs[[region]]
  is_ph <- identical(region, "philippines")
  grDevices::png(
    destination,
    width = spec$triptych_width,
    height = spec$triptych_height,
    res = 320,
    bg = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(1:4, nrow = 1), widths = spec$layout_widths)
  graphics::par(
    mar = spec$mar_triptych,
    oma = c(3.5, 0.8, 2.7, 0.2),
    bty = "n",
    tcl = -0.33,
    mgp = c(2.45, 0.8, 0),
    xaxs = "i",
    yaxs = "i"
  )
  for (id in names(periods)) panel_map(id, region, mode)
  graphics::par(mar = spec$legend_mar, xaxs = "i", yaxs = "i")
  shared_legend(mode)
  main_title <- paste0(
    SPECIES, " | ",
    if (is_ph) "Philippines map window (0-20 m depth)" else "Overall modelled range (all marine depths)",
    " | ", if (mode == "continuous") "Continuous suitability" else "Binary predicted distribution"
  )
  graphics::mtext(main_title, outer = TRUE, side = 3, line = 0.45, font = 2, cex = 1.30)
  note <- paste0(
    "Percent = (period area - present area) / present area x 100; ",
    "A: paleo-MARSPEC anomaly, B/C: Bio-ORACLE v3; ",
    "C: SSP5-8.5; assumed sea-level rise +",
    if (ADJUST_FUTURE_WATER_DEPTH) SEA_LEVEL_RISE_M else 0,
    " m (no coastal inundation model)"
  )
  graphics::mtext(note, outer = TRUE, side = 1, line = 0.8, cex = 0.81)
  invisible(destination)
}

save_single_map <- function(destination, id, region, mode) {
  spec <- plot_specs[[region]]
  is_ph <- identical(region, "philippines")
  grDevices::png(
    destination,
    width = spec$single_width,
    height = spec$single_height,
    res = 320,
    bg = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(1:2, nrow = 1), widths = c(1, 0.18))
  graphics::par(
    mar = spec$mar_single,
    oma = c(3.3, 0.7, 2.5, 0.2),
    bty = "o",
    tcl = -0.33,
    mgp = c(2.45, 0.8, 0),
    xaxs = "i",
    yaxs = "i"
  )
  panel_map(id, region, mode)
  graphics::par(mar = spec$legend_mar, xaxs = "i", yaxs = "i")
  shared_legend(mode)
  period_title <- paste0(periods[[id]]$label, " | ",
    if (is_ph) "Philippines map window (0-20 m depth)" else "Overall modelled range (all marine depths)")
  sub_title <- if (mode == "continuous") "Continuous suitability" else "Binary predicted distribution"
  graphics::mtext(paste0(SPECIES, " | ", period_title, " | ", sub_title),
    outer = TRUE, side = 3, line = 0.45, font = 2, cex = 1.22)
  note <- paste0(
    "Percent = (period area - present area) / present area x 100; ",
    pct_text(id, region), "; ",
    "C scenario = SSP5-8.5"
  )
  graphics::mtext(note, outer = TRUE, side = 1, line = 0.75, cex = 0.79)
  invisible(destination)
}

files_png <- c(
  overall_continuous = file.path(OUTPUT_DIR,
    "Overall_range_A_LGM_B_Present_C_Future_suitability.png"),
  overall_binary = file.path(OUTPUT_DIR,
    "Overall_range_A_LGM_B_Present_C_Future_binary.png"),
  philippines_continuous = file.path(OUTPUT_DIR,
    "Philippines_range_A_LGM_B_Present_C_Future_0_20m_suitability.png"),
  philippines_binary = file.path(OUTPUT_DIR,
    "Philippines_range_A_LGM_B_Present_C_Future_0_20m_binary.png")
)
save_triptych(files_png[["overall_continuous"]], "overall", "continuous")
save_triptych(files_png[["overall_binary"]], "overall", "binary")
save_triptych(files_png[["philippines_continuous"]], "philippines", "continuous")
save_triptych(files_png[["philippines_binary"]], "philippines", "binary")

single_png <- c(
  A_overall_continuous = file.path(OUTPUT_DIR,
    "A_LGM_overall_suitability.png"),
  B_overall_continuous = file.path(OUTPUT_DIR,
    "B_Present_overall_suitability.png"),
  C_overall_continuous = file.path(OUTPUT_DIR,
    "C_Future_overall_suitability.png"),
  A_overall_binary = file.path(OUTPUT_DIR,
    "A_LGM_overall_binary.png"),
  B_overall_binary = file.path(OUTPUT_DIR,
    "B_Present_overall_binary.png"),
  C_overall_binary = file.path(OUTPUT_DIR,
    "C_Future_overall_binary.png"),
  A_philippines_continuous = file.path(OUTPUT_DIR,
    "A_LGM_Philippines_0_20m_suitability.png"),
  B_philippines_continuous = file.path(OUTPUT_DIR,
    "B_Present_Philippines_0_20m_suitability.png"),
  C_philippines_continuous = file.path(OUTPUT_DIR,
    "C_Future_Philippines_0_20m_suitability.png"),
  A_philippines_binary = file.path(OUTPUT_DIR,
    "A_LGM_Philippines_0_20m_binary.png"),
  B_philippines_binary = file.path(OUTPUT_DIR,
    "B_Present_Philippines_0_20m_binary.png"),
  C_philippines_binary = file.path(OUTPUT_DIR,
    "C_Future_Philippines_0_20m_binary.png")
)
for (id in names(periods)) {
  save_single_map(single_png[[paste0(id, "_overall_continuous")]], id, "overall", "continuous")
  save_single_map(single_png[[paste0(id, "_overall_binary")]], id, "overall", "binary")
  save_single_map(single_png[[paste0(id, "_philippines_continuous")]], id, "philippines", "continuous")
  save_single_map(single_png[[paste0(id, "_philippines_binary")]], id, "philippines", "binary")
}

# -----------------------------------------------------------------------------
# 10. PROVENANCE AND FINAL CHECK
# -----------------------------------------------------------------------------
notes <- c(
  paste("Species:",SPECIES),
  "Script revision: v6: fixed Windows degree-symbol tick labels, reduced map/frame mismatch, improved output dimensions, and added individual period maps",
  "Predictors: sea surface temperature (C), sea surface salinity, signed bathymetry (m)",
  "Training: Bio-ORACLE v3 present 2010-2020; random background + binomial GLM",
  "LGM: MARSPEC LGM SST/SSS anomalies + Bio-ORACLE present baseline; LGM bathymetry",
  "Future: Bio-ORACLE v3 SSP5-8.5 2090-2100 SST/SSS",
  paste("Future assumed global SLR (m):",if(ADJUST_FUTURE_WATER_DEPTH) SEA_LEVEL_RISE_M else 0),
  "Fixed physical seabed, extra future water depth only, no newly inundated cells",
  "Area calculations: terra::cellSize(...,mask=TRUE) applied to presence-only mask",
  "Output: unique timestamped run directory avoids Windows file locks on old CSVs",
  "GLM glm.fit fitted probabilities ~0/1 warning may indicate separation, review coefficients and spatial validation",
  "Overall range: full MODEL_EXTENT; ALL modelled marine water depths",
  "Philippines: rectangular map window PH_EXTENT; shallow 0-20 m depth mask",
  "The Philippines window is not the Philippine exclusive economic zone",
  "Area percent: (LGM - present)/present*100, (future - present)/present*100",
  "Positive area change indicates expansion of MODELLED suitability, not actual colonization",
  "LGM exposed shelf is approximated from modern MARSPEC -120 to 0 m contours",
  "Modern coastline shown as land polygon for all scenarios; future inundation excluded",
  "Coarse Bio-ORACLE mean-bathymetry grid is NOT reef-scale shallow habitat resolution",
  "AUC from random holdout, not independent spatial cross-validation",
  "Cross-product delta-method paleo hindcast is exploratory; validate model transferability",
  sprintf("Random holdout AUC: %.4f",auc_test),
  sprintf("Training-derived Youden threshold: %.4f",threshold),
  paste("R version:",R.version.string),
  paste("biooracler:",as.character(utils::packageVersion("biooracler"))),
  paste("terra:",as.character(utils::packageVersion("terra"))),
  "Bio-ORACLE: https://www.bio-oracle.org/",
  "MARSPEC: https://lifewatch.github.io/sdmpredictors/",
  "IPCC AR6: https://www.ipcc.ch/report/ar6/wg1/"
)
writeLines(notes,log_path)
required_files <- c(files_png, single_png,
 file.path(OUTPUT_DIR,"A_B_C_suitable_area_summary.csv"),
 file.path(OUTPUT_DIR,"Present_vs_LGM_and_Future_percent_comparisons.csv"))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)>0) {
  stop("Missing output file(s): ",paste(missing_files,collapse="; "))
}
message("[10/10] Combined and individual maps plus both area CSVs verified")
message("\n==============================================")
message("COMPLETED: overall and Philippines combined + individual outputs (v6 plotting fix)")
message(sprintf("Holdout AUC = %.3f; threshold = %.3f",auc_test,threshold))
message("Overall and Philippines comparisons saved to: ", OUTPUT_DIR)
for (f in files_png) message("Combined figure: ", f)
for (f in single_png) message("Individual figure: ", f)
message("==============================================")
} # End run_acanthochromis_sdm

# Run as a SINGLE source() command in RStudio's Console.
run_acanthochromis_sdm()
###############################################################################
# END of script v6
###############################################################################
