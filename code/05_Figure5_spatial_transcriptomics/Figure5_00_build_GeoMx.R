#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GeomxTools)
  library(NanoStringNCTools)
  library(Biobase)
})

args <- commandArgs(trailingOnly = TRUE)

project_dir <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

raw_dir <- file.path(
  project_dir,
  "data",
  "05_Figure5_spatial_transcriptomics",
  "raw"
)

dcc_dir <- file.path(raw_dir, "dcc")

pkc_file <- file.path(
  raw_dir,
  "Mm_R_NGS_WTA_v1.0.pkc"
)

annotation_file <- file.path(
  raw_dir,
  "HSV_Mock_Mouse_annotations.xlsx"
)

output_dir <- file.path(
  project_dir,
  "data",
  "05_Figure5_spatial_transcriptomics",
  "processed"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(
  output_dir,
  "GeoMx_WTA_processed_rebuilt.rds"
)

reference_file <- "data/05_Figure5_spatial_transcriptomics/processed/GeoMx_WTA_processed_manuscript.rds"

dcc_files <- list.files(
  dcc_dir,
  pattern = "\\.dcc$",
  full.names = TRUE
)

stopifnot(length(dcc_files) > 0)
stopifnot(file.exists(pkc_file))
stopifnot(file.exists(annotation_file))

demoData <- readNanoStringGeoMxSet(
  dccFiles = dcc_files,
  pkcFiles = pkc_file,
  phenoDataFile = annotation_file,
  phenoDataSheet = NULL,
  phenoDataDccColName = "Sample_ID",
  protocolDataColNames = c("Aoi", "Roi"),
  experimentDataColNames = c("Panel")
)

demoData <- shiftCountsOne(
  demoData,
  useDALogic = TRUE
)

modules <- gsub(".pkc", "", annotation(demoData))

QC_params <- list(
  minSegmentReads = 1000,
  percentTrimmed = 80,
  percentStitched = 80,
  percentAligned = 75,
  percentSaturation = 50,
  minNegativeCount = 1,
  maxNTCCount = 9000,
  minNuclei = 20,
  minArea = 1000
)

demoData <- setSegmentQCFlags(
  demoData,
  qcCutoffs = QC_params
)

QCResults <- protocolData(demoData)[["QCFlags"]]

QCResults$QCStatus <- apply(
  QCResults,
  1,
  function(x) ifelse(sum(x) == 0, "PASS", "WARNING")
)

demoData <- demoData[, QCResults$QCStatus == "PASS"]

negativeGeoMeans <- esBy(
  negativeControlSubset(demoData),
  GROUP = "Module",
  FUN = function(x) {
    assayDataApply(
      x,
      MARGIN = 2,
      FUN = ngeoMean,
      elt = "exprs"
    )
  }
)

protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

negCols <- paste0("NegGeoMean_", modules)

pData(demoData)[, negCols] <-
  sData(demoData)[["NegGeoMean"]]

demoData <- setBioProbeQCFlags(
  demoData,
  qcCutoffs = list(
    minProbeRatio = 0.1,
    percentFailGrubbs = 20
  ),
  removeLocalOutliers = TRUE
)

demoData <- subset(
  demoData,
  fData(demoData)[["QCFlags"]][, "LowProbeRatio"] == FALSE &
    fData(demoData)[["QCFlags"]][, "GlobalGrubbsOutlier"] == FALSE
)

target_demoData <- aggregateCounts(demoData)

cutoff <- 2
minLOQ <- 2

LOQ <- data.frame(
  row.names = colnames(target_demoData)
)

for (module in modules) {

  vars <- paste0(
    c("NegGeoMean_", "NegGeoSD_"),
    module
  )

  if (all(vars %in% colnames(pData(target_demoData)))) {

    LOQ[, module] <- pmax(
      minLOQ,
      pData(target_demoData)[, vars[1]] *
        pData(target_demoData)[, vars[2]]^cutoff
    )
  }
}

pData(target_demoData)$LOQ <- LOQ

LOQ_Mat <- NULL

for (module in modules) {

  idx <- fData(target_demoData)$Module == module

  mat <- t(
    esApply(
      target_demoData[idx, ],
      MARGIN = 1,
      FUN = function(x) x > LOQ[, module]
    )
  )

  LOQ_Mat <- rbind(LOQ_Mat, mat)
}

LOQ_Mat <- LOQ_Mat[
  fData(target_demoData)$TargetName,
]

pData(target_demoData)$GenesDetected <-
  colSums(
    LOQ_Mat,
    na.rm = TRUE
  )

pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected /
  nrow(target_demoData)

target_demoData <-
  target_demoData[
    ,
    pData(target_demoData)$GeneDetectionRate >= 0.025
  ]

LOQ_Mat <- LOQ_Mat[
  ,
  colnames(target_demoData)
]

fData(target_demoData)$DetectedSegments <-
  rowSums(
    LOQ_Mat,
    na.rm = TRUE
  )

fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments /
  nrow(pData(target_demoData))

negativeProbefData <- subset(
  fData(target_demoData),
  CodeClass == "Negative"
)

neg_probes <- unique(
  negativeProbefData$TargetName
)

target_demoData <-
  target_demoData[
    fData(target_demoData)$DetectionRate >= 0.025 |
      fData(target_demoData)$TargetName %in% neg_probes,
  ]

target_demoData <- normalize(
  target_demoData,
  norm_method = "quant",
  desiredQuantile = 0.75,
  toElt = "q_norm"
)

target_demoData <- normalize(
  target_demoData,
  norm_method = "neg",
  fromElt = "exprs",
  toElt = "neg_norm"
)

saveRDS(
  target_demoData,
  output_file
)

if (file.exists(reference_file)) {

  reference <- readRDS(reference_file)

  cat("\nValidation\n")
  cat("Rows:", nrow(target_demoData), "\n")
  cat("Cols:", ncol(target_demoData), "\n")
  cat("Reference Rows:", nrow(reference), "\n")
  cat("Reference Cols:", ncol(reference), "\n")
  cat("Matching genes:", identical(rownames(target_demoData), rownames(reference)), "\n")
  cat("Matching ROIs:", identical(colnames(target_demoData), colnames(reference)), "\n")
}

cat("\nSaved:", output_file, "\n")
