## Build cells.tsv for CRC from the CRC_RAW directory.
##
## Filename pattern: GSM<id>_scTrioSeq2Met_<patient>_<location><lane>_<cell>.singleC.txt.gz
## Group column is 'location' (NC, PT, LN, ...). Patient is preserved.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--raw_dir",         type = "character"),
    make_option("--proto_patients",  type = "character", default = ""),
    make_option("--proto_locations", type = "character", default = ""),
    make_option("--cells_per_group", type = "integer",   default = 10L),
    make_option("--prototype",       type = "character", default = "true"),
    make_option("--out",             type = "character")
)))

prototype <- tolower(opt$prototype) %in% c("true", "1", "yes")
patients_keep  <- if (nchar(opt$proto_patients))  strsplit(opt$proto_patients,  ",")[[1]] else character()
locations_keep <- if (nchar(opt$proto_locations)) strsplit(opt$proto_locations, ",")[[1]] else character()

files <- list.files(opt$raw_dir, pattern = "\\.singleC\\.txt\\.gz$", full.names = TRUE)
if (!length(files)) stop("no singleC files in ", opt$raw_dir)

basenames <- basename(files)
m <- regmatches(basenames,
                regexec("^(GSM\\d+)_scTrioSeq2Met_(CRC\\d+)_([A-Z]+)\\d*_(\\d+)\\.singleC\\.txt\\.gz$", basenames))

valid <- vapply(m, length, integer(1)) == 5
if (!any(valid)) stop("no singleC filenames matched the scTrioSeq2Met pattern")
parts <- do.call(rbind, lapply(m[valid], function(x) x[2:5]))
colnames(parts) <- c("gsm", "patient", "location", "cell_idx")
dt <- data.table(parts)
dt[, path := normalizePath(files[valid])]
dt[, cell_id := paste(patient, location, gsm, cell_idx, sep = "_")]

if (prototype) {
    if (length(patients_keep))  dt <- dt[patient  %in% patients_keep]
    if (length(locations_keep)) dt <- dt[location %in% locations_keep]
    dt <- dt[order(patient, location)]
    dt <- dt[, head(.SD, opt$cells_per_group), by = .(patient, location)]
}

out <- dt[, .(cell_id, group = location, path, format = "bismark", patient, location)]
fwrite(out, opt$out, sep = "\t")
message(sprintf("[manifest] wrote %d cells across %d groups",
                nrow(out), uniqueN(out$group)))
