## Build cells.tsv for argelaguet from sample_metadata.txt + the local cell dir.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--metadata",   type = "character"),
    make_option("--cells_dir",  type = "character"),
    make_option("--group_col",  type = "character", default = "lineage10x"),
    make_option("--out",        type = "character")
)))

meta <- fread(opt$metadata, sep = "\t", header = TRUE)
meta <- meta[pass_metQC == TRUE & !is.na(id_met)]

files <- list.files(opt$cells_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
ids   <- sub("\\.tsv\\.gz$", "", basename(files))
have  <- data.table(cell_id = ids, path = normalizePath(files))

merged <- merge(have, meta, by.x = "cell_id", by.y = "id_met")

out <- merged[, .(
    cell_id,
    group  = get(opt$group_col),
    path,
    format = "scnmt",
    stage,
    embryo
)]
out <- out[!is.na(group)]

fwrite(out, opt$out, sep = "\t")
message(sprintf("[manifest] wrote %d cells across %d groups",
                nrow(out), uniqueN(out$group)))
