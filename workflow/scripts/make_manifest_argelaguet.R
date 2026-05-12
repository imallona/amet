## Build cells.tsv for argelaguet from sample_metadata.txt + the local cell dir.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--metadata",   type = "character"),
    make_option("--cells_dir",  type = "character"),
    make_option("--group_col",  type = "character", default = "lineage10x"),
    make_option("--proto_stages",   type = "character", default = ""),
    make_option("--proto_lineages", type = "character", default = ""),
    make_option("--prototype",  type = "character", default = "false"),
    make_option("--out",        type = "character")
)))

prototype <- tolower(opt$prototype) %in% c("true", "1", "yes")
stages_keep   <- if (nchar(opt$proto_stages))   strsplit(opt$proto_stages,   ",")[[1]] else character()
lineages_keep <- if (nchar(opt$proto_lineages)) strsplit(opt$proto_lineages, ",")[[1]] else character()

meta <- fread(opt$metadata, sep = "\t", header = TRUE)
meta <- meta[pass_metQC == TRUE & !is.na(id_met)]

if (prototype) {
    if (length(stages_keep))   meta <- meta[stage %in% stages_keep]
    if (length(lineages_keep)) meta <- meta[get(opt$group_col) %in% lineages_keep]
}

files <- list.files(opt$cells_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
ids   <- sub("\\.tsv\\.gz$", "", basename(files))
have  <- data.table(cell_id = ids, path = normalizePath(files))

merged <- merge(have, meta, by.x = "cell_id", by.y = "id_met")

extra_cols <- intersect(c("lineage10x", "lineage10x_2", "plate"),
                        colnames(merged))

## Per-cell coverage proxy: cpg_level tsv.gz is one line per observed CpG,
## so size on disk is monotonic in cell coverage. Written so the per-combo
## subset can pick the top-N highest-coverage cells per (stage, lineage)
## plate-balanced without re-stat-ing the filesystem.
merged[, size := file.size(path)]

out <- merged[, .(
    cell_id,
    group = get(opt$group_col),
    path,
    format = "scnmt",
    stage,
    embryo,
    size
)]
for (col in extra_cols)
    out[[col]] <- merged[[col]]

out <- out[!is.na(group)]

fwrite(out, opt$out, sep = "\t")
message(sprintf("[manifest] wrote %d cells across %d groups",
                nrow(out), uniqueN(out$group)))
