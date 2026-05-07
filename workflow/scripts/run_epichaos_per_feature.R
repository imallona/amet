## Run epiCHAOS once per genomic feature, returning a heterogeneity score per
## feature. This adapter is for the feature-variability sweep, where the natural
## comparison axis is per-feature dispersion (not per-group).
##
## Approach: for each feature in the BED, slice the cells x CpGs binary matrix
## to that feature's CpG range and call epiCHAOS::epiCHAOS with all cells
## assigned to a single dummy group. Output one (feature_id, eITH, eITH_adj) row
## per feature.

suppressPackageStartupMessages({
    library(optparse); library(data.table); library(epiCHAOS)
})

options <- list(
    make_option(c("--manifest"), type = "character"),
    make_option(c("--cpg_reference"), type = "character"),
    make_option(c("--bed"), type = "character"),
    make_option(c("--output"), type = "character"),
    make_option(c("--threads"), type = "integer", default = 4)
)
opt <- parse_args(OptionParser(option_list = options))

manifest <- fread(opt$manifest, header = TRUE, sep = "\t")
cpgs <- fread(opt$cpg_reference, header = FALSE, col.names = c("chrom", "pos"))
bed <- fread(opt$bed, header = FALSE,
             col.names = c("chrom", "start", "end", "name"))
n_cpg <- nrow(cpgs)
cell_keys <- paste(cpgs$chrom, cpgs$pos, sep = "_")
key_to_idx <- setNames(seq_len(n_cpg), cell_keys)

read_cell_bits <- function(path) {
    bits <- rep(NA_integer_, n_cpg)
    if (!file.exists(path)) return(bits)
    d <- fread(path, header = FALSE, sep = "\t",
               col.names = c("chrom", "pos1based", "strand", "context",
                             "mc", "cov", "methylated"))
    cpg_start_0based <- ifelse(d$strand == "-", d$pos1based - 2L, d$pos1based - 1L)
    keys <- paste(d$chrom, cpg_start_0based, sep = "_")
    idx <- key_to_idx[keys]
    valid <- !is.na(idx)
    bits[idx[valid]] <- as.integer(d$methylated[valid] > 0)
    bits
}

message(sprintf("[run_epichaos_per_feature] reading %d cells", nrow(manifest)))
mat_list <- parallel::mclapply(manifest$path, read_cell_bits, mc.cores = opt$threads)
mat <- do.call(cbind, mat_list)
colnames(mat) <- manifest$cell_id
rownames(mat) <- cell_keys
mat[is.na(mat)] <- 0L

per_feature <- data.frame(feature_id = character(), eITH = numeric(),
                          eITH_adj = numeric(), stringsAsFactors = FALSE)

## epiCHAOS rank-normalises het.raw across the groups in a single call, so a
## per-feature reading needs all features present as groups in one call. Each
## (cell, feature) becomes a pseudo-cell carrying that feature's CpG values;
## group label is the feature_id. The matrix has L rows (CpG positions within a
## feature) and n_cells * n_features columns, one column per pseudo-cell.
n_cells <- ncol(mat)
features <- bed$name
feat_ranges <- vector("list", nrow(bed))
for (i in seq_len(nrow(bed))) {
    feat_ranges[[i]] <- which(cpgs$chrom == bed$chrom[i] &
                              cpgs$pos >= bed$start[i] & cpgs$pos < bed$end[i])
}
feat_lens <- vapply(feat_ranges, length, integer(1))
if (length(unique(feat_lens)) != 1) {
    stop("per-feature stacking requires features of equal CpG count")
}
L <- feat_lens[1]

stack_cols <- vector("list", nrow(bed))
group_labels <- character()
pseudo_ids <- character()
for (i in seq_len(nrow(bed))) {
    sub_mat <- mat[feat_ranges[[i]], , drop = FALSE]
    stack_cols[[i]] <- sub_mat
    pseudo_ids <- c(pseudo_ids,
                    sprintf("%s_cell%03d", features[i], seq_len(n_cells)))
    group_labels <- c(group_labels, rep(features[i], n_cells))
}
combined <- do.call(cbind, stack_cols)
colnames(combined) <- pseudo_ids
rownames(combined) <- sprintf("pos%03d", seq_len(L))

stack_meta <- data.frame(cell_id = pseudo_ids, group = group_labels,
                         row.names = pseudo_ids, stringsAsFactors = FALSE)

message(sprintf("[run_epichaos_per_feature] stacked matrix: %d x %d, %d feature groups",
                nrow(combined), ncol(combined), length(features)))

res <- epiCHAOS::epiCHAOS(counts = combined, meta = stack_meta,
                          colname = "group", n = 100, plot = FALSE)
res$feature_id <- sub("^group-", "", res$state)
per_feature <- data.frame(feature_id = res$feature_id,
                          eITH = res$het.raw,
                          eITH_adj = res$het.adj,
                          stringsAsFactors = FALSE)

fwrite(per_feature, opt$output, sep = "\t")
message(sprintf("[run_epichaos_per_feature] wrote %s (%d features)",
                opt$output, nrow(per_feature)))
