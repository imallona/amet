## Run epiCHAOS per group on amet's simulated cells.
##
## epiCHAOS expects a binary matrix (features x cells) plus a metadata data.frame
## with a grouping column. We construct the matrix by re-reading each cell's allc
## file from the manifest and aligning calls to the CpG reference.

suppressPackageStartupMessages({
    library(optparse); library(data.table); library(epiCHAOS)
})

options <- list(
    make_option(c("--manifest"), type = "character"),
    make_option(c("--cpg_reference"), type = "character"),
    make_option(c("--output"), type = "character"),
    make_option(c("--threads"), type = "integer", default = 4)
)
opt <- parse_args(OptionParser(option_list = options))

manifest <- fread(opt$manifest, header = TRUE, sep = "\t")
cpgs <- fread(opt$cpg_reference, header = FALSE, col.names = c("chrom", "pos"))
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

message(sprintf("[run_epichaos] reading %d cells", nrow(manifest)))
mat_list <- parallel::mclapply(manifest$path, read_cell_bits, mc.cores = opt$threads)
mat <- do.call(cbind, mat_list)
colnames(mat) <- manifest$cell_id
rownames(mat) <- cell_keys

message(sprintf("[run_epichaos] matrix dim = %d x %d, NA fraction = %.3f",
                nrow(mat), ncol(mat), mean(is.na(mat))))

mat[is.na(mat)] <- 0L
meta <- data.frame(cell_id = manifest$cell_id, group = manifest$group,
                   row.names = manifest$cell_id, stringsAsFactors = FALSE)

message("[run_epichaos] computing eITH per group")
results <- epiCHAOS::epiCHAOS(counts = mat, meta = meta, colname = "group",
                              n = 100, plot = FALSE)
## epiCHAOS returns a data.frame with columns het.raw, het.adj, state.
## We expose both as eITH (raw) and eITH_adj (count-adjusted).
## epiCHAOS prefixes the state with the meta column name (e.g. "group-acvi01_p0.05");
## strip that prefix so the output joins cleanly against amet's `group` column.
out <- data.frame(group = sub("^group-", "", results$state),
                  eITH = results$het.raw,
                  eITH_adj = results$het.adj,
                  stringsAsFactors = FALSE)
fwrite(out, opt$output, sep = "\t")
message(sprintf("[run_epichaos] wrote %s (%d groups)", opt$output, nrow(out)))
