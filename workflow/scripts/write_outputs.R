## Helpers for writing synthetic per-cell methylation files in allc format and
## a corresponding CpG reference TSV.

write_allc <- function(path, chrom, pos_0based, bits, coverage = NULL) {
    if (is.null(coverage)) coverage <- rep(1L, length(bits))
    stopifnot(length(bits) == length(coverage))
    pos_1based <- pos_0based + 1L
    df <- data.frame(
        chrom = chrom,
        pos = pos_1based,
        strand = "+",
        context = "CGN",
        mc = bits,
        cov = coverage,
        methylated_flag = bits
    )
    con <- if (endsWith(path, ".gz")) gzfile(path, "w") else file(path, "w")
    on.exit(close(con))
    write.table(df, con, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
}

write_cpg_reference <- function(path, chrom, pos_0based) {
    df <- data.frame(chrom = chrom, pos = pos_0based)
    con <- if (endsWith(path, ".gz")) gzfile(path, "w") else file(path, "w")
    on.exit(close(con))
    write.table(df, con, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
}

write_bed_feature <- function(path, chrom, start_0based, end_0based, name) {
    df <- data.frame(chrom = chrom, start = start_0based, end = end_0based, name = name)
    write.table(df, path, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
}

write_manifest <- function(path, cell_ids, groups, paths, extras = NULL) {
    df <- data.frame(cell_id = cell_ids, group = groups, path = paths)
    if (!is.null(extras)) df <- cbind(df, extras)
    write.table(df, path, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
}
