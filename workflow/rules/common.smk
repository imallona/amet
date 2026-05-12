## Shared rules for the dataset analyses: amet binary build and chr19 FASTA
## fetches. amet derives the CpG reference internally from the FASTA via its
## --genome flag, so the workflow does not pre-build a CpG reference.

from glob import glob

REFS = op.join(RESULTS, "refs")

## Helper scripts every analytical Rmd sources via source() or via the
## AMET_RENDER_HELPERS env var. Declared as snakemake inputs so script edits
## invalidate stale HTMLs. Dataset-specific rules concatenate this with the
## per-Rmd extras (driver_utils.R, diff_testing.R, embedding_utils.R).
SCRIPTS_DIR = op.join(REPO_ROOT, "workflow", "scripts")
RMD_SHARED_SCRIPTS = [
    op.join(SCRIPTS_DIR, "render_logging.R"),
    op.join(SCRIPTS_DIR, "plot_theme.R"),
    op.join(SCRIPTS_DIR, "palettes.R"),
]
DRIVER_UTILS_R = op.join(SCRIPTS_DIR, "driver_utils.R")
EMBEDDING_UTILS_R = op.join(SCRIPTS_DIR, "embedding_utils.R")
DIFF_TESTING_R = op.join(SCRIPTS_DIR, "diff_testing.R")

METHOD = op.join(REPO_ROOT, "method")
## Cargo.lock is gitignored (binary build artifact), so it's not listed as
## an input. cargo regenerates it from Cargo.toml on each build.
AMET_SOURCES = (
    glob(op.join(METHOD, "src", "**", "*.rs"), recursive=True)
    + glob(op.join(METHOD, "tests", "*.rs"))
    + [op.join(METHOD, "Cargo.toml")]
)

## chr19 FASTA sources keyed by (assembly, chrom-naming convention). Picked
## per dataset so that FASTA chrom names match the cell data.
CHR19 = {
    "mm10_ucsc": "https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chr19.fa.gz",
    "mm10_ensembl": "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz",
    "hg19_ucsc": "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr19.fa.gz",
}

## Whole-genome FASTA sources. Used when the prototype chr19 restriction is
## lifted, so amet derives a genome-wide CpG reference. About 800 MB
## compressed each.
WHOLE_GENOME_FASTA = {
    "mm10_ucsc":    "https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz",
    "mm10_ensembl": "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz",
    "hg19_ucsc":    "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz",
}


rule build_amet:
    """Build the amet release binary. cargo handles incremental compilation."""
    conda:
        op.join("..", "envs", "rust.yml")
    input:
        sources = AMET_SOURCES,
    output:
        binary = AMET,
    log:
        op.join(RESULTS, "logs", "build_amet.log"),
    shell:
        """
        mkdir -p $(dirname {log})
        cd {METHOD} && cargo build --release &> {log}
        """


rule fetch_chr19_fasta:
    """Pull a single chr19 FASTA. Source decides chrom naming (UCSC = chrN, Ensembl = N)."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    wildcard_constraints:
        source = "|".join(CHR19),
    output:
        fa = op.join(REFS, "{source}", "chr19.fa"),
    params:
        url = lambda w: CHR19[w.source],
    log:
        op.join(REFS, "{source}", "logs", "fetch_chr19.log"),
    shell:
        """
        mkdir -p $(dirname {output.fa})
        echo "[fetch_chr19_fasta] source={wildcards.source} url={params.url}" > {log}
        curl -sSL {params.url} 2>> {log} | gunzip -c > {output.fa}
        echo "[fetch_chr19_fasta] wrote $(wc -l < {output.fa}) lines, $(wc -c < {output.fa}) bytes to {output.fa}" >> {log}
        """


rule fetch_whole_genome_fasta:
    """Pull a whole-genome FASTA. Used when the chr19 restriction is lifted.
    Heavy: ~800 MB compressed, ~2.5 GB uncompressed. amet derives the CpG
    reference from this on first use (cached as <fasta>.cpg)."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    wildcard_constraints:
        source = "|".join(WHOLE_GENOME_FASTA),
    output:
        fa = op.join(REFS, "{source}", "genome.fa"),
    params:
        url = lambda w: WHOLE_GENOME_FASTA[w.source],
    log:
        op.join(REFS, "{source}", "logs", "fetch_genome.log"),
    shell:
        """
        mkdir -p $(dirname {output.fa})
        echo "[fetch_whole_genome] source={wildcards.source} url={params.url}" > {log}
        curl -sSL {params.url} 2>> {log} | gunzip -c > {output.fa}
        echo "[fetch_whole_genome] wrote $(wc -c < {output.fa}) bytes to {output.fa}" >> {log}
        """


rule build_cpg_reference:
    """Pre-build <fasta>.cpg so scoring jobs don't race on first use."""
    conda:
        op.join("..", "envs", "rust.yml")
    input:
        fa = op.join(REFS, "{source}", "genome.fa"),
        amet = AMET,
    output:
        cpg = op.join(REFS, "{source}", "genome.fa.cpg"),
    log:
        op.join(REFS, "{source}", "logs", "build_cpg_reference.log"),
    shell:
        """
        {input.amet} --build-cpg-only --genome {input.fa} &> {log}
        """


rule whole_genome_sizes:
    """Chrom-sizes from a whole-genome FASTA, written as <chr>\\t<len>."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        fa = op.join(REFS, "{source}", "genome.fa"),
    output:
        sizes = op.join(REFS, "{source}", "genome.sizes"),
    shell:
        r"""
        awk 'BEGIN{{n=""; len=0}}
             /^>/{{if(n) print n"\t"len; n=substr($1,2); len=0; next}}
             {{len+=length($0)}}
             END{{if(n) print n"\t"len}}' {input.fa} > {output.sizes}
        """
