## Shared rules for the dataset analyses: amet binary build and chr19 FASTA
## fetches. amet derives the CpG reference internally from the FASTA via its
## --genome flag, so the workflow does not pre-build a CpG reference.

from glob import glob

REFS = op.join(RESULTS, "refs")

METHOD = op.join(REPO_ROOT, "method")
AMET_SOURCES = (
    glob(op.join(METHOD, "src", "**", "*.rs"), recursive=True)
    + glob(op.join(METHOD, "tests", "*.rs"))
    + [op.join(METHOD, "Cargo.toml"), op.join(METHOD, "Cargo.lock")]
)

## chr19 FASTA sources keyed by (assembly, chrom-naming convention). Picked
## per dataset so that FASTA chrom names match the cell data.
CHR19 = {
    "mm10_ucsc": "https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chr19.fa.gz",
    "mm10_ensembl": "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz",
    "hg19_ucsc": "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr19.fa.gz",
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
