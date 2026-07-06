import os
import csv

configfile: "config.yaml"

# ── paths ─────────────────────────────────────────────────────────────────────
RESULT      = config["result_dir"]
CTAT        = config["ref_ctat"]
SAMPLESHEET = config["samplesheet"]

FUSION_DIR = os.path.join(RESULT, "fusion_result")
LOG_DIR    = os.path.join(RESULT, "logs")
QC_DIR     = os.path.join(RESULT, "post_trim_fastqc")
MQC_DIR    = os.path.join(RESULT, "multiqc")

# ── load sample sheet ─────────────────────────────────────────────────────────
# expected columns (tab-separated, with header): sampleID  pathR1  pathR2
SAMPLES = {}   # { sampleID: {"r1": path, "r2": path} }

with open(SAMPLESHEET) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        sid = row["sampleID"].strip()
        r1  = row["pathR1"].strip()
        r2  = row["pathR2"].strip()

        # validate files exist at parse time — fail fast before any jobs launch
        if not os.path.isfile(r1):
            raise FileNotFoundError(f"[samplesheet] R1 not found for {sid}: {r1}")
        if not os.path.isfile(r2):
            raise FileNotFoundError(f"[samplesheet] R2 not found for {sid}: {r2}")

        SAMPLES[sid] = {"r1": r1, "r2": r2}

print(f"Loaded {len(SAMPLES)} samples from {SAMPLESHEET}")

SAMPLE_IDS = list(SAMPLES.keys())

# ── create output dirs once at DAG-build time ─────────────────────────────────
for d in [FUSION_DIR, LOG_DIR, QC_DIR, MQC_DIR]:
    os.makedirs(d, exist_ok=True)


# ── helper: get R1/R2 path for a given sample ────────────────────────────────
def get_r1(wildcards): return SAMPLES[wildcards.s]["r1"]
def get_r2(wildcards): return SAMPLES[wildcards.s]["r2"]


# =============================================================================
rule all:
    input:
        expand(
            os.path.join(FUSION_DIR, "{s}", "star-fusion.fusion_predictions.tsv"),
            s=SAMPLE_IDS,
        ),
        os.path.join(MQC_DIR, "multiqc_report.html"),


# =============================================================================
rule fastqc:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        r1_zip = os.path.join(QC_DIR, "{s}_1_val_1_fastqc.zip"),
        r2_zip = os.path.join(QC_DIR, "{s}_2_val_2_fastqc.zip"),
    threads: config["fastqc_threads"]
    log:
        os.path.join(LOG_DIR, "fastqc_{s}.log")
    shell:
        """
        fastqc -o {QC_DIR} -t {threads} {input.r1} {input.r2} > {log} 2>&1
        """


# =============================================================================
rule multiqc:
    input:
        expand(os.path.join(QC_DIR, "{s}_1_val_1_fastqc.zip"), s=SAMPLE_IDS),
    output:
        os.path.join(MQC_DIR, "multiqc_report.html"),
    log:
        os.path.join(LOG_DIR, "multiqc.log")
    shell:
        """
        multiqc -o {MQC_DIR} {QC_DIR} > {log} 2>&1
        """


# =============================================================================
# star_fusion — run 2 in parallel (128 GB RAM, 64 cores machine)
# serialise with: --resources star_fusion_slot=2  on the CLI
# =============================================================================
rule star_fusion:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        os.path.join(FUSION_DIR, "{s}", "star-fusion.fusion_predictions.tsv"),
    params:
        ctat   = CTAT,
        outdir = os.path.join(FUSION_DIR, "{s}"),
    threads: config["star_fusion_threads"]
    resources:
        star_fusion_slot = 1
    log:
        os.path.join(LOG_DIR, "star_fusion_{s}.log")
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}

        STAR-Fusion \
            --genome_lib_dir             {params.ctat}       \
            --left_fq                    {input.r1}          \
            --right_fq                   {input.r2}          \
            --CPU                        {threads}           \
            --output_dir                 {params.outdir}     \
            --min_junction_reads         2                   \
            --min_sum_frags              5                   \
            --min_FFPM                   0.1                 \
            --min_novel_junction_support 3                   \
            --examine_coding_effect                          \
            --FusionInspector            validate            \
            > {log} 2>&1
        """
