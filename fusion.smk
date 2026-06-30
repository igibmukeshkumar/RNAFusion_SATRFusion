import os, glob

configfile: "config.yaml"

# ── paths ─────────────────────────────────────────────────────────────────────
TRIMMED    = config["trimmed_dir"]
RESULT     = config["result_dir"]
CTAT       = config["ref_ctat"]
R1_SUF     = config["r1_suffix"]
R2_SUF     = config["r2_suffix"]
TMP        = config["tmp_dir"]

FUSION_DIR = os.path.join(RESULT, "fusion_result")
LOG_DIR    = os.path.join(RESULT, "logs")
QC_DIR     = os.path.join(RESULT, "post_trim_fastqc")
MQC_DIR    = os.path.join(RESULT, "multiqc")

# ── sample discovery ──────────────────────────────────────────────────────────
SAMPLES = [
    os.path.basename(f).replace(R1_SUF, "")
    for f in glob.glob(os.path.join(TRIMMED, f"*{R1_SUF}"))
]

# ── create dirs once at DAG-build time ───────────────────────────────────────
for d in [FUSION_DIR, LOG_DIR, QC_DIR, MQC_DIR, TMP]:
    os.makedirs(d, exist_ok=True)


# =============================================================================
rule all:
    input:
        expand(
            os.path.join(FUSION_DIR, "{s}", "star-fusion.fusion_predictions.tsv"),
            s=SAMPLES,
        ),
        os.path.join(MQC_DIR, "multiqc_report.html"),


# =============================================================================
rule fastqc:
    input:
        r1 = os.path.join(TRIMMED, "{s}" + R1_SUF),
        r2 = os.path.join(TRIMMED, "{s}" + R2_SUF),
    output:
        os.path.join(QC_DIR, "{s}_1_val_1_fastqc.zip"),
        os.path.join(QC_DIR, "{s}_2_val_2_fastqc.zip"),
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
        expand(os.path.join(QC_DIR, "{s}_1_val_1_fastqc.zip"), s=SAMPLES),
    output:
        os.path.join(MQC_DIR, "multiqc_report.html"),
    log:
        os.path.join(LOG_DIR, "multiqc.log")
    shell:
        """
        multiqc -o {MQC_DIR} {QC_DIR} > {log} 2>&1
        """


# =============================================================================
# star_fusion — serialised: threads=6 + star_fusion_slot=1
# guarantees exactly ONE job runs at a time regardless of --cores value
# =============================================================================
rule star_fusion:
    input:
        r1 = os.path.join(TRIMMED, "{s}" + R1_SUF),
        r2 = os.path.join(TRIMMED, "{s}" + R2_SUF),
    output:
        os.path.join(FUSION_DIR, "{s}", "star-fusion.fusion_predictions.tsv"),
    params:
        outdir     = os.path.join(FUSION_DIR, "{s}"),
        tmp_sample = os.path.join(TMP, "{s}"),
        ctat       = CTAT,
    threads: config["star_fusion_threads"]
    resources:
        star_fusion_slot = 1,
        mem_mb           = config.get("star_fusion_mem_mb", 40000),
    log:
        os.path.join(LOG_DIR, "star_fusion_{s}.log")
    shell:
        """
        set -euo pipefail

        mkdir -p {params.tmp_sample} {params.outdir}

        STAR-Fusion \
            --genome_lib_dir             {params.ctat}       \
            --left_fq                    {input.r1}          \
            --right_fq                   {input.r2}          \
            --CPU                        {threads}            \
            --output_dir                 {params.tmp_sample} \
            --min_junction_reads         2                   \
            --min_sum_frags              5                   \
            --min_FFPM                   0.1                 \
            --min_novel_junction_support 3                   \
            --examine_coding_effect                          \
            --FusionInspector            validate            \
            > {log} 2>&1

        rsync -rL {params.tmp_sample}/ {params.outdir}/
        rm -rf {params.tmp_sample}
        """
