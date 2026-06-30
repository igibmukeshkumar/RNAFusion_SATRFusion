#!/usr/bin/env python3
"""
filter_aml_fusions.py

Filters a merged STAR-Fusion master sheet (800 AML samples, no FusionInspector)
through 6 sequential stages: read evidence, large anchor, isoform consolidation,
blacklist removal, cohort recurrence cap (with AML driver whitelist), and
coding-effect / in-frame flagging.

Input columns expected (tab-separated, matches real STAR-Fusion abridged output):
    #FusionName, JunctionReadCount, SpanningFragCount, est_J, est_S, SpliceType,
    LeftGene, LeftBreakpoint, RightGene, RightBreakpoint, LargeAnchorSupport,
    FFPM, LeftBreakDinuc, LeftBreakEntropy, RightBreakDinuc, RightBreakEntropy,
    annots, path, sample_name

Usage:
    python filter_aml_fusions.py \
        --input combined_star-fusion.fusion_predictions.abridged.tsv \
        --output filtered_800_aml_fusions.tsv \
        --babiceanu babiceanu_normal_fusions.txt \
        --gtex gtex_fusion_panel.txt
"""

import argparse
import sys
import pandas as pd


# ── known AML driver fusions — never discarded by recurrence cap ──────────────
AML_WHITELIST = {
    "PML--RARA", "RUNX1--RUNX1T1", "CBFB--MYH11", "DEK--NUP214",
    "KMT2A--MLLT3", "KMT2A--AFF1", "KMT2A--MLLT1", "KMT2A--ELL",
    "KMT2A--MLLT10", "KMT2A--MLLT4", "NUP98--NSD1", "NUP98--KDM5A",
    "NUP98--HOXA9", "BCR--ABL1", "ETV6--RUNX1", "MYH11--CBFB",
    "RUNX1T1--RUNX1", "RARA--PML", "NPM1--MLF1", "ETV6--ABL1",
}


def load_blacklist(path):
    """Load a single-column blacklist file of fusion names. Returns empty set if path is None."""
    if not path:
        return set()
    try:
        bl = pd.read_csv(path, header=None, names=["fusion_name"])
        return set(bl["fusion_name"].astype(str).str.strip())
    except FileNotFoundError:
        print(f"WARNING: blacklist file not found: {path} — skipping", file=sys.stderr)
        return set()


def main():
    ap = argparse.ArgumentParser(description="Filter merged STAR-Fusion AML cohort TSV")
    ap.add_argument("--input",     required=True, help="merged combined TSV")
    ap.add_argument("--output",    required=True, help="filtered output TSV")
    ap.add_argument("--babiceanu", default=None,  help="Babiceanu normal-tissue fusion panel (1 col, no header)")
    ap.add_argument("--gtex",      default=None,  help="GTEx recurrent fusion panel (1 col, no header)")
    ap.add_argument("--min_junction_reads",  type=int,   default=3)
    ap.add_argument("--min_spanning_frags",  type=int,   default=5)
    ap.add_argument("--min_ffpm",            type=float, default=0.1)
    ap.add_argument("--isoform_frac",        type=float, default=0.10,
                     help="keep isoforms with >= this fraction of the dominant isoform's JunctionReadCount")
    ap.add_argument("--recurrence_cap",      type=float, default=0.10,
                     help="discard fusions seen in more than this fraction of samples (unless whitelisted)")
    args = ap.parse_args()

    # ── load ────────────────────────────────────────────────────────────────
    df = pd.read_csv(args.input, sep="\t")

    # normalise the leading '#' on the first column name if present
    df.columns = [c.lstrip("#") for c in df.columns]

    n_start = len(df)
    n_samples_start = df["sample_name"].nunique()
    print(f"Loaded {n_start} fusion calls across {n_samples_start} samples")

    # ── Stage 1 — read evidence floor ──────────────────────────────────────
    df = df[
        (df["JunctionReadCount"] >= args.min_junction_reads) &
        (df["SpanningFragCount"] >= args.min_spanning_frags) &
        (df["FFPM"] >= args.min_ffpm)
    ]
    print(f"After Stage 1 (read evidence): {len(df)} calls")

    # ── Stage 2 — large anchor mandatory ───────────────────────────────────
    # real values are YES_LDAS / NO_LDAS, not plain YES/NO
    df = df[df["LargeAnchorSupport"].astype(str).str.startswith("YES")]
    print(f"After Stage 2 (large anchor): {len(df)} calls")

    # ── Stage 3 — isoform consolidation ────────────────────────────────────
    # use transform to avoid pandas version differences in groupby().apply()
    # dropping the grouping columns (changed behaviour across pandas versions)
    max_support = df.groupby(["sample_name", "FusionName"])["JunctionReadCount"].transform("max")
    df = df[df["JunctionReadCount"] >= args.isoform_frac * max_support]
    print(f"After Stage 3 (isoform consolidation): {len(df)} calls")

    # ── Stage 4 — blacklist removal ────────────────────────────────────────
    blacklist = load_blacklist(args.babiceanu) | load_blacklist(args.gtex)
    if blacklist:
        before = len(df)
        df = df[~df["FusionName"].isin(blacklist)]
        print(f"After Stage 4 (blacklist, {len(blacklist)} entries loaded): "
              f"{len(df)} calls ({before - len(df)} removed)")
    else:
        print("Stage 4 (blacklist): skipped — no blacklist files provided")

    # ── Stage 5 — cohort recurrence cap, with AML driver whitelist ────────
    n_samples = df["sample_name"].nunique()
    recurrence = df.groupby("FusionName")["sample_name"].nunique()
    over_cap = recurrence[recurrence / n_samples > args.recurrence_cap].index
    artefacts = [f for f in over_cap if f not in AML_WHITELIST]

    if artefacts:
        print(f"Stage 5: {len(artefacts)} fusion(s) exceed "
              f"{args.recurrence_cap*100:.0f}% recurrence cap and are NOT whitelisted — removing:")
        for f in sorted(artefacts):
            pct = 100 * recurrence[f] / n_samples
            print(f"    {f}  ({recurrence[f]}/{n_samples} samples, {pct:.1f}%)")
    df = df[~df["FusionName"].isin(artefacts)]
    print(f"After Stage 5 (cohort recurrence): {len(df)} calls")

    # ── Stage 6 — coding effect flag (informational, not a filter) ────────
    # PROT_FUSION_TYPE is only present when --examine_coding_effect was used;
    # fall back gracefully if the column is absent in this run.
    if "PROT_FUSION_TYPE" in df.columns:
        df["is_in_frame"] = df["PROT_FUSION_TYPE"].astype(str).str.contains("INFRAME", na=False)
    else:
        df["is_in_frame"] = pd.NA
        print("NOTE: PROT_FUSION_TYPE column not found — "
              "in-frame flag set to NA (re-run STAR-Fusion with --examine_coding_effect to populate this)")

    # flag whitelisted AML drivers explicitly for easy downstream sorting
    df["is_aml_driver"] = df["FusionName"].isin(AML_WHITELIST)

    df = df.sort_values(
        ["is_aml_driver", "is_in_frame", "JunctionReadCount"],
        ascending=[False, False, False]
    )

    # ── write output ────────────────────────────────────────────────────
    df.to_csv(args.output, sep="\t", index=False)

    n_final_samples = df["sample_name"].nunique()
    print("\n" + "=" * 60)
    print(f"FINAL: {len(df)} fusion calls across {n_final_samples} samples")
    print(f"  (started with {n_start} calls across {n_samples_start} samples)")
    print(f"  AML driver fusions retained: {df['is_aml_driver'].sum()}")
    print(f"  Output written to: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
