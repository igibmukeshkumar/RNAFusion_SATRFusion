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
        --gtex gtex_fusion_panel.txt \
        --whitelist-file my_aml_drivers.txt

Whitelist file format (one entry per line, '#' comments allowed):
    PML--RARA          full fusion pair — matches this exact pair only
                        (either gene order: PML--RARA or RARA--PML both match)
    KMT2A               single gene — matches ANY fusion where KMT2A is
                        either partner (KMT2A--MLLT3, KMT2A--AFF1, etc.)
"""

import argparse
import sys
import pandas as pd


# ── built-in default AML driver whitelist — used only if --whitelist is not supplied ──
DEFAULT_AML_WHITELIST = {
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


def load_whitelist(path, enabled=True):
    """
    Load a whitelist of driver fusions and/or single genes (one entry per
    line, no header, '#' lines and blank lines ignored).

    Two entry formats are supported and may be mixed freely in the same file:
        GENE_A--GENE_B   exact fusion pair match (either gene order)
        GENE_A           single-gene match — protects ANY fusion where
                          GENE_A appears as either the left or right partner
                          (e.g. "KMT2A" matches KMT2A--MLLT3, KMT2A--AFF1, etc.)

    Falls back to the built-in DEFAULT_AML_WHITELIST when no path is supplied.
    Returns (pair_set, gene_set) — both empty when enabled=False.
    """
    if not enabled:
        print("Whitelist DISABLED (--no-whitelist) — no fusion is protected "
              "from the recurrence cap, including known AML drivers")
        return set(), set()

    if not path:
        print(f"No --whitelist file supplied — using built-in default "
              f"({len(DEFAULT_AML_WHITELIST)} known AML driver fusions)")
        raw_entries = set(DEFAULT_AML_WHITELIST)
    else:
        try:
            with open(path) as fh:
                raw_entries = {
                    line.strip() for line in fh
                    if line.strip() and not line.strip().startswith("#")
                }
            print(f"Loaded {len(raw_entries)} whitelist entr(y/ies) from: {path}")
        except FileNotFoundError:
            print(f"WARNING: whitelist file not found: {path} — "
                  f"falling back to built-in default whitelist", file=sys.stderr)
            raw_entries = set(DEFAULT_AML_WHITELIST)

    # split into full fusion-pair entries (contain "--") and single-gene entries
    pair_set = {e for e in raw_entries if "--" in e}
    gene_set = {e for e in raw_entries if "--" not in e}

    if gene_set:
        print(f"  → {len(pair_set)} full fusion-pair entries, "
              f"{len(gene_set)} single-gene entries (matches either partner): "
              f"{sorted(gene_set)}")

    return pair_set, gene_set


def fusion_matches_whitelist(fusion_name, pair_set, gene_set):
    """
    True if fusion_name (format 'GENE_A--GENE_B') matches the whitelist —
    either as an exact pair (checked both orientations) or because either
    individual partner gene is in gene_set.
    """
    if fusion_name in pair_set:
        return True
    parts = fusion_name.split("--")
    if len(parts) == 2:
        gene_a, gene_b = parts
        # exact pair, reversed orientation (B--A also counts as a match)
        if f"{gene_b}--{gene_a}" in pair_set:
            return True
        # single-gene match — either partner
        if gene_a in gene_set or gene_b in gene_set:
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description="Filter merged STAR-Fusion AML cohort TSV")
    ap.add_argument("--input",     required=True, help="merged combined TSV")
    ap.add_argument("--output",    required=True, help="filtered output TSV")
    ap.add_argument("--babiceanu", default=None,  help="Babiceanu normal-tissue fusion panel (1 col, no header)")
    ap.add_argument("--gtex",      default=None,  help="GTEx recurrent fusion panel (1 col, no header)")
    ap.add_argument("--whitelist-file", default=None, dest="whitelist_file",
                     help="custom driver fusion whitelist file, 1 fusion name per line "
                          "(e.g. PML--RARA). Overrides built-in default AML whitelist. "
                          "Ignored if whitelisting is disabled (see --no-whitelist).")
    ap.add_argument("--whitelist", dest="whitelist_enabled", action="store_true", default=True,
                     help="enable whitelist protection in Stage 5 (default: ON). "
                          "Uses --whitelist-file if supplied, otherwise the built-in "
                          "20-fusion AML driver list.")
    ap.add_argument("--no-whitelist", dest="whitelist_enabled", action="store_false",
                     help="disable whitelist protection entirely — Stage 5 recurrence cap "
                          "applies to ALL fusions including known AML drivers like PML--RARA. "
                          "Use this for an unbiased discovery run with no built-in assumptions.")
    ap.add_argument("--min_junction_reads",  type=int,   default=3)
    ap.add_argument("--min_spanning_frags",  type=int,   default=5)
    ap.add_argument("--min_ffpm",            type=float, default=0.1)
    ap.add_argument("--isoform_frac",        type=float, default=0.10,
                     help="keep isoforms with >= this fraction of the dominant isoform's JunctionReadCount")
    ap.add_argument("--recurrence_cap",      type=float, default=0.10,
                     help="discard fusions seen in more than this fraction of samples (unless whitelisted)")
    ap.add_argument("--recurrence_action", choices=["remove", "flag"], default="remove",
                     help="'remove' (default): drop non-whitelisted fusions exceeding the recurrence cap. "
                          "'flag': keep all rows, add an is_recurrent_artefact column "
                          "(True/False) instead of removing them.")
    args = ap.parse_args()

    wl_pairs, wl_genes = load_whitelist(args.whitelist_file, enabled=args.whitelist_enabled)

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
    artefacts = [f for f in over_cap if not fusion_matches_whitelist(f, wl_pairs, wl_genes)]

    if args.recurrence_action == "remove":
        if artefacts:
            print(f"Stage 5: {len(artefacts)} fusion(s) exceed "
                  f"{args.recurrence_cap*100:.0f}% recurrence cap and are NOT whitelisted — removing:")
            for f in sorted(artefacts):
                pct = 100 * recurrence[f] / n_samples
                print(f"    {f}  ({recurrence[f]}/{n_samples} samples, {pct:.1f}%)")
        df = df[~df["FusionName"].isin(artefacts)]
        print(f"After Stage 5 (cohort recurrence, removed): {len(df)} calls")
    else:
        # flag mode — keep every row, mark recurrence status instead of dropping
        df["is_recurrent_artefact"] = df["FusionName"].isin(artefacts)
        if artefacts:
            print(f"Stage 5: {len(artefacts)} fusion(s) exceed "
                  f"{args.recurrence_cap*100:.0f}% recurrence cap and are NOT whitelisted — "
                  f"FLAGGED as is_recurrent_artefact=True (rows retained):")
            for f in sorted(artefacts):
                pct = 100 * recurrence[f] / n_samples
                print(f"    {f}  ({recurrence[f]}/{n_samples} samples, {pct:.1f}%)")
        print(f"After Stage 5 (cohort recurrence, flagged not removed): {len(df)} calls "
              f"({df['is_recurrent_artefact'].sum()} flagged True)")

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
    df["is_aml_driver"] = df["FusionName"].apply(
        lambda f: fusion_matches_whitelist(f, wl_pairs, wl_genes)
    )

    sort_cols = ["is_aml_driver", "is_in_frame", "JunctionReadCount"]
    sort_asc  = [False, False, False]
    if args.recurrence_action == "flag":
        # push flagged recurrent artefacts to the bottom, drivers still float to top
        sort_cols.insert(1, "is_recurrent_artefact")
        sort_asc.insert(1, True)

    df = df.sort_values(sort_cols, ascending=sort_asc)

    # ── write output ────────────────────────────────────────────────────
    df.to_csv(args.output, sep="\t", index=False)

    n_final_samples = df["sample_name"].nunique()
    print("\n" + "=" * 60)
    print(f"FINAL: {len(df)} fusion calls across {n_final_samples} samples")
    print(f"  (started with {n_start} calls across {n_samples_start} samples)")
    print(f"  AML driver fusions retained: {df['is_aml_driver'].sum()}")
    if args.recurrence_action == "flag":
        print(f"  Flagged as recurrent artefact (kept, not removed): "
              f"{df['is_recurrent_artefact'].sum()}")
    print(f"  Output written to: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
