#!/usr/bin/env python3
"""
fusion_filter_fi.py
-------------------
FusionInspector-aware RNA fusion filtering for STAR-Fusion output.

FILTER MODE (--filter-mode):
  strict   (default) : ALL filters must pass  (A AND B AND C ...)
  any                : ANY ONE filter passes  (A OR  B OR  C ...)
  support            : (JunctionRead AND SpanningFrag AND FFPM) pass — other filters ignored
  custom             : specify exact conditions via --require and --any-of groups

CUSTOM MODE EXAMPLES:
  Require junction AND ffpm, but accept if EITHER spanning OR large_anchor passes:
    --filter-mode custom
    --require junction ffpm
    --any-of spanning large_anchor

  Require junction AND (spanning OR est_j):
    --filter-mode custom
    --require junction
    --any-of spanning est_j

WHITELIST:
  Known fusions bypass ALL filters when --whitelist-override (default ON).
  Use --no-whitelist-override to apply filters even to known fusions.
  Whitelist entries: full pairs (GENE_A--GENE_B) or single-gene (GENE_A matches any partner).
"""

import argparse
import re
import sys
import pandas as pd


# ── default cut-offs ──────────────────────────────────────────────────────────
DEFAULTS = {
    "min_junction":              3,
    "min_spanning":              2,
    "min_est_j":                 3,
    "min_est_s":                 1,
    "min_ffpm":                  0.10,
    "require_large_anchor":      False,
    "max_counter_left":          1000,
    "max_counter_right":         1000,
    "min_left_entropy":          1.0,
    "min_right_entropy":         1.0,
    "min_far_left":              0.0,
    "min_far_right":             0.0,
    "max_microh_brkpt_dist":     10000,
    "max_num_microh_near_brkpt": 10,
    "whitelist_override":        True,
}

# required columns enforced at load time
REQUIRED_COLS = {"JunctionReadCount", "SpanningFragCount", "FFPM"}

# valid filter keys for --require / --any-of in custom mode
FILTER_KEYS = [
    "junction", "spanning", "est_j", "est_s", "ffpm",
    "large_anchor", "counter_left", "counter_right",
    "left_entropy", "right_entropy",
    "far_left", "far_right",
    "microh_dist", "microh_count",
]


# ── helpers ───────────────────────────────────────────────────────────────────
def norm_fusion(x):
    """Normalise FusionName for reliable whitelist matching."""
    if pd.isna(x):
        return ""
    x = str(x).strip()
    x = re.sub(r"\s+", "", x)
    x = x.replace("::", "--").replace("__", "--").replace("/", "--")
    return x


def load_whitelist(path):
    """
    Load driver fusion whitelist.
    Supports full pairs (GENE_A--GENE_B) and single-gene entries (GENE_A).
    Raises FileNotFoundError clearly if file missing.
    Returns (pair_set, gene_set).
    """
    if not path:
        return set(), set()

    import os
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Whitelist file not found: {path}")

    raw_entries = set()

    try:
        df = pd.read_csv(path, sep=None, engine="python", comment="#", dtype=str)
        cols = {c.lower(): c for c in df.columns}
        for name in ("fusionname", "fusion_name", "fusion", "knownfusion"):
            if name in cols:
                raw_entries = {
                    norm_fusion(x)
                    for x in df[cols[name]].dropna()
                    if norm_fusion(x)
                }
                break
    except Exception:
        pass

    if not raw_entries:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                value = re.split(r"[\t,]", line)[0].strip()
                if value.lower() in {"fusionname", "fusion_name", "fusion", "knownfusion"}:
                    continue
                value = norm_fusion(value)
                if value:
                    raw_entries.add(value)

    pair_set = {e for e in raw_entries if "--" in e}
    gene_set  = {e for e in raw_entries if "--" not in e}
    print(f"Whitelist: {len(pair_set)} fusion pairs, {len(gene_set)} single-gene entries")
    return pair_set, gene_set


def is_whitelisted_fusion(fusion_name, pair_set, gene_set, direction="exact"):
    """
    Check if a fusion name matches the whitelist.

    direction modes:
      exact       X--Y matches ONLY X--Y as written in whitelist (conservative, default)
      reverse     X--Y matches ONLY Y--X in whitelist (flipped orientation)
      both_dir    X--Y matches X--Y OR Y--X in whitelist (direction-agnostic pair match)
      any         X--Y matches X--Y, Y--X, plus single-gene entries (broadest pair match)
      pair_1match X--Y matches if EITHER gene of ANY whitelist pair equals X or Y
                  e.g. whitelist has A--B → matches A--C, C--A, B--C, C--B

    Single-gene entries (gene_set) always match either partner for:
      both_dir, any, pair_1match
    For exact and reverse, gene_set is NOT used (strict pair-only matching).
    """
    parts = fusion_name.split("--")
    if len(parts) != 2:
        return False

    a, b    = parts[0], parts[1]
    flipped = f"{b}--{a}"

    if direction == "exact":
        return fusion_name in pair_set

    elif direction == "reverse":
        return flipped in pair_set

    elif direction == "both_dir":
        if fusion_name in pair_set or flipped in pair_set:
            return True
        # gene_set active for both_dir
        if a in gene_set or b in gene_set:
            return True
        return False

    elif direction == "any":
        # full pair either direction + single-gene
        if fusion_name in pair_set or flipped in pair_set:
            return True
        if a in gene_set or b in gene_set:
            return True
        return False

    elif direction == "pair_1match":
        # whitelist pair X--Y → extract all individual genes from pair_set
        # match if data fusion shares at least one gene with any whitelist pair
        # e.g. whitelist A--B matches data A--C or C--A or B--C or C--B
        wl_genes_from_pairs = set()
        for pair in pair_set:
            pp = pair.split("--")
            if len(pp) == 2:
                wl_genes_from_pairs.update(pp)
        if a in wl_genes_from_pairs or b in wl_genes_from_pairs:
            return True
        # also check gene_set entries
        if a in gene_set or b in gene_set:
            return True
        return False

    else:
        raise ValueError(f"Unknown match_direction: {direction}")


def num(df, col, default=0):
    """Safe numeric column."""
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def yes(df, col):
    """YES/YES_LDAS/TRUE/1/Y flag column."""
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].astype(str).str.upper().isin({"YES", "YES_LDAS", "TRUE", "1", "Y"})


def build_filter_masks(df, args):
    """
    Build one boolean Series per named filter key.
    Returns dict { key: Series } for use in AND/OR/custom logic.
    """
    return {
        "junction":      num(df, "JunctionReadCount")  >= args.min_junction,
        "spanning":      num(df, "SpanningFragCount")  >= args.min_spanning,
        "est_j":         num(df, "est_J")              >= args.min_est_j,
        "est_s":         num(df, "est_S")              >= args.min_est_s,
        "ffpm":          num(df, "FFPM")               >= args.min_ffpm,
        "large_anchor":  yes(df, "LargeAnchorSupport"),
        "counter_left":  num(df, "NumCounterFusionLeft")  <= args.max_counter_left,
        "counter_right": num(df, "NumCounterFusionRight") <= args.max_counter_right,
        "left_entropy":  num(df, "LeftBreakEntropy")   >= args.min_left_entropy,
        "right_entropy": num(df, "RightBreakEntropy")  >= args.min_right_entropy,
        "far_left":      num(df, "FAR_left")           >= args.min_far_left,
        "far_right":     num(df, "FAR_right")          >= args.min_far_right,
        "microh_dist":   num(df, "microh_brkpt_dist")  <= args.max_microh_brkpt_dist,
        "microh_count":  num(df, "num_microh_near_brkpt") <= args.max_num_microh_near_brkpt,
    }


def apply_filter_mode(masks, args):
    """
    Combine individual filter masks according to --filter-mode.

    strict   : ALL masks must be True  (AND of everything)
    any      : ANY mask is True        (OR of everything)
    support  : junction AND spanning AND ffpm only
    custom   : --require keys all True AND at least one --any-of key True
    """
    mode = args.filter_mode
    idx  = next(iter(masks.values())).index   # shared index

    if mode == "strict":
        mask = pd.Series(True, index=idx)
        for m in masks.values():
            mask &= m
        return mask

    elif mode == "any":
        mask = pd.Series(False, index=idx)
        for m in masks.values():
            mask |= m
        return mask

    elif mode == "support":
        return masks["junction"] & masks["spanning"] & masks["ffpm"]

    elif mode == "custom":
        require_keys = args.require or []
        anyof_keys   = args.any_of  or []

        # validate keys
        bad = [k for k in require_keys + anyof_keys if k not in masks]
        if bad:
            raise ValueError(
                f"Unknown filter key(s): {bad}\n"
                f"Valid keys: {FILTER_KEYS}"
            )

        mask = pd.Series(True, index=idx)

        # all --require keys must pass (AND)
        for k in require_keys:
            mask &= masks[k]

        # at least one --any-of key must pass (OR), if any specified
        if anyof_keys:
            anyof_mask = pd.Series(False, index=idx)
            for k in anyof_keys:
                anyof_mask |= masks[k]
            mask &= anyof_mask

        return mask

    else:
        raise ValueError(f"Unknown filter-mode: {mode}")


def print_filter_summary(masks, final_mask, wl_hits, whitelist_override, mode):
    """Print per-filter kept/dropped counts."""
    labels = {
        "junction":      "JunctionReadCount",
        "spanning":      "SpanningFragCount",
        "est_j":         "est_J",
        "est_s":         "est_S",
        "ffpm":          "FFPM",
        "large_anchor":  "LargeAnchorSupport",
        "counter_left":  "NumCounterFusionLeft",
        "counter_right": "NumCounterFusionRight",
        "left_entropy":  "LeftBreakEntropy",
        "right_entropy": "RightBreakEntropy",
        "far_left":      "FAR_left",
        "far_right":     "FAR_right",
        "microh_dist":   "microh_brkpt_dist",
        "microh_count":  "num_microh_near_brkpt",
    }
    n = len(final_mask)
    print(f"\nFilters (mode={mode}):")
    for key, m in masks.items():
        dropped = (~m & (~wl_hits if whitelist_override else pd.Series(True, index=m.index))).sum()
        print(f"  {labels[key]:45s} kept={m.sum():5d}  dropped={dropped:5d}")


def fusion_filter(df, pair_set, gene_set, args):
    """Apply filters and return filtered DataFrame."""
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"Input missing required columns: {missing}")

    fusion_col = "#FusionName" if "#FusionName" in df.columns else "FusionName"
    df     = df.copy()
    fusion = df[fusion_col].map(norm_fusion)

    # whitelist is active ONLY when --known-fusions was supplied (Bug 1)
    use_whitelist = bool(pair_set or gene_set)
    wl_hits = (
        fusion.apply(lambda f: is_whitelisted_fusion(f, pair_set, gene_set, args.match_direction))
        if use_whitelist
        else pd.Series(False, index=df.index)
    )

    # build individual masks
    masks = build_filter_masks(df, args)

    # print summary before combining
    combined = apply_filter_mode(masks, args)
    print_filter_summary(masks, combined, wl_hits, args.whitelist_override, args.filter_mode)

    # whitelist override only when --known-fusions was supplied (Bug 1)
    if use_whitelist and args.whitelist_override:
        combined |= wl_hits

    out = df.loc[combined].copy()

    # KNOWN_FUSION column only when --known-fusions was supplied (Bug 2)
    if use_whitelist:
        out["KNOWN_FUSION"] = wl_hits.loc[out.index].map({True: "YES", False: "NO"})

    if args.sample_name:
        out.insert(0, "sample_name", args.sample_name)

    return out


def print_quick_summary(args):
    """Print a clean pre-run summary of active settings."""
    W = 70
    D = DEFAULTS

    def changed(val, default): return " *" if val != default else ""

    # ── filter mode block ────────────────────────────────────────────────────
    lines = ["", "=" * W, " QUICK SUMMARY", "=" * W]

    # filter mode description
    lines.append(f" Filter mode : {args.filter_mode.upper()}")
    if args.filter_mode == "strict":
        lines.append("   ALL filters must pass (AND of every condition)")
    elif args.filter_mode == "any":
        lines.append("   ANY ONE filter must pass (OR of every condition)")
    elif args.filter_mode == "support":
        lines.append("   Junction AND Spanning AND FFPM only")
    elif args.filter_mode == "custom":
        req  = " AND ".join(args.require) if args.require else "(none)"
        anyf = " OR  ".join(args.any_of)  if args.any_of  else "(none)"
        lines.append(f"   --require : {req}")
        lines.append(f"   --any-of  : {anyf}")
        if args.require and args.any_of:
            lines.append(f"   Logic     : ({req}) AND ({anyf})")
        elif args.require:
            lines.append(f"   Logic     : {req}")
        else:
            lines.append(f"   Logic     : {anyf}")

    lines.append("-" * W)

    # cut-off table — mark user-changed values with *
    lines.append(" Cut-offs (* = changed from default):")
    lines.append(f"   {'Parameter':<35} {'Value':>10}  {'Default':>10}")
    lines.append(f"   {'-'*35} {'-'*10}  {'-'*10}")

    rows = [
        ("JunctionReadCount >=",     args.min_junction,              D["min_junction"]),
        ("SpanningFragCount >=",      args.min_spanning,              D["min_spanning"]),
        ("est_J >=",                  args.min_est_j,                 D["min_est_j"]),
        ("est_S >=",                  args.min_est_s,                 D["min_est_s"]),
        ("FFPM >=",                   args.min_ffpm,                  D["min_ffpm"]),
        ("LargeAnchorSupport",        "required" if args.require_large_anchor else "not required",
                                      "not required"),
        ("NumCounterFusionLeft <=",   args.max_counter_left,          D["max_counter_left"]),
        ("NumCounterFusionRight <=",  args.max_counter_right,         D["max_counter_right"]),
        ("LeftBreakEntropy >=",       args.min_left_entropy,          D["min_left_entropy"]),
        ("RightBreakEntropy >=",      args.min_right_entropy,         D["min_right_entropy"]),
        ("FAR_left >=",               args.min_far_left,              D["min_far_left"]),
        ("FAR_right >=",              args.min_far_right,             D["min_far_right"]),
        ("microh_brkpt_dist <=",      args.max_microh_brkpt_dist,     D["max_microh_brkpt_dist"]),
        ("num_microh_near_brkpt <=",  args.max_num_microh_near_brkpt, D["max_num_microh_near_brkpt"]),
    ]

    for label, val, default in rows:
        flag = " *" if val != default else "  "
        lines.append(f"   {label:<35} {str(val):>10}  {str(default):>10}{flag}")

    lines.append("-" * W)

    # whitelist block
    if args.known_fusions:
        lines.append(f" Whitelist    : {args.known_fusions}")
        lines.append(f" Match mode   : {args.match_direction}")
        override = "ON  — known fusions bypass all filters" if args.whitelist_override \
                   else "OFF — known fusions still filtered"
        lines.append(f" Override     : {override}")
        if getattr(args, "whitelist_only", False):
            lines.append(f" Output scope : WHITELIST ONLY (--whitelist-only active)")
    else:
        lines.append(" Whitelist    : not supplied")

    # sample filter block
    if args.sample_name:
        lines.append(f" Sample filter: '{args.sample_name}'")

    lines.append("=" * W)
    print("\n".join(lines))


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── I/O ──────────────────────────────────────────────────────────────────
    p.add_argument("-i", "--input",        required=True,
                   help="STAR-Fusion / FusionInspector TSV input")
    p.add_argument("-o", "--output",       required=True,
                   help="filtered output TSV")
    p.add_argument("--sample-name",        default=None,
                   help="prepend sample_name column (useful for multi-sample merging)")
    p.add_argument("--known-fusions", default=None,
                   help=(
                       "Whitelist file of driver fusions to protect from the recurrence cap.\n"
                       "Accepted formats:\n"
                       "  one entry per line, or TSV/CSV with a FusionName column.\n\n"
                       "Three entry types (can be mixed in the same file):\n"
                       "  X--Y   exact pair  — matches only X--Y as written (default behaviour)\n"
                       "  Y--X   reversed    — use --match-direction reverse or any to match X--Y in data\n"
                       "  X      single gene — matches X--<anything> OR <anything>--X regardless of direction\n\n"
                       "See --match-direction to control how pair entries are matched."
                   ))
    p.add_argument(
        "--match-direction",
        choices=["exact", "reverse", "both_dir", "any", "pair_1match"],
        default="exact",
        dest="match_direction",
        help=(
            "Controls how whitelist fusion pairs (X--Y) are matched against fusion calls in the data.\n"
            "Five modes:\n\n"
            "  exact      (default)\n"
            "             Whitelist X--Y matches ONLY X--Y in data.\n"
            "             Most conservative — use when whitelist and STAR-Fusion output\n"
            "             share the same gene orientation.\n\n"
            "  reverse\n"
            "             Whitelist X--Y matches ONLY Y--X in data.\n"
            "             Use when your whitelist orientation is opposite to STAR-Fusion output.\n\n"
            "  both_dir\n"
            "             Whitelist X--Y matches X--Y OR Y--X in data.\n"
            "             Direction-agnostic pair match. Also activates single-gene entries.\n"
            "             Example: whitelist PML--RARA matches PML--RARA and RARA--PML.\n\n"
            "  any\n"
            "             Whitelist X--Y matches X--Y or Y--X, PLUS single-gene entries\n"
            "             (gene_set) match any fusion containing that gene.\n"
            "             Broadest pair+gene matching mode.\n\n"
            "  pair_1match\n"
            "             A whitelist pair X--Y matches ANY fusion that contains X or Y\n"
            "             as either partner, regardless of the other gene.\n"
            "             Example: whitelist A--B matches A--C, C--A, B--C, C--B.\n"
            "             Use for gene-centric whitelisting without listing every partner.\n\n"
            "Single-gene entries (e.g. KMT2A without '--') always match any fusion\n"
            "containing that gene in both_dir, any, and pair_1match modes.\n"
            "In exact and reverse modes, only full pair entries are used."
        ),
    )

    # ── filter mode ───────────────────────────────────────────────────────────
    p.add_argument(
        "--filter-mode",
        choices=["strict", "any", "support", "custom"],
        default="strict",
        help=(
            "How to combine filter conditions:\n"
            "  strict  (default) ALL filters must pass (A AND B AND C ...)\n"
            "  any              ANY one filter passes  (A OR  B OR  C ...)\n"
            "  support          junction AND spanning AND ffpm only\n"
            "  custom           use --require and --any-of to define exact logic\n"
            "                   e.g. --require junction ffpm --any-of spanning large_anchor\n"
            "                   means: (junction AND ffpm) AND (spanning OR large_anchor)"
        ),
    )
    p.add_argument(
        "--require",
        nargs="+",
        metavar="KEY",
        help=(
            "[custom mode] filter keys that ALL must pass (AND logic). "
            f"Valid keys: {', '.join(FILTER_KEYS)}"
        ),
    )
    p.add_argument(
        "--any-of",
        nargs="+",
        metavar="KEY",
        dest="any_of",
        help=(
            "[custom mode] filter keys where AT LEAST ONE must pass (OR logic). "
            f"Valid keys: {', '.join(FILTER_KEYS)}"
        ),
    )

    # ── read support ──────────────────────────────────────────────────────────
    p.add_argument("--min-junction", type=float, default=DEFAULTS["min_junction"],
                   help="min reads crossing exact fusion breakpoint (default: %(default)s)")
    p.add_argument("--min-spanning", type=float, default=DEFAULTS["min_spanning"],
                   help="min spanning fragment count (default: %(default)s)")
    p.add_argument("--min-est-j",    type=float, default=DEFAULTS["min_est_j"],
                   help="min estimated junction reads est_J (default: %(default)s)")
    p.add_argument("--min-est-s",    type=float, default=DEFAULTS["min_est_s"],
                   help="min estimated spanning reads est_S (default: %(default)s)")

    # ── expression ────────────────────────────────────────────────────────────
    p.add_argument("--min-ffpm",     type=float, default=DEFAULTS["min_ffpm"],
                   help="min fusion fragments per million reads (default: %(default)s)")

    # ── anchor ────────────────────────────────────────────────────────────────
    p.add_argument("--require-large-anchor", action="store_true",
                   help="require LargeAnchorSupport=YES (mandatory when spanning=0)")

    # ── counter fusion ────────────────────────────────────────────────────────
    p.add_argument("--max-counter-left",  type=float, default=DEFAULTS["max_counter_left"],
                   help="max counter-fusion reads supporting left gene (default: %(default)s)")
    p.add_argument("--max-counter-right", type=float, default=DEFAULTS["max_counter_right"],
                   help="max counter-fusion reads supporting right gene (default: %(default)s)")

    # ── breakpoint entropy ────────────────────────────────────────────────────
    p.add_argument("--min-left-entropy",  type=float, default=DEFAULTS["min_left_entropy"],
                   help="min sequence entropy at left breakpoint — low = repetitive (default: %(default)s)")
    p.add_argument("--min-right-entropy", type=float, default=DEFAULTS["min_right_entropy"],
                   help="min sequence entropy at right breakpoint (default: %(default)s)")

    # ── FAR ───────────────────────────────────────────────────────────────────
    p.add_argument("--min-far-left",  type=float, default=DEFAULTS["min_far_left"],
                   help="min fraction of abnormal reads at left breakpoint (default: %(default)s)")
    p.add_argument("--min-far-right", type=float, default=DEFAULTS["min_far_right"],
                   help="min fraction of abnormal reads at right breakpoint (default: %(default)s)")

    # ── microhomology ─────────────────────────────────────────────────────────
    p.add_argument("--max-microh-brkpt-dist",     type=float, default=DEFAULTS["max_microh_brkpt_dist"],
                   help="max microhomology breakpoint distance (default: %(default)s)")
    p.add_argument("--max-num-microh-near-brkpt", type=float, default=DEFAULTS["max_num_microh_near_brkpt"],
                   help="max number of microhomology events near breakpoint (default: %(default)s)")

    # ── whitelist ─────────────────────────────────────────────────────────────
    p.add_argument("--no-whitelist-override", dest="whitelist_override",
                   action="store_false", default=True,
                   help="disable whitelist override — known fusions still filtered if thresholds fail")
    p.add_argument("--whitelist-only", action="store_true", default=False,
                   help=(
                       "output ONLY fusions that matched the whitelist (KNOWN_FUSION=YES).\n"
                       "Normal filter thresholds are still applied unless --no-whitelist-override is set.\n"
                       "Requires --known-fusions to be supplied."
                   ))

    args = p.parse_args()

    # validate custom mode args
    if args.filter_mode == "custom" and not args.require and not args.any_of:
        p.error("--filter-mode custom requires at least --require or --any-of")
    if args.filter_mode != "custom" and (args.require or args.any_of):
        p.error("--require and --any-of are only used with --filter-mode custom")

    # validate whitelist-only
    if args.whitelist_only and not args.known_fusions:
        p.error("--whitelist-only requires --known-fusions to be supplied")

    # print summary before any processing
    print_quick_summary(args)

    # load input
    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    print(f"\nInput  : {len(df)} fusion calls")

    # ── sample filter ─────────────────────────────────────────────────────────
    # --sample-name behaviour depends on whether column already exists:
    #   column EXISTS in file → filter rows to that sample (master file use-case)
    #   column ABSENT         → add it as a new annotation column
    if args.sample_name:
        if "sample_name" in df.columns:
            before = len(df)
            df = df[df["sample_name"] == args.sample_name].reset_index(drop=True)
            print(f"Sample filter : '{args.sample_name}' → {len(df)} of {before} rows")
            if len(df) == 0:
                avail = sorted(df["sample_name"].unique().tolist()) if before > 0 else []
                sys.exit(
                    f"ERROR: no rows found for sample '{args.sample_name}'.\n"
                    f"Available sample names: {avail}"
                )
        # column absent — fusion_filter() will insert it

    # load whitelist
    pair_set, gene_set = load_whitelist(args.known_fusions)

    # filter
    filtered = fusion_filter(df, pair_set, gene_set, args)

    # write output
    if args.whitelist_only:
        before_wl = len(filtered)
        filtered = filtered[filtered["KNOWN_FUSION"] == "YES"].copy()
        print(f"Whitelist-only: {len(filtered)} of {before_wl} rows retained (KNOWN_FUSION=YES)")

    filtered.to_csv(args.output, sep="\t", index=False)

    print(f"\nOutput : {len(filtered)} fusion calls → {args.output}")
    if args.known_fusions and "KNOWN_FUSION" in filtered.columns:
        n_known = (filtered["KNOWN_FUSION"] == "YES").sum()
        print(f"         {n_known} whitelisted fusions retained")


if __name__ == "__main__":
    main()
