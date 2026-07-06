#!/bin/bash
set -euo pipefail

# ── usage ──────────────────────────────────────────────────────────────────
# ./combine_fusions.sh <search_dir> <outfile> <filename_pattern>
#
# <search_dir>        directory to search under, e.g. /path/to/results
# <outfile>           path to the combined output TSV
# <filename_pattern>  glob pattern passed to `find -name`,
#                      e.g. "*-fusion.fusion_predictions.abridged.tsv"
# ─────────────────────────────────────────────────────────────────────────

if [ $# -lt 3 ]; then
    echo "Usage: $0 <search_dir> <outfile> <filename_pattern>" >&2
    exit 1
fi

search_dir="$1"
outfile="$2"
name_pattern="$3"

rm -f "$outfile"

search_dir_abs=$(realpath "$search_dir")

find "$search_dir" -type f -name "$name_pattern" | sort | while read -r f; do
    f_abs=$(realpath "$f")
    rel=${f_abs#"$search_dir_abs"/}
    sample=${rel%%/*}          # first directory under search_dir = sample name, regardless of depth
    path=$(dirname "$f")

    if [ ! -f "$outfile" ]; then
        # first file: keep header, append literal column names "path" and "sample_name"
        awk 'BEGIN{OFS="\t"} NR==1{print $0,"path","sample_name"}' "$f" > "$outfile"
        awk -v p="$path" -v s="$sample" 'BEGIN{OFS="\t"} NR>1{print $0,p,s}' "$f" >> "$outfile"
    else
        # subsequent files: skip header, append actual path/sample values
        awk -v p="$path" -v s="$sample" 'BEGIN{OFS="\t"} NR>1{print $0,p,s}' "$f" >> "$outfile"
    fi
done

if [ -f "$outfile" ]; then
    n_files=$(find "$search_dir" -type f -name "$name_pattern" | wc -l)
    n_rows=$(($(wc -l < "$outfile") - 1))
    echo "Combined output written to: $outfile  (${n_files} files, ${n_rows} data rows)"
else
    echo "No files matched pattern '$name_pattern' under '$search_dir' — nothing written." >&2
    exit 1
fi
