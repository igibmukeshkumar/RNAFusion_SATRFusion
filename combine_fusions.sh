#!/bin/bash

while getopts "i:o:p:" opt; do
    case $opt in
        i) studies+=("$OPTARG") ;;
        o) out="$OPTARG" ;;
        p) pattern="$OPTARG" ;;
        *) echo "Usage: $0 -i study_path [-i study_path ...] -o output.tsv -p 'pattern'"; exit 1 ;;
    esac
done

rm -f "$out"

first=1

for study in "${studies[@]}"; do
    study=$(realpath "$study")

    while read -r f; do

        f=$(realpath "$f")
        rel="${f#"$study"/}"
        sample="${rel%%/*}"

        if [ "$first" -eq 1 ]; then

            awk -v p="$f" -v s="$sample" \
                'BEGIN{OFS="\t"}
                 NR==1 {print $0,"path","sample_name"; next}
                 {print $0,p,s}' \
                "$f" > "$out"

            first=0

        else

            awk -v p="$f" -v s="$sample" \
                'NR>1 {print $0,p,s}' \
                "$f" >> "$out"

        fi

    done < <(find "$study" -type f -name "$pattern" | sort)

done

echo "Output: $(realpath "$out")"
