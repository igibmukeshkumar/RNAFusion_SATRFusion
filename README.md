# RNA Fusion Detection Pipeline

## Step 1: Run the Snakemake Pipeline

For systems without a job scheduler (e.g., workstation or standalone server), execute the pipeline in the background using `nohup`:

```bash
nohup snakemake -s fusion.smk \
    --configfile config.yaml \
    --cores 8 \
    --resources star_fusion_slot=2 \
    --until star_fusion \
    > v2.log 2>&1 &
```

Detach the job from the current shell:

```bash
disown
```

---

## Step 2: Combine FusionInspector Results

Merge all `FusionInspector` abridged fusion files into a single tab-delimited file.

```bash
./combine_fusions.sh ./ all_samples_abridged.tsv "*.FusionInspector.fusions.abridged.tsv"
```

This generates a combined TSV containing all FusionInspector fusion predictions with the corresponding sample name.

---

## Step 3: Filter High-Confidence Fusions

Filter the combined fusion calls using the provided whitelist of normal fusions and user-defined thresholds.

```bash
python3 filter_fusion.py \
    --input test1.tsv \
    --output filtered_discovery_800aml.tsv \
    --whitelist-file normal_fusions.txt \
    --min_junction_reads 2 \
    --min_spanning_frags 0 \
    --min_ffpm 0.1 \
    --recurrence_cap 0.4
```

### Filtering Parameters

| Parameter | Description | Value |
|-----------|-------------|------:|
| `--min_junction_reads` | Minimum number of junction-supporting reads | 2 |
| `--min_spanning_frags` | Minimum number of spanning fragments | 0 |
| `--min_ffpm` | Minimum fusion expression (FFPM) | 0.1 |
| `--recurrence_cap` | Maximum recurrence frequency allowed before filtering | 0.4 |
| `--whitelist-file` | List of known normal/recurrent fusions to remove | `known_normal_fusions.txt` |

The final output (`filtered_discovery_800aml.tsv`) contains high-confidence fusion candidates for downstream analysis.
