#========================= Run in Background on workstation or where no job scheduler =========================
`` nohup snakemake -s fusion_v2_final.smk --configfile config.yaml --cores 8 --resources star_fusion_slot=2 --until star_fusion > v2.log 2>&1 & ``
``disown``

#================================================== Run Filter ==================================================
``python3 ft-fusion.py --input  s800.tsv --output filtered_discovery_800aml.tsv --whitelist-file babiceanu_normal_fusions.txt --min_junction_reads 2 --min_spanning_frags 0 --min_ffpm 0.1 --recurrence_cap 0.4``
