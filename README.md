#========================= Run in Background on workstation or where no job scheduler =========================
``nohup snakemake -s fusion_v2_final.smk --configfile config.yaml --cores 8 --resources star_fusion_slot=2 --until star_fusion > v2.log 2>&1 &``
``disown``
