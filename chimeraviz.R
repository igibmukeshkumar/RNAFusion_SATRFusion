#############################################################
## Load packages
#############################################################

library(chimeraviz); library(EnsDb.Hsapiens.v86)
library(GenomicFeatures)
library(GenomicAlignments)
library(ggplot2)
setwd("/home/adminsb/disk2/meta-analysis/rnafusion/")

#############################################################
## Filters
#############################################################
library(dplyr)

fusion <- read.delim("all_samples_ft.tsv")

fusion.filtered <- fusion %>%
  filter(
    JunctionReadCount >= 10,
    SpanningFragCount >= 5,
    FFPM >= 0.1,
    LargeAnchorSupport == "YES",
    is_in_frame %in% c(TRUE, "True"),
    !grepl("readthrough|GTEx|BodyMap|Normal", annots, ignore.case = TRUE),
    !grepl("^IG[HKL]|^TRA|^TRB|^TRD|^TRG|^HLA|^RPL|^RPS|^MT-", LeftGene),
    !grepl("^IG[HKL]|^TRA|^TRB|^TRD|^TRG|^HLA|^RPL|^RPS|^MT-", RightGene)
  ) %>%
  distinct(sample_name, FusionName, LeftBreakpoint, RightBreakpoint, .keep_all = TRUE)

write.table(
  fusion.filtered,
  "all_samples_ft_filtered.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)





#############################################################
## chimeraviz v1.38.0 workflow
#############################################################

library(chimeraviz)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

## Import STAR-Fusion output
fusion <- import_starfusion(
  filename = "all_samples_ft_filtered.tsv",
  genome_version = "hg38",
  limit = Inf
)

## Select one fusion (e.g. first fusion)
fusion1 <- fusion[[1]]

library(methods)


# check for sample_name it the html report summary 
df <- read.delim("all_samples_ft_filtered.tsv", stringsAsFactors = FALSE)

for (i in seq_along(fusion)) {
  slot(fusion[[i]], "id") <- df$sample_name[i]
}
fusion[[1]]

create_fusion_report(fusion, "output.html")

# Plot 1: Fusion events overview
plot_circle(fusion)
create_fusion_report(fusion, "output.html")
# Plot 2: Fusion partner transcript structures
plot_transcripts(fusion = fusion1, edb = edb)
plot_transcripts(fusion1, edb = edb, bamfile = NULL, which_transcripts = "exonBoundary", non_ucsc = TRUE, ylim = c(0, 1000), reduce_transcripts = FALSE, bedgraphfile = NULL)

# Plot 3: Predicted fusion transcript model
plot_fusion_transcript(fusion = fusion1, edb = edb)

# Plot 4: Transcript plots for all fusions
for(i in seq_along(fusion)){pdf(paste0("Fusion_transcripts_",i,".pdf")); plot_transcripts(fusion[[i]], edb); dev.off()}

# Plot 5: Fusion transcript models for all fusions
for(i in seq_along(fusion)){pdf(paste0("Fusion_transcript_model_",i,".pdf")); plot_fusion_transcript(fusion[[i]], edb); dev.off()}

# Plot 6: Full fusion plot with BAM coverage (requires BAM)
# plot_fusion(fusion1, bamfile="sample.bam", edb=edb)

# Plot 7: Fusion supporting reads plot (requires BAM)
# plot_fusion_reads(fusion1, bamfile="sample.bam")

# Plot 8: Full fusion plots for all fusions with BAM (requires BAM)
# for(i in seq_along(fusion)){pdf(paste0("Fusion_full_",i,".pdf")); plot_fusion(fusion[[i]], bamfile="sample.bam", edb=edb); dev.off()}

# Plot 9: Fusion read evidence plots for all fusions with BAM (requires BAM)
# for(i in seq_along(fusion)){pdf(paste0("Fusion_reads_",i,".pdf")); plot_fusion_reads(fusion[[i]], bamfile="sample.bam"); dev.off()}
