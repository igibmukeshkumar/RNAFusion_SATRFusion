library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
pfam <- getBM( attributes = c( "ensembl_transcript_id", "ensembl_peptide_id", "pfam",
    "pfam_start", "pfam_end"), mart = ensembl)
write.table(pfam,  "pfam_annotations.tsv",  sep = "\t",  quote = FALSE,  row.names = FALSE)

## Prepare chimeraviz BED file
pfam <- read.delim("pfam_annotations.tsv", stringsAsFactors = FALSE)
colnames(pfam) <- c("Transcript_id","Protein_id","Pfam_id","Start","End")
pfam$Domain_name_abbreviation <- pfam$Pfam_id
pfam$Domain_name_full <- pfam$Pfam_id
pfam <- pfam[, c("Transcript_id","Pfam_id","Domain_name_abbreviation","Domain_name_full","Start","End")]
write.table(pfam, "pfam_domains.bed", sep="\t", quote=FALSE, row.names=FALSE)

## Plot protein domains
library(chimeraviz); library(EnsDb.Hsapiens.v86)
plot_fusion_transcript_with_protein_domain(fusion = fusion1,edb = EnsDb.Hsapiens.v86, bedfile = "pfam_domains.bed",
  gene_upstream_transcript = "ENST00000300305", gene_downstream_transcript = "ENST00000520724")
