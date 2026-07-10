library(dplyr); library(tidyr)
library(ggplot2); library("readr")
setwd("/home/adminsb/pc100_80_42_95/SF_analysis/")

#========================
# User options
#========================
known_only <- T      # TRUE = only known fusions; FALSE = all recurrent fusions
top_n <- 130             # Inf = plot all


fusion_original <- read_tsv("all_samples_ft_filtered.tsv", show_col_types = FALSE)

# Canonicalize observed fusions
fusion <- fusion_original %>%
  mutate(
    Fusion = sapply(strsplit(FusionName, "--"),
                    function(x) paste(sort(trimws(x)), collapse = "--"))
  )

# Known leukemia fusions
fusion_list <- c(
  #======================================
  # Reference: MedGenome and Literature #
  #======================================
  "BCR--ABL1","ETV6--ABL1","ZMIZ1--ABL1","ETV6--ABL2","RANBP2--ALK","CHIC2--ETV6","MEF2D--CSF1R","DEK--NUP214",
  "FUS--ERG","MNX1--ETV6","HERVK--FGFR1","TRIM24--FGFR1","CBFA2T3--GLIS2","HOXA11--BZW2","TRB--HOXA11","NUP98--HOXA13",
  "GATA2--HOXA9","HOXA9--ANGPT1","NIPBL--HOXA9","NUP98--HOXA9","NUP98--HOXC11","NUP98--HOXC13","BCR--JAK2","PCM1--JAK2",
  "CREBBP--KAT6A","KAT6A--EP300","NUP98--KDM5A","TPM4--KLF2","KMT2A--ELL","KMT2A--EP300","KMT2A--EPS15","KMT2A--FOXO3",
  "KMT2A--FOXO4","KMT2A--MLLT1","KMT2A--MLLT10","KMT2A--MLLT11","KMT2A--MLLT3","KMT2A--MLLT6","KMT2A--PICALM","KMT2A--SEPT6",
  "KMT2A--SEPT9","NUP98--KMT2A","ETV6--MECOM","GATA2--MECOM","MECOM--CDK6","RBM15--MKL1","PICALM--MLLT10","MYB--GATA1",
  "CBFB--MYH11","NF1--LRRC37B","NUP98--NSD1","ETV6--NTRK3","SQSTM1--NUP214","ERC1--PDGFRB","TRIP11--PDGFRB","ETV6--PRDM16",
  "ETV6--RARA","FIP1L1--RARA","GTF2I--RARA","NPM1--RARA","NUMA1--RARA","PML--RARA","STAT5B--RARA","ZBTB16--RARA",
  "MECOM--RUNX1","RUNX1--CBFA2T3","RUNX1--PRDM16","RUNX1--RUNX1T1","NUP98--TOP1",
  #======================================
  # Additional AML fusions (Literature) #
  #======================================
  "IRF2BP2--RARA","TBL1XR1--RARA","BCOR--RARA","PRKAR1A--RARA","FNDC3B--RARA","STAT3--RARA","FIP1L1--PDGFRA","ETV6--PDGFRB",
  "FGFR1OP--FGFR1","BCR--FGFR1","CNTRL--FGFR1","CEP110--FGFR1","ZMYM2--FGFR1","KMT2A--AFDN","KMT2A--ARHGEF12"
)

# Canonicalize the reference list (X--Y == Y--X)
fusion_list <- unique(sapply(strsplit(fusion_list, "--"),
         function(x) paste(sort(trimws(x)), collapse = "--")))

# Count unique samples per fusion
fusion_summary <- fusion %>%
  distinct(sample_name, Fusion) %>%
  count(Fusion, name = "Sample_Count", sort = TRUE)

# Optional: keep any fusion involving a known fusion gene
if (known_only) {
  # Extract all genes from the known fusion list
  known_genes <- unique(unlist(strsplit(fusion_list, "--")))
  fusion_summary <- fusion %>%
    separate(Fusion, into = c("Gene1", "Gene2"), sep = "--", remove = FALSE) %>%
    filter(Gene1 %in% known_genes | Gene2 %in% known_genes) %>%
    distinct(sample_name, Fusion) %>%
    count(Fusion, name = "Sample_Count", sort = TRUE)
}

# Save results
write.csv(fusion_summary, ifelse(known_only, "Known_Fusions_Sample_Count.csv", "All_Fusions_Sample_Count.csv"),
          row.names = FALSE)

# Plot
plot_data <- if (is.infinite(top_n)) fusion_summary else head(fusion_summary, top_n)

ggplot(plot_data,
       aes(x = reorder(Fusion, Sample_Count), y = Sample_Count)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 8) +
  labs(title = ifelse(known_only, "Known leukemia fusions", "All recurrent fusions"),
       x = "Fusion", y = "Number of samples")



#================================
# sample-fusion summary details #
#================================
if (known_only) {
  # Extract all genes from the known fusion list
  known_genes <- unique(unlist(strsplit(fusion_list, "--")))
  fusion_summary_1 <- fusion %>%
    separate(Fusion, into = c("Gene1", "Gene2"), sep = "--", remove = FALSE) %>%
    filter(Gene1 %in% known_genes | Gene2 %in% known_genes) %>%
    distinct(sample_name, Fusion)}

# fusion_summary_1 <- fusion # ---> Open for All Samples/Fusions
# Number of unique samples with known AML fusion genes
total_uniq_samples <- length(unique(fusion_summary_1$sample_name))
print(paste("Total Unique samples in Fusion Detection:", total_uniq_samples))
# Number of unique RNA fusions per sample
fusion_per_sample <- fusion_summary_1 %>%
  group_by(sample_name) %>%
  summarise(Fusion_Count = n(), .groups = "drop") %>%
  arrange(desc(Fusion_Count))

# Number of samples per fusion
samples_per_fusion<- fusion_summary_1 %>%
  group_by(Fusion) %>%
  summarise(Sample_Count = n(), .groups = "drop") %>%
  arrange(desc(Sample_Count))
