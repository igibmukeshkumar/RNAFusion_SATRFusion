library(dplyr)
library(tidyr)
library(readr)

url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
hgnc <- readr::read_tsv(url, show_col_types = FALSE)

gene_map <- hgnc %>%
  dplyr::filter(status == "Approved") %>%
  dplyr::select(approved = symbol, alias_symbol, prev_symbol)

gene_map <- dplyr::bind_rows(
  gene_map %>% dplyr::transmute(old = approved, approved),
  gene_map %>% tidyr::pivot_longer(c(prev_symbol), values_to = "old") %>%
    tidyr::separate_rows(old, sep = "\\|") %>%
    dplyr::filter(!is.na(old), old != "", !grepl("^[0-9]+(\\.[0-9]+)?$", old)) %>%
    dplyr::transmute(old, approved)
) %>%
  dplyr::distinct(old, approved) %>%
  dplyr::arrange(old)

readr::write_tsv(gene_map, "./HGNC_symbol_mapping.tsv")
