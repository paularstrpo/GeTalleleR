library(tidyverse)
library(devtools)
devtools::load_all()
data <- read_tsv('GeTallele/data/known_rs_013readcounts.tsv.zip')

test_vafs <- data %>%
    dplyr::select(CHROM, POS, REF, ALTA, AlignedReads, SNVCount, RefCount) %>%
    mutate(AlignedReads = substr(AlignedReads, 1, 29),
           VAF=get_vaf(ref=RefCount, alt=SNVCount)) %>%
    dplyr::select(-SNVCount, -RefCount) %>%
    pivot_wider(names_from=AlignedReads, values_from=VAF, values_fn=median) %>%
    data.frame()

rownames(test_vafs) <- paste(paste0('chr', test_vafs$CHROM), test_vafs$POS, test_vafs$REF, test_vafs$ALTA, sep='_')

layer_map <- colnames(test_vafs)[startsWith(colnames(test_vafs), 'X013')]
names(layer_map) <- substr(layer_map, 6, 8)

layer_map


test_reads <- data %>%
    dplyr::select(CHROM, POS, REF, ALTA, AlignedReads, SNVCount, RefCount) %>%
    mutate(AlignedReads = substr(AlignedReads, 1, 29),
           TotalReads=RefCount+SNVCount) %>%
    dplyr::select(-SNVCount, -RefCount) %>%
    pivot_wider(names_from=AlignedReads, values_from=TotalReads, values_fn=median) %>%
    data.frame()

rownames(test_reads) <- paste(paste0('chr', test_reads$CHROM), test_reads$POS, test_reads$REF, test_reads$ALTA, sep='_')

test_reads1 <- test_reads$X013_Ntr_BRCA_TCGA.BH.A0B8.11A
