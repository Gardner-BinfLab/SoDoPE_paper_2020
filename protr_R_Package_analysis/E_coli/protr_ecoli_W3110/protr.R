library("protr")
library("dplyr")


ecoli <- readFASTA("ecoli_W3110_tags.faa")
# Amino acid composition
AAC <- t(sapply(ecoli, extractAAC)) %>% as_tibble(rownames = "Accession") # Amino acid composition
write.table(AAC, file = "protr_AAC_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
DC <- t(sapply(ecoli, extractDC)) %>% as_tibble(rownames = "Accession") # Dipeptide composition
write.table(DC, file = "protr_DC_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
TC <- t(sapply(ecoli, extractTC)) %>% as_tibble(rownames = "Accession") # Tripeptide composition
write.table(TC, file = "protr_TC_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# Autocorrelation
MoreauBroto <- t(sapply(ecoli, extractMoreauBroto)) %>% as_tibble(rownames = "Accession") # Normalized Moreau-Broto autocorrelation
write.table(MoreauBroto, file = "protr_MoreauBroto_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
Moran <- t(sapply(ecoli, extractMoran)) %>% as_tibble(rownames = "Accession") # Moran autocorrelation
write.table(Moran, file = "protr_Moran_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
Geary <- t(sapply(ecoli, extractGeary)) %>% as_tibble(rownames = "Accession") # Geary autocorrelation
write.table(Geary, file = "protr_Geary_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# CTD descriptors
CTDC <- t(sapply(ecoli, extractCTDC)) %>% as_tibble(rownames = "Accession") # Composition
write.table(CTDC, file = "protr_CTDC_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
CTDT <- t(sapply(ecoli, extractCTDT)) %>% as_tibble(rownames = "Accession") # Transition
write.table(CTDT, file = "protr_CTDT_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
CTDD <- t(sapply(ecoli, extractCTDD)) %>% as_tibble(rownames = "Accession") # Distribution
write.table(CTDD, file = "protr_CTDD_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# Conjoint triad descriptors
CTriad <- t(sapply(ecoli, extractCTriad)) %>% as_tibble(rownames = "Accession") # Conjoint triad descriptors
write.table(CTriad, file = "protr_CTriad_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# Quasi-sequence-order descriptors
SOCN <- t(sapply(ecoli, extractSOCN)) %>% as_tibble(rownames = "Accession") # Sequence-order-coupling number
write.table(SOCN, file = "protr_SOCN_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
QSO <- t(sapply(ecoli, extractQSO)) %>% as_tibble(rownames = "Accession") # Quasi-sequence-order descriptors
write.table(QSO, file = "protr_QSO_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# Pseudo-amino acid composition
PAAC <- t(sapply(ecoli, extractPAAC)) %>% as_tibble(rownames = "Accession") # Pseudo-amino acid composition (PseAAC)
write.table(PAAC, file = "protr_PAAC_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
APAAC <- t(sapply(ecoli, extractAPAAC)) %>% as_tibble(rownames = "Accession") # Amphiphilic pseudo-amino acid composition (APseAAC)
write.table(APAAC, file = "protr_APAAC_ecoli_W3110.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)