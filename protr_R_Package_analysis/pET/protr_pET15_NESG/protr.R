library("protr")
library("dplyr")


pET15 <- readFASTA("pET15_NESG.faa")
# Amino acid composition
AAC <- t(sapply(pET15, extractAAC)) %>% as_tibble(rownames = "Accession") # Amino acid composition
write.table(AAC, file = "protr_AAC_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
DC <- t(sapply(pET15, extractDC)) %>% as_tibble(rownames = "Accession") # Dipeptide composition
write.table(DC, file = "protr_DC_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
TC <- t(sapply(pET15, extractTC)) %>% as_tibble(rownames = "Accession") # Tripeptide composition
write.table(TC, file = "protr_TC_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# Autocorrelation
MoreauBroto <- t(sapply(pET15, extractMoreauBroto)) %>% as_tibble(rownames = "Accession") # Normalized Moreau-Broto autocorrelation
write.table(MoreauBroto, file = "protr_MoreauBroto_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
Moran <- t(sapply(pET15, extractMoran)) %>% as_tibble(rownames = "Accession") # Moran autocorrelation
write.table(Moran, file = "protr_Moran_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
Geary <- t(sapply(pET15, extractGeary)) %>% as_tibble(rownames = "Accession") # Geary autocorrelation
write.table(Geary, file = "protr_Geary_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# CTD descriptors
CTDC <- t(sapply(pET15, extractCTDC)) %>% as_tibble(rownames = "Accession") # Composition
write.table(CTDC, file = "protr_CTDC_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
CTDT <- t(sapply(pET15, extractCTDT)) %>% as_tibble(rownames = "Accession") # Transition
write.table(CTDT, file = "protr_CTDT_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
CTDD <- t(sapply(pET15, extractCTDD)) %>% as_tibble(rownames = "Accession") # Distribution
write.table(CTDD, file = "protr_CTDD_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# Conjoint triad descriptors
CTriad <- t(sapply(pET15, extractCTriad)) %>% as_tibble(rownames = "Accession") # Conjoint triad descriptors
write.table(CTriad, file = "protr_CTriad_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# Quasi-sequence-order descriptors
SOCN <- t(sapply(pET15, extractSOCN)) %>% as_tibble(rownames = "Accession") # Sequence-order-coupling number
write.table(SOCN, file = "protr_SOCN_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
QSO <- t(sapply(pET15, extractQSO)) %>% as_tibble(rownames = "Accession") # Quasi-sequence-order descriptors
write.table(QSO, file = "protr_QSO_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
# Pseudo-amino acid composition
PAAC <- t(sapply(pET15, extractPAAC)) %>% as_tibble(rownames = "Accession") # Pseudo-amino acid composition (PseAAC)
write.table(PAAC, file = "protr_PAAC_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)
APAAC <- t(sapply(pET15, extractAPAAC)) %>% as_tibble(rownames = "Accession") # Amphiphilic pseudo-amino acid composition (APseAAC)
write.table(APAAC, file = "protr_APAAC_pET15_NESG.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)


