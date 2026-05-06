# ============================================
# ORGANISMOS SOPORTADOS
# ============================================

ORGANISMS <- list(
  # Mamíferos
  human       = list(taxid = 9606,   orgdb = "org.Hs.eg.db",      kegg = "hsa", keytype = "SYMBOL", nombre = "Homo sapiens"),
  mouse       = list(taxid = 10090,  orgdb = "org.Mm.eg.db",      kegg = "mmu", keytype = "SYMBOL", nombre = "Mus musculus"),
  rat         = list(taxid = 10116,  orgdb = "org.Rn.eg.db",      kegg = "rno", keytype = "SYMBOL", nombre = "Rattus norvegicus"),
  bovine      = list(taxid = 9913,   orgdb = "org.Bt.eg.db",      kegg = "bta", keytype = "SYMBOL", nombre = "Bos taurus"),
  canine      = list(taxid = 9615,   orgdb = "org.Cf.eg.db",      kegg = "cfa", keytype = "SYMBOL", nombre = "Canis lupus familiaris"),
  pig         = list(taxid = 9823,   orgdb = "org.Ss.eg.db",      kegg = "ssc", keytype = "SYMBOL", nombre = "Sus scrofa"),
  rhesus      = list(taxid = 9544,   orgdb = "org.Mmu.eg.db",     kegg = "mcc", keytype = "SYMBOL", nombre = "Macaca mulatta"),
  chimp       = list(taxid = 9598,   orgdb = "org.Pt.eg.db",      kegg = "ptr", keytype = "SYMBOL", nombre = "Pan troglodytes"),
  chicken     = list(taxid = 9031,   orgdb = "org.Gg.eg.db",      kegg = "gga", keytype = "SYMBOL", nombre = "Gallus gallus"),
  xenopus     = list(taxid = 8364,   orgdb = "org.Xl.eg.db",      kegg = "xla", keytype = "SYMBOL", nombre = "Xenopus laevis"),
  # Modelo
  zebrafish   = list(taxid = 7955,   orgdb = "org.Dr.eg.db",      kegg = "dre", keytype = "SYMBOL", nombre = "Danio rerio"),
  fly         = list(taxid = 7227,   orgdb = "org.Dm.eg.db",      kegg = "dme", keytype = "SYMBOL", nombre = "Drosophila melanogaster"),
  worm        = list(taxid = 6239,   orgdb = "org.Ce.eg.db",      kegg = "cel", keytype = "SYMBOL", nombre = "Caenorhabditis elegans"),
  yeast       = list(taxid = 4932,   orgdb = "org.Sc.sgd.db",     kegg = "sce", keytype = "ORF",    nombre = "Saccharomyces cerevisiae"),
  # Insecto
  mosquito    = list(taxid = 7165,   orgdb = "org.Ag.eg.db",      kegg = "aga", keytype = "SYMBOL", nombre = "Anopheles gambiae"),
  # Planta
  arabidopsis = list(taxid = 3702,   orgdb = "org.At.tair.db",    kegg = "ath", keytype = "TAIR",   nombre = "Arabidopsis thaliana"),
  # Bacteria
  ecoli_k12   = list(taxid = 83333,  orgdb = "org.EcK12.eg.db",   kegg = "eco", keytype = "SYMBOL", nombre = "Escherichia coli K12"),
  ecoli_sakai = list(taxid = 160492, orgdb = "org.EcSakai.eg.db", kegg = "ecs", keytype = "SYMBOL", nombre = "Escherichia coli Sakai"),
  # Parásito
  malaria     = list(taxid = 36329,  orgdb = "org.Pf.plasmo.db",  kegg = "pfa", keytype = "ORF",    nombre = "Plasmodium falciparum")
)

#' @importFrom grDevices rainbow
#' @importFrom stats complete.cases na.omit setNames
#' @importFrom utils head
NULL

