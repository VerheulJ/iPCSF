
# iPCSF <img src="man/figures/logo.png" align="right" height="139" alt="" />

> Improved Prize-Collecting Steiner Forest with automatic interactome and modern network visualization

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)


---

## What is iPCSF?

iPCSF is an R package that extends the original [PCSF](https://bioconductor.org/packages/PCSF/) algorithm with three major improvements:

| Feature | PCSF (original) | **iPCSF** |
|---|---|---|
| Interactome | old version only for human | **Automatic via STRING API** |
| Organisms supported | Human (limited) | **20+ organisms** |
| Visualization | Static | **Interactive HTML with clusters** |
| Enrichment | Basic | **GO (BP/MF/CC) + KEGG per cluster** |
| Multi-condition | No | **Yes — any number of conditions** |

---

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install iPCSF from GitHub
devtools::install_github("VerheulJ/iPCSF")
```

### Bioconductor dependencies

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("PCSF", "clusterProfiler", "AnnotationDbi"))
```

The organism annotation packages (e.g. `org.Rn.eg.db`, `org.Hs.eg.db`) are
installed automatically on first use — you don't need to install them manually.

---

## Quick start

## How it works

iPCSF()
│
├─ get_string_interactome()   # Query STRING API → from/to/cost table
│
├─ construir_red()            # PCSF algorithm → igraph subnet
│    ├─ terminals  = |−log10(p)|
│    └─ clustering = edge betweenness
│
├─ aplicar_enriquecimiento()  # GO + KEGG per cluster
│    ├─ enrichGO()  (BP, MF, CC)
│    └─ enrichKEGG()
│
└─ generar_html()             # Self-contained interactive HTML

---
### Minimum example

```r
library(iPCSF)

iPCSF(
  conditions = list(
    treatment = list(data = my_data, label = "Treatment", color = "#E41A1C")
  ),
  org        = "rat",
  gene_col   = "Gene.Symbol",
  log2fc_col = "log2FC",
  pval_col   = "neg_log10_pval"
)
```

### Two conditions

```r
iPCSF(
  conditions = list(
    females = list(data = df_females, label = "Females", color = "#E41A1C"),
    males   = list(data = df_males,   label = "Males",   color = "#377EB8")
  ),
  org             = "rat",
  gene_col        = "Gene.Symbol",
  log2fc_col      = "log2cociente",
  pval_col        = "log10pvalor",
  score_threshold = 400,
  cache_dir       = "~/.iPCSF_cache",
  output_file     = "my_network.html"
)
```

### Any number of conditions

```r
iPCSF(
  conditions = list(
    ctrl    = list(data = df_ctrl,    label = "Control",     color = "#999999"),
    dose_lo = list(data = df_low,     label = "Low dose",    color = "#4DAF4A"),
    dose_hi = list(data = df_high,    label = "High dose",   color = "#E41A1C"),
    rescue  = list(data = df_rescue,  label = "Rescue",      color = "#377EB8")
  ),
  org = "mouse"
)
```

---

## Input data format

`data` inside each condition must be a `data.frame` with at least three columns:

| Column | Description | Example name |
|---|---|---|
| Gene symbols | HGNC/MGI gene names | `Gene.Symbol`, `gene`, `symbol` |
| log2 fold-change | Numeric | `log2FC`, `log2cociente` |
| −log10(p-value) | Numeric, positive | `neg_log10_pval`, `log10pvalor` |


---

## Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `conditions` | named list | — | List of conditions (see above) |
| `org` | string | `"human"` | Organism key (see table below) |
| `gene_col` | string | `"gene"` | Column name for gene symbols |
| `log2fc_col` | string | `"log2FC"` | Column name for log2 fold-change |
| `pval_col` | string | `"pvalue"` | Column name for −log10(p-value) |
| `score_threshold` | integer | `400` | Minimum STRING interaction score (0–1000) |
| `cache_dir` | string | `NULL` | Folder to cache STRING interactome |
| `output_file` | string | `"iPCSF_network.html"` | Output HTML file name |
| `w` | numeric | `0.85` | PCSF omega parameter |
| `b` | numeric | `1` | PCSF beta parameter |
| `mu` | numeric | `0.00005` | PCSF mu parameter |

---

## Supported organisms

| Key | Species | OrgDb |
|---|---|---|
| `"human"` | *Homo sapiens* | org.Hs.eg.db |
| `"mouse"` | *Mus musculus* | org.Mm.eg.db |
| `"rat"` | *Rattus norvegicus* | org.Rn.eg.db |
| `"bovine"` | *Bos taurus* | org.Bt.eg.db |
| `"canine"` | *Canis lupus familiaris* | org.Cf.eg.db |
| `"pig"` | *Sus scrofa* | org.Ss.eg.db |
| `"rhesus"` | *Macaca mulatta* | org.Mmu.eg.db |
| `"chimp"` | *Pan troglodytes* | org.Pt.eg.db |
| `"chicken"` | *Gallus gallus* | org.Gg.eg.db |
| `"xenopus"` | *Xenopus laevis* | org.Xl.eg.db |
| `"zebrafish"` | *Danio rerio* | org.Dr.eg.db |
| `"fly"` | *Drosophila melanogaster* | org.Dm.eg.db |
| `"worm"` | *Caenorhabditis elegans* | org.Ce.eg.db |
| `"yeast"` | *Saccharomyces cerevisiae* | org.Sc.sgd.db |
| `"mosquito"` | *Anopheles gambiae* | org.Ag.eg.db |
| `"arabidopsis"` | *Arabidopsis thaliana* | org.At.tair.db |
| `"ecoli_k12"` | *Escherichia coli* K12 | org.EcK12.eg.db |
| `"ecoli_sakai"` | *Escherichia coli* Sakai | org.EcSakai.eg.db |
| `"malaria"` | *Plasmodium falciparum* | org.Pf.plasmo.db |

---

## Output

iPCSF generates a self-contained **interactive HTML file** that includes:

- **Network visualization** with nodes sized by −log10(p-value) and colored by cluster
- **Condition switcher** — one button per condition, dynamically generated
- **Cluster panel** — ranked by mean |log2FC|, click to highlight
- **Tooltips** — hover any node to see its cluster's GO and KEGG enrichment
- **Node shapes** — circles = terminal (differential), triangles = Steiner (connector)
- **Node borders** — red = up-regulated, blue = down-regulated, grey = Steiner

---


## Citation

If you use iPCSF in your research, please cite:
Verheul-Campos, J. (2026). iPCSF: Improved Prize-Collecting Steiner Forest
with automatic interactome and modern network visualization.
R package version 0.0.1. https://github.com/VerheulJ/iPCSF
And the original PCSF algorithm:
Akhmedov M, et al. (2017). PCSF: An R-package for network-based
interpretation of high-throughput data. PLOS Computational Biology.

---

## License

MIT © Julia Verheul-Campos
=======
# iPCSF
Improved Prize-Collecting Steiner Forest with automatic interactome and modern network visualization
>>>>>>> 11add726e6522accc0fe9cb15b715c3f742f01ee
