# ============================================
# ENRIQUECIMIENTO FUNCIONAL (GO + KEGG)
# ============================================

#' Carga el OrgDb del organismo, instalandolo si es necesario
#' @keywords internal
cargar_orgdb <- function(org_info) {
  pkg <- org_info$orgdb
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Instalando ", pkg, " desde Bioconductor...")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  return(get(pkg, envir = asNamespace(pkg)))
}

#' Enriquecimiento KEGG para un cluster
#' @keywords internal
get_kegg_summary <- function(genes, org_info, top_n = 3) {
  if (length(genes) < 3) return("Cluster pequeno para KEGG")

  tryCatch({
    orgdb    <- cargar_orgdb(org_info)
    entrez   <- AnnotationDbi::mapIds(orgdb, keys = genes,
                                      column = "ENTREZID", keytype = "SYMBOL",
                                      multiVals = "first")
    entrez   <- na.omit(entrez)
    if (length(entrez) < 3) return("No significant KEGG pathways")

    kegg <- clusterProfiler::enrichKEGG(
      gene          = entrez,
      organism      = org_info$kegg,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )

    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      top <- head(kegg@result, top_n)
      return(paste0("<b>KEGG:</b><br>",
                    paste0(top$Description,
                           " (p=", formatC(top$pvalue, format = "e", digits = 1), ")",
                           collapse = "<br>")))
    }
  }, error = function(e) {})

  return("No significant KEGG pathways")
}

#' Enriquecimiento GO para un cluster
#' @keywords internal
get_go_summary <- function(genes, org_info, top_n = 2) {
  if (length(genes) < 3) return("Cluster pequeno para GO")

  orgdb    <- cargar_orgdb(org_info)
  keytype  <- if (org_info$keytype %in% c("ORF", "TAIR")) org_info$keytype else "SYMBOL"

  get_ont <- function(ont) {
    tryCatch({
      go <- clusterProfiler::enrichGO(
        gene          = genes,
        OrgDb         = orgdb,
        keyType       = keytype,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05
      )
      if (!is.null(go) && nrow(go@result) > 0) {
        top   <- head(go@result, top_n)
        terms <- paste0(top$Description,
                        " (p=", formatC(top$pvalue, format = "e", digits = 1), ")",
                        collapse = "<br>")
        return(paste0("<b>", ont, ":</b><br>", terms))
      }
    }, error = function(e) {})
    return(paste0("<b>", ont, ":</b> No significant"))
  }

  paste(get_ont("BP"), "<br>", get_ont("MF"), "<br>", get_ont("CC"))
}

#' Aplica enriquecimiento a todos los clusters de una subred
#' @keywords internal
aplicar_enriquecimiento <- function(resultado, org_info) {
  if (is.null(resultado)) return(NULL)

  subnet        <- resultado$subnet
  cluster_means <- resultado$cluster_means
  tooltips      <- list()

  for (cluster_id in unique(igraph::V(subnet)$cluster)) {
    genes_cluster <- igraph::V(subnet)$name[igraph::V(subnet)$cluster == cluster_id]

    if (length(genes_cluster) >= 3) {
      kegg_txt <- get_kegg_summary(genes_cluster, org_info)
      go_txt   <- get_go_summary(genes_cluster, org_info)

      cluster_order <- order(cluster_means, decreasing = TRUE)
      rank_pos      <- which(cluster_order == cluster_id)
      mean_val      <- round(cluster_means[cluster_id], 3)

      tooltip_html <- paste0(
        "<div style='max-width:340px;padding:14px;background:#fff;color:#333;",
        "border-radius:12px;box-shadow:0 4px 20px rgba(0,0,0,0.25);'>",
        "<h4 style='margin:0 0 8px 0;color:#2E4057;border-bottom:2px solid #eee;",
        "padding-bottom:6px;'>Rank ", rank_pos, " -- Cluster ", cluster_id, "</h4>",
        "<p style='margin:0 0 8px 0;color:#666;font-size:11px;'>",
        "|log2FC| medio: <b>", mean_val, "</b></p>",
        "<b>Genes (", length(genes_cluster), "):</b><br>",
        "<span style='font-size:11px;'>", paste(genes_cluster, collapse = ", "), "</span>",
        "<br><br>", kegg_txt, "<br><br>", go_txt, "</div>"
      )
    } else {
      tooltip_html <- paste0(
        "<div style='padding:10px;background:#fff;border-radius:8px;'>",
        "<h4 style='color:#666;'>Cluster ", cluster_id, "</h4>",
        "<b>Genes:</b> ", paste(genes_cluster, collapse = ", "),
        "<br><i style='color:#999;'>Cluster pequeno</i></div>"
      )
    }

    tooltip_html <- gsub('"', '\\"', tooltip_html, fixed = TRUE)
    tooltips[[as.character(cluster_id)]] <- tooltip_html
  }

  resultado$tooltips <- tooltips
  return(resultado)
}
