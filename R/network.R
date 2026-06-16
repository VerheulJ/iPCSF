# ============================================
# CONSTRUCCION Y PROCESADO DE LA RED PCSF
# ============================================

#' Construye la subred PCSF a partir de una lista de genes desregulados
#'
#' @param desregulated_list Data.frame con genes, log2FC y p-valor.
#' @param org_info Lista de informacion del organismo (de \code{ORGANISMS}).
#' @param gene_col Columna con Gene Symbols.
#' @param log2fc_col Columna con log2 fold-change.
#' @param pval_col Columna con p-valor (-log10).
#' @param score_threshold Score minimo STRING (0-1000). Default 400.
#' @param cache_dir Carpeta de cache para el interactoma. NULL = sin cache.
#' @param w Parametro omega de PCSF. Default 0.85.
#' @param b Parametro beta de PCSF. Default 1.
#' @param mu Parametro mu de PCSF. Default 0.00005.
#'
#' @return Lista con \code{subnet}, \code{cluster_means}, \code{gene_scores}.
#' @noRd
construir_red <- function(desregulated_list,
                          org_info,
                          gene_col,
                          log2fc_col,
                          pval_col,
                          score_threshold = 400,
                          cache_dir       = NULL,
                          w  = 0.85,
                          b  = 1,
                          mu = 0.00005) {

  df <- desregulated_list

  # ── 1. Limpiar duplicados (misma logica que tu script original) ──
  if (gene_col %in% colnames(df)) {
    df <- df[!duplicated(df[[gene_col]]), ]
    df <- df[!is.na(df[[gene_col]]) & df[[gene_col]] != "", ]
  }

  # ── 2. Preparar terminals (prizes = |pval|) ──────────────────────
  genes     <- as.character(df[[gene_col]])
  terminals <- setNames(abs(as.numeric(df[[pval_col]])), genes)
  terminals <- terminals[!is.na(terminals) & terminals > 0]

  # ── 3. Preparar gene scores (|log2FC| para ranking clusters) ─────
  gene_scores <- setNames(abs(as.numeric(df[[log2fc_col]])), genes)
  gene_scores <- gene_scores[!is.na(gene_scores) & !is.na(names(gene_scores))]

  message("  Genes terminals : ", length(terminals))

  if (length(terminals) < 3) {
    stop("Se necesitan al menos 3 genes validos. ",
         "Revisa las columnas '", pval_col, "' y '", gene_col, "'.")
  }

  # ── 4. Interactoma desde STRING ───────────────────────────────────
  # Buscar el key del organismo en ORGANISMS para pasarlo a get_string_interactome
  org_key <- names(ORGANISMS)[
    sapply(ORGANISMS, function(x) identical(x$taxid, org_info$taxid))
  ]
  if (length(org_key) == 0) stop("No se pudo identificar el organismo en ORGANISMS.")

  interactome <- get_string_interactome(
    genes           = names(terminals),
    org             = org_key,
    score_threshold = score_threshold,
    cache_dir       = cache_dir
  )

  # ── 5. Intersectar genes con STRING ──────────────────────────────
  string_graph  <- igraph::graph_from_data_frame(
    interactome[, c("from", "to")], directed = FALSE
  )
  string_graph  <- igraph::simplify(string_graph)
  genes_comunes <- intersect(names(terminals), igraph::V(string_graph)$name)

  message("  Genes en STRING : ", length(genes_comunes), "/", length(terminals))

  if (length(genes_comunes) < 3) {
    stop(
      "Menos de 3 genes encontrados en STRING.\n",
      "Sugerencias:\n",
      "  - Verifica que los nombres de genes son Gene Symbols validos\n",
      "  - Prueba a bajar score_threshold (actual: ", score_threshold, ")\n",
      "  - Comprueba que el organismo es correcto: ", org_info$nombre
    )
  }

  terminals <- terminals[genes_comunes]

  # ── 6. Algoritmo PCSF ─────────────────────────────────────────────
  message("  Ejecutando PCSF...")
  ppi    <- PCSF::construct_interactome(interactome)
  subnet <- PCSF::PCSF(ppi = ppi, terminals = terminals,
                       w = w, b = b, mu = mu)

  if (igraph::vcount(subnet) < 2) {
    stop(
      "La subred resultante tiene menos de 2 nodos.\n",
      "Prueba a ajustar los parametros w, b o mu."
    )
  }

  # ── 7. Anotar nodos ───────────────────────────────────────────────

  # Tipo de nodo: Terminal vs Steiner
  igraph::V(subnet)$tipo_nodo <- ifelse(
    igraph::V(subnet)$name %in% names(terminals),
    "Terminal", "Steiner"
  )

  # Prize original
  igraph::V(subnet)$original_prize <- ifelse(
    igraph::V(subnet)$name %in% names(terminals),
    terminals[igraph::V(subnet)$name], 0
  )

  # Status de regulacion — generico: busca columna "Status" si existe
  if ("Status" %in% colnames(df)) {
    status_map <- setNames(as.character(df$Status), df[[gene_col]])
    igraph::V(subnet)$status <- ifelse(
      igraph::V(subnet)$name %in% names(status_map),
      status_map[igraph::V(subnet)$name],
      "Steiner"
    )
  } else {
    # Si no hay columna Status, inferir desde log2FC
    fc_map <- setNames(as.numeric(df[[log2fc_col]]), df[[gene_col]])
    igraph::V(subnet)$status <- sapply(igraph::V(subnet)$name, function(g) {
      if (!g %in% names(fc_map)) return("Steiner")
      fc <- fc_map[[g]]
      if (is.na(fc))  return("Steiner")
      if (fc > 0)     return("up-regulated")
      return("down-regulated")
    })
  }

  # ── 8. Clustering ─────────────────────────────────────────────────
  message("  Clustering...")
clusters <- igraph::cluster_louvain(subnet)
  igraph::V(subnet)$cluster <- clusters$membership

  genes_by_cluster <- split(igraph::V(subnet)$name, clusters$membership)
  cluster_means    <- sapply(genes_by_cluster, function(g) {
    scores <- gene_scores[names(gene_scores) %in% g]
    if (length(scores) > 0) mean(scores, na.rm = TRUE) else 0
  })

  message("  Red lista       : ",
          igraph::vcount(subnet), " nodos, ",
          igraph::ecount(subnet), " aristas, ",
          length(unique(igraph::V(subnet)$cluster)), " clusters")

  return(list(
    subnet        = subnet,
    cluster_means = cluster_means,
    gene_scores   = gene_scores
  ))
}


#' @importFrom grDevices rainbow
#' @importFrom stats complete.cases na.omit setNames
#' @importFrom utils head
NULL
