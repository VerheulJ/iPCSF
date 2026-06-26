# ============================================
# NETWORK CONSTRUCTION -- PCSF ALGORITHM
# ============================================
#
# The following two functions (construct_interactome and pcsf_algorithm)
# are adapted from the PCSF R package by Murodzhon Akhmedov et al. (2017),
# originally licensed under MIT License.
#
# Original package: https://github.com/IOR-Bioinformatics/PCSF
# Citation: Akhmedov M, et al. (2017). PCSF: An R-package for network-based
#   interpretation of high-throughput data. PLOS Computational Biology.
# License: MIT © 2017 Murodzhon Akhmedov
# ============================================

#' Build igraph interactome from edge data.frame
#' Adapted from PCSF::construct_interactome (Akhmedov et al. 2017, MIT License)
#' @noRd
construct_interactome <- function(ppi) {
  if (missing(ppi))
    stop("Need to specify a list of edges to construct an interaction network.")
  if (nrow(ppi) < 1 || ncol(ppi) != 3 || !is.data.frame(ppi))
    stop("Need to provide a data.frame with three columns: from, to, cost.")

  node_names  <- unique(c(as.character(ppi[,1]), as.character(ppi[,2])))
  ppi.graph   <- igraph::graph.data.frame(ppi[,1:2], vertices = node_names, directed = FALSE)
  igraph::E(ppi.graph)$weight <- as.numeric(ppi[,3])
  ppi.graph   <- igraph::simplify(ppi.graph)
  return(ppi.graph)
}

#' Prize-Collecting Steiner Forest algorithm
#' Adapted from PCSF::PCSF (Akhmedov et al. 2017, MIT License)
#' @noRd
pcsf_algorithm <- function(ppi, terminals, w = 2, b = 1, mu = 0.0005) {
  if (missing(ppi) || !inherits(ppi, "igraph"))
    stop("Need to specify an igraph interaction network.")
  if (missing(terminals) || is.null(names(terminals)))
    stop("terminals must be a named numeric vector.")

  terminal_names  <- names(terminals)
  terminal_values <- as.numeric(terminals)

  node_names <- igraph::V(ppi)$name
  node_prz   <- vector(mode = "numeric", length = length(node_names))
  index      <- match(terminal_names, node_names)
  percent    <- signif((length(index) - sum(is.na(index))) / length(index) * 100, 4)

  if (percent < 5)
    stop("Less than 5% of terminal nodes matched in the interactome. Check your gene symbols.")

  message("  ", percent, "% of terminal nodes matched in the interactome")

  terminal_names  <- terminal_names[!is.na(index)]
  terminal_values <- terminal_values[!is.na(index)]
  index           <- index[!is.na(index)]
  node_prz[index] <- terminal_values

  dummies <- terminal_names

  message("  Solving PCSF...")

  node_degrees    <- igraph::degree(ppi)
  hub_penalization <- -mu * node_degrees
  node_prizes     <- b * node_prz
  zero_idx        <- which(node_prizes == 0)
  node_prizes[zero_idx] <- hub_penalization[zero_idx]

  edges <- igraph::ends(ppi, es = igraph::E(ppi))
  from  <- c(rep("DUMMY", length(dummies)), edges[,1])
  to    <- c(dummies, edges[,2])
  cost  <- c(rep(w, length(dummies)), igraph::E(ppi)$weight)

  output <- call_sr(from, to, cost, node_names, node_prizes)

  if (length(output[[1]]) == 0)
    stop("No subnetwork identified. Try adjusting w, b, or mu parameters.")

  e        <- data.frame(output[[1]], output[[2]], output[[3]])
  names(e) <- c("from", "to", "weight")

  type     <- rep("Steiner", length(output[[4]]))
  idx      <- match(terminal_names, output[[4]])
  idx      <- idx[!is.na(idx)]
  type[idx] <- "Terminal"

  v        <- data.frame(output[[4]], output[[5]], type)
  names(v) <- c("terminals", "prize", "type")

  subnet <- igraph::graph.data.frame(e, vertices = v, directed = FALSE)
  igraph::E(subnet)$weight <- as.numeric(output[[3]])
  subnet <- igraph::delete_vertices(subnet, "DUMMY")
  subnet <- igraph::delete_vertices(subnet, names(which(igraph::degree(subnet) == 0)))

  return(subnet)
}

#' Build PCSF network from differential expression data
#' @noRd
construir_red <- function(desregulated_list,
                           org_info,
                           gene_col        = "gene",
                           log2fc_col      = "log2FC",
                           pval_col        = "pvalue",
                           score_threshold = 400,
                           cache_dir       = NULL,
                           w  = 2,
                           b  = 1,
                           mu = 0.0005) {

  df <- desregulated_list

  # ── 1. Validate columns ───────────────────────────────────────────
  for (col in c(gene_col, log2fc_col, pval_col)) {
    if (!col %in% colnames(df))
      stop("Column '", col, "' not found. Available: ",
           paste(colnames(df), collapse = ", "))
  }

  # ── 2. Clean data ─────────────────────────────────────────────────
  df <- df[!is.na(df[[gene_col]]) & df[[gene_col]] != "", ]
  df <- df[!duplicated(df[[gene_col]]), ]

  genes <- as.character(df[[gene_col]])

  # ── 3. Gene scores (prizes = |-log10(pval)|) ──────────────────────
  gene_scores <- setNames(abs(as.numeric(df[[pval_col]])), genes)
  gene_scores <- gene_scores[!is.na(gene_scores) & gene_scores > 0]

  if (length(gene_scores) < 5)
    stop("Less than 5 genes with valid prizes. Check pval_col values.")

  # ── 4. Get organism key ───────────────────────────────────────────
  org_key <- names(ORGANISMS)[
    sapply(ORGANISMS, function(x) identical(x$taxid, org_info$taxid))
  ]

  # ── 5. Download interactome ───────────────────────────────────────
  message("  Downloading interactome...")
  interactome <- get_string_interactome(
    genes           = genes,
    org             = org_key,
    score_threshold = score_threshold,
    cache_dir       = cache_dir
  )

  # ── 6. Match terminals to interactome ────────────────────────────
  genes_comunes <- intersect(names(gene_scores),
                             unique(c(interactome$from, interactome$to)))
  if (length(genes_comunes) < 3)
    stop("Less than 3 terminal genes found in the interactome. ",
         "Try lowering score_threshold.")

  terminals <- gene_scores[genes_comunes]

  # ── 7. PCSF algorithm ────────────────────────────────────────────
  ppi    <- construct_interactome(interactome)
  subnet <- pcsf_algorithm(ppi = ppi, terminals = terminals,
                            w = w, b = b, mu = mu)

  if (igraph::vcount(subnet) < 2)
    stop("Resulting subnetwork has fewer than 2 nodes. ",
         "Try adjusting w, b or mu parameters.")

  # ── 8. Annotate nodes ────────────────────────────────────────────
  # Node type: Terminal vs Steiner
  igraph::V(subnet)$tipo_nodo <- ifelse(
    igraph::V(subnet)$name %in% names(terminals),
    "Terminal", "Steiner"
  )

  # Original prize
  igraph::V(subnet)$original_prize <- ifelse(
    igraph::V(subnet)$name %in% names(terminals),
    terminals[igraph::V(subnet)$name], 0
  )

  # Regulation status inferred from log2FC sign
  fc_map <- setNames(as.numeric(df[[log2fc_col]]), genes)
  igraph::V(subnet)$status <- sapply(igraph::V(subnet)$name, function(g) {
    if (!g %in% names(fc_map)) return("Steiner")
    fc <- fc_map[[g]]
    if (is.na(fc))  return("Steiner")
    if (fc > 0)     return("up-regulated")
    return("down-regulated")
  })

  # ── 9. Clustering ─────────────────────────────────────────────────
  message("  Clustering...")
  clusters                  <- igraph::cluster_louvain(subnet)
  igraph::V(subnet)$cluster <- clusters$membership
  genes_by_cluster          <- split(igraph::V(subnet)$name, clusters$membership)
  cluster_means             <- sapply(genes_by_cluster, function(g) {
    scores <- gene_scores[names(gene_scores) %in% g]
    if (length(scores) > 0) mean(scores, na.rm = TRUE) else 0
  })

  message("  Network ready: ",
          igraph::vcount(subnet), " nodes, ",
          igraph::ecount(subnet), " edges, ",
          length(unique(igraph::V(subnet)$cluster)), " clusters")

  return(list(
    subnet        = subnet,
    cluster_means = cluster_means,
    gene_scores   = gene_scores
  ))
}
