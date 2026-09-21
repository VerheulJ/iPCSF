# ============================================
# FUNCTIONAL CLUSTERING
# GO SEMANTIC SIMILARITY
# ============================================

#' Load the organism annotation database
#' @noRd
cargar_orgdb <- function(org_info) {

  pkg <- org_info$orgdb

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Install it with: BiocManager::install('", pkg, "')"
    )
  }

  get(pkg, envir = asNamespace(pkg))
}


#' Calculate GO semantic similarity between proteins
#' @noRd
calcular_similitud_funcional <- function(
    genes,
    org_info,
    ontology = "BP",
    measure = "Wang",
    combine = "BMA"
) {

  if (length(genes) < 2) {
    stop("At least two proteins are required.")
  }

  if (!requireNamespace("GOSemSim", quietly = TRUE)) {
    stop(
      "Package 'GOSemSim' is required.\n",
      "Install it with: BiocManager::install('GOSemSim')"
    )
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(
      "Package 'AnnotationDbi' is required.\n",
      "Install it with: BiocManager::install('AnnotationDbi')"
    )
  }

  orgdb <- cargar_orgdb(org_info)

  # ------------------------------------------
  # Identifier type
  # ------------------------------------------

  keytype <- if (
    org_info$keytype %in% c("ORF", "TAIR")
  ) {
    org_info$keytype
  } else {
    "SYMBOL"
  }

  # ------------------------------------------
  # Map input identifiers to ENTREZID
  # ------------------------------------------

  mapping <- AnnotationDbi::mapIds(
    orgdb,
    keys = unique(genes),
    column = "ENTREZID",
    keytype = keytype,
    multiVals = "first"
  )

  mapping <- mapping[!is.na(mapping)]
  mapping <- mapping[!duplicated(mapping)]

  if (length(mapping) < 2) {
    stop(
      "Fewer than two proteins could be mapped to ENTREZID."
    )
  }

  entrez <- unname(mapping)

  # ------------------------------------------
  # GO semantic data
  # ------------------------------------------

  sem_data <- GOSemSim::godata(
    OrgDb = orgdb,
    ont = ontology,
    computeIC = TRUE
  )

  # ------------------------------------------
  # Gene-gene semantic similarity
  # ------------------------------------------

  similarity <- GOSemSim::mgeneSim(
    entrez,
    semData = sem_data,
    measure = measure,
    combine = combine
  )

  similarity <- as.matrix(similarity)

  # ------------------------------------------
  # Replace Entrez IDs by original gene names
  # ------------------------------------------

  gene_names <- names(mapping)

  row_names <- rownames(similarity)

  translated_names <- gene_names[
    match(row_names, unname(mapping))
  ]

  translated_names[
    is.na(translated_names)
  ] <- row_names[
    is.na(translated_names)
  ]

  rownames(similarity) <- translated_names
  colnames(similarity) <- translated_names

  # ------------------------------------------
  # Return
  # ------------------------------------------

  list(
    similarity = similarity,
    genes = translated_names,
    mapping = mapping,
    ontology = ontology,
    measure = measure,
    combine = combine
  )
}


#' Cluster proteins according to functional similarity
#' @noRd
cluster_funcional <- function(
    similarity,
    method = "average",
    n_clusters = NULL
) {

  similarity <- as.matrix(similarity)

  similarity[is.na(similarity)] <- 0

  similarity <- (
    similarity + t(similarity)
  ) / 2

  diag(similarity) <- 1

  # Functional distance
  distance <- 1 - similarity

  distance[distance < 0] <- 0
  diag(distance) <- 0

  distance <- stats::as.dist(distance)

  # Hierarchical clustering
  hc <- stats::hclust(
    distance,
    method = method
  )

  # Number of clusters
  if (is.null(n_clusters)) {

    # First version:
    # use a data-driven cut based on the dendrogram.
    #
    # We will refine this criterion after seeing
    # real iPCSF networks.

    n_clusters <- max(
      2,
      round(sqrt(nrow(similarity)))
    )

    n_clusters <- min(
      n_clusters,
      nrow(similarity)
    )
  }

  clusters <- stats::cutree(
    hc,
    k = n_clusters
  )

  list(
    clusters = clusters,
    hclust = hc,
    distance = distance
  )
}


#' Perform functional clustering of an iPCSF subnet
#' @noRd
hacer_clustering_funcional <- function(
    genes,
    org_info,
    ontology = "BP",
    measure = "Wang",
    combine = "BMA",
    cluster_method = "average",
    n_clusters = NULL
) {

  # ------------------------------------------
  # Functional similarity
  # ------------------------------------------

  functional <- calcular_similitud_funcional(
    genes = genes,
    org_info = org_info,
    ontology = ontology,
    measure = measure,
    combine = combine
  )

  # ------------------------------------------
  # Clustering
  # ------------------------------------------

  clustering <- cluster_funcional(
    similarity = functional$similarity,
    method = cluster_method,
    n_clusters = n_clusters
  )

  # ------------------------------------------
  # Return
  # ------------------------------------------

  list(
    similarity = functional$similarity,
    clusters = clustering$clusters,
    hclust = clustering$hclust,
    distance = clustering$distance,
    genes = functional$genes,
    mapping = functional$mapping,
    ontology = ontology,
    measure = measure,
    combine = combine
  )
}
