# ============================================
# DESCARGA DE INTERACTOMA DESDE STRING API
# ============================================
#' Descarga el interactoma desde STRING API
#'
#' @param genes Vector de Gene Symbols.
#' @param org Organismo (ej. "rat", "human", "mouse"). Ver \code{ORGANISMS}.
#' @param score_threshold Score minimo STRING (0-1000). Default 400.
#' @param cache_dir Carpeta para cachear el resultado. NULL = sin cache.
#'
#' @return Data.frame con columnas \code{from}, \code{to}, \code{cost}.
#' @export
get_string_interactome <- function(genes,
                                   org = "human",
                                   score_threshold = 400,
                                   cache_dir = NULL) {
  org_info <- ORGANISMS[[org]]
  if (is.null(org_info)) {
    stop("Organismo '", org, "' no reconocido.\nOpciones: ",
         paste(names(ORGANISMS), collapse = ", "))
  }

  # Cache
  if (!is.null(cache_dir)) {
    cache_file <- file.path(cache_dir,
                            paste0("STRING_", org, "_", score_threshold, ".rds"))
    if (file.exists(cache_file)) {
      message("Cargando interactoma desde cache: ", cache_file)
      return(readRDS(cache_file))
    }
  }

  message("Consultando STRING API para ", org_info$nombre, "...")

  # 1. Interacciones entre los genes de entrada
  resp_network <- httr::POST(
    "https://string-db.org/api/json/network",
    body = list(
      identifiers     = paste(genes, collapse = "%0d"),
      species         = as.character(org_info$taxid),
      required_score  = as.character(score_threshold),
      network_type    = "functional",
      caller_identity = "iPCSF"
    ),
    encode = "form"
  )
  httr::stop_for_status(resp_network)
  network <- jsonlite::fromJSON(
    httr::content(resp_network, "text", encoding = "UTF-8")
  )

  # 2. Vecinos externos (potenciales Steiner nodes)
  resp_partners <- httr::POST(
    "https://string-db.org/api/json/interaction_partners",
    body = list(
      identifiers     = paste(genes, collapse = "%0d"),
      species         = as.character(org_info$taxid),
      required_score  = as.character(score_threshold),
      limit           = "5",
      caller_identity = "iPCSF"
    ),
    encode = "form"
  )
  httr::stop_for_status(resp_partners)
  partners <- jsonlite::fromJSON(
    httr::content(resp_partners, "text", encoding = "UTF-8")
  )

  # 3. Combinar ambas redes
  if (!is.null(partners) && length(partners) > 0 && nrow(partners) > 0) {
    network <- rbind(network, partners)
    network <- network[!duplicated(paste(network$preferredName_A,
                                         network$preferredName_B)), ]
  }

  if (is.null(network) || length(network) == 0 || nrow(network) == 0) {
    stop(
      "STRING no devolvio interacciones.\n",
      "Sugerencias:\n",
      "  - Prueba a bajar score_threshold (actual: ", score_threshold, ")\n",
      "  - Verifica que los gene symbols son correctos para ", org_info$nombre
    )
  }

  # Contar genes mapeados y posibles Steiner
  genes_encontrados <- unique(c(network$preferredName_A, network$preferredName_B))
  genes_steiner     <- setdiff(genes_encontrados, genes)
  message("  ", length(intersect(genes_encontrados, genes)), "/", length(genes),
          " genes mapeados en STRING")
  message("  ", length(genes_steiner), " posibles Steiner nodes disponibles")

  # Convertir a formato iPCSF: from, to, cost
  interactome <- data.frame(
    from = network$preferredName_A,
    to   = network$preferredName_B,
    cost = 1 - (network$score / 1000),
    stringsAsFactors = FALSE
  )

  interactome <- interactome[complete.cases(interactome), ]
  message("  Interactoma listo: ", nrow(interactome), " interacciones")

  # Guardar cache
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(interactome, cache_file)
    message("  Cache guardado en: ", cache_file)
  }

  return(interactome)
}
