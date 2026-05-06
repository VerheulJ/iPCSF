# ============================================
# DESCARGA DE INTERACTOMA DESDE STRING API
# ============================================

#' Descarga el interactoma desde STRING API
#'
#' @param genes Vector de Gene Symbols.
#' @param org Organismo (ej. "rat", "human", "mouse"). Ver \code{ORGANISMS}.
#' @param score_threshold Score mínimo STRING (0-1000). Default 400.
#' @param cache_dir Carpeta para cachear el resultado. NULL = sin caché.
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

  # 1. Mapear genes a STRING IDs
  genes_str <- paste(genes, collapse = "%0d")
  url_map <- paste0(
    "https://string-db.org/api/json/get_string_ids?",
    "identifiers=", genes_str,
    "&species=", org_info$taxid,
    "&limit=1&echo_query=1"
  )

  resp_map <- httr::GET(url_map)
  httr::stop_for_status(resp_map)
  mapping <- jsonlite::fromJSON(
    httr::content(resp_map, "text", encoding = "UTF-8")
  )

  if (is.null(mapping) || nrow(mapping) == 0) {
    stop("No se encontraron genes en STRING para el organismo: ", org_info$nombre)
  }

  string_ids    <- unique(mapping$stringId)
  message("  ", length(string_ids), "/", length(genes), " genes mapeados en STRING")

  # 2. Obtener red de interacciones
  ids_str <- paste(string_ids, collapse = "%0d")
  url_net <- paste0(
    "https://string-db.org/api/json/network?",
    "identifiers=", ids_str,
    "&species=", org_info$taxid,
    "&required_score=", score_threshold,
    "&network_type=functional"
  )

  resp_net <- httr::GET(url_net)
  httr::stop_for_status(resp_net)
  network <- jsonlite::fromJSON(
    httr::content(resp_net, "text", encoding = "UTF-8")
  )

  if (is.null(network) || nrow(network) == 0) {
    stop("STRING no devolvió interacciones. Prueba a bajar score_threshold.")
  }

  # 3. Convertir a formato iPCSF
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
