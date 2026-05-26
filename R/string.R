# ============================================
# DOWNLOAD INTERACTOME FROM STRING API
# ============================================

#' Download the full STRING interactome for a given organism
#'
#' Downloads the complete protein interaction network from the STRING database,
#' filters by confidence score and returns only interactions involving
#' the input genes or their neighbors (potential Steiner nodes).
#'
#' @param genes Character vector of Gene Symbols.
#' @param org Organism code. One of: "human", "mouse", "rat", "bovine", "canine",
#'   "pig", "rhesus", "chimp", "chicken", "xenopus", "zebrafish", "fly", "worm",
#'   "yeast", "mosquito", "arabidopsis", "ecoli_k12", "ecoli_sakai", "malaria".
#' @param score_threshold Minimum STRING interaction score (0-1000). Default 400.
#'   \itemize{
#'     \item 900-1000: very high confidence
#'     \item 700-900: high confidence
#'     \item 400-700: medium confidence (recommended)
#'     \item 150-400: low confidence
#'   }
#' @param cache_dir Folder to cache the downloaded interactome. \code{NULL} = no cache.
#'
#' @return Data.frame with three columns:
#'   \itemize{
#'     \item \code{from} Gene symbol A
#'     \item \code{to} Gene symbol B
#'     \item \code{cost} Edge cost = 1 - (STRING score / 1000).
#'       Lower cost = more reliable interaction = preferred by PCSF.
#'       With default threshold 400, cost ranges from 0.0 to 0.6.
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' genes <- c("Actb", "Tp53", "Mapk1", "Akt1")
#' interactome <- get_string_interactome(
#'   genes           = genes,
#'   org             = "rat",
#'   score_threshold = 400
#' )
#' head(interactome)
#' }
get_string_interactome <- function(genes,
                                   org = "human",
                                   score_threshold = 400,
                                   cache_dir = NULL) {
  org_info <- ORGANISMS[[org]]
  if (is.null(org_info)) {
    stop("Organism '", org, "' not recognized.\nOptions: ",
         paste(names(ORGANISMS), collapse = ", "))
  }

  # Cache
  if (!is.null(cache_dir)) {
    cache_file <- file.path(cache_dir,
                            paste0("STRING_", org, "_", score_threshold, ".rds"))
    if (file.exists(cache_file)) {
      message("Loading interactome from cache: ", cache_file)
      return(readRDS(cache_file))
    }
  }

  message("Querying STRING API for ", org_info$nombre, "...")

  # 1. Interactions between input genes
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

  # 2. External neighbors (potential Steiner nodes)
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

  # 3. Combine both networks
  if (!is.null(partners) && length(partners) > 0 && nrow(partners) > 0) {
    network <- rbind(network, partners)
    network <- network[!duplicated(paste(network$preferredName_A,
                                         network$preferredName_B)), ]
  }

  if (is.null(network) || length(network) == 0 || nrow(network) == 0) {
    stop(
      "STRING returned no interactions.\n",
      "Suggestions:\n",
      "  - Try lowering score_threshold (current: ", score_threshold, ")\n",
      "  - Verify that gene symbols are correct for ", org_info$nombre
    )
  }

  # Count mapped genes and potential Steiner nodes
  genes_found  <- unique(c(network$preferredName_A, network$preferredName_B))
  steiner_candidates <- setdiff(genes_found, genes)
  message("  ", length(intersect(genes_found, genes)), "/", length(genes),
          " genes mapped in STRING")
  message("  ", length(steiner_candidates), " potential Steiner nodes available")

  # Convert to iPCSF format: from, to, cost
  # cost = 1 - (combined_score / 1000)
  # Score 1000 -> cost 0.0 (very reliable, PCSF prefers it)
  # Score 700  -> cost 0.3 (high confidence)
  # Score 400  -> cost 0.6 (medium confidence, default threshold)
  interactome <- data.frame(
    from = network$preferredName_A,
    to   = network$preferredName_B,
    cost = 1 - (network$score / 1000),
    stringsAsFactors = FALSE
  )

  interactome <- interactome[complete.cases(interactome), ]
  message("  Interactome ready: ", nrow(interactome), " interactions")

  # Save cache
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(interactome, cache_file)
    message("  Cache saved: ", cache_file)
  }

  return(interactome)
}
