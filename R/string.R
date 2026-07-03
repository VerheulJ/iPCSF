# ============================================
# DOWNLOAD INTERACTOME FROM STRING API
# ============================================

#' Download the full STRING interactome for a given organism
#'
#' Downloads the complete protein interaction network from the STRING database,
#' filters by confidence score and returns only interactions involving
#' the input genes or their neighbors (potential Steiner nodes).
#' The full interactome is cached to disk on first use — subsequent calls
#' load from cache instantly regardless of which genes are passed.
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
#'   On first run downloads the full organism interactome (~75MB for rat).
#'   Subsequent runs load from cache instantly.
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
#' \donttest{
#' genes <- c("Actb", "Tp53", "Mapk1", "Akt1")
#' interactome <- get_string_interactome(
#'   genes           = genes,
#'   org             = "rat",
#'   score_threshold = 400,
#'   cache_dir       = tempdir()
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

  # ── Cache: full interactome (not filtered by genes) ──────────────
  # Key includes org and score_threshold but NOT genes --
  # this ensures the same full interactome is reused for all conditions
  full_interactome <- NULL
  cache_file <- NULL

  if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    cache_file <- file.path(cache_dir,
                            paste0("STRING_", org, "_", score_threshold, "_FULL.rds"))
    if (file.exists(cache_file)) {
      message("Loading interactome from cache: ", cache_file)
      full_interactome <- readRDS(cache_file)
    }
  }

  # ── Download full interactome if not cached ───────────────────────
  if (is.null(full_interactome)) {
    message("Downloading STRING interactome for ", org_info$nombre, "...")
    message("(this may take a few minutes on first run)")

    url_links <- paste0(
      "https://stringdb-downloads.org/download/protein.links.v12.0/",
      org_info$taxid, ".protein.links.v12.0.txt.gz"
    )
    url_info <- paste0(
      "https://stringdb-downloads.org/download/protein.info.v12.0/",
      org_info$taxid, ".protein.info.v12.0.txt.gz"
    )

    tmp_links <- tempfile(fileext = ".txt.gz")
    tmp_info  <- tempfile(fileext = ".txt.gz")

    tryCatch({
      download.file(url_links, tmp_links, mode = "wb", quiet = FALSE)
      download.file(url_info,  tmp_info,  mode = "wb", quiet = FALSE)
    }, error = function(e) {
      stop("Error downloading STRING: ", e$message,
           "\nCheck your internet connection.")
    })

    message("Processing interactome...")

    links <- data.table::fread(tmp_links, sep = " ",  data.table = FALSE)
    info  <- data.table::fread(tmp_info,  sep = "\t", data.table = FALSE, quote = "")

    # Filter by score threshold
    links <- links[links$combined_score >= score_threshold, ]

    # Remove taxid prefix from IDs
    prefix        <- paste0(org_info$taxid, ".")
    info$clean_id <- sub(paste0("^", prefix), "", info[,1])
    links$clean1  <- sub(paste0("^", prefix), "", links$protein1)
    links$clean2  <- sub(paste0("^", prefix), "", links$protein2)

    # Map protein IDs to gene symbols
    id_to_symbol <- setNames(info$preferred_name, info$clean_id)
    links$gene1  <- id_to_symbol[links$clean1]
    links$gene2  <- id_to_symbol[links$clean2]

    # Build full interactome
    # cost = 1 - (combined_score / 1000)
    # Score 1000 -> cost 0.0 (very reliable, PCSF prefers it)
    # Score 700  -> cost 0.3 (high confidence)
    # Score 400  -> cost 0.6 (medium confidence, default threshold)
    full_interactome <- data.frame(
      from = links$gene1,
      to   = links$gene2,
      cost = 1 - (links$combined_score / 1000),
      stringsAsFactors = FALSE
    )
    full_interactome <- full_interactome[complete.cases(full_interactome), ]

    message("  Full interactome: ", nrow(full_interactome), " interactions")

    # Save full interactome to cache
    if (!is.null(cache_file)) {
      saveRDS(full_interactome, cache_file)
      message("  Cache saved: ", cache_file)
    }

    unlink(tmp_links)
    unlink(tmp_info)
  }

  # ── Filter by input genes + their neighbors ───────────────────────
  # Keep interactions where at least one end is a terminal gene
  # The other end becomes a potential Steiner node
  idx      <- full_interactome$from %in% genes | full_interactome$to %in% genes
  filtered <- full_interactome[idx, ]

  if (nrow(filtered) == 0) {
    stop(
      "No interactions found for the input genes.\n",
      "Suggestions:\n",
      "  - Try lowering score_threshold (current: ", score_threshold, ")\n",
      "  - Verify that gene symbols are correct for ", org_info$nombre
    )
  }

  genes_found        <- unique(c(filtered$from, filtered$to))
  steiner_candidates <- setdiff(genes_found, genes)

  message("  ", length(intersect(genes_found, genes)), "/", length(genes),
          " genes mapped in STRING")
  message("  ", length(steiner_candidates), " potential Steiner nodes available")
  message("  Interactome ready: ", nrow(filtered), " interactions")

  return(filtered)
}
