# ============================================
# MAIN FUNCTION iPCSF
# ============================================
#' @useDynLib iPCSF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
# Null-coalesce operator (internal)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @importFrom grDevices rainbow
#' @importFrom stats complete.cases na.omit setNames sd
#' @importFrom utils head
NULL

#' Calculate optimal PCSF parameters w, b, mu from data
#' @noRd
calcular_parametros <- function(terminals, interactome) {
  prize_mean   <- mean(abs(terminals), na.rm = TRUE)
  n_terminals  <- length(terminals)
  n_genes      <- length(unique(c(interactome$from, interactome$to)))
  n_edges      <- nrow(interactome)
  mean_degree  <- (2 * n_edges) / n_genes

  w_opt  <- round(sqrt(n_terminals) * 0.4, 2)
  b_opt  <- 1
  mu_opt <- round(prize_mean / (mean_degree * 100), 5)

  message("  Parameters calculated automatically:")
  message("    w  = ", w_opt,  " (sqrt(", n_terminals, ") * 0.4)")
  message("    b  = ", b_opt)
  message("    mu = ", mu_opt, " (prize_mean / (mean_degree * 100))")

  list(w = w_opt, b = b_opt, mu = mu_opt)
}

#' Collect condition data interactively
#' @noRd
recopilar_condiciones <- function(n) {
  conditions      <- list()
  default_colors  <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                        "#FF7F00", "#A65628", "#F781BF", "#999999")

  for (i in seq_len(n)) {
    cat(sprintf("\n--- Condition %d of %d ---\n", i, n))

    id    <- readline(prompt = "  Short ID (e.g. females): ")
    label <- readline(prompt = "  HTML label (e.g. Females): ")

    color_default <- default_colors[((i - 1) %% length(default_colors)) + 1]
    color_input   <- readline(prompt = sprintf("  Hex color [default: %s]: ", color_default))
    color         <- if (nchar(trimws(color_input)) == 0) color_default else trimws(color_input)

    path <- readline(prompt = "  Path to .xlsx file: ")
    while (!file.exists(path)) {
      cat("  File not found. Please try again.\n")
      path <- readline(prompt = "  Path to .xlsx file: ")
    }

    df <- readxl::read_xlsx(path)
    cat("  Available columns:", paste(colnames(df), collapse = ", "), "\n")

    gene_col   <- readline(prompt = "  Gene column (e.g. Gene.Symbol): ")
    log2fc_col <- readline(prompt = "  log2FC column: ")
    pval_col   <- readline(prompt = "  p-value column (-log10): ")

    conditions[[id]] <- list(
      data       = df,
      label      = if (nchar(trimws(label)) == 0) id else label,
      color      = color,
      gene_col   = gene_col,
      log2fc_col = log2fc_col,
      pval_col   = pval_col
    )
  }
  return(conditions)
}

#' Run a complete PCSF analysis with automatic interactome download
#'
#' @details
#' Two ways to use iPCSF:
#'
#' \strong{1. List mode} -- pass your data directly:
#' \preformatted{
#' iPCSF(
#'   conditions = list(case1 = list(data = df, label = "Case 1")),
#'   org = "rat", gene_col = "Gene.Symbol"
#' )
#' }
#'
#' @param conditions Named list of conditions, or an integer N to collect
#'   N conditions interactively. Each list element must contain:
#'   \itemize{
#'     \item \code{data} Data.frame with desregulated genes (required)
#'     \item \code{label} Label shown in the HTML (optional)
#'     \item \code{color} Hex color for the condition (optional)
#'     \item \code{gene_col} Gene column for this condition (optional, overrides global)
#'     \item \code{log2fc_col} log2FC column for this condition (optional, overrides global)
#'     \item \code{pval_col} p-value column for this condition (optional, overrides global)
#'   }
#' @param org Organism code. One of: "human", "mouse", "rat", "bovine", "canine",
#'   "pig", "rhesus", "chimp", "chicken", "xenopus", "zebrafish", "fly", "worm",
#'   "yeast", "mosquito", "arabidopsis", "ecoli_k12", "ecoli_sakai", "malaria".
#' @param gene_col Column name for Gene Symbols (global default). Default \code{"gene"}.
#' @param log2fc_col Column name for log2 fold-change (global default). Default \code{"log2FC"}.
#' @param pval_col Column name for p-value (-log10) (global default). Default \code{"pvalue"}.
#' @param score_threshold Minimum STRING interaction score (0-1000). Default 400.
#' @param cache_dir Folder to cache the downloaded interactome. \code{NULL} = no cache.
#' @param output_file Name of the generated HTML file. Default \code{"iPCSF_network.html"}.
#' @param w PCSF omega parameter. \code{NULL} = calculated automatically from data.
#' @param b PCSF beta parameter. \code{NULL} = calculated automatically from data.
#' @param mu PCSF mu parameter. \code{NULL} = calculated automatically from data.
#'
#' @return Invisible list with \code{resultados}, \code{html} and \code{org}.
#' @export
#' @examples
#' \donttest{
#' # One condition
#' iPCSF(
#'   conditions = list(
#'     alcohol = list(data = df, label = "Alcohol", color = "#E41A1C")
#'   ),
#'   org         = "rat",
#'   gene_col    = "Gene.Symbol",
#'   log2fc_col  = "log2cociente_alcohol",
#'   pval_col    = "log10pvalor_alcohol",
#'   cache_dir   = tempdir(),
#'   output_file = file.path(tempdir(), "network.html")
#' )
#'
#' # Two conditions with different column names per condition
#' iPCSF(
#'   conditions = list(
#'     females = list(data = df_f, label = "Females", color = "#E41A1C",
#'                    log2fc_col = "log2cociente_hembras",
#'                    pval_col   = "log10pvalor_hembras"),
#'     males   = list(data = df_m, label = "Males",   color = "#377EB8",
#'                    log2fc_col = "log2cociente_machos",
#'                    pval_col   = "log10pvalor_machos")
#'   ),
#'   org         = "rat",
#'   gene_col    = "Gene.Symbol",
#'   cache_dir   = tempdir(),
#'   output_file = file.path(tempdir(), "network.html")
#' )
#' }


iPCSF <- function(conditions,
                  org             = "human",
                  gene_col        = "gene",
                  log2fc_col      = "log2FC",
                  pval_col        = "pvalue",
                  score_threshold = 400,
                  cache_dir       = NULL,
                  output_file     = "iPCSF_network.html",
                  w  = NULL,
                  b  = NULL,
                  mu = NULL) {

  # Validate organism
  org_info <- ORGANISMS[[org]]
  if (is.null(org_info)) {
    stop(
      "Organism '", org, "' not supported.\n",
      "Available options: ", paste(names(ORGANISMS), collapse = ", ")
    )
  }

  # Interactive mode if conditions is an integer
  if (is.numeric(conditions) && length(conditions) == 1) {
    n <- as.integer(conditions)
    cat(sprintf("\n[iPCSF] Interactive mode: %d condition(s)\n", n))
    conditions <- recopilar_condiciones(n)
  }

  # Validate conditions list
  if (!is.list(conditions) || is.null(names(conditions)) ||
      any(names(conditions) == "")) {
    stop(
      "'conditions' must be a named list or an integer N.\n",
      "Examples:\n",
      "  conditions = 2  (interactive mode)\n",
      "  conditions = list(case1 = list(data = df1, ...))\n"
    )
  }

  message("\n========================================")
  message("  iPCSF -- ", org_info$nombre)
  message("  ", length(conditions), " condition(s): ",
          paste(names(conditions), collapse = ", "))
  if (is.null(w)) message("  w, b, mu: calculated automatically")
  else            message("  w=", w, " b=", b, " mu=", mu)
  message("========================================\n")

  resultados <- list()

  for (cond_id in names(conditions)) {
    cond <- conditions[[cond_id]]

    message("[", cond_id, "] Processing...")

    # Validate data field
    if (is.null(cond$data) || !is.data.frame(cond$data)) {
      stop("Condition '", cond_id, "': field 'data' must be a data.frame.")
    }

    # Use per-condition columns if provided, otherwise use globals
    cond_gene_col   <- cond$gene_col   %||% gene_col
    cond_log2fc_col <- cond$log2fc_col %||% log2fc_col
    cond_pval_col   <- cond$pval_col   %||% pval_col

    message("  gene_col   : ", cond_gene_col)
    message("  log2fc_col : ", cond_log2fc_col)
    message("  pval_col   : ", cond_pval_col)

    # Validate columns
    for (col in c(cond_gene_col, cond_log2fc_col, cond_pval_col)) {
      if (!col %in% colnames(cond$data)) {
        stop(
          "Condition '", cond_id, "': column '", col, "' not found.\n",
          "Available columns: ", paste(colnames(cond$data), collapse = ", ")
        )
      }
    }

    # Get organism key
    org_key <- names(ORGANISMS)[
      sapply(ORGANISMS, function(x) identical(x$taxid, org_info$taxid))
    ]
    genes_cond <- as.character(cond$data[[cond_gene_col]])
    genes_cond <- genes_cond[!is.na(genes_cond) & genes_cond != ""]

    # Auto-calculate parameters if not provided
    w_use  <- w
    b_use  <- b
    mu_use <- mu

    if (is.null(w) || is.null(b) || is.null(mu)) {
      message("  Downloading interactome to calculate parameters...")
      interactome_tmp <- tryCatch(
        get_string_interactome(
          genes           = genes_cond,
          org             = org_key,
          score_threshold = score_threshold,
          cache_dir       = cache_dir
        ),
        error = function(e) NULL
      )

      if (!is.null(interactome_tmp)) {
        terminals_tmp <- setNames(
          abs(as.numeric(cond$data[[cond_pval_col]])),
          genes_cond
        )
        terminals_tmp <- terminals_tmp[!is.na(terminals_tmp) & terminals_tmp > 0]

        params <- calcular_parametros(terminals_tmp, interactome_tmp)
        w_use  <- w  %||% params$w
        b_use  <- b  %||% params$b
        mu_use <- mu %||% params$mu
      } else {
        w_use  <- w  %||% 2
        b_use  <- b  %||% 1
        mu_use <- mu %||% 0.0005
        message("  Could not calculate parameters, using defaults: w=", w_use,
                " b=", b_use, " mu=", mu_use)
      }
    }

    # Build network
    res <- tryCatch(
      construir_red(
        desregulated_list = cond$data,
        org_info          = org_info,
        gene_col          = cond_gene_col,
        log2fc_col        = cond_log2fc_col,
        pval_col          = cond_pval_col,
        score_threshold   = score_threshold,
        cache_dir         = cache_dir,
        w = w_use, b = b_use, mu = mu_use
      ),
      error = function(e) {
        message("  [", cond_id, "] ERROR in network: ", e$message)
        return(NULL)
      }
    )

    if (is.null(res)) {
      message("  [", cond_id, "] Skipping condition due to error.")
      resultados[[cond_id]] <- NULL
      next
    }

    # Condition metadata
    res$label <- cond$label %||% cond_id
    res$color <- cond$color %||% "#377EB8"

    # Functional enrichment
    message("[", cond_id, "] Running functional enrichment...")
    res <- tryCatch(
      aplicar_enriquecimiento(res, org_info),
      error = function(e) {
        message("  [", cond_id, "] WARNING enrichment: ", e$message)
        res$tooltips <- list()
        return(res)
      }
    )

    resultados[[cond_id]] <- res
    message("[", cond_id, "] OK -- ",
            igraph::vcount(res$subnet), " nodes, ",
            length(unique(igraph::V(res$subnet)$cluster)), " clusters\n")
  }

  # Check at least one condition succeeded
  validos <- Filter(Negate(is.null), resultados)
  if (length(validos) == 0) {
    stop("No conditions produced results. Check your input data.")
  }

  # Generate HTML
  message("[iPCSF] Generating HTML visualization...")
  generar_html(
    resultados  = resultados,
    output_file = output_file,
    titulo      = paste0("iPCSF \u2014 ", org_info$nombre)
  )

  message("[iPCSF] Done.")
  message("[iPCSF] HTML: ", output_file, "\n")

  return(invisible(list(
    resultados = resultados,
    html       = output_file,
    org        = org_info
  )))
}
