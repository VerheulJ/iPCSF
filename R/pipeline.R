# ============================================
# FUNCIoN PRINCIPAL iPCSF
# ============================================

# Null-coalesce operator (internal)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @importFrom grDevices rainbow
#' @importFrom stats complete.cases na.omit setNames sd
#' @importFrom utils head readline
NULL

#' Calcula parametros optimos w, b, mu a partir de los datos
#' @keywords internal
calcular_parametros <- function(terminals, interactome) {
  prize_medio     <- mean(abs(terminals), na.rm = TRUE)
  n_terminals     <- length(terminals)
  n_genes_total   <- length(unique(c(interactome$from, interactome$to)))
  n_interacciones <- nrow(interactome)
  grado_medio     <- (2 * n_interacciones) / n_genes_total

  w_opt  <- round(sqrt(n_terminals) * 0.4, 2)
  b_opt  <- 1
  mu_opt <- round(prize_medio / (grado_medio * 100), 5)

  message("  Parametros calculados automaticamente:")
  message("    w  = ", w_opt,  " (sqrt(", n_terminals, ") * 0.4)")
  message("    b  = ", b_opt)
  message("    mu = ", mu_opt, " (prize_medio / (grado_medio * 100))")

  list(w = w_opt, b = b_opt, mu = mu_opt)
}

#' Recopila interactivamente los datos de N condiciones
#' @keywords internal
recopilar_condiciones <- function(n) {
  conditions <- list()
  colores_default <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#A65628", "#F781BF", "#999999")

  for (i in seq_len(n)) {
    cat(sprintf("\n--- Condicion %d de %d ---\n", i, n))

    id    <- readline(prompt = "  ID corto (ej: females): ")
    label <- readline(prompt = "  Etiqueta HTML (ej: Females): ")

    color_default <- colores_default[((i - 1) %% length(colores_default)) + 1]
    color_input   <- readline(prompt = sprintf("  Color hex [default: %s]: ", color_default))
    color         <- if (nchar(trimws(color_input)) == 0) color_default else trimws(color_input)

    ruta <- readline(prompt = "  Ruta al archivo .xlsx: ")
    while (!file.exists(ruta)) {
      cat("  Archivo no encontrado. Intenta de nuevo.\n")
      ruta <- readline(prompt = "  Ruta al archivo .xlsx: ")
    }

    df        <- readxl::read_xlsx(ruta)
    cat("  Columnas disponibles:", paste(colnames(df), collapse = ", "), "\n")

    gene_col   <- readline(prompt = "  Columna de genes (ej: Gene.Symbol): ")
    log2fc_col <- readline(prompt = "  Columna log2FC: ")
    pval_col   <- readline(prompt = "  Columna pvalor (-log10): ")

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

#' Analisis PCSF completo con interactoma automatico
#'
#' @param conditions Lista nombrada de condiciones, o un entero N para
#'   recopilar N condiciones interactivamente. Cada elemento de la lista debe contener:
#'   \itemize{
#'     \item \code{data} Data.frame con los genes desregulados (obligatorio)
#'     \item \code{label} Etiqueta para mostrar en el HTML (opcional)
#'     \item \code{color} Color hex para la condicion (opcional)
#'     \item \code{gene_col} Columna de genes (opcional, sobreescribe el global)
#'     \item \code{log2fc_col} Columna log2FC (opcional, sobreescribe el global)
#'     \item \code{pval_col} Columna pvalor (opcional, sobreescribe el global)
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
#'
#' @examples
#' \dontrun{
#' # Interactive mode -- iPCSF asks for data of each condition
#' iPCSF(conditions = 2, org = "rat", gene_col = "Gene.Symbol")
#'
#' # One condition, automatic parameters
#' iPCSF(
#'   conditions = list(
#'     alcohol = list(data = df, label = "Alcohol", color = "#E41A1C")
#'   ),
#'   org        = "rat",
#'   gene_col   = "Gene.Symbol",
#'   log2fc_col = "log2cociente_alcohol",
#'   pval_col   = "log10pvalor_alcohol"
#' )
#'
#' # Two conditions with different columns, manual parameters
#' iPCSF(
#'   conditions = list(
#'     females = list(data = df_f, label = "Females", color = "#E41A1C",
#'                    log2fc_col = "log2cociente_hembras",
#'                    pval_col   = "log10pvalor_hembras"),
#'     males   = list(data = df_m, label = "Males",   color = "#377EB8",
#'                    log2fc_col = "log2cociente_machos",
#'                    pval_col   = "log10pvalor_machos")
#'   ),
#'   org      = "rat",
#'   gene_col = "Gene.Symbol",
#'   w = 5, b = 1, mu = 0.005
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

  # Si conditions es un entero, modo interactivo
  if (is.numeric(conditions) && length(conditions) == 1) {
    n <- as.integer(conditions)
    cat(sprintf("\n[iPCSF] Modo interactivo: %d condicion(es)\n", n))
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
  if (is.null(w)) message("  w, b, mu: calculados automaticamente")
  else            message("  w=", w, " b=", b, " mu=", mu)
  message("========================================\n")

  resultados <- list()

  for (cond_id in names(conditions)) {
    cond <- conditions[[cond_id]]

    message("[", cond_id, "] Processing...")

    if (is.null(cond$data) || !is.data.frame(cond$data)) {
      stop("Condition '", cond_id, "': field 'data' must be a data.frame.")
    }

    cond_gene_col   <- cond$gene_col   %||% gene_col
    cond_log2fc_col <- cond$log2fc_col %||% log2fc_col
    cond_pval_col   <- cond$pval_col   %||% pval_col

    message("  gene_col   : ", cond_gene_col)
    message("  log2fc_col : ", cond_log2fc_col)
    message("  pval_col   : ", cond_pval_col)

    for (col in c(cond_gene_col, cond_log2fc_col, cond_pval_col)) {
      if (!col %in% colnames(cond$data)) {
        stop(
          "Condition '", cond_id, "': column '", col, "' not found.\n",
          "Available columns: ", paste(colnames(cond$data), collapse = ", ")
        )
      }
    }

    # Obtener interactoma para calcular parametros si son NULL
    org_key <- names(ORGANISMS)[
      sapply(ORGANISMS, function(x) identical(x$taxid, org_info$taxid))
    ]
    genes_cond <- as.character(cond$data[[cond_gene_col]])
    genes_cond <- genes_cond[!is.na(genes_cond) & genes_cond != ""]

    # Calcular parametros automaticamente si no se especifican
    w_use  <- w
    b_use  <- b
    mu_use <- mu

    if (is.null(w) || is.null(b) || is.null(mu)) {
      message("  Descargando interactoma para calcular parametros...")
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
        message("  No se pudo calcular parametros, usando defaults: w=", w_use,
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

    res$label <- cond$label %||% cond_id
    res$color <- cond$color %||% "#377EB8"

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

  validos <- Filter(Negate(is.null), resultados)
  if (length(validos) == 0) {
    stop("No conditions produced results. Check your input data.")
  }

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
