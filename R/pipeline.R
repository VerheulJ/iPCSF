# ============================================
# FUNCIÓN PRINCIPAL iPCSF
# ============================================

#' Operador null-coalesce interno
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Análisis PCSF completo con interactoma automático
#'
#' @param conditions Lista nombrada de condiciones. Cada elemento debe contener:
#'   \itemize{
#'     \item \code{data} Data.frame con los genes desregulados (obligatorio)
#'     \item \code{label} Etiqueta para mostrar en el HTML (opcional, default = nombre de la condición)
#'     \item \code{color} Color hex para la condición (opcional, default = "#377EB8")
#'   }
#' @param org Organismo. Ej: "rat", "human", "mouse". Ver \code{ORGANISMS}.
#' @param gene_col Nombre de la columna con Gene Symbols. Default \code{"gene"}.
#' @param log2fc_col Nombre de la columna con log2 fold-change. Default \code{"log2FC"}.
#' @param pval_col Nombre de la columna con p-valor (-log10). Default \code{"pvalue"}.
#' @param score_threshold Score mínimo STRING (0-1000). Default 400.
#' @param cache_dir Carpeta para cachear el interactoma descargado. \code{NULL} = sin caché.
#' @param output_file Nombre del archivo HTML generado. Default \code{"iPCSF_network.html"}.
#' @param w Parámetro omega de PCSF. Default 0.85.
#' @param b Parámetro beta de PCSF. Default 1.
#' @param mu Parámetro mu de PCSF. Default 0.00005.
#'
#' @return Lista invisible con \code{resultados} (lista por condición),
#'   \code{html} (ruta al HTML) y \code{org} (info del organismo).
#' @export
#'
#' @examples
#' \dontrun{
#' # Ejemplo con una condición
#' iPCSF(
#'   conditions = list(
#'     case1 = list(data = df_case1, label = "Treatment A", color = "#E41A1C")
#'   ),
#'   org        = "rat",
#'   gene_col   = "Gene.Symbol",
#'   log2fc_col = "log2cociente",
#'   pval_col   = "log10pvalor"
#' )
#'
#' # Ejemplo con dos condiciones
#' iPCSF(
#'   conditions = list(
#'     case1 = list(data = df_case1, label = "Treatment A", color = "#E41A1C"),
#'     case2 = list(data = df_case2, label = "Treatment B", color = "#377EB8")
#'   ),
#'   org             = "human",
#'   gene_col        = "gene",
#'   log2fc_col      = "log2FC",
#'   pval_col        = "pvalue",
#'   score_threshold = 700,
#'   cache_dir       = "~/.iPCSF_cache",
#'   output_file     = "my_network.html"
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
                  w  = 0.85,
                  b  = 1,
                  mu = 0.00005) {

  # ── Validar organismo ──────────────────────────────────────────
  org_info <- ORGANISMS[[org]]
  if (is.null(org_info)) {
    stop(
      "Organismo '", org, "' no soportado.\n",
      "Opciones disponibles: ", paste(names(ORGANISMS), collapse = ", ")
    )
  }

  # ── Validar conditions ─────────────────────────────────────────
  if (!is.list(conditions) || is.null(names(conditions)) ||
      any(names(conditions) == "")) {
    stop(
      "'conditions' debe ser una lista nombrada.\n",
      "Ejemplo:\n",
      "  list(\n",
      "    case1 = list(data = df1, label = 'Group A', color = '#E41A1C'),\n",
      "    case2 = list(data = df2, label = 'Group B', color = '#377EB8')\n",
      "  )"
    )
  }

  message("\n========================================")
  message("  iPCSF — ", org_info$nombre)
  message("  ", length(conditions), " condición(es): ",
          paste(names(conditions), collapse = ", "))
  message("========================================\n")

  resultados <- list()

  for (cond_id in names(conditions)) {
    cond <- conditions[[cond_id]]

    message("[", cond_id, "] Procesando...")

    # Validar campo data
    if (is.null(cond$data) || !is.data.frame(cond$data)) {
      stop("Condición '", cond_id, "': el campo 'data' debe ser un data.frame.")
    }

    # Validar columnas
    for (col in c(gene_col, log2fc_col, pval_col)) {
      if (!col %in% colnames(cond$data)) {
        stop(
          "Condición '", cond_id, "': columna '", col, "' no encontrada.\n",
          "Columnas disponibles: ", paste(colnames(cond$data), collapse = ", ")
        )
      }
    }

    # Construir red
    res <- tryCatch(
      construir_red(
        desregulated_list = cond$data,
        org_info          = org_info,
        gene_col          = gene_col,
        log2fc_col        = log2fc_col,
        pval_col          = pval_col,
        score_threshold   = score_threshold,
        cache_dir         = cache_dir,
        w = w, b = b, mu = mu
      ),
      error = function(e) {
        message("  [", cond_id, "] ERROR en red: ", e$message)
        return(NULL)
      }
    )

    if (is.null(res)) {
      message("  [", cond_id, "] Omitiendo condición por error.")
      resultados[[cond_id]] <- NULL
      next
    }

    # Metadatos de la condición
    res$label <- cond$label %||% cond_id
    res$color <- cond$color %||% "#377EB8"

    # Enriquecimiento
    message("[", cond_id, "] Calculando enriquecimiento funcional...")
    res <- tryCatch(
      aplicar_enriquecimiento(res, org_info),
      error = function(e) {
        message("  [", cond_id, "] WARNING enriquecimiento: ", e$message)
        res$tooltips <- list()
        return(res)
      }
    )

    resultados[[cond_id]] <- res
    message("[", cond_id, "] OK — ",
            igraph::vcount(res$subnet), " nodos, ",
            length(unique(igraph::V(res$subnet)$cluster)), " clusters\n")
  }

  # Verificar que al menos una condición tuvo éxito
  validos <- Filter(Negate(is.null), resultados)
  if (length(validos) == 0) {
    stop("Ninguna condición produjo resultados. Revisa los datos de entrada.")
  }

  # Generar HTML
  message("[iPCSF] Generando visualización HTML...")
  generar_html(
    resultados  = resultados,
    output_file = output_file,
    titulo      = paste0("iPCSF \u2014 ", org_info$nombre)
  )

  message("[iPCSF] Completado.")
  message("[iPCSF] HTML: ", output_file, "\n")

  return(invisible(list(
    resultados = resultados,
    html       = output_file,
    org        = org_info
  )))
}
