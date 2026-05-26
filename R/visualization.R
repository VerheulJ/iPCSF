# ============================================
# VISUALIZACIoN -- generar_html()
# ============================================

#' Genera paleta de colores para los clusters
#' @noRd
generar_paleta_clusters <- function(n_clusters) {
  colores_base <- c(
    "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
    "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
    "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#B3B3B3",
    "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E"
  )
  if (n_clusters <= length(colores_base)) return(colores_base[1:n_clusters])
  c(colores_base, rainbow(n_clusters - length(colores_base), s = 0.7, v = 0.9))
}

#' Prepara los datos de una condicion para JSON
#' @noRd
preparar_datos_condicion <- function(resultado) {

  subnet        <- resultado$subnet
  cluster_means <- resultado$cluster_means
  tooltips_map  <- resultado$tooltips %||% list()

  clusters_unicos   <- sort(unique(igraph::V(subnet)$cluster))
  paleta            <- generar_paleta_clusters(length(clusters_unicos))
  cluster_color_map <- setNames(paleta, as.character(clusters_unicos))

  # Layout
  set.seed(42)
  layout_coords <- igraph::layout_with_fr(subnet, niter = 1000) * 5
  rango_x <- max(layout_coords[,1]) - min(layout_coords[,1]) + 0.001
  rango_y <- max(layout_coords[,2]) - min(layout_coords[,2]) + 0.001
  layout_coords[,1] <- (layout_coords[,1] - min(layout_coords[,1])) / rango_x * 3000 + 100
  layout_coords[,2] <- (layout_coords[,2] - min(layout_coords[,2])) / rango_y * 2500 + 100

  # Tamanos
  prize <- abs(igraph::V(subnet)$original_prize)
  sizes <- if (max(prize) > min(prize)) {
    ((prize - min(prize)) / (max(prize) - min(prize))) * 50 + 15
  } else rep(20, length(prize))

  # Bordes segun regulacion
  border_colors <- sapply(igraph::V(subnet)$status, function(s) {
    if (is.na(s) || s %in% c("Steiner", "unknown")) "#888888"
    else if (grepl("up",   s, ignore.case = TRUE))   "#FF0000"
    else if (grepl("down", s, ignore.case = TRUE))   "#0000FF"
    else "#888888"
  })

  # Tooltips
  node_titles <- sapply(igraph::V(subnet)$cluster, function(cl) {
    tt <- tooltips_map[[as.character(cl)]]
    if (is.null(tt)) "" else tt
  })

  # Nodos
  nodos <- data.frame(
    id      = as.character(igraph::V(subnet)$name),
    label   = as.character(igraph::V(subnet)$name),
    x       = as.numeric(layout_coords[,1]),
    y       = as.numeric(layout_coords[,2]),
    size    = as.numeric(sizes),
    color   = as.character(sapply(igraph::V(subnet)$cluster,
                                  function(cl) cluster_color_map[[as.character(cl)]])),
    border  = as.character(border_colors),
    shape   = as.character(ifelse(igraph::V(subnet)$tipo_nodo == "Terminal",
                                  "dot", "triangle")),
    cluster = as.integer(igraph::V(subnet)$cluster),
    title   = as.character(node_titles),
    status  = as.character(igraph::V(subnet)$status),
    stringsAsFactors = FALSE
  )

  # Aristas
  edge_list <- igraph::as_edgelist(subnet)
  aristas <- data.frame(
    from   = as.character(edge_list[,1]),
    to     = as.character(edge_list[,2]),
    weight = as.numeric(
      if (!is.null(igraph::E(subnet)$weight)) igraph::E(subnet)$weight
      else rep(1, igraph::ecount(subnet))
    ),
    stringsAsFactors = FALSE
  )

  # Clusters ordenados por |log2FC| medio
  cluster_order <- order(cluster_means, decreasing = TRUE)
  cluster_info  <- data.frame(
    cluster = as.integer(clusters_unicos[cluster_order]),
    color   = as.character(paleta[cluster_order]),
    count   = as.integer(sapply(clusters_unicos[cluster_order],
                                function(cl) sum(igraph::V(subnet)$cluster == cl))),
    mean_fc = as.numeric(round(cluster_means[cluster_order], 3)),
    rank    = as.integer(seq_along(cluster_order)),
    stringsAsFactors = FALSE
  )

  list(
    nodos    = nodos,
    aristas  = aristas,
    clusters = cluster_info,
    stats    = list(
      nodos     = as.integer(igraph::vcount(subnet)),
      aristas   = as.integer(igraph::ecount(subnet)),
      clusters  = as.integer(length(clusters_unicos)),
      terminals = as.integer(sum(igraph::V(subnet)$tipo_nodo == "Terminal"))
    )
  )
}

#' Genera el HTML interactivo con la red iPCSF
#'
#' @param resultados Lista nombrada de resultados (una entrada por condicion).
#'   Cada elemento debe ser la salida de \code{construir_red()} + \code{aplicar_enriquecimiento()},
#'   con campos adicionales \code{label} y \code{color}.
#' @param output_file Ruta del archivo HTML de salida.
#' @param titulo Titulo que aparece en el header del HTML.
#'
#' @return Ruta del archivo HTML generado (invisible).
#' @export
generar_html <- function(resultados,
                         output_file = "iPCSF_network.html",
                         titulo      = "iPCSF Network") {

  # Filtrar condiciones con datos
  validos <- Filter(Negate(is.null), resultados)
  if (length(validos) == 0) stop("No hay resultados validos para visualizar.")

  # Preparar datos JSON para cada condicion
  datos_completos <- list()
  for (cond_id in names(validos)) {
    datos_completos[[cond_id]] <- preparar_datos_condicion(validos[[cond_id]])
  }

  datos_json <- jsonlite::toJSON(
    datos_completos,
    auto_unbox = TRUE, na = "null", force = TRUE
  )

  # Metadatos de condiciones para los botones (label + color)
  meta_condiciones <- lapply(names(validos), function(cond_id) {
    list(
      id    = cond_id,
      label = validos[[cond_id]]$label %||% cond_id,
      color = validos[[cond_id]]$color %||% "#377EB8"
    )
  })
  meta_json <- jsonlite::toJSON(meta_condiciones, auto_unbox = TRUE)

  # Primera condicion activa por defecto
  primera_condicion <- names(validos)[1]

  html <- paste0(
    '<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>', titulo, '</title>
  <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap" rel="stylesheet">
  <style>
    :root {
      --bg:#0a0e1a; --surface:#111827; --surface2:#1a2234;
      --border:rgba(255,255,255,0.07); --text:#e2e8f0; --muted:#64748b;
      --accent:#38bdf8; --up:#f87171; --down:#60a5fa;
    }
    *{margin:0;padding:0;box-sizing:border-box}
    body{font-family:"IBM Plex Sans",sans-serif;background:var(--bg);
         color:var(--text);height:100vh;display:flex;flex-direction:column;overflow:hidden}

    /* Header */
    .header{padding:14px 28px;border-bottom:1px solid var(--border);
            display:flex;align-items:center;justify-content:space-between;
            background:linear-gradient(90deg,var(--surface) 0%,var(--bg) 100%)}
    .header-left h1{font-size:1.1em;font-weight:600;color:#fff;letter-spacing:-0.3px}
    .header-left p{font-size:0.72em;color:var(--muted);margin-top:2px;
                   font-family:"IBM Plex Mono",monospace}
    .badge{background:rgba(56,189,248,0.1);border:1px solid rgba(56,189,248,0.25);
           color:var(--accent);padding:4px 12px;border-radius:20px;
           font-size:0.72em;font-family:"IBM Plex Mono",monospace}

    /* Controls -- generados dinamicamente */
    .controls{padding:10px 28px;border-bottom:1px solid var(--border);
              display:flex;align-items:center;gap:10px;flex-wrap:wrap;
              background:var(--surface)}
    .ctrl-label{font-size:0.65em;color:var(--muted);text-transform:uppercase;
                letter-spacing:1.2px;margin-right:4px}
    .cond-btn{background:var(--surface2);border:2px solid var(--c,#444);
              color:#fff;padding:7px 20px;border-radius:8px;cursor:pointer;
              font-family:"IBM Plex Sans",sans-serif;font-size:0.85em;font-weight:500;
              transition:all 0.2s ease}
    .cond-btn:hover{filter:brightness(1.25)}
    .cond-btn.active{
      background:color-mix(in srgb,var(--c) 18%,transparent);
      box-shadow:0 0 0 1px var(--c),0 4px 16px color-mix(in srgb,var(--c) 25%,transparent)
    }

    /* Stats */
    .stats{display:flex;border-bottom:1px solid var(--border)}
    .stat{flex:1;padding:9px 16px;text-align:center;border-right:1px solid var(--border)}
    .stat:last-child{border-right:none}
    .stat-val{font-family:"IBM Plex Mono",monospace;font-size:1.35em;
              font-weight:600;color:var(--accent)}
    .stat-lbl{font-size:0.65em;color:var(--muted);text-transform:uppercase;
              letter-spacing:1.2px;margin-top:1px}

    /* Workspace */
    .workspace{display:flex;flex:1;overflow:hidden}
    #network-wrap{flex:1;background:#f8fafc}
    #network{width:100%;height:100%}

    /* Sidebar */
    .sidebar{width:250px;background:var(--surface);border-left:1px solid var(--border);
             display:flex;flex-direction:column;overflow:hidden}
    .sidebar-section{padding:14px 14px 10px;border-bottom:1px solid var(--border)}
    .sidebar-title{font-size:0.65em;font-weight:600;color:var(--muted);
                   text-transform:uppercase;letter-spacing:1.2px;margin-bottom:10px}
    .cluster-scroll{flex:1;overflow-y:auto;padding:8px}
    .cluster-scroll::-webkit-scrollbar{width:4px}
    .cluster-scroll::-webkit-scrollbar-thumb{background:var(--border);border-radius:2px}

    /* Cluster items */
    .cluster-item{display:flex;align-items:center;gap:8px;padding:7px 8px;
                  border-radius:7px;cursor:pointer;transition:background 0.15s;
                  border:1px solid transparent}
    .cluster-item:hover{background:rgba(255,255,255,0.04);border-color:var(--border)}
    .cluster-item.active{background:rgba(56,189,248,0.07);
                         border-color:rgba(56,189,248,0.2)}
    .cluster-item.dimmed{opacity:0.18}
    .swatch{width:10px;height:10px;border-radius:2px;flex-shrink:0}
    .cl-info{flex:1;min-width:0}
    .cl-name{font-size:0.8em;font-weight:500}
    .cl-meta{font-size:0.67em;color:var(--muted);font-family:"IBM Plex Mono",monospace}
    .cl-rank{font-size:0.6em;font-weight:700;color:var(--muted);
             background:rgba(255,255,255,0.06);padding:2px 5px;border-radius:3px;flex-shrink:0}

    /* Legend */
    .legend-row{display:flex;align-items:center;gap:8px;
                font-size:0.75em;color:var(--muted);margin-bottom:6px}
    .dot{width:10px;height:10px;border-radius:50%;flex-shrink:0}
    .tri{width:0;height:0;border-left:6px solid transparent;
         border-right:6px solid transparent;border-bottom:10px solid;flex-shrink:0}
    .ring{width:10px;height:10px;border-radius:50%;border:3px solid;flex-shrink:0}

    /* Tooltip */
    #tooltip{position:fixed;display:none;z-index:9999;pointer-events:none;max-width:340px;font-size:0.8em}
  </style>
</head>
<body>

<div class="header">
  <div class="header-left">
    <h1>', titulo, '</h1>
    <p>Prize-Collecting Steiner Forest &nbsp;&middot;&nbsp; STRING PPI &nbsp;&middot;&nbsp; iPCSF</p>
  </div>
  <div class="badge">iPCSF v0.0.1</div>
</div>

<!-- Botones generados dinamicamente desde R -->
<div class="controls" id="controls">
  <span class="ctrl-label">Condition</span>
</div>

<div class="stats">
  <div class="stat"><div class="stat-val" id="s-nodos">--</div><div class="stat-lbl">Nodes</div></div>
  <div class="stat"><div class="stat-val" id="s-edges">--</div><div class="stat-lbl">Edges</div></div>
  <div class="stat"><div class="stat-val" id="s-clusters">--</div><div class="stat-lbl">Clusters</div></div>
  <div class="stat"><div class="stat-val" id="s-terminals">--</div><div class="stat-lbl">Terminals</div></div>
</div>

<div class="workspace">
  <div id="network-wrap"><div id="network"></div></div>
  <div class="sidebar">
    <div class="sidebar-section">
      <div class="sidebar-title">Clusters -- ranked by |log2FC|</div>
    </div>
    <div class="cluster-scroll" id="cluster-list"></div>
    <div class="sidebar-section">
      <div class="sidebar-title">Node type</div>
      <div class="legend-row"><div class="dot" style="background:#818cf8"></div>Terminal (differential)</div>
      <div class="legend-row"><div class="tri" style="border-bottom-color:#818cf8"></div>Steiner (connector)</div>
    </div>
    <div class="sidebar-section">
      <div class="sidebar-title">Regulation (border)</div>
      <div class="legend-row"><div class="ring" style="border-color:var(--up)"></div>Up-regulated</div>
      <div class="legend-row"><div class="ring" style="border-color:var(--down)"></div>Down-regulated</div>
      <div class="legend-row"><div class="ring" style="border-color:#888"></div>Steiner / unknown</div>
    </div>
    <div class="sidebar-section" style="font-size:0.67em;color:var(--muted);
         font-family:\'IBM Plex Mono\',monospace;line-height:1.6;border-bottom:none">
      Click cluster to highlight<br>Scroll to zoom &middot; Drag to pan
    </div>
  </div>
</div>
<div id="tooltip"></div>

<script>
const DATA = ', datos_json, ';
const META = ', meta_json, ';
let activeCond = "', primera_condicion, '";
let net = null;
let selCluster = null;

// Generar botones dinamicamente desde META
const controls = document.getElementById("controls");
META.forEach(m => {
  const btn = document.createElement("button");
  btn.className = "cond-btn" + (m.id === activeCond ? " active" : "");
  btn.dataset.cond = m.id;
  btn.textContent = m.label;
  btn.style.setProperty("--c", m.color);
  btn.addEventListener("click", function() {
    document.querySelectorAll(".cond-btn").forEach(b => b.classList.remove("active"));
    this.classList.add("active");
    activeCond = this.dataset.cond;
    render();
  });
  controls.appendChild(btn);
});

function render() {
  selCluster = null;
  const d = DATA[activeCond];

  if (!d || !d.nodos || d.nodos.length === 0) {
    ["s-nodos","s-edges","s-clusters","s-terminals"].forEach(id =>
      document.getElementById(id).textContent = "0");
    if (net) net.destroy();
    document.getElementById("network").innerHTML =
      \'<div style="display:flex;align-items:center;justify-content:center;\' +
      \'height:100%;color:#94a3b8;font-size:1em;">No data for this condition</div>\';
    document.getElementById("cluster-list").innerHTML = "";
    return;
  }

  // Stats
  document.getElementById("s-nodos").textContent    = d.stats.nodos;
  document.getElementById("s-edges").textContent    = d.stats.aristas;
  document.getElementById("s-clusters").textContent = d.stats.clusters;
  document.getElementById("s-terminals").textContent = d.stats.terminals;

  // Red
  const nodes = new vis.DataSet(d.nodos.map(n => ({
    id: n.id, label: n.label, x: n.x, y: n.y, size: n.size,
    color: { background: n.color, border: n.border,
             highlight: { background: n.color, border: n.border } },
    borderWidth: 3, shape: n.shape,
    font: { size: 9, color: "#1e293b", face: "IBM Plex Sans" },
    cluster: n.cluster, title: n.title
  })));

  const edges = new vis.DataSet(d.aristas.map((e, i) => ({
    id: i, from: e.from, to: e.to,
    color: { color: "#cbd5e1", opacity: 0.6 }, width: 1.2
  })));

  if (net) net.destroy();
  net = new vis.Network(
    document.getElementById("network"),
    { nodes, edges },
    { physics: false,
      interaction: { hover: true, zoomView: true, dragView: true, dragNodes: true } }
  );

  // Tooltip
  const tooltip = document.getElementById("tooltip");
  net.on("hoverNode", p => {
    const n = nodes.get(p.node);
    if (n && n.title) { tooltip.innerHTML = n.title; tooltip.style.display = "block"; }
  });
  net.on("blurNode", () => { tooltip.style.display = "none"; });
  document.getElementById("network-wrap").addEventListener("mousemove", e => {
    if (tooltip.style.display === "block") {
      const x = e.clientX + 14, y = e.clientY + 14;
      const tw = tooltip.offsetWidth, th = tooltip.offsetHeight;
      tooltip.style.left = (x + tw > window.innerWidth  ? e.clientX - tw - 10 : x) + "px";
      tooltip.style.top  = (y + th > window.innerHeight ? e.clientY - th - 10 : y) + "px";
    }
  });

  // Leyenda clusters
  const list = document.getElementById("cluster-list");
  list.innerHTML = "";
  d.clusters.forEach(c => {
    const el = document.createElement("div");
    el.className = "cluster-item";
    el.dataset.cluster = c.cluster;
    el.innerHTML =
      `<div class="swatch" style="background:${c.color}"></div>
       <div class="cl-info">
         <div class="cl-name">Cluster ${c.cluster}</div>
         <div class="cl-meta">${c.count} genes &middot; |FC|=${c.mean_fc}</div>
       </div>
       <div class="cl-rank">#${c.rank}</div>`;
    el.onclick = () => toggleCluster(c.cluster, nodes, d);
    list.appendChild(el);
  });
}

function toggleCluster(cid, nodes, d) {
  if (selCluster === cid) {
    selCluster = null;
    nodes.update(d.nodos.map(n =>
      ({ id: n.id, color: { background: n.color, border: n.border } })));
    document.querySelectorAll(".cluster-item").forEach(x =>
      x.classList.remove("dimmed", "active"));
  } else {
    selCluster = cid;
    nodes.update(d.nodos.map(n => ({
      id: n.id,
      color: n.cluster === cid
        ? { background: n.color, border: n.border }
        : { background: "#e2e8f0", border: "#cbd5e1" }
    })));
    document.querySelectorAll(".cluster-item").forEach(x => {
      const same = parseInt(x.dataset.cluster) === cid;
      x.classList.toggle("dimmed", !same);
      x.classList.toggle("active",  same);
    });
  }
}

render();
</script>
</body>
</html>')

  con <- file(output_file, "w", encoding = "UTF-8")
  writeLines(html, con)
  close(con)

  message("[viz] HTML generado: ", output_file,
          " (", round(nchar(html) / 1024, 1), " KB)")
  return(invisible(output_file))
}


#' @importFrom grDevices rainbow
#' @importFrom stats complete.cases na.omit setNames
#' @importFrom utils head
NULL
