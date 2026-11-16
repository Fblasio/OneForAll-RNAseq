# run_app.R: dependency check + launcher for One for All.
# Run from the application folder (heatmap_app/), e.g. via lanzar_app.bat.

# 1. usar siempre una librería de usuario dentro de la carpeta de la app
app_dir  <- normalizePath(".", winslash = "/", mustWork = FALSE)
user_lib <- file.path(app_dir, "rlib")

if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

# 2. configurar mirror de CRAN
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 3. asegurar BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib)
}

# 4. paquetes necesarios
pkgs <- c(
  "shiny", "colourpicker",
  "circlize", "readxl", "dplyr", "ggplot2", "DT",
  "digest", "jsonlite",
  "ComplexHeatmap", "InteractiveComplexHeatmap",
  "org.At.tair.db", "AnnotationDbi",
  "clusterProfiler", "DOSE", "enrichplot", "shinyjs",
  "shinycssloaders", "pathview", "base64enc",
  "httr", "DESeq2", "openxlsx", "ggrepel",
  "later", "ggVennDiagram", "sf", "sp", "tximport",
  "ggvenn", "stringr", "zip", "ashr",
  "tibble", "scales", "cluster", "limma", "cicerone"
)

# org.Slycopersicum.eg.db (tomato annotation) is not on CRAN or Bioconductor,
# so it cannot be installed automatically. It is only needed for the tomato
# workflow; see README.md ("Tomato annotation package").
if (!requireNamespace("org.Slycopersicum.eg.db", quietly = TRUE)) {
  message("Note: 'org.Slycopersicum.eg.db' is not installed. The tomato workflow ",
          "will not work until it is copied into: ", user_lib,
          " (see README.md, 'Tomato annotation package').")
}

installed <- rownames(installed.packages())
new <- pkgs[!pkgs %in% installed]

if (length(new) > 0) {
  message("Instalando paquetes faltantes en: ", user_lib)
  BiocManager::install(new, ask = FALSE, update = FALSE, lib = user_lib)
} else {
  message("Todos los paquetes ya están instalados ✔")
}

# 5. cargar paquetes
for (p in pkgs) {
  ok <- require(p, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
  if (!ok) stop("No se pudo cargar el paquete: ", p)
}


# 6. lanzar la app (app.R está en esta misma carpeta)
message("Launching app...")
Sys.setenv(CHESKI_AUTOCLOSE = "1")
shiny::runApp("app.R", launch.browser = TRUE)

# runApp() only returns once stopApp() has fired (real browser close,
# not a refresh). Exit cleanly with status 0 so lanzar_app.bat can tell
# this apart from a crash.
message("App stopped (browser closed).")
quit(save = "no", status = 0)