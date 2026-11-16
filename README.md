# One for All

An interactive Shiny application for RNA-seq differential expression, functional
enrichment and publication-quality figures — with a reproducibility bundle that records
**every** analysis run in a session, not just the last one.

Built for bench biologists who want the gold-standard methods (DESeq2, clusterProfiler)
without writing R, and for bioinformaticians who need to know exactly what was run.

> **Requires R 4.5.3.** One for All is developed and tested against **R 4.5.3 with
> Bioconductor 3.22**. The Windows installer ships this exact R version together with every
> package, so nothing else needs to be installed. If you run it from the source code, use
> R 4.5.3 (see [Option B](#option-b--from-this-repository)).

📖 **User manual:** [`heatmap_app/USER_MANUAL.docx`](heatmap_app/USER_MANUAL.docx)

---

## Contents

- [What it does](#what-it-does)
- [Installation](#installation)
  - [Option A — Windows installer (recommended)](#option-a--windows-installer-recommended)
  - [Option B — From this repository](#option-b--from-this-repository)
  - [Tomato annotation package](#tomato-annotation-package)
- [Folder layout](#folder-layout)
- [Quick start](#quick-start)
- [Input file format](#input-file-format)
- [Updating](#updating)
  - [Adding R packages that were not bundled](#adding-r-packages-that-were-not-bundled)
  - [Adding a new species](#adding-a-new-species)
- [The reproducibility bundle](#the-reproducibility-bundle)
- [Troubleshooting](#troubleshooting)
- [Licence and third-party components](#licence-and-third-party-components)

---

## What it does

| Area | Capability |
|---|---|
| **Import** | Excel / CSV count matrices, or Salmon `quant.sf` via tximport. Validation with plain-language warnings. |
| **Metadata** | Auto-parsed from sample names, built manually in-app, or uploaded as CSV. |
| **Differential expression** | DESeq2 with an arbitrary design formula, including interaction terms. Reference levels, count pre-filter and LFC shrinkage are all under user control and all recorded. |
| **Quality control** | PCA, MDS, sample-correlation heatmap, PC-vs-factor variance, batch correction (refused when batch and biology are confounded). |
| **Enrichment** | GO (BP / MF / CC) and KEGG over-representation, plus GSEA. Three selectable background sets. KEGG pathway diagrams via pathview. |
| **Figures** | Interactive heatmaps (ComplexHeatmap), volcano plots, Venn diagrams up to 5 sets, GO/KEGG dot, lollipop and bar plots — all exportable as PNG and PDF. |
| **Reproducibility** | One zip containing a manifest, session info, a session log, runnable re-run scripts, and the frozen gene list **and background** for every enrichment run. |

Species currently shipped: *Arabidopsis thaliana* and *Solanum lycopersicum*.
See [Adding a new species](#adding-a-new-species).

---

## Installation

Two supported routes. Pick **A** if you just want to use the application, **B** if you
want the source, or want to modify or extend it.

### Option A — Windows installer (recommended)

1. Go to the [Releases](../../releases) page of this repository.
2. Download the latest **`OneForAll_Setup_vX.Y.exe`**.
3. Run it and follow the wizard.

**What the installer does**

- Installs **per-user**, so **no administrator rights and no UAC prompt** are needed.
- Default location `%LOCALAPPDATA%\OneForAll` (changeable in the wizard).
- Bundles **R 4.5.3 (portable)** and **every required R package** (including the tomato
  annotation), so **no internet connection is needed to install or to run** — the
  exceptions are KEGG enrichment and the STRING network, which query those servers live.
- You do **not** need R installed on your computer, and an existing R installation is
  neither used nor modified.
- Optionally creates Start Menu and desktop shortcuts.

**Launching** — use the Start Menu or desktop shortcut. A console window opens (it is the
R process; leave it open) and the application appears in your default browser. Closing the
browser tab shuts the application down after a short grace period.

**Uninstalling** — *Uninstall One for All* in the Start Menu, or Windows
Settings ▸ Apps. This removes the whole install directory, including the bundled R.

> **Windows SmartScreen.** The installer is not code-signed, so Windows may show
> *"Windows protected your PC"*. Click **More info ▸ Run anyway**. If you would rather not,
> use Option B instead.

### Option B — From this repository

**Requirements**

The repository contains the application source code only. **R itself and the R packages
are not stored in the repository** (they are several GB); you install R 4.5.3 yourself and
the packages are downloaded on the first launch.

**Requirements**

- Windows 10 or later (the launcher is a `.bat`; on Linux/macOS see the note below)
- **R 4.5.3** — the version the application is developed and tested with. Other versions
  may install different package versions and give different results or errors.
- About 4–6 GB of free disk space (R + packages)
- An internet connection **for the first launch only**, so packages can be fetched

**Steps**

1. Download the repository — *Code ▸ Download ZIP*, or:

```bash
git clone https://github.com/Fblasio/OneForAll-RNAseq.git
```

2. **Install R 4.5.3 inside the repository folder.** Download the Windows installer
   [`R-4.5.3-win.exe`](https://cran.r-project.org/bin/windows/base/old/4.5.3/R-4.5.3-win.exe)
   from CRAN and, when it asks for the destination folder, choose
   `OneForAll-RNAseq\R-4.5.3`. To leave any R you already have untouched, untick
   *Save version number in registry* and *Associate R with .RData files* in the
   *Additional tasks* step. The result must match [Folder layout](#folder-layout): the
   launcher looks for R at `..\R-4.5.3\bin\Rscript.exe`, one level **above** the
   application folder.

3. *(Tomato only)* install the tomato annotation package — see
   [Tomato annotation package](#tomato-annotation-package).

4. Launch, from the `heatmap_app` folder:

```bash
lanzar_app.bat
```

   or, from any R 4.5.3 session started in `heatmap_app`:

```bash
Rscript run_app.R
```

**First launch** takes several minutes (up to ~30 min on a slow connection). `run_app.R`
checks the package list, installs anything missing into the local `rlib/` folder via
`BiocManager` (Bioconductor 3.22), then starts the application. Later launches take about
10–15 seconds.

> **Using an R 4.5.3 installed elsewhere.** Edit the path in `lanzar_app.bat`, or run
> `Rscript run_app.R` with that R. Note that `run_app.R` always prepends the app's own
> `rlib/` to `.libPaths()`, so packages install there rather than into your system library.

> **Linux / macOS.** The application itself is platform-independent R, but only the Windows
> launcher is provided and these platforms are not tested. Run `Rscript run_app.R` from the
> `heatmap_app` directory with R 4.5.3.

### Tomato annotation package

The tomato workflow needs **`org.Slycopersicum.eg.db`**, a third-party annotation package
(Laura Rodriguez, UPM; Artistic-2.0) that is **not available from CRAN or Bioconductor**, so
`run_app.R` cannot install it automatically. The Windows installer already includes it.

For Option B, download `org.Slycopersicum.eg.db.zip` from the
[Releases](../../releases) page and unzip it into `heatmap_app/rlib/`, so that the file
`heatmap_app/rlib/org.Slycopersicum.eg.db/DESCRIPTION` exists. Arabidopsis works without it.

---

## Folder layout

The launcher depends on this structure — R must sit **beside** the application folder,
not inside it. The repository already has this layout; only `R-4.5.3/` and `rlib/` are
missing until you install R and launch the app for the first time.

```
OneForAll-RNAseq/            (or the installer's folder, %LOCALAPPDATA%\OneForAll)
├── R-4.5.3/                 R 4.5.3 — bundled by the installer, or installed by you (Option B)
│   └── bin/Rscript.exe
└── heatmap_app/
    ├── app.R                UI + top-level server
    ├── global.R             package loading, shared helpers, module sourcing
    ├── run_app.R            dependency check + launcher
    ├── lanzar_app.bat       Windows entry point
    ├── modules/             Shiny modules (DE, GO, KEGG, GSEA, Venn, repro, …)
    ├── utils/               enrichment, normalisation, validation, species config
    ├── rlib/                private R library — created on first launch (not in git)
    ├── data/                side ID maps (e.g. tomato KEGG/STRING mapping)
    ├── images/              icons and UI assets
    └── USER_MANUAL.docx     user manual
```

---

## Quick start

1. **Choose the species** in the sidebar, then click *Switch species & reload*.
2. Pick the **RNAseq** feature.
3. **Upload your count matrix** (step 1 in the sidebar) — see the format below.
4. **Describe the experiment.** Either let the app parse the sample names, build the table
   manually, or upload a metadata CSV.
5. **Set the design formula** (for example `~ genotype * treatment`), choose the reference
   level for each factor and the minimum count filter, then **Build the DESeq2 object**.
6. **Run a contrast.** Choose the factor and the two levels to compare — or select
   *Interaction term* if the design contains one. Set the adjusted p-value and fold-change
   cut-offs, and decide whether to apply LFC shrinkage.
7. **Explore**: results table, volcano plot, PCA, heatmaps, GO and KEGG enrichment on the
   significant genes, Venn diagrams against earlier gene lists.
8. **Download the reproducibility bundle** before you close the session.

---

## Input file format

A gene-level count matrix: **gene identifiers in the first column**, one column per sample.

| Gene ID | Sample_A_R1 | Sample_A_R2 | Sample_B_R1 | Sample_B_R2 |
|---|---|---|---|---|
| AT1G01010 | 44 | 67 | 26 | 35 |
| AT1G01020 | 87 | 105 | 127 | 139 |

- Accepted formats: `.xlsx`, `.xls`, `.csv`
- **Raw, un-normalised integer counts** for DESeq2. The application will warn you if the
  matrix looks already normalised — DESeq2 models counting noise and needs raw counts.
- Gene IDs must match the selected species (`AT1G01010` for Arabidopsis,
  `Solyc01g005000.3` for tomato).
- Sample names containing `_` parse automatically into factors and a replicate
  (`Col_warm_R1` → `Col` / `warm` / `R1`). If your names do not follow that pattern, use
  *Create manual table* or upload a metadata CSV.

A metadata CSV needs a `sample` column matching the count-matrix column names exactly, plus
one column per experimental factor:

```
sample,genotype,treatment
Col_cont_R1,Col,CONT
Col_warm_R1,Col,HOT
```

---

## Updating

### Updating the application

- **Installed via the .exe** — download the newer `OneForAll_Setup_vX.Y.exe` from Releases
  and run it. It installs over the previous version. Your `rlib/` and any packages you added
  are replaced, so re-apply local modifications afterwards.
- **Installed from the repository** — `git pull`, or download the new ZIP. Keep your existing
  `rlib/` folder to avoid re-downloading every package.

Analysis outputs are never written inside the install directory, so updating cannot destroy
your results.

### Adding R packages that were not bundled

The bundled `rlib/` carries what the shipped features need. If a new release adds a
dependency, or you want an optional package, install it **with the bundled R and into the
app's own library** — not into a system library:

```bash
"..\R-4.5.3\bin\Rscript.exe" -e ".libPaths('rlib'); install.packages('<pkg>', lib='rlib')"
```

For Bioconductor packages:

```bash
"..\R-4.5.3\bin\Rscript.exe" -e ".libPaths('rlib'); BiocManager::install('<pkg>', lib='rlib', update=FALSE, ask=FALSE)"
```

`run_app.R` also does this automatically for anything in its own `pkgs` list that is
missing at launch, so simply starting the application after an update usually suffices.

### Adding a new species

This needs three steps: install an annotation package, register the species, and add it to
the selector. Steps 2 and 3 are small source edits — there is no configuration file yet.

**1. Install an OrgDb annotation package into `rlib/`.**

If Bioconductor publishes one for your organism:

```bash
"..\R-4.5.3\bin\Rscript.exe" -e ".libPaths('rlib'); BiocManager::install('org.Dm.eg.db', lib='rlib', update=FALSE, ask=FALSE)"
```

If it does not (common for non-model plants), build one with
[`AnnotationForge::makeOrgPackage()`](https://bioconductor.org/packages/AnnotationForge/)
and install the resulting package into `rlib/`. Include an `ENTREZID` column if you want
KEGG enrichment to work, since KEGG keys on NCBI Gene IDs.

**2. Register the species** — add a branch to `get_species_config()` in
`utils/species_config.R`, following the two that are already there:

```r
} else if (species == "fruitfly") {
  library(org.Dm.eg.db)
  ORG_DB              <- org.Dm.eg.db
  ORG_DB_PKG          <- "org.Dm.eg.db"
  SPECIES_LABEL       <- "Fruit fly"
  SPECIES_COMMON      <- "Drosophila melanogaster"
  KEY_TYPE            <- "FLYBASE"     # keytype your IDs use
  SYMBOL_COL          <- "SYMBOL"      # or NA_character_ if absent
  ENTREZ_COL          <- "ENTREZID"    # or NA_character_ if absent
  KEGG_CODE           <- "dme"         # KEGG 3-letter organism code
  GENE_IDS_UPPERCASE  <- FALSE         # TRUE only if IDs are case-normalised upward
  UNIPROT_TAXID       <- "7227"
  GENE_DB_LABEL       <- "FlyBase"
  GENE_DB_URL         <- "https://flybase.org/reports/"
  GENE_DETAIL_URL     <- "https://flybase.org/reports/"
  GENE_ID_PLACEHOLDER <- "FBgn0000008\nFBgn0000015\nFBgn0000017"
  canonicalize_tomato_ids <- NULL
  collapse_duplicate_rows <- NULL
}
```

**3. Add it to the selector** — extend the `choices` of `species_switch` in `app.R`
(around line 143):

```r
choices = c("Arabidopsis thaliana"    = "arabidopsis",
            "Solanum lycopersicum"    = "tomato",
            "Drosophila melanogaster" = "fruitfly")
```

Restart the application. If KEGG returns nothing, the usual cause is that your IDs do not
map to Entrez — the application reports the mapped fraction, so check that number first.

---

## The reproducibility bundle

*Download analysis bundle* produces `OneForAll_analysis_<timestamp>.zip` containing:

| File | Contents |
|---|---|
| `manifest.json` | Every parameter in scope — heatmap config, DE contrasts, GO, KEGG, GSEA |
| `session_info.txt` | R, package and annotation versions, plus full `sessionInfo()` |
| `session_log.txt` | Timeline of uploads, warnings and errors |
| `rerun_DE.R` | Reproduces **every** DE contrast run in the session, pre-filter included |
| `rerun_GO.R`, `rerun_KEGG.R`, `rerun_GSEA.R` | Reproduce every enrichment run |
| `gene_lists/` | The exact gene vector **and background** frozen for each run |

Because the count pre-filter is stored on the DESeq2 object and written into `rerun_DE.R`
as executable code, the re-run script reproduces its own adjusted p-values — the pre-filter
determines independent filtering, which determines every `padj`.

---

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `there is no package called '…'` | You are running a different R. Use the bundled `..\R-4.5.3\bin\Rscript.exe`, or make sure `rlib/` is on `.libPaths()`. |
| The console window opens and closes immediately | An error before startup. Run `lanzar_app.bat` from a command prompt to keep the message on screen. |
| Browser does not open | Go to the `http://127.0.0.1:<port>` address printed in the console. |
| KEGG enrichment returns nothing | Needs internet — it queries the KEGG servers live. If you are online, check the reported ID-mapping fraction; KEGG keys on Entrez IDs. |
| STRING network is empty or fails | Needs internet — it queries string-db.org live. |
| Tomato selected: `there is no package called 'org.Slycopersicum.eg.db'` | Option B only: install the [tomato annotation package](#tomato-annotation-package). |
| *Windows protected your PC* | The installer is unsigned. **More info ▸ Run anyway**, or install from source. |
| GO/KEGG plot omits the most significant term | Plots rank terms by the chosen x-axis metric, not by adjusted p-value. Switch the X-axis to *Enrichment Factor*, or read the exported table, which is ranked by significance. |
| Two browser tabs disagree about the species | The active species is process-wide. Use one tab per session. |

---

## Licence and third-party components

**One for All** is released under the **GNU General Public License v3.0** — see
[`LICENSE`](LICENSE). Third-party R packages remain under their own licences;
see [`third_part_licenses.txt`](heatmap_app/third_part_licenses.txt). The Windows installer
redistributes R 4.5.3 (GPL-2 | GPL-3), whose source code is available from
[CRAN](https://cran.r-project.org/src/base/R-4/).

Two restrictions do **not** come from the software licence and are worth knowing:

- **KEGG.** Pathway enrichment queries the KEGG REST API and pathway diagrams are drawn
  with `pathview`. KEGG data is free for academic use; **commercial use requires a licence
  from Kanehisa Laboratories** (<https://www.kegg.jp/kegg/legal.html>).
- **`org.Slycopersicum.eg.db`**, the tomato annotation, is a third-party package
  (Artistic-2.0) by Laura Rodriguez (UPM). It is required only for the tomato workflow.

## Citation

If you use One for All in published work, please cite the repository
(<https://github.com/Fblasio/OneForAll-RNAseq>) and the version you used. *(A full citation
and DOI will be added on publication.)*
