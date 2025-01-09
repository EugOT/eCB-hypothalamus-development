---
title: "Metabolic and opioid system expression analysis of Pubertal and Adult Hypothalamus with focus on PVN"
authors:
  - name: Evgenii O. Tretiakov, PhD
    affiliation: 
      - Department of Molecular Neurosciences, Center for Brain Research, Medical University of Vienna, Vienna A-1090, Austria
    roles: 
      - Conceptualization
      - Formal analysis
      - Data curation
      - Investigation
      - Methodology
      - Software
      - Visualization
      - Resources
      - Writing
    corresponding: true
    orcid: 0000-0001-5920-2190
    email: evgenii.tretiakov@meduniwien.ac.at
  - name: Tibor Harkany, PhD
    affiliation: 
      - Department of Molecular Neurosciences, Center for Brain Research, Medical University of Vienna, Vienna A-1090, Austria
      - Division of Molecular and Cellular Neuroendocrinology, Department of Neuroscience, Biomedicum 7D, Karolinska Institutet, Solna SE-17165, Sweden
    roles: 
      - Conceptualization
      - Funding acquisition
      - Resources
      - Supervision
      - Project administration
      - Writing
    corresponding: true
    orcid: 0000-0002-6637-5900
    email: tibor.harkany@meduniwien.ac.at
date: "2025-01-09"
format:
  html: 
    default: true
    toc: true
    df-print: paged
    code-fold: true
    fig-width: 9
    fig-height: 12
    fig-format: retina
    fig-responsive: true
    fig-dpi: 300
    embed-resources: true
  pdf: 
    default: true
    colorlinks: true
    fontsize: 12pt
execute:
  keep-md: true
  echo: true
  error: false
  message: false
  warning: false
  debug: false
knitr:
  opts_chunk:
    autodep: true
    fig.align: center
    fig.retina: 2
    fig.width: 14
    fig.height: 12
bibliography: references.bib
---







## Setup parameters




::: {.cell layout-align="center"}

```{.r .cell-code}
# Load tidyverse infrastructure packages
suppressPackageStartupMessages({
  library(future)
  library(here)
  library(tidyverse)
  library(magrittr)
  library(stringr)
  library(skimr)
  library(RColorBrewer)
  library(viridis)
})


# Load packages for scRNA-seq analysis and visualisation
suppressPackageStartupMessages({
  library(UpSetR)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(ggstatsplot)
  library(Seurat)
  library(SeuratWrappers)
  library(scCustomize)
})
```
:::




### Set paths




::: {.cell layout-align="center"}

```{.r .cell-code}
src_dir <- here("code")
data_dir <- here("data")
output_dir <- here("output")
plots_dir <- here(output_dir, "figures/")
tables_dir <- here(output_dir, "tables/")
```
:::




### Load helper functions and gene-sets




::: {.cell layout-align="center"}

```{.r .cell-code}
source(here(src_dir, "genes.R"))
source(here(src_dir, "functions.R"))
```
:::




### Set fixed variables




::: {.cell layout-align="center"}

```{.r .cell-code}
# set seed
reseed <- 42
set.seed(seed = reseed)

# Parameters for parallel execution
n_cores <- parallelly::availableCores()/2 - 1
plan("multisession", workers = n_cores)
options(
  future.globals.maxSize = 100000 * 1024^2,
  future.rng.onMisuse = "ignore"
)
plan()
```

::: {.cell-output .cell-output-stdout}

```
multisession:
- args: function (..., workers = c(cgroups2.cpu.max = 5), envir = parent.frame())
- tweaked: TRUE
- call: plan("multisession", workers = n_cores)
```


:::

```{.r .cell-code}
# ggplot2 theme
theme_set(ggmin::theme_powerpoint())
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
bioproject <- "PRJNA547712"
project <- "kim2020_Hypoth-dev"
cb_fpr <- 0.001
low_cutoff_gene <- 500
high_cutoff_gene <- NULL
high_cutoff_gene <- 5000
low_cutoff_umis <- NULL
low_cutoff_umis <- -Inf
high_cutoff_umis <- 25000
high_cutoff_pc_mt <- 15
high_cutoff_pc_ribo <- 20
high_cutoff_pc_hb <- 0.1
high_cutoff_doublet_score <- 0.33
high_cutoff_complexity <- 0.85
connectivity_model <- "min_tree"
k <- 10
metric <- "euclidean"
signature <- 100
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
opioid_system_genes <- c(
  # Classical opioid receptors
  "Oprd1", # Delta opioid receptor
  "Oprk1", # Kappa opioid receptor
  "Oprl1", # Nociceptin/orphanin FQ receptor
  "Oprm1", # Mu opioid receptor

  # Processing enzymes
  "Pcsk1", # Proprotein convertase 1
  "Pcsk2", # Proprotein convertase 2

  # Endogenous opioid precursors
  "Pdyn", # Prodynorphin
  "Penk", # Proenkephalin
  # "Pomc",   # Proopiomelanocortin
  "Pnoc" # Prepronociceptin
)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
metabolic_signaling_genes <- c(
  # Receptor tyrosine kinases and ligands
  "Alk", # Anaplastic lymphoma kinase - neural development, metabolism
  "Fam150a", # ALK ligand 1/Augmentor-β - ALK receptor activator
  "Fam150b", # ALK ligand 2/Augmentor-α - ALK receptor activator

  # Melanocortin system
  "Mc3r", # Melanocortin 3 receptor - energy homeostasis, inflammation
  "Mc4r", # Melanocortin 4 receptor - appetite control, energy balance

  # Metabolic hormone receptors
  # "Lepr", # Leptin receptor - energy balance, satiety
  # "Insr", # Insulin receptor - glucose homeostasis
  # "Igf1r",    # Insulin-like growth factor 1 receptor - growth, development

  # Signaling adaptors/regulators
  "Lmo4" # LIM domain only 4 - transcriptional regulation, metabolism
  # "Irs1", # Insulin receptor substrate 1 - insulin signaling
  # "Irs4" # Insulin receptor substrate 4 - insulin/leptin signaling
)
```
:::




## Load data from Xu et al (2020)




::: {.cell layout-align="center"}

```{.r .cell-code}
expr_counts <- read_csv(here(data_dir, "GSE148568_compiled_data.csv"))
ercc_counts <- read_csv(here(data_dir, "GSE148568_compiled_data_erccs.csv"))
meta_data <- read_csv(here(data_dir, "GSE148568_cell_metadata_after_qc.csv")) |> dplyr::select(-1)

rownames(meta_data) <- meta_data$sample_id

gene.names <- expr_counts |> pull(1)
expr_counts <- expr_counts |> dplyr::select(-1) |> as.matrix()
rownames(expr_counts) <- gene.names
# expr_counts <- expr_counts[ , meta_data$sample_id]

ercc.names <- ercc_counts |> pull(1)
ercc_counts <- ercc_counts |> dplyr::select(-1) |> as.matrix()
rownames(ercc_counts) <- ercc.names
# ercc_counts <- ercc_counts[ , meta_data$sample_id]

library(SingleCellExperiment)
sce.xu <- SingleCellExperiment(
  assays = list(counts = expr_counts),
  altExps = list(ERCC = SummarizedExperiment(assays = list(counts = ercc_counts))))#,
  # colData = metadata)
# altExp(sce.xu, "ERCC") <- SummarizedExperiment(ercc_counts)

# rm(expr_counts, ercc_counts, gene.names, ercc.names)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
library(scater)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
ensdb_genes <- genes(EnsDb.Mmusculus.v79)
MT_names <- ensdb_genes[seqnames(ensdb_genes) == "MT"]$symbol
is_mito <- rownames(sce.xu) %in% MT_names
table(is_mito)
```

::: {.cell-output .cell-output-stdout}

```
is_mito
FALSE  TRUE 
48360    37 
```


:::

```{.r .cell-code}
sce.xu_cell <- perCellQCMetrics(sce.xu, subsets=list(Mito=is_mito), use.altexps = "ERCC")
sce.xu_feature <- perFeatureQCMetrics(sce.xu)
sce.xu <- addPerCellQC(sce.xu, subsets=list(Mito=is_mito), use.altexps = TRUE)
sce.xu <- addPerFeatureQC(sce.xu)

glimpse(sce.xu_cell)
```

::: {.cell-output .cell-output-stdout}

```
Formal class 'DFrame' [package "S4Vectors"] with 6 slots
  ..@ rownames       : chr [1:960] "Q5_A01" "Q5_A02" "Q5_A03" "Q5_A04" ...
  ..@ nrows          : int 960
  ..@ elementType    : chr "ANY"
  ..@ elementMetadata: NULL
  ..@ metadata       : list()
  ..@ listData       :List of 9
  .. ..$ sum                  : num [1:960] 216068 272288 360886 168088 144742 ...
  .. ..$ detected             : num [1:960] 4771 5589 5857 4169 4246 ...
  .. ..$ subsets_Mito_sum     : num [1:960] 8483 9700 14176 5376 6246 ...
  .. ..$ subsets_Mito_detected: num [1:960] 10 11 12 9 10 10 12 12 13 9 ...
  .. ..$ subsets_Mito_percent : num [1:960] 3.93 3.56 3.93 3.2 4.32 ...
  .. ..$ altexps_ERCC_sum     : num [1:960] 34733 25722 34362 36404 36779 ...
  .. ..$ altexps_ERCC_detected: num [1:960] 30 34 33 27 35 29 33 30 29 26 ...
  .. ..$ altexps_ERCC_percent : num [1:960] 13.85 8.63 8.69 17.8 20.26 ...
  .. ..$ total                : num [1:960] 250801 298010 395248 204492 181521 ...
```


:::

```{.r .cell-code}
glimpse(sce.xu_feature)
```

::: {.cell-output .cell-output-stdout}

```
Formal class 'DFrame' [package "S4Vectors"] with 6 slots
  ..@ rownames       : chr [1:48397] "4933401J01Rik" "Gm26206" "Xkr4" "Gm18956" ...
  ..@ nrows          : int 48397
  ..@ elementType    : chr "ANY"
  ..@ elementMetadata: NULL
  ..@ metadata       : list()
  ..@ listData       :List of 2
  .. ..$ mean    : num [1:48397] 0 0 0.939 0 0.629 ...
  .. ..$ detected: num [1:48397] 0 0 6.67 0 6.46 ...
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
qc.lib2 <- isOutlier(sce.xu_cell$sum, log=TRUE, type="lower")
attr(qc.lib2, "thresholds")
```

::: {.cell-output .cell-output-stdout}

```
  lower  higher 
5102.72     Inf 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
qc.nexprs2 <- isOutlier(sce.xu_cell$detected, log=TRUE, type="lower")
attr(qc.nexprs2, "thresholds")
```

::: {.cell-output .cell-output-stdout}

```
   lower   higher 
2061.251      Inf 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
qc.spike2 <- isOutlier(sce.xu_cell$altexps_ERCC_percent, type="higher")
attr(qc.spike2, "thresholds")
```

::: {.cell-output .cell-output-stdout}

```
   lower   higher 
    -Inf 28.53717 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
qc.mito2 <- isOutlier(sce.xu_cell$subsets_Mito_percent, type="higher")
attr(qc.mito2, "thresholds")
```

::: {.cell-output .cell-output-stdout}

```
   lower   higher 
    -Inf 5.981342 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2), SpikeProp=sum(qc.spike2, na.rm = T), MitoProp=sum(qc.mito2, na.rm = T), Total=sum(discard2))
```

::: {.cell-output .cell-output-stdout}

```
DataFrame with 1 row and 5 columns
    LibSize    NExprs SpikeProp  MitoProp     Total
  <integer> <integer> <integer> <integer> <integer>
1        93       109        84        20       130
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
reasons <- quickPerCellQC(sce.xu_cell, sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))
```

::: {.cell-output .cell-output-stdout}

```
             low_lib_size            low_n_features high_subsets_Mito_percent 
                       93                       109                        NA 
high_altexps_ERCC_percent                   discard 
                       NA                       130 
```


:::

```{.r .cell-code}
sce.xu$discard <- reasons$discard
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
plotColData(sce.xu, x="sum", y="subsets_Mito_percent", colour_by="discard")
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-qc-discard-mito-1.png){#fig-xu2020-PVN-qc-discard-mito fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
plotColData(sce.xu, x="sum", y="detected", colour_by="discard")
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-qc-discard-detected-1.png){#fig-xu2020-PVN-qc-discard-detected fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
plotColData(sce.xu, x="altexps_ERCC_percent", y="subsets_Mito_percent", colour_by="discard")
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-qc-discard-ERCC-mito-1.png){#fig-xu2020-PVN-qc-discard-ERCC-mito fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sce.xu <- sce.xu[ , meta_data$sample_id]
colData(sce.xu) %<>% cbind(meta_data)
library(scales)
plotColData(sce.xu, x="sum", y="detected", colour_by="discard", other_fields = "core_label") + 
  facet_wrap(~core_label) + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-qc-discard-core_label-1.png){#fig-xu2020-PVN-qc-discard-core_label fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
plotColData(sce.xu, x="sum", y="detected", colour_by="discard", other_fields = "batch_label") + 
  facet_wrap(~batch_label)  + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-qc-discard-batch_label-1.png){#fig-xu2020-PVN-qc-discard-batch_label fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
plotHighestExprs(sce.xu, exprs_values = "counts",  colour_cells_by="detected")
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-qc-highest-expression-1.png){#fig-xu2020-PVN-qc-highest-expression fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sce.xu.keep_feature <- nexprs(
  sce.xu,
  byrow = TRUE,
  detection_limit = 1) >= 2
rowData(sce.xu)$discard <- ! sce.xu.keep_feature
table(rowData(sce.xu)$discard)
```

::: {.cell-output .cell-output-stdout}

```

FALSE  TRUE 
19700 28697 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sce.xu <- sce.xu[!rowData(sce.xu)$discard, 
                 !colData(sce.xu)$discard]

library(scran)
sce.xu <- computeSumFactors(sce.xu)
sce.xu <- logNormCounts(sce.xu)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sce.xu <- runPCA(
  sce.xu,
  ncomponents = 30)
sce.xu$cluster1_label %<>% as_factor()
names(sce.xu$cluster1_color) <- sce.xu$cluster1_label
sce.xu$cluster1_color %<>% as_factor()
sce.xu$cluster1_color <- fct_reorder(sce.xu$cluster1_color, sce.xu$cluster1_id)
plotPCA(sce.xu, ncomponents = 4, colour_by = "cluster1_label",size_by = "detected", shape_by = "batch_label") * scale_colour_discrete(type = levels(sce.xu$cluster1_color))

sce.xu <- runUMAP(sce.xu, pca = 10)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-qc-pca-1.png){#fig-xu2020-PVN-qc-pca fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
dec <- modelGeneVar(sce.xu)
sce.xu.hvg <- getTopHVGs(sce.xu, prop=0.1)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
genes_stk_vln_repr <- unique(c(
  "Sim1",
  "Slc17a6",
  "Slc32a1",
  "Fos",
  "Gad2",
  "Oxt",
  "Pdyn",
  "Avp",
  "Npy1r",
  "Sst",
  "Trh",
  "Mc3r",
  "Mc4r",
  "Crh",
  "Scgn",
  "Fam150b",
  "Alk",
  "Ntng1",
  "Reln",
  "Penk",
  opioid_system_genes
))
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
plotGroupedHeatmap(
  sce.xu,
  features=genes_stk_vln_repr,
  group=c("cluster1_label"), 
  scale = F,
  center=F)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-gr-heatmap-1.png){#fig-xu2020-PVN-gr-heatmap fig-align='center' width=1800}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
col_mat <- levels(sce.xu$cluster1_color)
names(col_mat) <- 1:15
plotHeatmap(
  sce.xu,
  features = genes_stk_vln_repr,
  colour_columns_by = "cluster1_label",
  column_annotation_colours = list(cluster1_label = col_mat),
  order_columns_by = c("cluster1_label", "Slc17a6", "Gad2")
  #order_columns_by = c("cluster1_label", "Crh", "Trh", "Fam150b", "Alk", "Scgn")
  )
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-cluster-heatmap-1.png){#fig-xu2020-PVN-cluster-heatmap fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
col_mat <- levels(sce.xu$cluster1_color)
names(col_mat) <- 1:15
plotHeatmap(
  sce.xu,
  features = genes_stk_vln_repr,
  colour_columns_by = "cluster1_label",
  column_annotation_colours = list(cluster1_label = col_mat),
  order_columns_by = c("Crh", "Fam150b", "Alk", "Scgn", "Trh", "Sst", "Avp", "Oxt")
  )
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-sort-heatmap-1.png){#fig-xu2020-PVN-sort-heatmap fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
plotExpression(
  sce.xu,
  features = genes_stk_vln_repr,
  x = "cluster1_label",
  colour_by = "cluster1_label", ncol = 4) * 
  scale_colour_discrete(type = levels(sce.xu$cluster1_color))
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-vln-1.png){#fig-xu2020-PVN-vln fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.xu <- Seurat::as.Seurat(
  sce.xu,
  counts = "counts",
  data = NULL)
srt.xu <- RenameAssays(srt.xu, assay.name = "originalexp", new.assay.name = "RNA")
DefaultAssay(srt.xu) <- "RNA"
srt.xu <- Store_Palette_Seurat(seurat_object = srt.xu, palette = rev(brewer.pal(n = 11, name = "Spectral")), palette_name = "expr_Colour_Pal")
srt.xu <- Store_Palette_Seurat(seurat_object = srt.xu, palette = levels(sce.xu$cluster1_color), palette_name = "cluster_Colour_Pal")
names(srt.xu@misc$cluster_Colour_Pal) <- levels(sce.xu$cluster1_label)

srt.xu$cluster1_label %<>% as.numeric()
srt.xu$cluster1_label %<>% as_factor()

srt.xu.npcs <- 30
Idents(srt.xu) <- "cluster1_label"
srt.xu@meta.data |> glimpse()
```

::: {.cell-output .cell-output-stdout}

```
Rows: 693
Columns: 45
$ orig.ident            <fct> Q5, Q5, Q5, Q5, Q5, Q5, Q5, Q5, Q5, Q5, Q5, Q5, …
$ nCount_originalexp    <dbl> 215968, 272258, 360826, 168068, 192132, 374445, …
$ nFeature_originalexp  <int> 4749, 5561, 5819, 4152, 4362, 5991, 5889, 6900, …
$ sum                   <dbl> 216068, 272288, 360886, 168088, 192206, 374496, …
$ detected              <dbl> 4771, 5589, 5857, 4169, 4384, 6034, 5944, 6967, …
$ subsets_Mito_sum      <dbl> 8483, 9700, 14176, 5376, 8272, 13611, 15283, 175…
$ subsets_Mito_detected <dbl> 10, 11, 12, 9, 10, 12, 12, 13, 9, 14, 10, 9, 8, …
$ subsets_Mito_percent  <dbl> 3.926079, 3.562405, 3.928110, 3.198325, 4.303716…
$ altexps_ERCC_sum      <dbl> 34733, 25722, 34362, 36404, 32207, 32061, 32897,…
$ altexps_ERCC_detected <dbl> 30, 34, 33, 27, 29, 33, 30, 29, 26, 27, 32, 32, …
$ altexps_ERCC_percent  <dbl> 13.848828, 8.631254, 8.693782, 17.802163, 14.351…
$ total                 <dbl> 250801, 298010, 395248, 204492, 224413, 406557, …
$ discard               <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,…
$ sample_id             <chr> "Q5_A01", "Q5_A02", "Q5_A03", "Q5_A04", "Q5_A06"…
$ plate_id              <dbl> 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, …
$ plate_label           <chr> "Q5", "Q5", "Q5", "Q5", "Q5", "Q5", "Q5", "Q5", …
$ plate_suffix          <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, …
$ plate_color           <chr> "#BFA64A", "#BFA64A", "#BFA64A", "#BFA64A", "#BF…
$ totalUMIs_label       <dbl> 213736, 269821, 358377, 165957, 189989, 371913, …
$ totalUMIs_id          <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1…
$ totalUMIs_color       <chr> "#1772E7", "#1E8FFD", "#2471BB", "#1258D2", "#15…
$ totalgenes_label      <dbl> 4700, 5519, 5796, 4082, 4306, 5967, 5883, 6915, …
$ totalgenes_id         <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1…
$ totalgenes_color      <chr> "#56362A", "#A43C16", "#BE3F10", "#303D49", "#32…
$ totalerccs_label      <dbl> 34733, 25722, 34362, 36404, 32207, 32061, 32897,…
$ totalerccs_id         <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1…
$ totalerccs_color      <chr> "#D2410B", "#3F342F", "#CC400C", "#ED4304", "#A9…
$ cluster1_label        <fct> 2, 1, 2, 8, 8, 1, 1, 1, 15, 6, 1, 1, 14, 14, 2, …
$ cluster1_id           <dbl> 2, 1, 2, 8, 8, 1, 1, 1, 15, 6, 1, 1, 14, 14, 2, …
$ cluster1_color        <fct> #6ED24F, #CF4C30, #6ED24F, #5D3C25, #5D3C25, #CF…
$ batch_label           <chr> "Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q"…
$ batch_id              <dbl> 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, …
$ batch_color           <chr> "#7D9DA6", "#7D9DA6", "#7D9DA6", "#7D9DA6", "#7D…
$ memb_label            <dbl> 10.0000000, 10.0000000, 10.0000000, 1.3946248, 4…
$ memb_id               <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1…
$ memb_color            <chr> "#FF0000", "#FF0000", "#FF0000", "#1050CB", "#2D…
$ core_label            <dbl> 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, …
$ core_id               <dbl> 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, …
$ core_color            <chr> "#FF0000", "#FF0000", "#FF0000", "#00008B", "#FF…
$ wellposition_id       <chr> "A01", "A02", "A03", "A04", "A06", "A07", "A08",…
$ sizeFactor            <dbl> 1.0092536, 1.2989806, 1.6779947, 0.8693471, 0.98…
$ nCount_ERCC           <dbl> 34733, 25722, 34362, 36404, 32207, 32061, 32897,…
$ nFeature_ERCC         <int> 30, 34, 33, 27, 29, 33, 30, 29, 26, 27, 32, 32, …
$ nCount_RNA            <dbl> 215968, 272258, 360826, 168068, 192132, 374445, …
$ nFeature_RNA          <int> 4749, 5561, 5819, 4152, 4362, 5991, 5889, 6900, …
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.xu <- FindVariableFeatures(srt.xu, selection.method = "vst", nfeatures = 1000)

srt.xu.hvg <- VariableFeatures(srt.xu)
var_regex <- "^Hla-|^Ig[hjkl]|^Rna|^mt-|^Rp[sl]|^Hb[^(p)]|^Gm"
srt.xu.hvg <- srt.xu.hvg[str_detect(pattern = var_regex, string = srt.xu.hvg, negate = TRUE)]

srt.xu.keep_genes <-
  c(gene_int, srt.xu.hvg) %>%
  unique() %>%
  .[!. %in% housekeeping_mouse] %>%
  .[!. %in% sex_genes] %>%
  .[!. %in% c(stress_genes, ieg_gene_list$Mus_musculus_IEG)]
glimpse(srt.xu.keep_genes)
```

::: {.cell-output .cell-output-stdout}

```
 chr [1:917] "Adcyap1r1" "Avpr1a" "Calcr" "Calcrl" "Cckar" "Cckbr" "Cntfr" ...
```


:::

```{.r .cell-code}
srt.xu.hvg <- srt.xu.hvg[srt.xu.hvg %in% srt.xu.keep_genes]


srt.xu.top20 <- head(srt.xu.hvg, 20)
srt.xu.top20 
```

::: {.cell-output .cell-output-stdout}

```
 [1] "Sst"      "Trh"      "Cartpt"   "Penk"     "Tac1"     "Spink8"  
 [7] "Gal"      "Avp"      "Sncg"     "Crh"      "Th"       "Prdm8"   
[13] "Prph"     "Cbln2"    "Arhgap36" "Tcf7l2"   "Unc13c"   "Hspa1a"  
[19] "Gad2"     "Lhx1os"  
```


:::

```{.r .cell-code}
plot1 <- VariableFeaturePlot(srt.xu)
LabelPoints(plot = plot1, points = srt.xu.top20, repel = TRUE, xnudge = 0, ynudge = 0)


srt.xu.all.genes <- rownames(srt.xu)
srt.xu <- ScaleData(srt.xu, features = srt.xu.all.genes)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-srt-init-vst-features-1.png){#fig-xu2020-PVN-srt-init-vst-features fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
library(gprofiler2)
mmus_s <- gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m <- gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

srt.xu <- CellCycleScoring(srt.xu, s.features = mmus_s, g2m.features = mmus_g2m)
table(srt.xu[[]]$Phase)
```

::: {.cell-output .cell-output-stdout}

```

 G1 G2M   S 
372  44 277 
```


:::

```{.r .cell-code}
VlnPlot(srt.xu, features = c("S.Score","G2M.Score"), cols = srt.xu@misc$cluster_Colour_Pal) & 
  theme(plot.title = element_text(size=10))
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-vln-cell-cycle-1.png){#fig-xu2020-PVN-vln-cell-cycle fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.xu <- SCTransform(
  srt.xu,
  variable.features.n = 2000,
  method = "glmGamPoi",
  ncells = ncol(srt.xu),
  return.only.var.genes = FALSE,
  vars.to.regress = c(
    "subsets_Mito_percent",
    "altexps_ERCC_percent",
    "S.Score",
    "G2M.Score"),
  verbose = F)

srt.xu.hvg <- head(VariableFeatures(srt.xu), 1000)
srt.xu.hvg <- srt.xu.hvg[str_detect(pattern = var_regex, string = srt.xu.hvg, negate = TRUE)]

srt.xu.keep_genes <-
  c(gene_int, srt.xu.hvg) %>%
  unique() %>%
  .[!. %in% housekeeping_mouse] %>%
  .[!. %in% sex_genes] %>%
  .[!. %in% c(stress_genes, ieg_gene_list$Mus_musculus_IEG)]
glimpse(srt.xu.keep_genes)
```

::: {.cell-output .cell-output-stdout}

```
 chr [1:957] "Adcyap1r1" "Avpr1a" "Calcr" "Calcrl" "Cckar" "Cckbr" "Cntfr" ...
```


:::

```{.r .cell-code}
srt.xu.hvg <- srt.xu.hvg[srt.xu.hvg %in% srt.xu.keep_genes]

srt.xu <- srt.xu %>%
  RunPCA(features = srt.xu.keep_genes, npcs = srt.xu.npcs, seed.use = reseed, verbose = FALSE)

srt.xu <- BuildClusterTree(
  object = srt.xu,
  features = srt.xu.keep_genes,
  reorder = TRUE,
  verbose = TRUE)
srt.xu@tools$BuildClusterTree$tip.label %<>% str_remove("g")

srt.xu
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
37765 features across 693 samples within 3 assays 
Active assay: SCT (17973 features, 2000 variable features)
 3 layers present: counts, data, scale.data
 2 other assays present: ERCC, RNA
 3 dimensional reductions calculated: PCA, UMAP, pca
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
PlotClusterTree(srt.xu)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-cluster-tree-1.png){#fig-xu2020-PVN-cluster-tree fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
ElbowPlot(srt.xu, ndims = srt.xu.npcs) # 21
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-PVN-pca-elbow-1.png){#fig-xu2020-PVN-pca-elbow fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
invisible(gc())
set.seed(reseed)

selected_pcs <-
  seq_len(21)

if (!file.exists(here(output_dir, "xu2020-metabolic-and-opioids-in-neurons-init-umap-search-ref.Rds"))) {
  
  source(here(src_dir, "scDEED.R"))
  library(furrr)

  permuted.srt.xu <- Permuted(srt.xu, K = 21)

  invisible(gc())
  set.seed(reseed)
  
  umap_example <- scDEED(
    input_data = srt.xu,
    K = 21,
    n_neighbors = seq(from = 5, to = 35, by = 10),
    min.dist = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.8),
    reduction.method = "umap",
    rerun = FALSE,
    permuted = permuted.srt.xu,
    default_assay = "SCT"
  )

  readr::write_rds(
    x = umap_example,
    file = here(output_dir, "xu2020-metabolic-and-opioids-in-neurons-init-umap-search-ref.Rds")
  )
} else {
  umap_example <-
    read_rds(here(output_dir, "xu2020-metabolic-and-opioids-in-neurons-init-umap-search-ref.Rds"))
}

umap_example$num_dubious
```

::: {.cell-output-display}

`````{=html}
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["n_neighbors"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["min.dist"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["number_dubious_cells"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"5","2":"0.01","3":"27"},{"1":"15","2":"0.01","3":"74"},{"1":"25","2":"0.01","3":"25"},{"1":"35","2":"0.01","3":"27"},{"1":"5","2":"0.05","3":"23"},{"1":"15","2":"0.05","3":"52"},{"1":"25","2":"0.05","3":"18"},{"1":"35","2":"0.05","3":"44"},{"1":"5","2":"0.10","3":"43"},{"1":"15","2":"0.10","3":"18"},{"1":"25","2":"0.10","3":"20"},{"1":"35","2":"0.10","3":"19"},{"1":"5","2":"0.25","3":"26"},{"1":"15","2":"0.25","3":"15"},{"1":"25","2":"0.25","3":"15"},{"1":"35","2":"0.25","3":"67"},{"1":"5","2":"0.50","3":"35"},{"1":"15","2":"0.50","3":"55"},{"1":"25","2":"0.50","3":"59"},{"1":"35","2":"0.50","3":"53"},{"1":"5","2":"0.80","3":"18"},{"1":"15","2":"0.80","3":"10"},{"1":"25","2":"0.80","3":"13"},{"1":"35","2":"0.80","3":"16"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
`````

:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
invisible(gc())
set.seed(seed = reseed)

srt.xu <-
  srt.xu |>
  FindNeighbors(
    dims = selected_pcs,
    k.param = umap_example$num_dubious |>
      dplyr::slice_min(
        order_by = c(number_dubious_cells),
        n = 1
      ) |>
      dplyr::slice_min(
        order_by = c(min.dist),
        n = 1
      ) |>
      pull(n_neighbors),
    annoy.metric = "cosine",
    n.trees = 100,
    verbose = FALSE
  )

srt.xu <-
  srt.xu |>
  RunUMAP(
    dims = selected_pcs,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    return.model = TRUE,
    n.epochs = 1000L,
    n.neighbors = umap_example$num_dubious |>
      dplyr::slice_min(
        order_by = c(number_dubious_cells),
        n = 1
      ) |>
      dplyr::slice_min(
        order_by = c(min.dist),
        n = 1
      ) |>
      pull(n_neighbors),
    min.dist = umap_example$num_dubious |>
      dplyr::slice_min(
        order_by = c(number_dubious_cells),
        n = 1
      ) |>
      dplyr::slice_min(
        order_by = c(min.dist),
        n = 1
      ) |>
      pull(min.dist),
    seed.use = reseed,
    verbose = FALSE
  )
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
ColorDimSplit(srt.xu, node = 25)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-umap-pvn-split-25-1.png){#fig-xu2020-umap-pvn-split-25 fig-align='center' width=1500}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
ColorDimSplit(srt.xu, node = 27)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-umap-pvn-split-27-1.png){#fig-xu2020-umap-pvn-split-27 fig-align='center' width=1500}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
ColorDimSplit(srt.xu, node = 29)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-umap-pvn-split-29-1.png){#fig-xu2020-umap-pvn-split-29 fig-align='center' width=1500}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
DimPlot_scCustom(srt.xu, reduction = "umap", alpha = 0.7, label = T, colors_use = srt.xu@misc$cluster_Colour_Pal)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-umap-pvn-types-1.png){#fig-xu2020-umap-pvn-types fig-align='center' width=1500}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.xu.markers <- FindAllMarkers(
  srt.xu,
  only.pos = TRUE,
  test.use = "LR",
  random.seed = reseed,
  return.thresh = 0.0005,
  verbose = F) %>%
    Add_Pct_Diff()

srt.xu.markers %<>%
  pull(gene) %>%
  gprofiler2::gconvert(
  .,
  organism = "mmusculus",
  target = "MGI",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE
) %>%
  dplyr::select(name, description) %>%
  right_join(srt.xu.markers, by = c("name" = "gene"))

srt.xu.markers %<>% 
  dplyr::rename(gene = "name")

xu.top_5_unique <- Extract_Top_Markers(
  marker_dataframe = srt.xu.markers,
  num_genes = 5,
  rank_by = "avg_log2FC",
  named_vector = FALSE,
  make_unique = TRUE)

xu.top_5_unique
```

::: {.cell-output .cell-output-stdout}

```
 [1] "Tmem37"        "Dct"           "Fzd7"          "Inmt"         
 [5] "Defb18"        "Dennd2c"       "Gm37465"       "5430420F09Rik"
 [9] "Npff"          "Plekha2"       "Sgk3"          "Igfbp2"       
[13] "Calcrl"        "Kcnj2"         "Colgalt2"      "Gm43486"      
[17] "Gm26891"       "Prrg1"         "Gm43321"       "Mfge8"        
[21] "H2-Ab1"        "Gm4847"        "C1qtnf7"       "Tmem45a"      
[25] "Pmch"          "Tmem81"        "Csn1s2b"       "Dnase2a"      
[29] "Pabpc1l"       "Dkk1"          "Gm42615"       "Inhbb"        
[33] "AC151602.15"   "Acp5"          "Gm16551"       "Onecut3"      
[37] "Onecut1"       "Six3os1"       "Npas1"         "BC049352"     
[41] "5430416O09Rik" "Sim2"          "Defb25"        "9530026P05Rik"
[45] "D130058E05Rik" "Igfbp7"        "Gbx2"          "Tcf7l2"       
[49] "Lhx9"          "Six6"          "Nr5a2"         "Isl1"         
[53] "Dlx5"          "Sox6"          "Pirt"          "Adamtsl2"     
[57] "Prrxl1"        "Gm37064"       "Ccna1"         "F2rl2"        
[61] "Gm44536"       "Gas1"          "Tle6"          "BC030499"     
[65] "Gm5067"        "Dnah11"        "Gnmt"          "Rho"          
[69] "Scgn"          "Gm43632"       "Glyat"         "4732419C18Rik"
[73] "Brs3"          "Glp2r"        
```


:::

```{.r .cell-code}
xu.pvn.to.plot <- unique(c(
    metabolic_signaling_genes, 
    opioid_system_genes[opioid_system_genes %in% srt.xu@assays$SCT@data@Dimnames[[1]]],
    "Sim1",
    "Fos",
    "Slc17a6",
    "Slc32a1",
    "Gad1",
    "Gad2",
    "Scgn",
    "Crh",
    "Trh",
    "Oxt",
    "Avp",
    "Sst",
    "Cartpt", 
    "Adcyap1",
    "Bdnf",
    "Cck",
    "Gal",
    "Ghrh",
    "Grp",
    "Nmb",
    "Nts",
    "Pmch",
    "Reln",
    "Tac1",
    "Pomc",
    "Npy1r",
    "Ntng1"))

xu.pvn.to.plot.sorted <- Extract_Top_Markers(
  marker_dataframe = srt.xu.markers |> dplyr::filter(gene %in% xu.pvn.to.plot),
  num_genes = 20,
  rank_by = "avg_log2FC",
  named_vector = FALSE,
  make_unique = TRUE)

xu.pvn.to.plot[!xu.pvn.to.plot %in% xu.pvn.to.plot.sorted]
```

::: {.cell-output .cell-output-stdout}

```
[1] "Fam150a" "Mc3r"    "Oprm1"   "Cck"     "Ghrh"    "Nmb"    
```


:::

```{.r .cell-code}
inserting_elements <- c("Nmb", "Mc3r", "Oprm1", "Ghrh")


xu.pvn.to.plot.sorted <-
  c(
    append(xu.pvn.to.plot.sorted[1:which(xu.pvn.to.plot.sorted == "Pcsk2")], inserting_elements[1]),
    append(xu.pvn.to.plot.sorted[which(xu.pvn.to.plot.sorted == "Pmch"):which(xu.pvn.to.plot.sorted == "Npy1r")] , inserting_elements[2]),
    append(xu.pvn.to.plot.sorted[which(xu.pvn.to.plot.sorted == "Mc4r"):which(xu.pvn.to.plot.sorted == "Alk")], inserting_elements[3]),
    append(xu.pvn.to.plot.sorted[which(xu.pvn.to.plot.sorted == "Cartpt"):which(xu.pvn.to.plot.sorted == "Fos")], inserting_elements[4]),
    xu.pvn.to.plot.sorted[which(xu.pvn.to.plot.sorted == "Reln"):length(xu.pvn.to.plot.sorted)]
  )

xu.pvn.to.plot <- xu.pvn.to.plot.sorted

# Print the updated vector
print(xu.pvn.to.plot)
```

::: {.cell-output .cell-output-stdout}

```
 [1] "Oxt"     "Oprk1"   "Pcsk1"   "Pomc"    "Gal"     "Sim1"    "Avp"    
 [8] "Tac1"    "Pdyn"    "Pcsk2"   "Nmb"     "Pmch"    "Nts"     "Sst"    
[15] "Npy1r"   "Mc3r"    "Mc4r"    "Grp"     "Alk"     "Oprm1"   "Cartpt" 
[22] "Lmo4"    "Fos"     "Ghrh"    "Reln"    "Ntng1"   "Oprl1"   "Gad2"   
[29] "Bdnf"    "Pnoc"    "Trh"     "Adcyap1" "Slc17a6" "Gad1"    "Slc32a1"
[36] "Fam150b" "Penk"    "Crh"     "Scgn"   
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.xu,
  reduction = "umap",
  features = xu.pvn.to.plot,
  layer = "data",
  alpha_exp = 0.6,
  max.cut = "q97",
  label = F,
  num_columns = 4
) #* NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-xu2020-pvn-feature-1.png){#fig-xu2020-pvn-feature fig-align='center' width=3600}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
Stacked_VlnPlot(srt.xu, features = xu.pvn.to.plot, colors_use = srt.xu@misc$cluster_Colour_Pal)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-stack-vln-xu2020-pvn-features-1.png){#fig-stack-vln-xu2020-pvn-features fig-align='center' width=3600}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
DoHeatmap(
  object = srt.xu,
  features = xu.pvn.to.plot,
  cells = scCustomize::Random_Cells_Downsample(
    seurat_object = srt.xu,
    num_cells = 150,
    allow_lower = T),
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 8,
  hjust = 0.5,
  vjust = -2,
  angle = 0,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.07,
  group.colors = srt.xu@misc$cluster_Colour_Pal) +
  scale_fill_viridis(
    option = "viridis",
    direction = 1) &
  NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-heatmap-xu2020-pvn-features-1.png){#fig-heatmap-xu2020-pvn-features fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
DoHeatmap(
  object = srt.xu,
  features = xu.top_5_unique,
  cells = scCustomize::Random_Cells_Downsample(
    seurat_object = srt.xu,
    num_cells = 150,
    allow_lower = T),
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 8,
  hjust = 0.5,
  vjust = -2.4,
  angle = 0,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.035,
  group.colors = srt.xu@misc$cluster_Colour_Pal) +
  scale_fill_viridis(
    option = "viridis",
    direction = 1) &
  NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-heatmap-xu2020-pvn-top-features-1.png){#fig-heatmap-xu2020-pvn-top-features fig-align='center' width=4200}
:::
:::





## Load data from Lopez JP et al (2021)




::: {.cell layout-align="center"}

```{.r .cell-code}
srt.lopez <-  schard::h5ad2seurat(here(
  "/data/1_heteroAstrocytes/PRJNA679294/data",
  "class_cello/PRJNA679294-whole_dataset-0.001-cello_annotation.h5ad"
), use.raw = TRUE)


Idents(srt.lopez) <- "ora_celltype"
table(Idents(srt.lopez))
```

::: {.cell-output .cell-output-stdout}

```

               Endothelial cells                          Neurons 
                             493                             3309 
Oligodendrocyte progenitor cells                 Oligodendrocytes 
                             383                             1275 
                       Microglia                       Astrocytes 
                             463                             2838 
                       Pericytes                  Ependymal cells 
                             257                              472 
                     Macrophages                  Mesangial cells 
                              30                               52 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.lopez,
  reduction = "Xpacmap_",
  features = xu.pvn.to.plot,
  layer = "data",
  alpha_exp = 0.3,
  max.cut = "q95",
  label = F,
  num_columns = 4
) #* NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-lopez2021-pvn-feature-all-cells-1.png){#fig-lopez2021-pvn-feature-all-cells fig-align='center' width=3600}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
DimPlot_scCustom(srt.lopez, reduction = "Xpacmap_", group.by = "ora_celltype", split.by = "libname", alpha = 0.3, label = F)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-lopez2021-pacmap-pvn-types-all-cells-1.png){#fig-lopez2021-pacmap-pvn-types-all-cells fig-align='center' width=3000}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.lopez <- subset(srt.lopez, idents = c("Neurons"), subset = libname == "control")

Idents(srt.lopez) <- "k_tree"

print(srt.lopez)
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
22835 features across 1677 samples within 1 assay 
Active assay: RNA (22835 features, 0 variable features)
 2 layers present: counts, data
 5 dimensional reductions calculated: Xpacmap_, Xpca_, Xumap_, oraestimate_, orapvals_
```


:::

```{.r .cell-code}
invisible(gc())
table(Idents(srt.lopez))
```

::: {.cell-output .cell-output-stdout}

```

  3   9  13   2  19  15  16  17 
390 219 150 446  89 180  92 111 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
DimPlot_scCustom(srt.lopez, reduction = "Xpacmap_", alpha = 0.3, label = F)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-lopez2021-pacmap-pvn-types-subset-neurons-1.png){#fig-lopez2021-pacmap-pvn-types-subset-neurons fig-align='center' width=1500}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
DimPlot_scCustom(srt.lopez, reduction = "Xumap_", alpha = 0.3, label = F)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-lopez2021-umap-pvn-types-subset-neurons-1.png){#fig-lopez2021-umap-pvn-types-subset-neurons fig-align='center' width=1500}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.lopez <- SCTransform(
  srt.lopez,
  variable.features.n = 2000,
  method = "glmGamPoi",
  ncells = ncol(srt.lopez),
  return.only.var.genes = FALSE,
  vars.to.regress = c(
    "log10GenesPerUMI",
    "percent_mito_ribo",
    "S.Score",
    "G2M.Score"),
  verbose = F)

srt.lopez.hvg <- head(VariableFeatures(srt.lopez), 1000)
srt.lopez.hvg <- srt.lopez.hvg[str_detect(pattern = var_regex, string = srt.lopez.hvg, negate = TRUE)]

srt.lopez.keep_genes <-
  c(gene_int, srt.lopez.hvg) %>%
  unique() %>%
  .[!. %in% housekeeping_mouse] %>%
  .[!. %in% sex_genes] %>%
  .[!. %in% c(stress_genes, ieg_gene_list$Mus_musculus_IEG)]
glimpse(srt.lopez.keep_genes)
```

::: {.cell-output .cell-output-stdout}

```
 chr [1:1039] "Adcyap1r1" "Avpr1a" "Calcr" "Calcrl" "Cckar" "Cckbr" "Cntfr" ...
```


:::

```{.r .cell-code}
srt.lopez.hvg <- srt.lopez.hvg[srt.lopez.hvg %in% srt.lopez.keep_genes]

srt.lopez <- srt.lopez %>%
  RunPCA(features = srt.lopez.keep_genes, npcs = srt.xu.npcs, seed.use = reseed, verbose = FALSE)

srt.lopez
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
41230 features across 1677 samples within 2 assays 
Active assay: SCT (18395 features, 2000 variable features)
 3 layers present: counts, data, scale.data
 1 other assay present: RNA
 6 dimensional reductions calculated: Xpacmap_, Xpca_, Xumap_, oraestimate_, orapvals_, pca
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
# Find transfer anchors
anchors <- FindTransferAnchors(
  reference = srt.xu,
  query = srt.lopez,
  dims = selected_pcs,
  reference.reduction = "pca"
)

# Map the query data onto the reference UMAP and transfer cell type annotations
srt.lopez <- MapQuery(
  anchorset = anchors,
  reference = srt.xu,
  query = srt.lopez,
  refdata = list(cluster1_label = "cluster1_label"), # Transfer cell type labels
  reference.reduction = "pca",
  reduction.model = "umap"
)

srt.lopez <- Store_Palette_Seurat(seurat_object = srt.lopez, palette = levels(sce.xu$cluster1_color), palette_name = "cluster_Colour_Pal")
names(srt.lopez@misc$cluster_Colour_Pal) <- levels(sce.xu$cluster1_label)

srt.lopez$predicted.cluster1_label %<>% fct(levels = levels(srt.xu)[levels(srt.xu) %in% unique(srt.lopez$predicted.cluster1_label)])
# The projected UMAP coordinates are in srt.lopez[["ref.umap"]]
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
# Plot the reference UMAP colored by cell types
p1 <- DimPlot_scCustom(
  seurat_object = srt.xu,
  reduction = "umap",
  group.by = "cluster1_label",
  pt.size = 1,
  colors_use = srt.lopez@misc$cluster_Colour_Pal,
  shuffle = TRUE,
  seed = reseed,
  alpha = 0.5,
  repel = TRUE,
  label = TRUE,
  label.size = 5
) + ggtitle("Reference: Cell Type Annotations")

# Plot the query cells projected onto the reference UMAP, colored by the predicted cell types
p2 <- DimPlot_scCustom(
  seurat_object = srt.lopez,
  reduction = "ref.umap",
  group.by = "predicted.cluster1_label",
  pt.size = 1,
  colors_use = srt.lopez@misc$cluster_Colour_Pal,
  shuffle = TRUE,
  seed = reseed,
  alpha = 0.3,
  repel = TRUE,
  label = TRUE,
  label.size = 5
) + NoLegend() + ggtitle("Query: Transferred Cell Type Labels")

# Combine the plots
p1 + p2 + plot_layout(guides = "collect")
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-types-all-cells-with-query,-1.png){#fig-types-all-cells-with-query, fig-align='center' width=3000}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
lopez.pvn.to.plot <- unique(c(
    "Oprd1",
    "Cck", 
    xu.pvn.to.plot))

FeaturePlot_scCustom(
  srt.lopez,
  reduction = "ref.umap",
  features = lopez.pvn.to.plot,
  layer = "data",
  alpha_exp = 0.3,
  max.cut = "q99",
  label = F,
  num_columns = 4
) #* NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-lopez2021-ref-umap-pvn-features-1.png){#fig-lopez2021-ref-umap-pvn-features fig-align='center' width=3600}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
Stacked_VlnPlot(
  srt.lopez,
  features = lopez.pvn.to.plot,
  group.by = "predicted.cluster1_label",
  colors_use = srt.lopez@misc$cluster_Colour_Pal)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-stack-vln-lopez2021-pvn-features-1.png){#fig-stack-vln-lopez2021-pvn-features fig-align='center' width=3600}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
DoHeatmap(
  object = srt.lopez,
  features = lopez.pvn.to.plot,
  cells = scCustomize::Random_Cells_Downsample(
    seurat_object = srt.lopez,
    num_cells = 25,
    allow_lower = T),
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  group.by = "predicted.cluster1_label",
  assay = NULL,
  label = TRUE,
  size = 8,
  hjust = 0.5,
  vjust = -2,
  angle = 0,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.07,
  group.colors = srt.lopez@misc$cluster_Colour_Pal) +
  scale_fill_viridis(
    option = "viridis",
    direction = 1) &
  NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-heatmap-lopez2021-pvn-features-1.png){#fig-heatmap-lopez2021-pvn-features fig-align='center' width=4200}
:::
:::





## Load Kim DW et al., 2020[@kim2020]




::: {.cell layout-align="center"}

```{.r .cell-code}
srt.kim <-  schard::h5ad2seurat(here(
  data_dir,
  "kim2020_combined.h5ad"
),use.raw = TRUE)

X_umap <- srt.kim@meta.data |> dplyr::select(X, Y, Z) |> as.matrix()
colnames(X_umap) <- c("UMAP_1", "UMAP_2", "UMAP_3")
rownames(X_umap) <- colnames(srt.kim)
srt.kim[["umap"]] <- CreateDimReducObject(embeddings = X_umap, key = "umap_", assay = DefaultAssay(srt.kim))
srt.kim$Age %<>% forcats::fct(levels = c(
  "E10", "E11", "E12", "E13", "E14", 
  "E15", "E16", "E17", "E18", "P0", 
  "P2", "P4", "P8", "P10", "P14", "P23", "P45"))

Idents(srt.kim) <- "Age"
srt.kim <- Store_Palette_Seurat(seurat_object = srt.kim, palette = rev(brewer.pal(n = 11, name = "Spectral")), palette_name = "expr_Colour_Pal")
```
:::




## Load Romanov et al., 2020[@romanov2020]




::: {.cell layout-align="center"}

```{.r .cell-code}
print(srt.kim)
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
27998 features across 128006 samples within 1 assay 
Active assay: RNA (27998 features, 0 variable features)
 2 layers present: counts, data
 1 dimensional reduction calculated: umap
```


:::

```{.r .cell-code}
srt.romanov.pub <- readRDS("/data/1_heteroAstrocytes/PRJNA548917/old/oldCCA_nae_srt.rds")
srt.romanov.pub <- UpdateSeuratObject(srt.romanov.pub)
Idents(srt.romanov.pub) <-
  factor(srt.romanov.pub$wtree,
    ordered = TRUE
  )

# Consistent colours and clusters names
colours_wtree <- setNames(read_lines(here(data_dir, "colours_wtree.tsv")), 1:45)

srt.romanov.pub$age <-
  Cells(srt.romanov.pub) |>
  str_split(pattern = ":", simplify = T) %>%
  .[, 1] %>%
  str_split_fixed(pattern = "_", n = 3) %>%
  .[, 3]
print(srt.romanov.pub)
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
24340 features across 51199 samples within 1 assay 
Active assay: RNA (24340 features, 3500 variable features)
 3 layers present: counts, data, scale.data
 3 dimensional reductions calculated: pca, tsne, umap
```


:::

```{.r .cell-code}
glimpse(srt.romanov.pub@meta.data)
```

::: {.cell-output .cell-output-stdout}

```
Rows: 51,199
Columns: 20
$ nGene            <int> 1652, 782, 447, 1706, 1106, 894, 727, 734, 669, 617, …
$ nUMI             <dbl> 2787, 1090, 544, 2709, 1817, 1220, 995, 1036, 920, 86…
$ orig.ident       <fct> Hypothalamus, Hypothalamus, Hypothalamus, Hypothalamu…
$ res.0.2          <chr> "23", "23", "23", "23", "23", "23", "23", "23", "23",…
$ res.0.4          <chr> "34", "34", "34", "34", "34", "34", "34", "34", "34",…
$ res.0.8          <chr> "42", "42", "42", "42", "42", "42", "42", "42", "42",…
$ res.1.2          <chr> "47", "47", "47", "47", "47", "47", "47", "47", "47",…
$ res.2            <chr> "54", "54", "54", "54", "54", "54", "54", "54", "54",…
$ tree.ident       <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
$ pro_Inter        <chr> "41", "41", "41", "41", "41", "41", "41", "41", "41",…
$ pro_Enter        <chr> "41", "41", "41", "41", "41", "41", "41", "41", "41",…
$ tree_final       <fct> 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 1…
$ subtree          <fct> 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 4…
$ prim_walktrap    <fct> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 3…
$ umi_per_gene     <dbl> 1.687046, 1.393862, 1.217002, 1.587925, 1.642857, 1.3…
$ log_umi_per_gene <dbl> 0.22712693, 0.14421974, 0.08529138, 0.20082998, 0.215…
$ nCount_RNA       <dbl> 2787, 1090, 544, 2709, 1817, 1220, 995, 1036, 920, 86…
$ nFeature_RNA     <int> 1652, 782, 447, 1706, 1106, 894, 727, 734, 669, 617, …
$ wtree            <fct> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 3…
$ age              <chr> "P23", "3P2", "3P2", "P2", "P2", "P2", "P2", "P2", "P…
```


:::

```{.r .cell-code}
table(Idents(srt.romanov.pub))
```

::: {.cell-output .cell-output-stdout}

```

    1     2     3     4     5     6     7     8     9    10    11    12    13 
 2344  8146   395   402  3234   712   552   374   259   952 13727  1615   765 
   14    15    16    17    18    19    20    21    22    23    24    25    26 
  832  1244   792   590   808  2486  1683   628  1039  1750   292   394   547 
   27    28    29    30    31    32    33    34    35    36    37    38    39 
  391   407   507    93    81   402   143   701   222   353   324    73    78 
   40    41    42    43    44    45 
  328   190    73    37   179    55 
```


:::

```{.r .cell-code}
srt.romanov.pub %<>% RenameIdents(object = ., `43` = "mneOXY")
srt.romanov.pub %<>% RenameIdents(object = ., `26` = "mneVAS")
srt.romanov.pub %<>% RenameIdents(object = ., `31` = "pneSS")
srt.romanov.pub %<>% RenameIdents(object = ., `24` = "pneCRH")
srt.romanov.pub %<>% RenameIdents(object = ., `15` = "pneTRH")
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.romanov.pub$stage <-
  srt.romanov.pub$age %>%
  forcats::fct_collapse(
    Embryonic = c("E15", "E17"),
    Neonatal = c("P0", "P2", "3P2"),
    Pubertal = c("1P10", "P10"),
    Adult = c("P23")
  )
srt.romanov.pub$stage %<>% factor(levels = c("Embryonic", "Neonatal", "Pubertal", "Adult"), ordered = TRUE)
srt.romanov.pub$stage %>% forcats::fct_count()
```

::: {.cell-output-display}

`````{=html}
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["f"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n"],"name":[2],"type":["int"],"align":["right"]}],"data":[{"1":"Embryonic","2":"19503"},{"1":"Neonatal","2":"20316"},{"1":"Pubertal","2":"8965"},{"1":"Adult","2":"2415"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
`````

:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.romanov.pub$age <-
  plyr::mapvalues(
    x = srt.romanov.pub$age,
    from = c("E15", "E17", "P0", "P2", "3P2", "1P10", "P10", "P23"),
    to = c("E15", "E17", "P00", "P02", "P02", "P10", "P10", "P23")
  )



srt.romanov.pub$age %>% forcats::fct_count()
```

::: {.cell-output-display}

`````{=html}
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["f"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n"],"name":[2],"type":["int"],"align":["right"]}],"data":[{"1":"E15","2":"8290"},{"1":"E17","2":"11213"},{"1":"P00","2":"7492"},{"1":"P02","2":"12824"},{"1":"P10","2":"8965"},{"1":"P23","2":"2415"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
`````

:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot(
  srt.romanov.pub,
  features = c(neurotrans, metabolic_signaling_genes, "Crh", "Trh", "Oxt"),
  label = F,
  blend = F,
  order = TRUE,
  pt.size = 1.2,
  raster.dpi = c(1024, 1024),
  alpha = 0.5,
  split.by = "age"
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-feature-metabolic-romanov2020-1.png){#fig-feature-metabolic-romanov2020 fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  srt.romanov.pub@assays$RNA@data %>%
  as.data.frame() %>%
  t()
rownames(sbs_mtx) <- colnames(srt.romanov.pub)

# Filter features
filt_low_genes <-
  colSums(sbs_mtx) %>%
  .[. > quantile(., 0.4)] %>%
  names()
sbs_mtx %<>% .[, filt_low_genes]

min_filt_vector2 <-
  sbs_mtx %>%
  as_tibble() %>%
  dplyr::select(all_of(filt_low_genes)) %>%
  summarise(across(.fns = ~ quantile(.x, .005))) %>%
  as.list() %>%
  map(as.double) %>%
  simplify() %>%
  .[filt_low_genes]

# Prepare table of intersection sets analysis
content_sbs_mtx <-
  (sbs_mtx > min_filt_vector2) %>%
  as_tibble() %>%
  mutate_all(as.numeric)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(content_sbs_mtx),
  order.by = "freq",
  group.by = "sets",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 30,
  sets = c(metabolic_signaling_genes, "Crh", "Trh", "Oxt") %>%
    .[. %in% colnames(content_sbs_mtx)],
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-all-romanov2020-1.png){#fig-upset-group-metabolic-all-romanov2020 fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(content_sbs_mtx),
  order.by = "freq",
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 15,
  sets = c(metabolic_signaling_genes, "Crh", "Trh", "Oxt") %>%
    .[. %in% colnames(content_sbs_mtx)],
  nintersects = 20,
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-not-grouped-metabolic-all-romanov2020-1.png){#fig-upset-not-grouped-metabolic-all-romanov2020 fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx_full <- content_sbs_mtx |>
  dplyr::select(any_of(c(neurotrans, metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
  dplyr::bind_cols(srt.romanov.pub@meta.data)

sbs_mtx_full |> glimpse()
```

::: {.cell-output .cell-output-stdout}

```
Rows: 51,199
Columns: 35
$ Slc17a6          <dbl> 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Slc17a8          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Slc1a1           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Slc1a2           <dbl> 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Slc1a6           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Gad1             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,…
$ Slc32a1          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Slc6a1           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Alk              <dbl> 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,…
$ Mc4r             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Lmo4             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,…
$ Crh              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Trh              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Oxt              <dbl> 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ nGene            <int> 1652, 782, 447, 1706, 1106, 894, 727, 734, 669, 617, …
$ nUMI             <dbl> 2787, 1090, 544, 2709, 1817, 1220, 995, 1036, 920, 86…
$ orig.ident       <fct> Hypothalamus, Hypothalamus, Hypothalamus, Hypothalamu…
$ res.0.2          <chr> "23", "23", "23", "23", "23", "23", "23", "23", "23",…
$ res.0.4          <chr> "34", "34", "34", "34", "34", "34", "34", "34", "34",…
$ res.0.8          <chr> "42", "42", "42", "42", "42", "42", "42", "42", "42",…
$ res.1.2          <chr> "47", "47", "47", "47", "47", "47", "47", "47", "47",…
$ res.2            <chr> "54", "54", "54", "54", "54", "54", "54", "54", "54",…
$ tree.ident       <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
$ pro_Inter        <chr> "41", "41", "41", "41", "41", "41", "41", "41", "41",…
$ pro_Enter        <chr> "41", "41", "41", "41", "41", "41", "41", "41", "41",…
$ tree_final       <fct> 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 1…
$ subtree          <fct> 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 4…
$ prim_walktrap    <fct> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 3…
$ umi_per_gene     <dbl> 1.687046, 1.393862, 1.217002, 1.587925, 1.642857, 1.3…
$ log_umi_per_gene <dbl> 0.22712693, 0.14421974, 0.08529138, 0.20082998, 0.215…
$ nCount_RNA       <dbl> 2787, 1090, 544, 2709, 1817, 1220, 995, 1036, 920, 86…
$ nFeature_RNA     <int> 1652, 782, 447, 1706, 1106, 894, 727, 734, 669, 617, …
$ wtree            <fct> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 3…
$ age              <chr> "P23", "P02", "P02", "P02", "P02", "P02", "P02", "P02…
$ stage            <ord> Adult, Neonatal, Neonatal, Neonatal, Neonatal, Neonat…
```


:::
:::




## Prepare query mapping between datasets




::: {.cell layout-align="center"}

```{.r .cell-code}
srt.kim <- NormalizeData(srt.kim)
srt.kim <- FindVariableFeatures(srt.kim, selection.method = "vst", nfeatures = 3000)
# all.genes <- rownames(srt.kim)
# srt.kim <- ScaleData(srt.kim, features = all.genes)
srt.kim <- ScaleData(srt.kim)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
hypoth.anchors <- FindTransferAnchors(
  reference = srt.romanov.pub, query = srt.kim, dims = 1:30,
  reference.reduction = "pca"
)
predictions <- TransferData(anchorset = hypoth.anchors, refdata = srt.romanov.pub$wtree, dims = 1:30)
srt.kim <- AddMetaData(srt.kim, metadata = predictions)
table(srt.kim$predicted.id)
```

::: {.cell-output .cell-output-stdout}

```

     1     10     11     12     13     15     16     17     19      2     20 
  3700      1   2198    225    545     58    218     19 100308   9886    924 
    22     23     26     27     28     29      3     38      4      5      6 
   163     11    115     24     76    481    951     49    156   2895   4756 
     7      8      9 
   199     47      1 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.romanov.pub <- RunUMAP(srt.romanov.pub, dims = 1:30, reduction = "pca", return.model = TRUE)
srt.kim <- IntegrateEmbeddings(
  anchorset = hypoth.anchors, reference = srt.romanov.pub, query = srt.kim,
  new.reduction.name = "ref.pca"
)
srt.kim <- ProjectUMAP(
  query = srt.kim, query.reduction = "ref.pca", reference = srt.romanov.pub,
  reference.reduction = "pca", reduction.model = "umap"
)
Idents(srt.kim) <- srt.kim$Cluster
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
all.genes <- rownames(srt.kim)
gene.scale <- c(
  cnbn,
  opioid_system_genes,
  metabolic_signaling_genes,
  np,
  npr,
  nmr,
  neurotrans
) |>
  unique() %>%
  .[. %in% all.genes]
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
p1 <- DimPlot(srt.romanov.pub,
  reduction = "umap", group.by = "wtree", label = F
) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(srt.kim,
  reduction = "ref.umap", group.by = "Age", label = F
) + NoLegend() + ggtitle("Query transferred Embedding (more ages)")
p1 + p2
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-reference-umap-transfered-1.png){#fig-reference-umap-transfered fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
p1 <- FeaturePlot_scCustom(
  srt.romanov.pub,
  reduction = "umap",
  features = c(
    "Oxt",
    "Avp",
    "Sst",
    "Crh",
    "Trh"
  ),
  label = F,
  num_columns = 5,
  min.cutoff = "q05",
  na_cutoff = 2
) * NoLegend()
p2 <- FeaturePlot_scCustom(
  srt.kim,
  reduction = "ref.umap",
  features = c(
    "Oxt",
    "Avp",
    "Sst",
    "Crh",
    "Trh"
  ),
  label = FALSE,
  num_columns = 5,
  min.cutoff = "q05",
  na_cutoff = 2
) * NoLegend()
(p1 / p2)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-reference-umap-transfered-genes-1.png){#fig-reference-umap-transfered-genes fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
p1 <- FeaturePlot_scCustom(
  srt.romanov.pub,
  reduction = "umap",
  features = c(
    "Oxt",
    "Avp"
  ),
  label = F,
  num_columns = 2,
  min.cutoff = "q05",
  na_cutoff = 5
) * NoLegend()
p2 <- FeaturePlot_scCustom(
  srt.kim,
  reduction = "ref.umap",
  features = c(
    "Oxt",
    "Avp"
  ),
  label = FALSE,
  num_columns = 2,
  min.cutoff = "q05",
  na_cutoff = 5
) * NoLegend()
(p1 / p2)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-reference-umap-transfered-genes-Avp-Oxt-1.png){#fig-reference-umap-transfered-genes-Avp-Oxt fig-align='center' width=2160}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.kim$stage <-
  srt.kim$Age %>%
  forcats::fct_collapse(
    Embryonic = c(
      "E10", "E11", "E12", "E13",
      "E14", "E15", "E16", "E18"
    ),
    Neonatal = c("P4", "P8"),
    Pubertal = c("P14"),
    Adult = c("P45")
  )
srt.kim$stage %<>% factor(levels = c("Embryonic", "Neonatal", "Pubertal", "Adult"), ordered = TRUE)
FeaturePlot_scCustom(
  srt.kim,
  reduction = "ref.umap",
  features = c(
    "Oxt",
    "Avp",
    "Sst",
    "Crh",
    "Trh"
  ),
  split.by = "stage",
  min.cutoff = "q05",
  na_cutoff = 2,
  label = F,
  num_columns = 4
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-ref-embedding-split-stage-kim2020-np-1.png){#fig-ref-embedding-split-stage-kim2020-np fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.kim,
  reduction = "ref.umap",
  features = c(
    "Fam150a",
    "Fam150b",
    "Alk",
    "Scgn",
    "Crh"
  ),
  split.by = "stage",
  max.cutoff = "q95",
  label = F,
  num_columns = 4
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-ref-embedding-split-stage-kim2020-crh-alk-1.png){#fig-ref-embedding-split-stage-kim2020-crh-alk fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.romanov.pub,
  reduction = "umap",
  features = c(
    "Oxt",
    "Avp",
    "Sst",
    "Crh",
    "Trh"
  ),
  split.by = "stage",
  min.cutoff = "q05",
  na_cutoff = 2,
  label = F,
  num_columns = 4
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-ref-embedding-split-stage-romanov2020-np-1.png){#fig-ref-embedding-split-stage-romanov2020-np fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.romanov.pub,
  reduction = "umap",
  features = c(
    "Fam150a",
    "Fam150b",
    "Alk",
    "Scgn",
    "Crh"
  ),
  split.by = "stage",
  max.cutoff = "q97.5",
  label = F,
  num_columns = 4
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-ref-embedding-split-stage-romanov2020-crh-alk-1.png){#fig-ref-embedding-split-stage-romanov2020-crh-alk fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
if (!file.exists(here(data_dir, "kim2020_pvn_neurons.txt"))) {
  plot <- DimPlot(object = srt.kim, reduction = "ref.umap")
  srt.kim <- CellSelector(plot = plot, object = srt.kim, ident = "SelectedCells")

  selected_cells <- Cells(subset(srt.kim, idents = "SelectedCells"))
  write_lines(selected_cells, file = here(data_dir, "kim2020_pvn_neurons.txt"))
}
selected_cells <- read_lines(here(data_dir, "kim2020_pvn_neurons.txt"))
srt.kim <- subset(srt.kim, cells = c(selected_cells, WhichCells(srt.kim, expression = (
    Crh > 0 & (Scgn > 0 | Alk > 0 | Fam150b > 0 | Fam150a > 0)
))))
# srt.kim <- subset(srt.kim, subset = refUMAP_1 > 4 & refUMAP_2 > -1)

srt.kim@meta.data <- srt.kim@meta.data |> dplyr::rename(wtree = predicted.id, age = Age)

srt.kim <- subset(srt.kim, subset = stage %in% c("Pubertal", "Adult"))

srt.kim
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
27998 features across 954 samples within 1 assay 
Active assay: RNA (27998 features, 3000 variable features)
 3 layers present: counts, data, scale.data
 3 dimensional reductions calculated: umap, ref.pca, ref.umap
```


:::
:::




We subset Kim et al., 2020 dataset to only Pubertal and Adult stages.

## Intersection sets analysis

### PVN Neurons from Kim et al. 2020, Nature Communications




::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.kim,
  reduction = "ref.umap",
  features = c(
    "Oxt",
    "Avp",
    "Sst",
    "Crh",
    "Trh"
  ),
  na_cutoff = 2,
  label = F,
  num_columns = 5
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-kim2020-pvn-feature-np-split-by-stages-1.png){#fig-kim2020-pvn-feature-np-split-by-stages fig-align='center' width=6300}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.kim,
  reduction = "ref.umap",
  features = c(
    metabolic_signaling_genes,
    opioid_system_genes,
    "Slc17a6", "Gad1", "Gad2", "Crh", "Trh", "Oxt", "Sst"
  ),
  label = F,
  num_columns = 4
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-kim2020-pvn-feature-metabopioid-split-by-stages-1.png){#fig-kim2020-pvn-feature-metabopioid-split-by-stages fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  srt.kim@assays$RNA@data %>%
  as.data.frame() %>%
  t()

rownames(sbs_mtx) <- colnames(srt.kim)
colnames(sbs_mtx) <- rownames(srt.kim)

# Filter features
filt_low_genes <-
  colSums(sbs_mtx) %>%
  .[. > quantile(., 0.4)] %>%
  names()
sbs_mtx %<>% .[, filt_low_genes]

min_filt_vector2 <-
  sbs_mtx %>%
  as_tibble() %>%
  dplyr::select(all_of(filt_low_genes)) %>%
  summarise(across(.fns = ~ quantile(.x, .005))) %>%
  as.list() %>%
  map(as.double) %>%
  simplify() %>%
  .[filt_low_genes]

# Prepare table of intersection sets analysis
content_sbs_mtx_kim <-
  (sbs_mtx > min_filt_vector2) %>%
  as_tibble() %>%
  mutate_all(as.numeric) %>%
  bind_cols(
    srt.kim@meta.data |> dplyr::select(wtree, age, stage)
  )
```
:::




#### All




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_kim |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 30,
  sets = colnames(data),
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-kim2020-pvn-1.png){#fig-upset-group-metabolic-kim2020-pvn fig-align='center' width=4200}
:::
:::




#### Pubertal




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_kim |>
      dplyr::filter(stage == "Pubertal") |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 30,
  sets = colnames(data),
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-kim2020-pvn-Pubertal-1.png){#fig-upset-group-metabolic-kim2020-pvn-Pubertal fig-align='center' width=4200}
:::
:::




#### Adult




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_kim |>
      dplyr::filter(stage == "Adult") |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 30,
  sets = colnames(data),
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-kim2020-pvn-Adult-1.png){#fig-upset-group-metabolic-kim2020-pvn-Adult fig-align='center' width=4200}
:::
:::




### PVN Neurons from Romanov et al. 2020, Nature




::: {.cell layout-align="center"}

```{.r .cell-code}
srt.romanov.pvn <-
    subset(
        x = srt.romanov.pub,
        cells = unique(c(
            WhichCells(srt.romanov.pub,
                       idents = c(
                           "mneOXY", "mneVAS",
                           "pneSS", "pneCRH", "pneTRH"
                       )), 
            WhichCells(
                srt.romanov.pub,
                expression = (Crh > 0 & (Scgn > 0 | Alk > 0 | Fam150b > 0 | Fam150a > 0)))
        )),
        invert = FALSE
    )

table(srt.romanov.pvn$age)
```

::: {.cell-output .cell-output-stdout}

```

E15 E17 P00 P02 P10 P23 
332 615 374 389 537  67 
```


:::

```{.r .cell-code}
table(srt.romanov.pvn$stage)
```

::: {.cell-output .cell-output-stdout}

```

Embryonic  Neonatal  Pubertal     Adult 
      947       763       537        67 
```


:::

```{.r .cell-code}
srt.romanov.pvn <- subset(srt.romanov.pvn, subset = stage %in% c("Pubertal", "Adult"))
```
:::




We subset Romanov et al., 2020 dataset to only Pubertal and Adult
stages.




::: {.cell layout-align="center"}

```{.r .cell-code}
table(srt.romanov.pvn$age)
```

::: {.cell-output .cell-output-stdout}

```

P10 P23 
537  67 
```


:::

```{.r .cell-code}
table(srt.romanov.pvn$stage)
```

::: {.cell-output .cell-output-stdout}

```

Embryonic  Neonatal  Pubertal     Adult 
        0         0       537        67 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.romanov.pvn,
  reduction = "umap",
  features = c(
    "Oxt",
    "Avp",
    "Sst",
    "Crh",
    "Trh"
  ),
  na_cutoff = 2,
  label = F,
  num_columns = 5
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-romanov2020-pvn-feature-np-1.png){#fig-romanov2020-pvn-feature-np fig-align='center' width=6300}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.romanov.pvn,
  reduction = "umap",
  features = c(
    "Fam150a",
    "Fam150b",
    "Alk",
    "Scgn",
    "Crh"
  ),
  max.cutoff = "q97.5",
  label = F,
  num_columns = 5
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-romanov2020-pvn-feature-crh-alk-1.png){#fig-romanov2020-pvn-feature-crh-alk fig-align='center' width=6300}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt.romanov.pvn,
  reduction = "umap",
  features = c(
    metabolic_signaling_genes,
    opioid_system_genes,
    "Slc17a6", "Gad1", "Gad2", "Crh", "Trh", "Oxt", "Sst"
  ),
  label = F,
  num_columns = 4
) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-romanov2020-pvn-feature-metabopioid-1.png){#fig-romanov2020-pvn-feature-metabopioid fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  srt.romanov.pvn@assays$RNA@data %>%
  as.data.frame() %>%
  t()

# Filter features
filt_low_genes <-
  colSums(sbs_mtx) %>%
  .[. > quantile(., 0.4)] %>%
  names()
sbs_mtx %<>% .[, filt_low_genes]

min_filt_vector2 <-
  sbs_mtx %>%
  as_tibble() %>%
  dplyr::select(all_of(filt_low_genes)) %>%
  summarise(across(.fns = ~ quantile(.x, .005))) %>%
  as.list() %>%
  map(as.double) %>%
  simplify() %>%
  .[filt_low_genes]

# Prepare table of intersection sets analysis
content_sbs_mtx_romanov <-
  (sbs_mtx > min_filt_vector2) %>%
  as_tibble() %>%
  mutate_all(as.numeric) %>%
  bind_cols(
    srt.romanov.pvn@meta.data |> dplyr::select(wtree, age, stage)
  )
```
:::




#### All




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_romanov |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-romanov2020-pvn-1.png){#fig-upset-group-metabolic-romanov2020-pvn fig-align='center' width=10800}
:::
:::




#### Pubertal




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_romanov |>
      dplyr::filter(stage == "Pubertal") |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-romanov2020-pvn-Pubertal-1.png){#fig-upset-group-metabolic-romanov2020-pvn-Pubertal fig-align='center' width=10800}
:::
:::




#### Adult




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_romanov |>
      dplyr::filter(stage == "Adult") |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-romanov2020-pvn-Adult-1.png){#fig-upset-group-metabolic-romanov2020-pvn-Adult fig-align='center' width=7200}
:::
:::




### PVN Neurons from both datasets joined




::: {.cell layout-align="center"}

```{.r .cell-code}
# Prepare table of intersection sets analysis
to_select <-
  c(gene.scale, "wtree", "age", "stage") %>%
  .[. %in% colnames(content_sbs_mtx_kim)] %>%
  .[. %in% colnames(content_sbs_mtx_romanov)]

content_sbs_mtx <-
  bind_rows(
    content_sbs_mtx_kim |> dplyr::select(all_of(to_select)),
    content_sbs_mtx_romanov |> dplyr::select(all_of(to_select))
  )
```
:::




#### All




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-pvn-1.png){#fig-upset-group-metabolic-pvn fig-align='center' width=10800}
:::
:::




#### Pubertal




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      dplyr::filter(stage == "Pubertal") |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)

skim(as.data.frame(
  content_sbs_mtx |>
    dplyr::filter(stage == "Pubertal") |>
    dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt")))
))
```

::: {#fig-upset-group-metabolic-pvn-Pubertal-f2-1 .cell-output-display}

Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |753                |
|Number of columns        |7                  |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |7                  |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.37| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Mc3r          |         0|             1| 0.01| 0.12|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc4r          |         0|             1| 0.03| 0.18|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lmo4          |         0|             1| 0.28| 0.45|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Crh           |         0|             1| 0.16| 0.37|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Trh           |         0|             1| 0.33| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oxt           |         0|             1| 0.39| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |


:::

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-pvn-Pubertal-f2-1.png){#fig-upset-group-metabolic-pvn-Pubertal-f2-2 fig-align='center' width=10800}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      dplyr::filter(stage == "Pubertal") |>
      dplyr::select(any_of(c(
        opioid_system_genes, "Crh", "Trh", "Oxt"
      ))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)

skim(as.data.frame(
  content_sbs_mtx |>
    dplyr::filter(stage == "Pubertal") |>
    dplyr::select(any_of(c(
      opioid_system_genes, "Crh", "Trh", "Oxt"
    ))) |>
    dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
))
```

::: {#fig-upset-group-e-opioid-pvn-Pubertal-f3-1 .cell-output-display}

Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |753                |
|Number of columns        |12                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |12                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Oprd1         |         0|             1| 0.02| 0.13|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.15| 0.35|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Oprl1         |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Oprm1         |         0|             1| 0.19| 0.40|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pcsk1         |         0|             1| 0.17| 0.38|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pcsk2         |         0|             1| 0.34| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Pdyn          |         0|             1| 0.19| 0.39|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Penk          |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pnoc          |         0|             1| 0.16| 0.37|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Crh           |         0|             1| 0.16| 0.37|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Trh           |         0|             1| 0.33| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oxt           |         0|             1| 0.39| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |


:::

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-e-opioid-pvn-Pubertal-f3-1.png){#fig-upset-group-e-opioid-pvn-Pubertal-f3-2 fig-align='center' width=10800}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      dplyr::filter(stage == "Pubertal") |>
      dplyr::select(any_of(c(
        metabolic_signaling_genes,
        opioid_system_genes, "Crh", "Trh", "Oxt"
      ))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)

skim(as.data.frame(
  content_sbs_mtx |>
    dplyr::filter(stage == "Pubertal") |>
    dplyr::select(any_of(c(
      metabolic_signaling_genes,
      opioid_system_genes, "Crh", "Trh", "Oxt"
    ))) |>
    dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
))
```

::: {#fig-upset-group-e-metabopioid-pvn-Pubertal-f4-1 .cell-output-display}

Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |753                |
|Number of columns        |16                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |16                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.37| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Mc3r          |         0|             1| 0.01| 0.12|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc4r          |         0|             1| 0.03| 0.18|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lmo4          |         0|             1| 0.28| 0.45|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oprd1         |         0|             1| 0.02| 0.13|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.15| 0.35|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Oprl1         |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Oprm1         |         0|             1| 0.19| 0.40|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pcsk1         |         0|             1| 0.17| 0.38|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pcsk2         |         0|             1| 0.34| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Pdyn          |         0|             1| 0.19| 0.39|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Penk          |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pnoc          |         0|             1| 0.16| 0.37|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Crh           |         0|             1| 0.16| 0.37|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Trh           |         0|             1| 0.33| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oxt           |         0|             1| 0.39| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |


:::

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-e-metabopioid-pvn-Pubertal-f4-1.png){#fig-upset-group-e-metabopioid-pvn-Pubertal-f4-2 fig-align='center' width=10800}
:::
:::




#### Adult




::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      dplyr::filter(stage == "Adult") |>
      dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)

skim(as.data.frame(
  content_sbs_mtx |>
    dplyr::filter(stage == "Adult") |>
    dplyr::select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
    dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
))
```

::: {#fig-upset-group-metabolic-pvn-Adult-f2-1 .cell-output-display}

Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |805                |
|Number of columns        |7                  |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |7                  |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.04| 0.20|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc3r          |         0|             1| 0.01| 0.09|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc4r          |         0|             1| 0.01| 0.08|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lmo4          |         0|             1| 0.09| 0.28|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Crh           |         0|             1| 0.01| 0.11|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oxt           |         0|             1| 0.40| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |


:::

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-metabolic-pvn-Adult-f2-1.png){#fig-upset-group-metabolic-pvn-Adult-f2-2 fig-align='center' width=10800}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      dplyr::filter(stage == "Adult") |>
      dplyr::select(any_of(c(opioid_system_genes, "Crh", "Trh", "Oxt"))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 70,
  sets = colnames(data),
  empty.intersections = NULL
)

skim(as.data.frame(
  content_sbs_mtx |>
    dplyr::filter(stage == "Adult") |>
    dplyr::select(any_of(c(opioid_system_genes, "Crh", "Trh", "Oxt"))) |>
    dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
))
```

::: {#fig-upset-group-e-opioid-pvn-Adult-f3-1 .cell-output-display}

Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |805                |
|Number of columns        |12                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |12                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Oprd1         |         0|             1| 0.00| 0.05|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.04| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprl1         |         0|             1| 0.11| 0.31|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprm1         |         0|             1| 0.02| 0.14|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pcsk1         |         0|             1| 0.09| 0.29|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pcsk2         |         0|             1| 0.29| 0.46|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Pdyn          |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Penk          |         0|             1| 0.04| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pnoc          |         0|             1| 0.02| 0.14|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Crh           |         0|             1| 0.01| 0.11|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oxt           |         0|             1| 0.40| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |


:::

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-e-opioid-pvn-Adult-f3-1.png){#fig-upset-group-e-opioid-pvn-Adult-f3-2 fig-align='center' width=10800}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      dplyr::filter(stage == "Adult") |>
      dplyr::select(any_of(c(
        metabolic_signaling_genes,
        opioid_system_genes, "Crh", "Trh", "Oxt"
      ))) |>
      dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 30,
  sets = colnames(data),
  empty.intersections = NULL
)

skim(as.data.frame(
  content_sbs_mtx |>
    dplyr::filter(stage == "Adult") |>
    dplyr::select(any_of(c(
      metabolic_signaling_genes,
      opioid_system_genes, "Crh", "Trh", "Oxt"
    ))) |>
    dplyr::select(where(~ is.numeric(.x) && sum(.x) > 0))
))
```

::: {#fig-upset-group-e-metabopioid-pvn-Adult-f4-1 .cell-output-display}

Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |805                |
|Number of columns        |16                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |16                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.04| 0.20|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc3r          |         0|             1| 0.01| 0.09|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc4r          |         0|             1| 0.01| 0.08|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lmo4          |         0|             1| 0.09| 0.28|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprd1         |         0|             1| 0.00| 0.05|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.04| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprl1         |         0|             1| 0.11| 0.31|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprm1         |         0|             1| 0.02| 0.14|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pcsk1         |         0|             1| 0.09| 0.29|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pcsk2         |         0|             1| 0.29| 0.46|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Pdyn          |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Penk          |         0|             1| 0.04| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pnoc          |         0|             1| 0.02| 0.14|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Crh           |         0|             1| 0.01| 0.11|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oxt           |         0|             1| 0.40| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |


:::

::: {.cell-output-display}
![](02-endo-metabopioids_files/figure-html/fig-upset-group-e-metabopioid-pvn-Adult-f4-1.png){#fig-upset-group-e-metabopioid-pvn-Adult-f4-2 fig-align='center' width=10800}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sessioninfo::session_info()
```

::: {.cell-output .cell-output-stdout}

```
─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.2 (2024-10-31)
 os       Ubuntu 22.04.5 LTS
 system   x86_64, linux-gnu
 ui       X11
 language en_US:en
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Etc/UTC
 date     2025-01-09
 pandoc   3.1.11.1 @ /home/etretiakov/micromamba/bin/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────
 package              * version     date (UTC) lib source
 abind                  1.4-8       2024-09-12 [2] RSPM
 AnnotationDbi        * 1.68.0      2024-10-29 [2] RSPM (R 4.4.2)
 AnnotationFilter     * 1.30.0      2024-10-29 [2] RSPM (R 4.4.2)
 ape                    5.8         2024-04-11 [2] RSPM (R 4.4.0)
 base64enc              0.1-3       2015-07-28 [2] RSPM (R 4.4.0)
 bayestestR             0.15.0      2024-10-17 [2] RSPM (R 4.4.0)
 beachmat               2.22.0      2024-10-29 [2] RSPM (R 4.4.2)
 beeswarm               0.4.0       2021-06-01 [2] RSPM (R 4.4.0)
 Biobase              * 2.66.0      2024-10-29 [2] RSPM (R 4.4.2)
 BiocGenerics         * 0.52.0      2024-10-29 [2] RSPM (R 4.4.2)
 BiocIO                 1.16.0      2024-10-29 [2] RSPM (R 4.4.2)
 BiocManager            1.30.25     2024-08-28 [2] RSPM (R 4.4.0)
 BiocNeighbors          2.0.1       2024-11-28 [2] RSPM (R 4.4.2)
 BiocParallel           1.40.0      2024-10-29 [2] RSPM (R 4.4.2)
 BiocSingular           1.22.0      2024-10-29 [2] RSPM (R 4.4.2)
 Biostrings             2.74.0      2024-10-29 [2] RSPM (R 4.4.2)
 bit                    4.5.0       2024-09-20 [2] RSPM
 bit64                  4.5.2       2024-09-22 [2] RSPM
 bitops                 1.0-9       2024-10-03 [2] RSPM
 blob                   1.2.4       2023-03-17 [2] RSPM
 bluster                1.16.0      2024-10-29 [2] RSPM (R 4.4.2)
 cachem                 1.1.0       2024-05-16 [2] RSPM (R 4.4.0)
 callr                  3.7.6       2024-03-25 [2] RSPM (R 4.4.0)
 circlize               0.4.16      2024-12-04 [2] Github (jokergoo/circlize@9b21578)
 cli                    3.6.3       2024-06-21 [2] RSPM (R 4.4.0)
 cluster                2.1.6       2023-12-01 [2] CRAN (R 4.4.2)
 coda                   0.19-4.1    2024-01-31 [2] RSPM
 codetools              0.2-20      2024-03-31 [2] CRAN (R 4.4.2)
 colorspace             2.1-1       2024-07-26 [2] RSPM (R 4.4.0)
 correlation            0.8.6       2024-10-26 [2] RSPM (R 4.4.0)
 cowplot              * 1.1.3       2024-01-22 [2] RSPM
 crayon                 1.5.3       2024-06-20 [2] RSPM (R 4.4.0)
 curl                   6.0.1       2024-11-14 [2] RSPM (R 4.4.0)
 data.table             1.16.2      2024-10-10 [2] RSPM
 datawizard             0.13.0.17   2024-12-04 [2] Github (easystats/datawizard@25f8ec4)
 DBI                    1.2.3       2024-06-02 [2] RSPM
 DelayedArray           0.32.0      2024-10-29 [2] RSPM (R 4.4.2)
 deldir                 2.0-4       2024-02-28 [2] RSPM (R 4.4.0)
 digest                 0.6.37      2024-08-19 [2] RSPM (R 4.4.0)
 dotCall64              1.2         2024-10-04 [2] RSPM
 dplyr                * 1.1.4       2023-11-17 [2] RSPM (R 4.4.0)
 dqrng                  0.4.1       2024-05-28 [2] RSPM (R 4.4.0)
 edgeR                  4.4.1       2024-12-02 [2] RSPM (R 4.4.2)
 effectsize             0.8.9       2024-07-03 [2] RSPM (R 4.4.0)
 emmeans                1.10.5      2024-10-14 [2] RSPM
 EnsDb.Mmusculus.v79  * 2.99.0      2024-12-04 [2] RSPM (R 4.4.2)
 ensembldb            * 2.30.0      2024-10-29 [2] RSPM (R 4.4.2)
 estimability           1.5.1       2024-05-12 [2] RSPM (R 4.4.0)
 evaluate               1.0.1       2024-10-10 [2] RSPM (R 4.4.0)
 fansi                  1.0.6       2023-12-08 [2] RSPM (R 4.4.0)
 farver                 2.1.2       2024-05-13 [2] RSPM (R 4.4.0)
 fastDummies            1.7.4       2024-08-16 [2] RSPM
 fastmap                1.2.0       2024-05-15 [2] RSPM (R 4.4.0)
 fitdistrplus           1.2-1       2024-07-12 [2] RSPM (R 4.4.0)
 FNN                    1.1.4.1     2024-09-22 [2] RSPM (R 4.4.0)
 forcats              * 1.0.0       2023-01-29 [2] RSPM
 fs                     1.6.5       2024-10-30 [2] RSPM (R 4.4.0)
 future               * 1.34.0      2024-07-29 [2] RSPM
 future.apply           1.11.3      2024-10-27 [2] RSPM
 generics               0.1.3       2022-07-05 [2] RSPM (R 4.4.0)
 GenomeInfoDb         * 1.42.1      2024-11-28 [2] RSPM (R 4.4.2)
 GenomeInfoDbData       1.2.13      2024-12-04 [2] RSPM (R 4.4.2)
 GenomicAlignments      1.42.0      2024-10-29 [2] RSPM (R 4.4.2)
 GenomicFeatures      * 1.58.0      2024-10-29 [2] RSPM (R 4.4.2)
 GenomicRanges        * 1.58.0      2024-10-29 [2] RSPM (R 4.4.2)
 getPass                0.2-4       2023-12-10 [2] RSPM
 ggbeeswarm             0.7.2       2024-12-04 [2] Github (eclarke/ggbeeswarm@14ef76c)
 ggmin                  0.0.0.9000  2024-12-04 [2] Github (sjessa/ggmin@8ada274)
 ggplot2              * 3.5.1       2024-04-23 [2] RSPM (R 4.4.0)
 ggprism                1.0.5       2024-12-04 [2] Github (csdaw/ggprism@b6e6c0e)
 ggrastr                1.0.2       2024-12-04 [2] Github (VPetukhov/ggrastr@50ca3e0)
 ggrepel                0.9.6.9999  2024-12-04 [2] Github (slowkow/ggrepel@e72a66d)
 ggridges               0.5.6       2024-01-23 [2] RSPM
 ggstatsplot          * 0.12.5.9000 2024-12-04 [2] Github (IndrajeetPatil/ggstatsplot@d312b9f)
 git2r                  0.35.0      2024-10-20 [2] RSPM
 glmGamPoi              1.18.0      2024-10-29 [2] RSPM (R 4.4.2)
 GlobalOptions          0.1.2       2020-06-10 [2] RSPM (R 4.4.0)
 globals                0.16.3      2024-03-08 [2] RSPM
 glue                   1.8.0       2024-09-30 [2] RSPM (R 4.4.0)
 goftest                1.2-3       2021-10-07 [2] RSPM
 gprofiler2           * 0.2.3       2024-02-23 [2] RSPM (R 4.4.0)
 gridExtra              2.3         2017-09-09 [2] RSPM
 gtable                 0.3.6       2024-10-25 [2] RSPM (R 4.4.0)
 here                 * 1.0.1       2020-12-13 [2] RSPM
 hms                    1.1.3       2023-03-21 [2] RSPM
 htmltools              0.5.8.1     2024-04-04 [2] RSPM (R 4.4.0)
 htmlwidgets            1.6.4       2023-12-06 [2] RSPM (R 4.4.0)
 httpuv                 1.6.15      2024-03-26 [2] RSPM (R 4.4.0)
 httr                   1.4.7       2023-08-15 [2] RSPM (R 4.4.0)
 ica                    1.0-3       2022-07-08 [2] RSPM
 igraph                 2.1.1       2024-10-19 [2] RSPM (R 4.4.0)
 insight                1.0.0.2     2024-12-04 [2] Github (easystats/insight@8e78b12)
 IRanges              * 2.40.1      2024-12-05 [2] RSPM (R 4.4.2)
 irlba                  2.3.5.1     2022-10-03 [2] RSPM
 janitor                2.2.0.9000  2024-12-04 [2] Github (sfirke/janitor@6ee7919)
 jsonlite               1.8.9       2024-09-20 [2] RSPM (R 4.4.0)
 KEGGREST               1.46.0      2024-10-29 [2] RSPM (R 4.4.2)
 KernSmooth             2.23-24     2024-05-17 [2] CRAN (R 4.4.2)
 knitr                  1.49        2024-11-08 [2] RSPM
 labeling               0.4.3       2023-08-29 [2] RSPM (R 4.4.0)
 later                  1.4.1       2024-11-27 [2] RSPM (R 4.4.0)
 lattice                0.22-6      2024-03-20 [2] CRAN (R 4.4.2)
 lazyeval               0.2.2       2019-03-15 [2] RSPM (R 4.4.0)
 leiden                 0.4.3.1     2023-11-17 [2] RSPM
 lifecycle              1.0.4       2023-11-07 [2] RSPM (R 4.4.0)
 limma                  3.62.1      2024-11-03 [2] RSPM (R 4.4.2)
 listenv                0.9.1       2024-01-29 [2] RSPM
 lmtest                 0.9-40      2022-03-21 [2] RSPM (R 4.4.0)
 locfit                 1.5-9.10    2024-06-24 [2] RSPM (R 4.4.0)
 lubridate            * 1.9.3       2023-09-27 [2] RSPM
 magrittr             * 2.0.3       2022-03-30 [2] RSPM (R 4.4.0)
 MASS                   7.3-61      2024-06-13 [2] CRAN (R 4.4.2)
 Matrix                 1.7-1       2024-10-18 [2] CRAN (R 4.4.2)
 MatrixGenerics       * 1.18.0      2024-10-29 [2] RSPM (R 4.4.2)
 matrixStats          * 1.4.1       2024-09-08 [2] RSPM (R 4.4.0)
 memoise                2.0.1       2021-11-26 [2] RSPM (R 4.4.0)
 metapod                1.14.0      2024-10-29 [2] RSPM (R 4.4.2)
 mime                   0.12        2021-09-28 [2] RSPM (R 4.4.0)
 miniUI                 0.1.1.1     2018-05-18 [2] RSPM
 multcomp               1.4-26      2024-07-18 [2] RSPM
 munsell                0.5.1       2024-04-01 [2] RSPM (R 4.4.0)
 mvtnorm                1.3-2       2024-11-04 [2] RSPM
 nlme                   3.1-166     2024-08-14 [2] CRAN (R 4.4.2)
 org.Mm.eg.db         * 3.20.0      2024-12-04 [2] RSPM (R 4.4.2)
 paletteer              1.6.0       2024-01-21 [2] RSPM
 parallelly             1.39.0      2024-11-07 [2] RSPM
 parameters             0.24.0.3    2024-12-04 [2] Github (easystats/parameters@eff54e5)
 patchwork            * 1.3.0.9000  2024-12-04 [2] Github (thomasp85/patchwork@2695a9f)
 pbapply                1.7-2       2023-06-27 [2] RSPM
 pheatmap               1.0.12      2019-01-04 [2] RSPM (R 4.4.0)
 pillar                 1.9.0       2023-03-22 [2] RSPM (R 4.4.0)
 pkgconfig              2.0.3       2019-09-22 [2] RSPM (R 4.4.0)
 plotly                 4.10.4      2024-01-13 [2] RSPM
 plyr                   1.8.9       2023-10-02 [2] RSPM
 png                    0.1-8       2022-11-29 [2] RSPM
 polyclip               1.10-7      2024-07-23 [2] RSPM (R 4.4.0)
 prismatic              1.1.2       2024-04-10 [2] RSPM
 processx               3.8.4       2024-03-16 [2] RSPM
 progressr              0.15.1      2024-11-22 [2] RSPM
 promises               1.3.2       2024-11-28 [2] RSPM (R 4.4.0)
 ProtGenerics           1.38.0      2024-10-29 [2] RSPM (R 4.4.2)
 ps                     1.8.1       2024-10-28 [2] RSPM (R 4.4.0)
 purrr                * 1.0.2       2023-08-10 [2] RSPM (R 4.4.0)
 R.methodsS3            1.8.2       2022-06-13 [2] RSPM (R 4.4.0)
 R.oo                   1.27.0      2024-11-01 [2] RSPM (R 4.4.0)
 R.utils                2.12.3      2023-11-18 [2] RSPM (R 4.4.0)
 R6                     2.5.1       2021-08-19 [2] RSPM (R 4.4.0)
 RANN                   2.6.2       2024-08-25 [2] RSPM (R 4.4.0)
 RColorBrewer         * 1.1-3       2022-04-03 [2] RSPM
 Rcpp                   1.0.13-1    2024-11-02 [2] RSPM (R 4.4.0)
 RcppAnnoy              0.0.22      2024-01-23 [2] RSPM
 RcppHNSW               0.6.0       2024-02-04 [2] RSPM
 RCurl                  1.98-1.16   2024-07-11 [2] RSPM
 readr                * 2.1.5       2024-01-10 [2] RSPM
 rematch2               2.1.2       2020-05-01 [2] RSPM
 remotes                2.5.0       2024-03-17 [2] RSPM
 repr                   1.1.7       2024-03-22 [2] RSPM
 reshape2               1.4.4       2020-04-09 [2] RSPM
 restfulr               0.0.15      2022-06-16 [2] RSPM (R 4.4.2)
 reticulate             1.40.0.9000 2024-12-04 [2] Github (rstudio/reticulate@61f0fa4)
 rhdf5                  2.50.0      2024-10-29 [2] RSPM (R 4.4.2)
 rhdf5filters           1.18.0      2024-10-29 [2] RSPM (R 4.4.2)
 Rhdf5lib               1.28.0      2024-10-29 [2] RSPM (R 4.4.2)
 rjson                  0.2.23      2024-09-16 [2] RSPM (R 4.4.0)
 rlang                  1.1.4       2024-06-04 [2] RSPM (R 4.4.0)
 rmarkdown              2.29        2024-11-04 [2] RSPM
 ROCR                   1.0-11      2020-05-02 [2] RSPM (R 4.4.0)
 rprojroot              2.0.4       2023-11-05 [2] RSPM (R 4.4.0)
 Rsamtools              2.22.0      2024-10-29 [2] RSPM (R 4.4.2)
 RSpectra               0.16-2      2024-07-18 [2] RSPM
 RSQLite                2.3.8       2024-11-17 [2] RSPM (R 4.4.0)
 rstudioapi             0.17.1      2024-10-22 [2] RSPM
 rsvd                   1.0.5       2021-04-16 [2] RSPM (R 4.4.0)
 rtracklayer            1.66.0      2024-10-29 [2] RSPM (R 4.4.2)
 Rtsne                  0.17        2023-12-07 [2] RSPM (R 4.4.0)
 S4Arrays               1.6.0       2024-10-29 [2] RSPM (R 4.4.2)
 S4Vectors            * 0.44.0      2024-10-29 [2] RSPM (R 4.4.2)
 sandwich               3.1-1       2024-09-15 [2] RSPM
 ScaledMatrix           1.14.0      2024-10-29 [2] RSPM (R 4.4.2)
 scales               * 1.3.0       2023-11-28 [2] RSPM (R 4.4.0)
 scater               * 1.34.0      2024-10-29 [2] RSPM (R 4.4.2)
 scattermore            1.2         2023-06-12 [2] RSPM
 scCustomize          * 3.0.1       2025-01-09 [2] Github (samuel-marsh/scCustomize@3299b95)
 schard                 0.0.1       2024-12-04 [2] Github (cellgeni/schard@c22b46d)
 scran                * 1.34.0      2024-10-29 [2] RSPM (R 4.4.2)
 sctransform            0.4.1       2023-10-19 [2] RSPM
 scuttle              * 1.16.0      2024-10-29 [2] RSPM (R 4.4.2)
 sessioninfo            1.2.2       2021-12-06 [2] RSPM
 Seurat               * 5.1.0       2024-12-04 [2] Github (satijalab/seurat@1549dcb)
 SeuratObject         * 5.0.99.9001 2024-12-04 [2] Github (satijalab/seurat-object@42e53ba)
 SeuratWrappers       * 0.4.0       2024-12-04 [2] Github (satijalab/seurat-wrappers@a1eb0d8)
 shape                  1.4.6.1     2024-02-23 [2] RSPM
 shiny                  1.9.1       2024-08-01 [2] RSPM (R 4.4.0)
 SingleCellExperiment * 1.28.1      2024-11-10 [2] RSPM (R 4.4.2)
 skimr                * 2.1.5       2024-12-04 [2] Github (ropensci/skimr@d5126aa)
 snakecase              0.11.1      2023-08-27 [2] RSPM (R 4.4.0)
 sp                   * 2.1-4       2024-04-30 [2] RSPM
 spam                   2.11-0      2024-10-03 [2] RSPM
 SparseArray            1.6.0       2024-10-29 [2] RSPM (R 4.4.2)
 spatstat.data          3.1-4       2024-11-15 [2] RSPM (R 4.4.0)
 spatstat.explore       3.3-3       2024-10-22 [2] RSPM
 spatstat.geom          3.3-4       2024-11-18 [2] RSPM (R 4.4.0)
 spatstat.random        3.3-2       2024-09-18 [2] RSPM (R 4.4.0)
 spatstat.sparse        3.1-0       2024-06-21 [2] RSPM
 spatstat.univar        3.1-1       2024-11-05 [2] RSPM (R 4.4.0)
 spatstat.utils         3.1-1       2024-11-03 [2] RSPM (R 4.4.0)
 statmod                1.5.0       2023-01-06 [2] RSPM (R 4.4.0)
 statsExpressions       1.6.1       2024-10-31 [2] RSPM
 stringi                1.8.4       2024-05-06 [2] RSPM (R 4.4.0)
 stringr              * 1.5.1       2023-11-14 [2] RSPM (R 4.4.0)
 SummarizedExperiment * 1.36.0      2024-10-29 [2] RSPM (R 4.4.2)
 survival               3.7-0       2024-06-05 [2] CRAN (R 4.4.2)
 tensor                 1.5         2012-05-05 [2] RSPM
 TH.data                1.1-2       2023-04-17 [2] RSPM (R 4.4.0)
 tibble               * 3.2.1       2023-03-20 [2] RSPM (R 4.4.0)
 tidyr                * 1.3.1       2024-01-24 [2] RSPM (R 4.4.0)
 tidyselect             1.2.1       2024-03-11 [2] RSPM (R 4.4.0)
 tidyverse            * 2.0.0.9000  2024-12-04 [2] Github (tidyverse/tidyverse@c06a3c9)
 timechange             0.3.0       2024-01-18 [2] RSPM
 tzdb                   0.4.0       2023-05-12 [2] RSPM
 UCSC.utils             1.2.0       2024-10-29 [2] RSPM (R 4.4.2)
 UpSetR               * 1.4.0       2024-12-04 [2] Github (hms-dbmi/UpSetR@b14854a)
 utf8                   1.2.4       2023-10-22 [2] RSPM (R 4.4.0)
 uwot                   0.2.2       2024-04-21 [2] RSPM (R 4.4.0)
 vctrs                  0.6.5       2023-12-01 [2] RSPM (R 4.4.0)
 vipor                  0.4.7       2023-12-18 [2] RSPM (R 4.4.0)
 viridis              * 0.6.5       2024-01-29 [2] RSPM
 viridisLite          * 0.4.2       2023-05-02 [2] RSPM (R 4.4.0)
 vroom                  1.6.5       2023-12-05 [2] RSPM
 whisker                0.4.1       2022-12-05 [2] RSPM
 withr                  3.0.2       2024-10-28 [2] RSPM (R 4.4.0)
 workflowr            * 1.7.1       2023-08-23 [2] RSPM
 xfun                   0.49        2024-10-31 [2] RSPM (R 4.4.0)
 XML                    3.99-0.17   2024-06-25 [2] RSPM (R 4.4.0)
 xtable                 1.8-4       2019-04-21 [2] RSPM (R 4.4.0)
 XVector                0.46.0      2024-10-29 [2] RSPM (R 4.4.2)
 yaml                   2.3.10      2024-07-26 [2] RSPM
 zeallot                0.1.0       2018-01-28 [2] RSPM
 zlibbioc               1.52.0      2024-10-29 [2] RSPM (R 4.4.2)
 zoo                    1.8-12      2023-04-13 [2] RSPM (R 4.4.0)

 [1] /home/etretiakov/R/x86_64-pc-linux-gnu-library/4.4
 [2] /opt/R/4.4.2/lib/R/library

──────────────────────────────────────────────────────────────────────────────
```


:::
:::
