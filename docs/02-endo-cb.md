---
title: "eCB expression analysis of Hypothalamus development with focus on PVN"
author: "Evgenii O. Tretiakov"
date: "2024-11-27"
format:
  html:
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
  library(anndata)
  library(sceasy)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(scCustomize)
})

sc <- import("scanpy", convert = FALSE)
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
n_cores <- 8
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
- args: function (..., workers = 8, envir = parent.frame())
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
    "Oprd1",  # Delta opioid receptor
    "Oprk1",  # Kappa opioid receptor 
    "Oprl1",  # Nociceptin/orphanin FQ receptor
    "Oprm1",  # Mu opioid receptor
    
    # Processing enzymes
    "Pcsk1",  # Proprotein convertase 1
    "Pcsk2",   # Proprotein convertase 2
    
    # Endogenous opioid precursors
    "Pdyn",   # Prodynorphin
    "Penk",   # Proenkephalin
    #"Pomc",   # Proopiomelanocortin
    "Pnoc"   # Prepronociceptin
)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
metabolic_signaling_genes <- c(
    # Receptor tyrosine kinases and ligands
    "Alk",      # Anaplastic lymphoma kinase - neural development, metabolism
    "Fam150a",  # ALK ligand 1/Augmentor-α - ALK receptor activator
    "Fam150b",  # ALK ligand 2/Augmentor-β - ALK receptor activator
    
    # Melanocortin system
    "Mc3r",     # Melanocortin 3 receptor - energy homeostasis, inflammation
    "Mc4r",     # Melanocortin 4 receptor - appetite control, energy balance
    
    # Metabolic hormone receptors
    "Lepr",     # Leptin receptor - energy balance, satiety
    "Insr",     # Insulin receptor - glucose homeostasis
    #"Igf1r",    # Insulin-like growth factor 1 receptor - growth, development
    
    # Signaling adaptors/regulators
    "Lmo4",     # LIM domain only 4 - transcriptional regulation, metabolism
    "Irs1",     # Insulin receptor substrate 1 - insulin signaling
    "Irs4"      # Insulin receptor substrate 4 - insulin/leptin signaling
)
```
:::


## Load Kim DW et al 2020


::: {.cell layout-align="center"}

```{.r .cell-code}
anndata <- sc$read(here(
  "kim2020_combined.h5ad"
))
```
:::


### Convert adata object to R AnnDataR6 object.

::: {.cell layout-align="center"}

```{.r .cell-code}
adata <- py_to_r(anndata)
class(adata)
```

::: {.cell-output .cell-output-stdout}
```
[1] "AnnDataR6" "R6"       
```
:::

```{.r .cell-code}
class(adata$X)
```

::: {.cell-output .cell-output-stdout}
```
[1] "dgRMatrix"
attr(,"package")
[1] "Matrix"
```
:::

```{.r .cell-code}
adata
```

::: {.cell-output .cell-output-stdout}
```
AnnData object with n_obs × n_vars = 128006 × 27998
    obs: 'bc_name', 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'X1', 'X', 'Y', 'Z', 'Age', 'Cluster'
    var: '_index', 'features'
```
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt_path <- here(
  "data",
  sprintf("%s-whole_dataset-simple.h5Seurat", bioproject)
)

expr_mtx <- t(as.matrix(adata$X))
colnames(expr_mtx) <- rownames(adata$X)
rownames(expr_mtx) <- adata$var$features
srt <- CreateSeuratObject(
  expr_mtx,
  assay = "RNA",
  project = "kim2020_Hypoth_dev",
  meta.data = as.data.frame(adata$obs)
)

X_umap <- adata$obs |>
  select(X, Y) |>
  as.matrix()
colnames(X_umap) <- c("UMAP_1", "UMAP_2")
rownames(X_umap) <- colnames(expr_mtx)
srt[["umap"]] <- CreateDimReducObject(embeddings = X_umap, key = "umap_", assay = DefaultAssay(srt))

Idents(srt) <- "age"
srt <- Store_Palette_Seurat(seurat_object = srt, palette = rev(brewer.pal(n = 11, name = "Spectral")), palette_name = "expr_Colour_Pal")
```
:::


## Load Romanov et al 2020


::: {.cell layout-align="center"}

```{.r .cell-code}
print(srt)
```

::: {.cell-output .cell-output-stdout}
```
An object of class Seurat 
27998 features across 128006 samples within 1 assay 
Active assay: RNA (27998 features, 0 variable features)
 1 layer present: counts
 1 dimensional reduction calculated: umap
```
:::

```{.r .cell-code}
rar2020.srt.pub <- readRDS("/data/1_heteroAstrocytes/PRJNA548917/old/oldCCA_nae_srt.rds")
rar2020.srt.pub <- UpdateSeuratObject(rar2020.srt.pub)
Idents(rar2020.srt.pub) <-
  factor(rar2020.srt.pub$wtree,
    ordered = TRUE
  )

# Consistent colours and clusters names
colours_wtree <- setNames(read_lines(here(data_dir, "colours_wtree.tsv")), 1:45)

rar2020.srt.pub$age <-
  Cells(rar2020.srt.pub) |>
  str_split(pattern = ":", simplify = T) %>%
  .[, 1] %>%
  str_split_fixed(pattern = "_", n = 3) %>%
  .[, 3]
print(rar2020.srt.pub)
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
glimpse(rar2020.srt.pub@meta.data)
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
table(Idents(rar2020.srt.pub))
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
rar2020.srt.pub %<>% RenameIdents(object = ., `43` = "mneOXY")
rar2020.srt.pub %<>% RenameIdents(object = ., `26` = "mneVAS")
rar2020.srt.pub %<>% RenameIdents(object = ., `31` = "pneSS")
rar2020.srt.pub %<>% RenameIdents(object = ., `24` = "pneCRH")
rar2020.srt.pub %<>% RenameIdents(object = ., `15` = "pneTRH")
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
rar2020.srt.pub$stage <-
  rar2020.srt.pub$age %>%
  forcats::fct_collapse(
    Embryonic = c("E15", "E17"),
    Neonatal = c("P0", "P2", "3P2"),
    Pubertal = c("1P10", "P10"),
    Adult = c("P23")
  )
rar2020.srt.pub$stage %<>% factor(levels = c("Embryonic", "Neonatal", "Pubertal", "Adult"), ordered = TRUE)
rar2020.srt.pub$stage %>% forcats::fct_count()
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
rar2020.srt.pub$age <-
  plyr::mapvalues(
    x = rar2020.srt.pub$age,
    from = c("E15", "E17", "P0", "P2", "3P2", "1P10", "P10", "P23"),
    to = c("E15", "E17", "P00", "P02", "P02", "P10", "P10", "P23")
  )



rar2020.srt.pub$age %>% forcats::fct_count()
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
  rar2020.srt.pub,
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
![](02-endo-cb_files/figure-html/plot-feature-cb-romanov2020-1.png){fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  rar2020.srt.pub@assays$RNA@data %>%
  as.data.frame() %>%
  t()
rownames(sbs_mtx) <- colnames(rar2020.srt.pub)

# Filter features
filt_low_genes <-
  colSums(sbs_mtx) %>%
  .[. > quantile(., 0.4)] %>%
  names()
sbs_mtx %<>% .[, filt_low_genes]

min_filt_vector2 <-
  sbs_mtx %>%
  as_tibble() %>%
  select(all_of(filt_low_genes)) %>%
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

::: {.cell layout-align="center" fig.asp='1.214'}

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
![](02-endo-cb_files/figure-html/upset-group-e-cb-all-romanov2020-1.png){fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center" fig.asp='1.214'}

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
![](02-endo-cb_files/figure-html/upset-not-grouped-e-cb-all-romanov2020-1.png){fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx_full <- content_sbs_mtx |>
  select(any_of(c(neurotrans, metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |>
  dplyr::bind_cols(rar2020.srt.pub@meta.data)

sbs_mtx_full |> glimpse()
```

::: {.cell-output .cell-output-stdout}
```
Rows: 51,199
Columns: 38
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
$ Lepr             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Lmo4             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,…
$ Irs1             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ Irs4             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
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
srt <- NormalizeData(srt)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 3000)
# all.genes <- rownames(srt)
# srt <- ScaleData(srt, features = all.genes)
srt <- ScaleData(srt)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
hypoth.anchors <- FindTransferAnchors(
  reference = rar2020.srt.pub, query = srt, dims = 1:30,
  reference.reduction = "pca"
)
predictions <- TransferData(anchorset = hypoth.anchors, refdata = rar2020.srt.pub$wtree, dims = 1:30)
srt <- AddMetaData(srt, metadata = predictions)
table(srt$predicted.id)
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
rar2020.srt.pub <- RunUMAP(rar2020.srt.pub, dims = 1:30, reduction = "pca", return.model = TRUE)
srt <- IntegrateEmbeddings(
  anchorset = hypoth.anchors, reference = rar2020.srt.pub, query = srt,
  new.reduction.name = "ref.pca"
)
srt <- ProjectUMAP(
  query = srt, query.reduction = "ref.pca", reference = rar2020.srt.pub,
  reference.reduction = "pca", reduction.model = "umap"
)
Idents(srt) <- srt$Cluster
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
all.genes <- rownames(srt)
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
p1 <- DimPlot(rar2020.srt.pub,
  reduction = "umap", group.by = "wtree", label = F
) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(srt,
  reduction = "ref.umap", group.by = "Age", label = F
) + NoLegend() + ggtitle("Query transferred Embedding (more ages)")
p1 + p2
```

::: {.cell-output-display}
![](02-endo-cb_files/figure-html/plot-reference-umap-transfered-1.png){fig-align='center' width=4200}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
p1 <- FeaturePlot_scCustom(
  rar2020.srt.pub,
  reduction = "umap",
  features = c(
    "Oxt",
    "Avp",
    "Sst",
    "Crh",
    "Trh"),
  label = F,
  num_columns = 5,
  min.cutoff = 'q05',
  na_cutoff = 2
) * NoLegend()
p2 <- FeaturePlot_scCustom(
  srt,
  reduction = "ref.umap",
  features = c(
    "Oxt",
    "Avp",
    "Sst",
    "Crh",
    "Trh"),
  label = FALSE,
  num_columns = 5,
  min.cutoff = 'q05',
  na_cutoff = 2
) * NoLegend()
(p1 / p2)
```

::: {.cell-output-display}
![](02-endo-cb_files/figure-html/plot-reference-umap-transfered-genes-1.png){fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
p1 <- FeaturePlot_scCustom(
  rar2020.srt.pub,
  reduction = "umap",
  features = c(
    "Oxt",
    "Avp"),
  label = F,
  num_columns = 2,
  min.cutoff = 'q05',
  na_cutoff = 5
) * NoLegend()
p2 <- FeaturePlot_scCustom(
  srt,
  reduction = "ref.umap",
  features = c(
    "Oxt",
    "Avp"),
  label = FALSE,
  num_columns = 2,
  min.cutoff = 'q05',
  na_cutoff = 5
) * NoLegend()
(p1 / p2)
```

::: {.cell-output-display}
![](02-endo-cb_files/figure-html/plot-reference-umap-transfered-genes-Avp-Oxt-1.png){fig-align='center' width=2160}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt$stage <-
  srt$Age %>%
  forcats::fct_collapse(
    Embryonic = c(
      "E10", "E11", "E12", "E13",
      "E14", "E15", "E16", "E18"
    ),
    Neonatal = c("P4", "P8"),
    Pubertal = c("P14"),
    Adult = c("P45")
  )
srt$stage %<>% factor(levels = c("Embryonic", "Neonatal", "Pubertal", "Adult"), ordered = TRUE)
FeaturePlot_scCustom(
  srt,
  reduction = "ref.umap",
  features = c(
    "Oxt", 
    "Avp", 
    "Sst", 
    "Crh", 
    "Trh"), 
  split.by = "stage", 
  min.cutoff = 'q05',
  na_cutoff = 2, 
  label = F,
  num_columns = 4
  ) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-cb_files/figure-html/ref-embedding-split-stage-kim2020-np-1.png){fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  rar2020.srt.pub,
  reduction = "umap",
  features = c(
    "Oxt", 
    "Avp", 
    "Sst", 
    "Crh", 
    "Trh"), 
  split.by = "stage", 
  min.cutoff = 'q05',
  na_cutoff = 2, 
  label = F,
  num_columns = 4
  ) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-cb_files/figure-html/ref-embedding-split-stage-romanov2020-np-1.png){fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
if (!file.exists(here(data_dir, "kim2020_pvn_neurons.txt"))) {
  plot <- DimPlot(object = srt, reduction = "ref.umap")
  srt <- CellSelector(plot = plot, object = srt, ident = "SelectedCells")

  selected_cells <- Cells(subset(srt, idents = "SelectedCells"))
  write_lines(selected_cells, file = here(data_dir, "kim2020_pvn_neurons.txt"))
}
selected_cells <- read_lines(here(data_dir, "kim2020_pvn_neurons.txt"))
srt <- subset(srt, cells = selected_cells)
srt <- subset(srt, subset = refUMAP_1 > 4 & refUMAP_2 > -1)

srt@meta.data <- srt@meta.data |> rename(wtree = predicted.id, age = Age)

srt <- subset(srt, subset = stage %in% c("Pubertal", "Adult"))

srt
```

::: {.cell-output .cell-output-stdout}
```
An object of class Seurat 
27998 features across 950 samples within 1 assay 
Active assay: RNA (27998 features, 3000 variable features)
 3 layers present: counts, data, scale.data
 3 dimensional reductions calculated: umap, ref.pca, ref.umap
```
:::
:::


## Intersection sets analysis

### PVN Neurons from Kim et al. 2020, Nature Communications


::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt,
  reduction = "ref.umap",
  features = c(
    "Oxt", 
    "Avp", 
    "Sst", 
    "Crh", 
    "Trh"), 
  na_cutoff = 2, 
  label = F,
  num_columns = 5
  ) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-cb_files/figure-html/plot-kim2020-pvn-feature-np-split-by-stages-1.png){fig-align='center' width=6300}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt, 
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
![](02-endo-cb_files/figure-html/plot-kim2020-pvn-feature-metabopioid-split-by-stages-1.png){fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  srt@assays$RNA@layers$data %>%
  as.data.frame() %>%
  t()

rownames(sbs_mtx) <- colnames(srt)
colnames(sbs_mtx) <- rownames(srt)

# Filter features
filt_low_genes <-
  colSums(sbs_mtx) %>%
  .[. > quantile(., 0.4)] %>%
  names()
sbs_mtx %<>% .[, filt_low_genes]

min_filt_vector2 <-
  sbs_mtx %>%
  as_tibble() %>%
  select(all_of(filt_low_genes)) %>%
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
    srt@meta.data |> select(wtree, age, stage)
  )
```
:::


#### All


::: {.cell layout-align="center" fig.asp='1.214'}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_kim |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-kim2020-pvn-1.png){fig-align='center' width=4200}
:::
:::



#### Pubertal


::: {.cell layout-align="center" fig.asp='1.214'}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_kim |>
    filter(stage == "Pubertal") |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-kim2020-pvn-Pubertal-1.png){fig-align='center' width=4200}
:::
:::


#### Adult


::: {.cell layout-align="center" fig.asp='1.214'}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_kim |>
    filter(stage == "Adult") |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-kim2020-pvn-Adult-1.png){fig-align='center' width=4200}
:::
:::


### PVN Neurons from Romanov et al. 2020, Nature


::: {.cell layout-align="center"}

```{.r .cell-code}
rar2020.srt.pvn <-
  subset(
    x = rar2020.srt.pub,
    idents = c(
      "mneOXY", "mneVAS",
      "pneSS", "pneCRH", "pneTRH"
    ),
    invert = FALSE
  )

rar2020.srt.pvn <- subset(rar2020.srt.pvn, subset = umap_1 > 4 & umap_2 > -1)

rar2020.srt.pvn$age <-
  plyr::mapvalues(
    x = rar2020.srt.pvn$age,
    from = c("E15", "E17", "P0", "P2", "3P2", "1P10", "P10", "P23"),
    to = c("E15", "E17", "P00", "P02", "P02", "P10", "P10", "P23")
  )

rar2020.srt.pvn <- subset(rar2020.srt.pvn, subset = stage %in% c("Pubertal", "Adult"))
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  rar2020.srt.pvn,
  reduction = "umap",
  features = c(
    "Oxt", 
    "Avp", 
    "Sst", 
    "Crh", 
    "Trh"),
  na_cutoff = 2, 
  label = F,
  num_columns = 5
  ) * NoLegend()
```

::: {.cell-output-display}
![](02-endo-cb_files/figure-html/plot-romanov2020-pvn-feature-np-split-by-stages-1.png){fig-align='center' width=6300}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  rar2020.srt.pvn, 
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
![](02-endo-cb_files/figure-html/plot-romanov2020-pvn-feature-cb-split-by-stages-1.png){fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  rar2020.srt.pvn@assays$RNA@data %>%
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
  select(all_of(filt_low_genes)) %>%
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
    rar2020.srt.pvn@meta.data |> select(wtree, age, stage)
  )
```
:::


#### All


::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_romanov |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-romanov2020-pvn-1.png){fig-align='center' width=10800}
:::
:::


#### Pubertal


::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_romanov |>
    filter(stage == "Pubertal") |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-romanov2020-pvn-Pubertal-1.png){fig-align='center' width=10800}
:::
:::


#### Adult


::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx_romanov |>
    filter(stage == "Adult") |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-romanov2020-pvn-Adult-1.png){fig-align='center' width=7200}
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
    content_sbs_mtx_kim |> select(all_of(to_select)),
    content_sbs_mtx_romanov |> select(all_of(to_select))
  )
```
:::


#### All


::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-pvn-1.png){fig-align='center' width=10800}
:::
:::


#### Pubertal


::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Pubertal") |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-pvn-Pubertal-f2-1.png){fig-align='center' width=10800}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Pubertal") |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt")))
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |641                |
|Number of columns        |10                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |10                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.29| 0.45|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Mc3r          |         0|             1| 0.01| 0.11|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc4r          |         0|             1| 0.03| 0.17|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lepr          |         0|             1| 0.02| 0.16|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lmo4          |         0|             1| 0.27| 0.45|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Irs1          |         0|             1| 0.07| 0.26|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Irs4          |         0|             1| 0.25| 0.43|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Crh           |         0|             1| 0.05| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Oxt           |         0|             1| 0.39| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Pubertal") |>
    select(any_of(c(
      opioid_system_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-pvn-Pubertal-f3-1.png){fig-align='center' width=10800}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Pubertal") |>
    select(any_of(c(
      opioid_system_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |641                |
|Number of columns        |12                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |12                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Oprd1         |         0|             1| 0.01| 0.10|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.14| 0.34|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprl1         |         0|             1| 0.32| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oprm1         |         0|             1| 0.16| 0.36|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pcsk1         |         0|             1| 0.16| 0.37|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pcsk2         |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Pdyn          |         0|             1| 0.19| 0.39|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Penk          |         0|             1| 0.07| 0.26|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pnoc          |         0|             1| 0.13| 0.34|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Crh           |         0|             1| 0.05| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Oxt           |         0|             1| 0.39| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Pubertal") |>
    select(any_of(c(
      metabolic_signaling_genes,
      opioid_system_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-pvn-Pubertal-f4-1.png){fig-align='center' width=10800}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Pubertal") |>
    select(any_of(c(
      metabolic_signaling_genes,
      opioid_system_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |641                |
|Number of columns        |19                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |19                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.29| 0.45|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Mc3r          |         0|             1| 0.01| 0.11|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc4r          |         0|             1| 0.03| 0.17|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lepr          |         0|             1| 0.02| 0.16|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lmo4          |         0|             1| 0.27| 0.45|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Irs1          |         0|             1| 0.07| 0.26|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Irs4          |         0|             1| 0.25| 0.43|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Oprd1         |         0|             1| 0.01| 0.10|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.14| 0.34|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprl1         |         0|             1| 0.32| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oprm1         |         0|             1| 0.16| 0.36|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pcsk1         |         0|             1| 0.16| 0.37|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pcsk2         |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Pdyn          |         0|             1| 0.19| 0.39|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Penk          |         0|             1| 0.07| 0.26|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pnoc          |         0|             1| 0.13| 0.34|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Crh           |         0|             1| 0.05| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Oxt           |         0|             1| 0.39| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |
:::
:::


#### Adult


::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Adult") |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-pvn-Adult-f2-1.png){fig-align='center' width=10800}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Adult") |>
    select(any_of(c(metabolic_signaling_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |801                |
|Number of columns        |10                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |10                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.04| 0.20|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc3r          |         0|             1| 0.01| 0.09|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc4r          |         0|             1| 0.01| 0.08|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lepr          |         0|             1| 0.01| 0.09|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lmo4          |         0|             1| 0.09| 0.28|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Irs1          |         0|             1| 0.01| 0.12|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Irs4          |         0|             1| 0.04| 0.19|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Crh           |         0|             1| 0.01| 0.10|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oxt           |         0|             1| 0.40| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Adult") |>
    select(any_of(c(opioid_system_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-pvn-Adult-f3-1.png){fig-align='center' width=10800}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Adult") |>
    select(any_of(c(opioid_system_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |801                |
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
|Crh           |         0|             1| 0.01| 0.10|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oxt           |         0|             1| 0.40| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Adult") |>
    select(any_of(c(
      metabolic_signaling_genes,
      opioid_system_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
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
![](02-endo-cb_files/figure-html/upset-group-e-cb-pvn-Adult-f4-1.png){fig-align='center' width=10800}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
    filter(stage == "Adult") |>
    select(any_of(c(
      metabolic_signaling_genes,
      opioid_system_genes, "Crh", "Trh", "Oxt"))) |> 
    select(where(~ is.numeric(.x) && sum(.x) > 0))
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |801                |
|Number of columns        |19                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |19                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.04| 0.20|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc3r          |         0|             1| 0.01| 0.09|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Mc4r          |         0|             1| 0.01| 0.08|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lepr          |         0|             1| 0.01| 0.09|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Lmo4          |         0|             1| 0.09| 0.28|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Irs1          |         0|             1| 0.01| 0.12|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Irs4          |         0|             1| 0.04| 0.19|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprd1         |         0|             1| 0.00| 0.05|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.04| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprl1         |         0|             1| 0.11| 0.31|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprm1         |         0|             1| 0.02| 0.14|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pcsk1         |         0|             1| 0.09| 0.29|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pcsk2         |         0|             1| 0.29| 0.46|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Pdyn          |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Penk          |         0|             1| 0.04| 0.21|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Pnoc          |         0|             1| 0.02| 0.14|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Crh           |         0|             1| 0.01| 0.10|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oxt           |         0|             1| 0.40| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▅ |
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
 version  R version 4.4.0 (2024-04-24)
 os       Ubuntu 22.04.4 LTS
 system   x86_64, linux-gnu
 ui       X11
 language en_US:en
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Etc/UTC
 date     2024-11-27
 pandoc   3.2 @ /opt/python/3.8.8/bin/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────
 package          * version     date (UTC) lib source
 abind              1.4-8       2024-09-12 [2] RSPM (R 4.4.0)
 anndata          * 0.7.5.6     2023-03-17 [2] RSPM (R 4.4.0)
 assertthat         0.2.1       2019-03-21 [2] RSPM (R 4.4.0)
 base64enc          0.1-3       2015-07-28 [2] RSPM (R 4.4.0)
 bayestestR         0.13.2      2024-02-12 [2] RSPM (R 4.4.0)
 beeswarm           0.4.0       2021-06-01 [2] RSPM (R 4.4.0)
 BiocManager        1.30.23     2024-05-04 [2] RSPM (R 4.4.0)
 bit                4.0.5       2022-11-15 [2] RSPM (R 4.4.0)
 bit64              4.0.5       2020-08-30 [2] RSPM (R 4.4.0)
 circlize           0.4.16      2024-06-19 [2] Github (jokergoo/circlize@9b21578)
 cli                3.6.3       2024-06-21 [2] RSPM (R 4.4.0)
 cluster            2.1.6       2023-12-01 [2] RSPM (R 4.4.0)
 coda               0.19-4.1    2024-01-31 [2] RSPM (R 4.4.0)
 codetools          0.2-20      2024-03-31 [2] RSPM (R 4.4.0)
 colorspace         2.1-1       2024-07-26 [2] RSPM (R 4.4.0)
 correlation        0.8.5       2024-06-16 [2] RSPM (R 4.4.0)
 cowplot          * 1.1.3       2024-01-22 [2] RSPM (R 4.4.0)
 crayon             1.5.3       2024-06-20 [2] RSPM (R 4.4.0)
 data.table         1.16.0      2024-08-27 [2] RSPM (R 4.4.0)
 datawizard         0.11.0      2024-06-05 [2] RSPM (R 4.4.0)
 deldir             2.0-4       2024-02-28 [2] RSPM (R 4.4.0)
 digest             0.6.37      2024-08-19 [2] RSPM (R 4.4.0)
 dotCall64          1.1-1       2023-11-28 [2] RSPM (R 4.4.0)
 dplyr            * 1.1.4       2023-11-17 [2] RSPM (R 4.4.0)
 effectsize         0.8.8       2024-05-12 [2] RSPM (R 4.4.0)
 emmeans            1.10.2      2024-05-20 [2] RSPM (R 4.4.0)
 estimability       1.5.1       2024-05-12 [2] RSPM (R 4.4.0)
 evaluate           1.0.0       2024-09-17 [2] RSPM (R 4.4.0)
 fansi              1.0.6       2023-12-08 [2] RSPM (R 4.4.0)
 farver             2.1.2       2024-05-13 [2] RSPM (R 4.4.0)
 fastDummies        1.7.4       2024-08-16 [2] RSPM (R 4.4.0)
 fastmap            1.2.0       2024-05-15 [2] RSPM (R 4.4.0)
 fitdistrplus       1.2-1       2024-07-12 [2] RSPM (R 4.4.0)
 forcats          * 1.0.0       2023-01-29 [2] RSPM (R 4.4.0)
 future           * 1.34.0      2024-07-29 [2] RSPM (R 4.4.0)
 future.apply       1.11.2      2024-03-28 [2] RSPM (R 4.4.0)
 generics           0.1.3       2022-07-05 [2] RSPM (R 4.4.0)
 ggbeeswarm         0.7.2       2024-06-19 [2] Github (eclarke/ggbeeswarm@ce2da8a)
 ggmin              0.0.0.9000  2024-06-19 [2] Github (sjessa/ggmin@8ada274)
 ggplot2          * 3.5.1       2024-04-23 [2] RSPM (R 4.4.0)
 ggprism            1.0.5       2024-06-19 [2] Github (csdaw/ggprism@b6e6c0e)
 ggrastr            1.0.2       2024-06-19 [2] Github (VPetukhov/ggrastr@50ca3e0)
 ggrepel            0.9.6       2024-09-20 [2] Github (slowkow/ggrepel@e94776b)
 ggridges           0.5.6       2024-01-23 [2] RSPM (R 4.4.0)
 ggstatsplot      * 0.12.3.9000 2024-06-19 [2] Github (IndrajeetPatil/ggstatsplot@d55f86a)
 GlobalOptions      0.1.2       2020-06-10 [2] RSPM (R 4.4.0)
 globals            0.16.3      2024-03-08 [2] RSPM (R 4.4.0)
 glue               1.7.0       2024-01-09 [2] RSPM (R 4.4.0)
 goftest            1.2-3       2021-10-07 [2] RSPM (R 4.4.0)
 gridExtra          2.3         2017-09-09 [2] RSPM (R 4.4.0)
 gtable             0.3.5       2024-04-22 [2] RSPM (R 4.4.0)
 hdf5r              1.3.10      2024-03-02 [2] RSPM (R 4.4.0)
 here             * 1.0.1       2020-12-13 [2] RSPM (R 4.4.0)
 hms                1.1.3       2023-03-21 [2] RSPM (R 4.4.0)
 htmltools          0.5.8.1     2024-04-04 [2] RSPM (R 4.4.0)
 htmlwidgets        1.6.4       2023-12-06 [2] RSPM (R 4.4.0)
 httpuv             1.6.15      2024-03-26 [2] RSPM (R 4.4.0)
 httr               1.4.7       2023-08-15 [2] RSPM (R 4.4.0)
 ica                1.0-3       2022-07-08 [2] RSPM (R 4.4.0)
 igraph             2.0.3       2024-03-13 [2] RSPM (R 4.4.0)
 insight            0.20.1      2024-06-11 [2] RSPM (R 4.4.0)
 irlba              2.3.5.1     2022-10-03 [2] RSPM (R 4.4.0)
 janitor            2.2.0.9000  2024-06-19 [2] Github (sfirke/janitor@80cd1eb)
 jsonlite           1.8.8       2023-12-04 [2] RSPM (R 4.4.0)
 KernSmooth         2.23-24     2024-05-17 [2] RSPM (R 4.4.0)
 knitr              1.48        2024-07-07 [2] RSPM (R 4.4.0)
 labeling           0.4.3       2023-08-29 [2] RSPM (R 4.4.0)
 later              1.3.2       2023-12-06 [2] RSPM (R 4.4.0)
 lattice            0.22-6      2024-03-20 [2] RSPM (R 4.4.0)
 lazyeval           0.2.2       2019-03-15 [2] RSPM (R 4.4.0)
 leiden             0.4.3.1     2023-11-17 [2] RSPM (R 4.4.0)
 lifecycle          1.0.4       2023-11-07 [2] RSPM (R 4.4.0)
 listenv            0.9.1       2024-01-29 [2] RSPM (R 4.4.0)
 lmtest             0.9-40      2022-03-21 [2] RSPM (R 4.4.0)
 lubridate        * 1.9.3       2023-09-27 [2] RSPM (R 4.4.0)
 magrittr         * 2.0.3       2022-03-30 [2] RSPM (R 4.4.0)
 MASS               7.3-61      2024-06-13 [2] RSPM (R 4.4.0)
 Matrix             1.7-0       2024-04-26 [2] RSPM (R 4.4.0)
 matrixStats        1.4.1       2024-09-08 [2] RSPM (R 4.4.0)
 mime               0.12        2021-09-28 [2] RSPM (R 4.4.0)
 miniUI             0.1.1.1     2018-05-18 [2] RSPM (R 4.4.0)
 multcomp           1.4-25      2023-06-20 [2] RSPM (R 4.4.0)
 munsell            0.5.1       2024-04-01 [2] RSPM (R 4.4.0)
 mvtnorm            1.2-5       2024-05-21 [2] RSPM (R 4.4.0)
 nlme               3.1-165     2024-06-06 [2] RSPM (R 4.4.0)
 paletteer          1.6.0       2024-01-21 [2] RSPM (R 4.4.0)
 parallelly         1.38.0      2024-07-27 [2] RSPM (R 4.4.0)
 parameters         0.21.7      2024-05-14 [2] RSPM (R 4.4.0)
 patchwork        * 1.3.0.9000  2024-09-20 [2] Github (thomasp85/patchwork@2695a9f)
 pbapply            1.7-2       2023-06-27 [2] RSPM (R 4.4.0)
 pillar             1.9.0       2023-03-22 [2] RSPM (R 4.4.0)
 pkgconfig          2.0.3       2019-09-22 [2] RSPM (R 4.4.0)
 plotly             4.10.4      2024-01-13 [2] RSPM (R 4.4.0)
 plyr               1.8.9       2023-10-02 [2] RSPM (R 4.4.0)
 png                0.1-8       2022-11-29 [2] RSPM (R 4.4.0)
 polyclip           1.10-7      2024-07-23 [2] RSPM (R 4.4.0)
 progressr          0.14.0      2023-08-10 [2] RSPM (R 4.4.0)
 promises           1.3.0       2024-04-05 [2] RSPM (R 4.4.0)
 purrr            * 1.0.2       2023-08-10 [2] RSPM (R 4.4.0)
 R.methodsS3        1.8.2       2022-06-13 [2] RSPM (R 4.4.0)
 R.oo               1.26.0      2024-01-24 [2] RSPM (R 4.4.0)
 R.utils            2.12.3      2023-11-18 [2] RSPM (R 4.4.0)
 R6                 2.5.1       2021-08-19 [2] RSPM (R 4.4.0)
 RANN               2.6.2       2024-08-25 [2] RSPM (R 4.4.0)
 RColorBrewer     * 1.1-3       2022-04-03 [2] RSPM (R 4.4.0)
 Rcpp               1.0.13      2024-07-17 [2] RSPM (R 4.4.0)
 RcppAnnoy          0.0.22      2024-01-23 [2] RSPM (R 4.4.0)
 RcppHNSW           0.6.0       2024-02-04 [2] RSPM (R 4.4.0)
 readr            * 2.1.5       2024-01-10 [2] RSPM (R 4.4.0)
 rematch2           2.1.2       2020-05-01 [2] RSPM (R 4.4.0)
 remotes            2.5.0       2024-03-17 [2] RSPM (R 4.4.0)
 repr               1.1.7       2024-03-22 [2] RSPM (R 4.4.0)
 reshape2           1.4.4       2020-04-09 [2] RSPM (R 4.4.0)
 reticulate       * 1.39.0      2024-09-05 [2] RSPM (R 4.4.0)
 rlang              1.1.4       2024-06-04 [2] RSPM (R 4.4.0)
 rmarkdown          2.28        2024-08-17 [2] RSPM (R 4.4.0)
 ROCR               1.0-11      2020-05-02 [2] RSPM (R 4.4.0)
 rprojroot          2.0.4       2023-11-05 [2] RSPM (R 4.4.0)
 RSpectra           0.16-2      2024-07-18 [2] RSPM (R 4.4.0)
 rstudioapi         0.16.0      2024-03-24 [2] RSPM (R 4.4.0)
 rsvd               1.0.5       2021-04-16 [2] RSPM (R 4.4.0)
 Rtsne              0.17        2023-12-07 [2] RSPM (R 4.4.0)
 sandwich           3.1-0       2023-12-11 [2] RSPM (R 4.4.0)
 scales             1.3.0       2023-11-28 [2] RSPM (R 4.4.0)
 scattermore        1.2         2023-06-12 [2] RSPM (R 4.4.0)
 scCustomize      * 2.1.2       2024-06-19 [2] Github (samuel-marsh/scCustomize@fc7a282)
 sceasy           * 0.0.7       2024-06-19 [2] Github (cellgeni/sceasy@c1c0bf9)
 sctransform        0.4.1       2023-10-19 [2] RSPM (R 4.4.0)
 sessioninfo        1.2.2       2021-12-06 [2] RSPM (R 4.4.0)
 Seurat           * 5.1.0.9005  2024-09-20 [2] Github (satijalab/seurat@95de9dc)
 SeuratDisk       * 0.0.0.9021  2024-06-19 [2] Github (mojaveazure/seurat-disk@877d4e1)
 SeuratObject     * 5.0.99.9001 2024-09-20 [2] Github (satijalab/seurat-object@1a140c7)
 SeuratWrappers   * 0.3.5       2024-06-19 [2] Github (satijalab/seurat-wrappers@8d46d6c)
 shape              1.4.6.1     2024-02-23 [2] RSPM (R 4.4.0)
 shiny              1.9.1       2024-08-01 [2] RSPM (R 4.4.0)
 skimr            * 2.1.5       2024-06-19 [2] Github (ropensci/skimr@d5126aa)
 snakecase          0.11.1      2023-08-27 [2] RSPM (R 4.4.0)
 sp               * 2.1-4       2024-04-30 [2] RSPM (R 4.4.0)
 spam               2.10-0      2023-10-23 [2] RSPM (R 4.4.0)
 spatstat.data      3.1-2       2024-06-21 [2] RSPM (R 4.4.0)
 spatstat.explore   3.3-2       2024-08-21 [2] RSPM (R 4.4.0)
 spatstat.geom      3.3-3       2024-09-18 [2] RSPM (R 4.4.0)
 spatstat.random    3.3-2       2024-09-18 [2] RSPM (R 4.4.0)
 spatstat.sparse    3.1-0       2024-06-21 [2] RSPM (R 4.4.0)
 spatstat.univar    3.0-1       2024-09-05 [2] RSPM (R 4.4.0)
 spatstat.utils     3.1-0       2024-08-17 [2] RSPM (R 4.4.0)
 statsExpressions   1.5.4       2024-03-20 [2] RSPM (R 4.4.0)
 stringi            1.8.4       2024-05-06 [2] RSPM (R 4.4.0)
 stringr          * 1.5.1       2023-11-14 [2] RSPM (R 4.4.0)
 survival           3.7-0       2024-06-05 [2] RSPM (R 4.4.0)
 tensor             1.5         2012-05-05 [2] RSPM (R 4.4.0)
 TH.data            1.1-2       2023-04-17 [2] RSPM (R 4.4.0)
 tibble           * 3.2.1       2023-03-20 [2] RSPM (R 4.4.0)
 tidyr            * 1.3.1       2024-01-24 [2] RSPM (R 4.4.0)
 tidyselect         1.2.1       2024-03-11 [2] RSPM (R 4.4.0)
 tidyverse        * 2.0.0.9000  2024-06-19 [2] Github (tidyverse/tidyverse@62f32d4)
 timechange         0.3.0       2024-01-18 [2] RSPM (R 4.4.0)
 tzdb               0.4.0       2023-05-12 [2] RSPM (R 4.4.0)
 UpSetR           * 1.4.0       2024-06-19 [2] Github (hms-dbmi/UpSetR@b14854a)
 utf8               1.2.4       2023-10-22 [2] RSPM (R 4.4.0)
 uwot               0.2.2       2024-04-21 [2] RSPM (R 4.4.0)
 vctrs              0.6.5       2023-12-01 [2] RSPM (R 4.4.0)
 vipor              0.4.7       2023-12-18 [2] RSPM (R 4.4.0)
 viridis          * 0.6.5       2024-01-29 [2] RSPM (R 4.4.0)
 viridisLite      * 0.4.2       2023-05-02 [2] RSPM (R 4.4.0)
 vroom              1.6.5       2023-12-05 [2] RSPM (R 4.4.0)
 withr              3.0.1       2024-07-31 [2] RSPM (R 4.4.0)
 xfun               0.47        2024-08-17 [2] RSPM (R 4.4.0)
 xtable             1.8-4       2019-04-21 [2] RSPM (R 4.4.0)
 yaml               2.3.10      2024-07-26 [2] RSPM (R 4.4.0)
 zeallot            0.1.0       2018-01-28 [2] RSPM (R 4.4.0)
 zoo                1.8-12      2023-04-13 [2] RSPM (R 4.4.0)

 [1] /home/etretiakov/R/x86_64-pc-linux-gnu-library/4.4
 [2] /opt/R/4.4.0/lib/R/library

─ Python configuration ───────────────────────────────────────────────────────
 python:         /opt/python/3.8.8/bin/python
 libpython:      /opt/python/3.8.8/lib/libpython3.8.so
 pythonhome:     /opt/python/3.8.8:/opt/python/3.8.8
 version:        3.8.8 | packaged by conda-forge | (default, Feb 20 2021, 16:22:27)  [GCC 9.3.0]
 numpy:          /opt/python/3.8.8/lib/python3.8/site-packages/numpy
 numpy_version:  1.23.5
 scanpy:         /home/etretiakov/.local/lib/python3.8/site-packages/scanpy
 
 NOTE: Python version was forced by RETICULATE_PYTHON

──────────────────────────────────────────────────────────────────────────────
```
:::
:::
