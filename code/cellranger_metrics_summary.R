library(tidyverse)
metrics_summary_SRR9214090 <- read_csv("/data/PRJNA547712/cellranger/SRR9214090/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214090")
metrics_summary_SRR9214091 <- read_csv("/data/PRJNA547712/cellranger/SRR9214091/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214091")
metrics_summary_SRR9214092 <- read_csv("/data/PRJNA547712/cellranger/SRR9214092/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214092")
metrics_summary_SRR9214093 <- read_csv("/data/PRJNA547712/cellranger/SRR9214093/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214093")
metrics_summary_SRR9214094 <- read_csv("/data/PRJNA547712/cellranger/SRR9214094/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214094")
metrics_summary_SRR9214095 <- read_csv("/data/PRJNA547712/cellranger/SRR9214095/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214095")
metrics_summary_SRR9214097 <- read_csv("/data/PRJNA547712/cellranger/SRR9214097/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214097")
metrics_summary_SRR9214099 <- read_csv("/data/PRJNA547712/cellranger/SRR9214099/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214099")
metrics_summary_SRR9214100 <- read_csv("/data/PRJNA547712/cellranger/SRR9214100/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214100")
metrics_summary_SRR9214101 <- read_csv("/data/PRJNA547712/cellranger/SRR9214101/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214101")
metrics_summary_SRR9214102 <- read_csv("/data/PRJNA547712/cellranger/SRR9214102/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214102")
metrics_summary_SRR9214103 <- read_csv("/data/PRJNA547712/cellranger/SRR9214103/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214103")
metrics_summary_SRR9214104 <- read_csv("/data/PRJNA547712/cellranger/SRR9214104/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214104")
metrics_summary_SRR9214105 <- read_csv("/data/PRJNA547712/cellranger/SRR9214105/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214105")
metrics_summary_SRR9214106 <- read_csv("/data/PRJNA547712/cellranger/SRR9214106/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214106")
metrics_summary_SRR9214107 <- read_csv("/data/PRJNA547712/cellranger/SRR9214107/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214107")
metrics_summary_SRR9214108 <- read_csv("/data/PRJNA547712/cellranger/SRR9214108/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214108")
metrics_summary_SRR9214109 <- read_csv("/data/PRJNA547712/cellranger/SRR9214109/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214109")
metrics_summary_SRR9214110 <- read_csv("/data/PRJNA547712/cellranger/SRR9214110/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214110")
metrics_summary_SRR9214111 <- read_csv("/data/PRJNA547712/cellranger/SRR9214111/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214111")
metrics_summary_SRR9214112 <- read_csv("/data/PRJNA547712/cellranger/SRR9214112/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214112")
metrics_summary_SRR9214113 <- read_csv("/data/PRJNA547712/cellranger/SRR9214113/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214113")
metrics_summary_SRR9214114 <- read_csv("/data/PRJNA547712/cellranger/SRR9214114/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214114")
metrics_summary_SRR9214115 <- read_csv("/data/PRJNA547712/cellranger/SRR9214115/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214115")
metrics_summary_SRR9214116 <- read_csv("/data/PRJNA547712/cellranger/SRR9214116/outs/metrics_summary.csv") %>% mutate(Run = "SRR9214116")
metrics_summary <-
    bind_rows(
        metrics_summary_SRR9214090,
        metrics_summary_SRR9214091,
        metrics_summary_SRR9214092,
        metrics_summary_SRR9214093,
        metrics_summary_SRR9214094,
        metrics_summary_SRR9214095,
        metrics_summary_SRR9214097,
        metrics_summary_SRR9214099,
        metrics_summary_SRR9214100,
        metrics_summary_SRR9214101,
        metrics_summary_SRR9214102,
        metrics_summary_SRR9214103,
        metrics_summary_SRR9214104,
        metrics_summary_SRR9214105,
        metrics_summary_SRR9214106,
        metrics_summary_SRR9214107,
        metrics_summary_SRR9214108,
        metrics_summary_SRR9214109,
        metrics_summary_SRR9214110,
        metrics_summary_SRR9214111,
        metrics_summary_SRR9214112,
        metrics_summary_SRR9214113,
        metrics_summary_SRR9214114,
        metrics_summary_SRR9214115,
        metrics_summary_SRR9214116)

metrics_summary |>
    select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, "/data/PRJNA547712/metrics_summary.tsv")

