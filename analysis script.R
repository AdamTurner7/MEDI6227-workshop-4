###MEDI6227 workshop 4
CRAN_packages <- c("tidyverse", "cli", "BiocManager", "viridis",
                   "pheatmap", "cowplot", "factoextra", "ggnewscale", "ggpubr", "ggrepel")
Bioc_packages <- c("edgeR", "DESeq2", "EnhancedVolcano", "apeglm", "vsn", "clusterProfiler",
                   "org.Hs.eg.db", "HPO.db", "enrichplot", "GOSemSim", "pathview")
# Install missing packages
installed_CRAN_packages <- CRAN_packages %in% rownames(installed.packages())
installed_Bioc_packages <- Bioc_packages %in% rownames(installed.packages())
current_packages <- c(installed_Bioc_packages, installed_CRAN_packages)

if (any(current_packages == FALSE)) {
  install.packages(CRAN_packages[!installed_CRAN_packages])
  BiocManager::install(Bioc_packages[!installed_Bioc_packages])
}


# Load required packages
invisible(lapply(c(CRAN_packages, Bioc_packages), library, character.only = TRUE))
