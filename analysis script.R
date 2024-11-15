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

#load datasets
data(geneList, package = "DOSE")
str(geneList)

gene.df <- bitr(names(geneList), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

#ORA
UR_genes <- names(geneList)[geneList > 2]
head(UR_genes)
#or to look at the equivalent gene Symbols, which are typically easier to interpret
gene.df$SYMBOL[gene.df$ENTREZID %in%  head(UR_genes)]

ego <- enrichGO(gene          = UR_genes,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "bonferroni",
                readable      = TRUE)
names(ego@result)
#“ID” - The GO term ID
#“Description” - The name of the GO term
#“GeneRatio” - The ratio of submitted genes in the GO term
#“BgRatio” - The ratio of gene “universe” in the GO term
#“pvalue” - The hypergeometric test nominal p-value
#“p.adjust” - The adjusted p value (in this case bonferroni adjusted)
#“qvalue” - The q value associated with False discovery rate correction
#“geneID” - The gene ids associated with each GO term
#“Count” - The number of submitted genes in each GO term

dotplot(ego, showCategory=10, label_format = 30) +
  ggtitle("GO-BP ORA") + 
  theme(axis.text.y = element_text(size = 7))

ego <- pairwise_termsim(ego)
emapplot(ego, showCategory = 20, cex.params = list(category_label = 0.5, cex_line = 0.1))

#ego_simple <- simplify(ego, cutoff=0.01, by="p.adjust", select_fun=min)
#emapplot(ego_simple, showCategory = 20, cex.params = list(category_label = 0.5, cex_line = 0.1))

cnetplot(ego, categorySize="pvalue",
         color.params = list(foldChange = geneList),
         cex.params = list(gene_label = 0.5, category_label = 0.7),
         layout = "kk")

x <- enrichKEGG(gene=UR_genes, universe = names(geneList),
                organism     = 'hsa',
                pvalueCutoff = 0.05)

head(x@result)

dotplot(x, showCategory=10, label_format = 30) +
  ggtitle("Kegg Pathway ORA") + 
  theme(axis.text.y = element_text(size = 7))

hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

DR_genes <- names(geneList)[geneList < -2]
ck <- compareCluster(geneCluster = list(UR = UR_genes, DR = DR_genes), fun = "enrichKEGG",
                     organism = "hsa", universe = names(geneList))
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(ck)

#run below in console instead of script
ck <- pairwise_termsim(ck)
emapplot(ck, showCategory = 10, 
         cex.params = list(category_label = 1, cex_line = 0.2, pie2axis = 0.1),
         pie.params = list(legend_n = 2))

cnetplot(ck, 
         cex.params = list(gene_label = 0.75, group_label = 1),
         pie.params = list(legend_n = 2)) + #legend_loc_x and legend_n edit parameters associated with the node size scatterpie legend.
  theme_void(base_size = 7) 

###GSEA
plot(geneList)

y <- gseGO(geneList, OrgDb = org.Hs.eg.db,
           pvalueCutoff = 0.2,
           pAdjustMethod = "BH", 
           eps = 1e-100,
           verbose = FALSE)
y <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')

dotplot(y, showCategory=10, label_format = 30) +
  ggtitle("dotplot for ORA") + 
  theme(axis.text.y = element_text(size = 7))

y@result %>%
  mutate(label = ifelse(rank(qvalue) <= 5, Description, "")) %>%
  ggplot(aes(x = NES, y = -log10(qvalue), label = label)) +
  geom_point(aes(size = setSize, colour = p.adjust < 0.01)) +
  scale_size(range = c(0.5,2)) +
  scale_color_manual(values = c("grey70", "red")) +
  ggrepel::geom_text_repel(size = 2.5, min.segment.length = 0, force = 10) +
  theme_classic() 

head(y@result)

gs2p <- 1

gseaplot2(y, geneSetID = gs2p, title = y$Description[gs2p],
          pvalue_table = T)

gs2p <- 1

plot <- gseaplot2(y, geneSetID = gs2p,
                  title = paste(y$Description[gs2p]))
plot[[1]] <- plot[[1]] + labs(subtitle = paste0("(NES=", signif(y@result$NES[gs2p],3), 
                                                "; FDR q=", signif(y@result$qvalue[gs2p],3), ")"))

plot

##Part 2

load("C:/Users/ahgod/Documents/MEDI6227/MEDI6227-workshop-4/Session4_practical(1).Rdata")

#convert gene symbols to ENTREZID ID - create in new DF
gene.df1 <- bitr(DE.blocked.treat$gene_name, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

DE.blocked.treat <- left_join(DE.blocked.treat, gene.df1, by = c("gene_name" = "SYMBOL"))

DE.blocked.treat$DE_cat <- 
  ifelse(DE.blocked.treat$logFC > 0 &  DE.blocked.treat$FDR < 0.01, "M2", 
         ifelse(DE.blocked.treat$logFC <0 &  DE.blocked.treat$FDR < 0.01, "M1", NA))

pA <- DE.blocked.treat %>%
  ggplot(aes(x = logFC, y = -log10(FDR), colour = DE_cat)) +
  geom_point(size = 0.1) +
  theme_pubr(base_size = 7) + 
  theme(legend.title = element_blank(), legend.position = "right", legend.key.size = unit(5, "pt")) +
  scale_color_discrete()

cols <- ggpubr::get_palette(k = 2)

ck <- compareCluster(ENTREZID ~ DE_cat, data = DE.blocked.treat,
                     fun = enrichGO,
                     OrgDb = org.Hs.eg.db, 
                     universe = DE.blocked.treat$ENTREZID,
                     ont = "BP", readable = T)

pB <- dotplot(ck, label_format = 50) +
  theme_pubr(base_size = 7) + 
  scale_size(range = c(0.5,2)) +
  theme(legend.position = "right", legend.key.size = unit(5, "pt"))

M1_ego <- enrichGO(gene          = DE.blocked.treat$ENTREZID[DE.blocked.treat$DE_cat == "M1"],
                   universe      = DE.blocked.treat$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   readable      = TRUE)

M1_ego <- pairwise_termsim(M1_ego)
pC <- cnetplot(M1_ego, showCategory = c("cytokine-mediated signaling pathway", "response to type II interferon"),
               cex.params = list(gene_label = 0.5, category_label = 0.5),
               color.params = list(category = cols[1])) +
  theme_void(base_size = 7)

M2_ego <- enrichGO(gene          = DE.blocked.treat$ENTREZID[DE.blocked.treat$DE_cat == "M2"],
                   universe      = DE.blocked.treat$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   readable      = TRUE)

M2_ego <- pairwise_termsim(M2_ego)

pD <- cnetplot(M2_ego, showCategory = c("response to toxic substance", "positive regulation of wound healing"),
               cex.params = list(gene_label = 0.5, category_label = 0.5),
               color.params = list(category = cols[2])) +
  theme_void(base_size = 7) 

pABCD <- 
  ggarrange(ggarrange(pA, pB, ncol = 2, nrow = 1, labels = "AUTO", widths = c(1,2), align = "h"),
            pC,pD,
            ncol = 1, nrow = 3, labels = c("", "C", "D") )
pABCD
