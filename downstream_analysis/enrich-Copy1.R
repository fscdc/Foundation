library(stringr)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(biomaRt)

hsa <- org.Hs.eg.db


folder_path <- "/home/foundation/program/foundation/record/pvalue-0.05"

# get all .txt files in the folder
txt_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

# iterate over all .txt files
for (file in txt_files) {
    file_name <- basename(file)

    pattern <- "^(.*)-(.*)-(.*)-(.*)-(.*)-DEG.txt$"
    matches <- str_match(file_name, pattern)

    dataset <- matches[1,2]
    genelist_name <- matches[1,3]
    cluster_method <- matches[1,4]
    resolution <- matches[1,5]
    embedding_method <- matches[1,6]

    cat("Doing enrichment for :", file_name, "\n")

    allmarker <- read.table(file, header = TRUE, sep = "\t")

    # GO & KEGG enrichment analysis
    for (i in unique(allmarker$group)) {
        cluster=subset(allmarker,group==i)

        cluster_sorted <- cluster %>% arrange(pvals_adj)
        
        cluster_top15 <- cluster_sorted[1:30, ]
        
        gene_ensembl <- unique(cluster_top15$names)
        
        genelist <- bitr(gene_ensembl, fromType = "SYMBOL", toType = c("ENTREZID", "GENENAME"), OrgDb = hsa)

        print(genelist)
        
        # GO-------------------------
        ego <- enrichGO(
        gene=genelist$ENTREZID,
        OrgDb = hsa, 
        ont='ALL',
        pAdjustMethod = 'BH',
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        readable = TRUE
        )
        
        ego <- setReadable(ego, OrgDb = hsa, keyType = "ENTREZID")
    #     go_data <- ego@result %>% as.data.frame()
    #     print("GO analysis result:")
    #     print(go_data)
        

        dotplot_figure <- dotplot(ego, split = "ONTOLOGY", showCategory = 4)
        dotplot_path <- paste0("/home/foundation/program/Foundation/record/figures/enrich/",dataset, "-", genelist_name, "-", cluster_method, "-", resolution, "-", embedding_method, "-GO_dotplot_cluster", i, ".pdf")
        ggsave(dotplot_path, plot = dotplot_figure, width = 8, height = 6)
        
        barplot_figure <- barplot(ego, split = "ONTOLOGY", showCategory = 4)
        barplot_path <- paste0("/home/foundation/program/Foundation/record/figures/enrich/",dataset, "-", genelist_name, "-", cluster_method, "-", resolution, "-", embedding_method, "-GO_barplot_cluster", i, ".pdf")
        ggsave(barplot_path, plot = barplot_figure, width = 8, height = 6)
        
        
        # KEGG-------------------------
        enrich.kegg <- enrichKEGG(
        gene = genelist$ENTREZID,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.05,
        use_internal_data = FALSE
        )
        
        enrich.kegg <- setReadable(enrich.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    #     kegg_data <- enrich.kegg@result %>% as.data.frame()
    #     print("KEGG analysis result:")
    #     print(kegg_data)
        

        dotplot_figure <- dotplot(enrich.kegg, split = "Description", showCategory=4)
        dotplot_path <- paste0("/home/foundation/program/Foundation/record/figures/enrich/",dataset, "-", genelist_name, "-", cluster_method, "-", resolution, "-", embedding_method, "-KEGG_dotplot_cluster", i, ".pdf")
        ggsave(dotplot_path, plot = dotplot_figure, width = 8, height = 6)
        

        barplot_figure <- barplot(enrich.kegg, split = "Description", showCategory=4)
        barplot_path <- paste0("/home/foundation/program/Foundation/record/figures/enrich/",dataset, "-", genelist_name, "-", cluster_method, "-", resolution, "-", embedding_method, "-KEGG_barplot_cluster", i, ".pdf")
        ggsave(barplot_path, plot = barplot_figure, width = 8, height = 6)
        
    }

    # GSEA enrichment analysis

    gene_symbols <- rownames(allmarker)
    genelist <- setNames(allmarker$logfoldchanges, allmarker$names)
    genelist <- sort(genelist, decreasing = TRUE)

    genelist_entrez <- bitr(names(genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = hsa)

    matched_indices <- match(names(genelist), genelist_entrez$SYMBOL)
    genelist <- genelist[!is.na(matched_indices)]
    genelist_entrez <- genelist_entrez[matched_indices[!is.na(matched_indices)], "ENTREZID"]

    names(genelist) <- genelist_entrez
    genelist <- sort(genelist, decreasing = TRUE)


    gseaGO <- gseGO(
    geneList = genelist,
    OrgDb = hsa,
    keyType = "ENTREZID",
    ont = "BP",
    nPerm = 1000,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
    )

    library(enrichplot)

    nes_results <- gseaGO@result

    sorted_nes <- nes_results[order(nes_results$NES), ]

    top_n <- 10
    top_geneSetIDs <- sorted_nes$ID[1:top_n]

    gseaplot <- gseaplot2(gseaGO, geneSetID = top_geneSetIDs, pvalue_table = TRUE, title = "Top10 NES Pathways", ES_geom = "line")
    gseaplot_path <- paste0("/home/foundation/program/Foundation/record/figures/enrich/",dataset, "-", genelist_name, "-", cluster_method, "-", resolution, "-", embedding_method, "-GSEA-NES.pdf")
    ggsave(gseaplot_path, plot = gseaplot, width = 30, height = 20)

    sorted_nes <- nes_results[order(nes_results$p.adjust), ]

    top_n <- 10
    top_geneSetIDs <- sorted_nes$ID[1:top_n]

    gseaplot <- gseaplot2(gseaGO, geneSetID = top_geneSetIDs, pvalue_table = TRUE, title = "Top10 P_adjust Pathways", ES_geom = "line")
    gseaplot_path <- paste0("/home/foundation/program/Foundation/record/figures/enrich/",dataset, "-", genelist_name, "-", cluster_method, "-", resolution, "-", embedding_method, "-GSEA-Pval.pdf")
    ggsave(gseaplot_path, plot = gseaplot, width = 30, height = 20)

}