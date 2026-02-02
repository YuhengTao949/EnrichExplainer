CleanInputData <- function(enrich_res, pathway_num){
  if (!is.data.frame(enrich_res)){
    enrich_data <- as.data.frame(enrich_res)
  }else{
    enrich_data <- enrich_res
  }

  enrich_data <- enrich_data %>%
    dplyr::filter(p.adjust < 0.05)

  if("FoldEnrichment" %in% colnames(enrich_data)){
    enrich_data <- enrich_data %>% dplyr::arrange(desc(FoldEnrichment))
  }else{
    enrich_data <- enrich_data %>% dplyr::arrange(desc(NES))
  }

  if (!is.null(pathway_num)){
    enrich_data <- enrich_data %>% dplyr::slice_head(n = pathway_num)
  }

  cols <- intersect(c("ID", "Description", "GeneRatio", "FoldEnrichment",  "enrichmentScore", "NES", "p.adjust", "geneID", "core_enrichment"), colnames(enrich_data))
  enrich_data <- enrich_data[, cols]

  return(enrich_data)
}

process_compareClusterResult <- function(enrich_res, pathway_num){

  enrich_data <- as.data.frame(enrich_res)
  enrich_list <- split(enrich_data, enrich_data$Cluster)
  processed_data <- lapply(names(enrich_list), function(cluster){
    data_use <- enrich_list[[cluster]]
    if (is.data.frame(data_use) && nrow(data_use) == 0){
      res <- NULL
    }else{
      res <- CleanInputData(data_use, pathway_num)
    }
    return(res)
  })
  names(processed_data) <- names(enrich_list)
  return(processed_data)
}

process_enrich_results <- function(enrich_res, diff_genes, pathway_num){

  if(!inherits(enrich_res, "enrichResult") && !inherits(enrich_res, "gseaResult") &&
     !inherits(enrich_res, "compareClusterResult") && !is.data.frame(x)){
    stop("Enrichment result must be an enrichResult, gseaResult, or compareClusterResult object.")
  }

  if(inherits(enrich_res, "compareClusterResult")){

    processed_data <- process_compareClusterResult(enrich_res, pathway_num)

    res <- lapply(names(processed_data), function(x){

      dat <- processed_data[[x]]

      if (is.null(dat)){
        r <- list(enrich_data_processed = dat)
      }else{
        genes_include <- strsplit(dat[,ncol(dat)], split = "/") %>% unlist() %>% unique()
        logfc <- diff_genes[intersect(names(diff_genes), genes_include)]
        logfc <- sort(logfc, decreasing = TRUE)
        # logfc <- head(logfc, 100)

        r <- list(
          enrich_data_processed = dat,
          gene_logFC = logfc,
          genes = names(logfc)
        )
      }

    })
    names(res) <- names(processed_data)
  }else{
    processed_data <- CleanInputData(enrich_res, pathway_num)

    genes_include <- strsplit(processed_data[, ncol(processed_data)], split = "/") %>% unlist() %>% unique()
    logfc <- diff_genes[intersect(names(diff_genes), genes_include)]
    logfc <- sort(logfc, decreasing = TRUE)
    # logfc <- head(logfc, 100)

    res <- list(
      enrich_data_processed = processed_data,
      gene_logFC = logfc,
      genes = names(logfc)
    )
  }

  return(res)
}
