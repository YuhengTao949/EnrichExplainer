#' Interpret Enrichment Results using Large Language Models
#'
#' Processes enrichment results and optionally generates pathway interpretations
#' using an LLM API. Includes PPI network construction and context integration.
#'
#' @param enrich_res Enrichment analysis result object. Accepts standard
#'        \code{clusterProfiler} enrichment objects (e.g., \code{enrichResult})
#'        and \code{compareClusterResult} objects for multi-group analyses.
#' @param database Character string specifying the annotation database used
#'        (e.g., "KEGG", "GO_BP", "Reactome"). Must match the database used in
#'        \code{enrich_res}.
#' @param diff_genes Named numeric vector of differentially expressed genes.
#'        Genes must be sorted by log2FC in descending order. Format:
#'        \preformatted{
#'        c(GeneA = 3.2, GeneB = 2.5, GeneC = -1.8)
#'        }
#'        Typically generated with:\cr
#'        \code{sort(de_results$log2FoldChange, decreasing = TRUE)}
#' @param context Optional biological context string (e.g., "breast cancer",
#'        "liver tissue") to guide LLM interpretation. Default = \code{NULL}.
#' @param contact_LLM Logical indicating whether to contact LLM API.
#'        When \code{FALSE}, returns constructed prompts. Default = \code{TRUE}.
#' @param pathway_num Maximum number of significant pathways to include.
#'        If \code{NULL}, uses all pathways in \code{enrich_res}. Default = \code{NULL}.
#' @param ppi Logical indicating whether to incorporate protein-protein
#'        interaction data. Default = \code{TRUE}.
#' @param ppi_limit Maximum number of genes to include in PPI network.
#'        Applied when \code{ppi = TRUE}. If \code{NULL}, uses all genes in
#'        \code{names(diff_genes)}. Default = \code{NULL}.
#' @param organism Organism identifier in standard format (e.g., \code{"hsa"}
#'        for human, \code{"mmu"} for mouse). Default = \code{"hsa"}.
#' @param base_url LLM API endpoint URL. Required when \code{contact_LLM = TRUE}.
#' @param api_key API key for LLM authentication. Required when
#'        \code{contact_LLM = TRUE}.
#' @param model LLM model identifier (e.g., \code{"gpt-4-turbo"}). Required when
#'        \code{contact_LLM = TRUE}.
#'
#' @return
#' If \code{contact_LLM = TRUE}:
#' \itemize{
#'   \item Single enrichment input: Character string containing LLM interpretation
#'   \item Multi-group input (\code{compareClusterResult}): Named list of
#'         interpretation strings (names = cluster IDs)
#' }
#' If \code{contact_LLM = FALSE}: Returns prompt(s) as character strings
#' ready for API submission.

#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Single enrichment interpretation
#' interpret_tool(
#'   enrich_res = kegg_result,
#'   database = "KEGG",
#'   diff_genes = c(CD4 = 2.1, IL6 = 1.8, ...),
#'   base_url = "https://api.openai.com/v1",
#'   api_key = "your_api_key",
#'   model = "gpt-4-turbo"
#' )
#'
#' # Multi-group interpretation
#' interpret_tool(
#'   enrich_res = compareCluster_result,
#'   database = "GO",
#'   diff_genes = c(TP53 = 3.5, MYC = 2.8, ...),
#'   pathway_num = 30,
#'   context = "Lung adenocarcinoma scRNA-seq",
#'   base_url = "https://api.anthropic.com/v1",
#'   api_key = "claude_key",
#'   model = "claude-3-opus"
#' )}
interpret_tool <- function(enrich_res,
                           database,
                           diff_genes,
                           context = NULL,
                           contact_LLM = FALSE,
                           pathway_num = NULL,
                           ppi = TRUE,
                           ppi_limit = NULL,
                           organism = "hsa",
                           base_url,
                           api_key,
                           model){

  cat("\n---------- Process enrich result ----------\n\n")
  processed_res <- process_enrich_results(enrich_res,
                                          diff_genes,
                                          pathway_num)

  cat("---------- Construct prompt ----------\n\n")
  if(inherits(enrich_res, "compareClusterResult")){
    prompt_list <- lapply(names(processed_res), function(name){
      per_res <- processed_res[[name]]

      if(is.null(per_res$enrich_data_processed)){
        return(NULL)
      }else{
        prompt <- paste0("\n\n========== Cluster: ", name, "========== \n\n")
        prompt <- paste0(prompt,
                         construct_interpretation_prompt(database, per_res, context, diff_genes, ppi, organism, ppi_limit))
        return(prompt)
      }

    })
    names(prompt_list) <- names(processed_res)
  }else{

    prompt <- construct_interpretation_prompt(database, processed_res, context, diff_genes, ppi, organism, ppi_limit)
  }

  if(contact_LLM){
    cat("\n\n---------- Contact with LLM ----------\n\n")
    if(inherits(enrich_res, "compareClusterResult")){
      reslist <- list()
      for(i in names(prompt_list)){
        cat("\n\n========== Cluster: ", i, "==========\n\n")
        prompt <- prompt_list[[i]]
        if(is.null(prompt)){
          cat("Cluster:", i, "No significant pathways enriched and no marker genes available for interpretation. Skip this.")
          reslist[[i]] <- NULL
        }else{
          reslist[[i]] <- contact_with_LLM(prompt, api_key, base_url, model)
        }

      }
      return(reslist)
    }else{
      res <- contact_with_LLM(prompt, api_key, base_url, model)
      return(res)
    }
  }else{
    if(inherits(enrich_res, "compareClusterResult")){
      for (i in names(prompt_list)){
        return(prompt_list)
      }
    }else{
      return(prompt)
    }
  }

}
