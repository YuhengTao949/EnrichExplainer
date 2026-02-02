construct_interpretation_prompt <- function(database, processed_res, context, diff_genes, ppi, organism, ppi_limit){

  base_prompt <- "You are an expert biologist and bioinformatics researcher."

  if(!is.null(database)){
    base_prompt <- paste(base_prompt, "I have performed a functional enrichment analysis use", database, "database.")
  }else{
    base_prompt <- paste(base_prompt, "I have performed functional enrichment analyses using multiple databases (e.g., KEGG, Reactome, GO, ChEA/Transcription Factors, Disease Ontologies).")
  }

  if (!is.null(context)) {
    base_prompt <- paste0(base_prompt, "\n\n--Experimental Context:\n\n", context)
  }

  item <- processed_res$enrich_data_processed

  if(is.null(item)){
    base_prompt <- NULL
    return(base_prompt)
  }else{
    pathway_text <- paste(
      sapply(1:nrow(item), function(i) {
        row <- item[i, ]
        paste0(i, ". ", paste(names(row), row, sep=": ", collapse="\n"))
      }),
      collapse = "\n\n"
    )
    base_prompt <- paste0(base_prompt, "\n\nTop Enriched Terms:\n\n", pathway_text)

    if(!is.null(diff_genes)){
      fc <- processed_res$gene_logFC
      fc_text <- paste0(names(fc), "=", sprintf("%.2f", round(fc, 2)), collapse = "\n")
      base_prompt <- paste0(base_prompt, "\n\nTop Regulated Genes (Log2 Fold Change):\n\n", fc_text)
    }

    if(ppi){
      ppi_network <- get_ppi_context_text(processed_res$genes, organism, ppi_limit)
      if(!is.null(ppi_network)){base_prompt <- paste0(base_prompt, "\n\nPPI Network (Edge List from STRING):\n\n", ppi_network)}
    }

    base_prompt <- paste0(base_prompt,
                          "\n\n## INSTRUCTIONS FOR ANALYSIS:
Please use a **Chain-of-Thought** approach to analyze these results before generating the final report. Follow these reasoning steps:
1. **Source Deconvolution**: Identify the nature of the enriched terms. Distinguish between:
    - **Biological Processes/Pathways** (e.g., 'Cell Cycle', 'TCR Signaling') -> WHAT is happening.
    - **Upstream Regulators/TFs** (e.g., 'E2F1', 'NFKB1 target genes') -> WHO is driving it.
    - **Phenotypes/Diseases** (e.g., 'Inflammation') -> WHAT is the outcome.
2. **Gene-Level Analysis**: Identify shared key genes and unique functional modules across the pathways.
3. **Causal Integration**: Construct a regulatory narrative connecting the Regulators (TFs) to the Processes (Pathways). For example: 'Enrichment of E2F1 targets explains the observed upregulation of Cell Cycle pathways'.
4. **Contextual Mapping**: Map these themes to the provided experimental context (e.g., tissue, treatment, timepoint) to distinguish expected vs. unexpected findings.
5. **Network Analysis**: Use the provided PPI connections to identify **functional modules** (e.g., a receptor-ligand pair or a protein complex). Identify key hubs that drive the biological process.
6. **Network Refinement**: Prune the provided PPI network to retain only the most biologically relevant edges that explain the context.

Based on this deep analysis, provide a comprehensive biological interpretation in the EXACT following sections in order. Use clear headings (##) for each section:

### I. Integrative Overview
[Concise summary synthesizing the major biological themes]

### II. Regulatory Drivers
[Identify key transcription factors or master regulators found in the enrichment list (or inferred from hubs) and explain their potential role in driving the observed pathways.]

### III. Core Mechanisms
[Explain the underlying biological mechanisms including how pathways interact functionally, grouping pathways into major themes.]

### IV. Pathway Crosstalk
[Discuss potential interactions and regulatory networks between these pathways.]

### V. Biological Hypothesis
[Formulate a coherent biological hypothesis connecting the pathways ('What') to the biological meaning ('So What').Format: regulators → mechanisms → outcomes]

### VI. Narrative Summary
[Write a cohesive paragraph suitable for the 'Results' or 'Discussion' section of a scientific paper.]

### VII. Refined Network
[Based on PPI Network result, provide a list of objects representing the refined core regulatory network. Each object should have 'source', 'target', 'interaction' (e.g., binding, activation, inhibition), and 'reason' (why this edge is relevant).
 The format is:
 List of refined core regulatory network edges in format:
  - Source: [node1], Target: [node2], Interaction: [type], Reason: [justification]
  (Repeat for each relevant edge)]

### VIII. Network Evidence
[Based on PPI Network result, describe specific protein complexes or signaling modules found in the network that support your conclusion.]

## CRITICAL GUIDELINES:

1. **GROUND EVERY CLAIM**: Each biological interpretation MUST reference specific evidence from the enrichment results.
   - Format: '...process (supported by: [pathway names])'
   - Format: '...regulated by (supported by: [TF names])'

2. **PRIORITIZE BY EVIDENCE**: Focus interpretation on terms with strongest statistical support (highest enrichment).

3. **DISTINGUISH OBSERVATION FROM SPECULATION**: Clearly separate:
   - Direct observations from the data
   - Established biological knowledge
   - Reasonable inferences
   - Speculative hypotheses

4. **ACKNOWLEDGE LIMITATIONS**: Note if results are:
   - Consistent with expected biology
   - Potentially artifactual
   - Limited by database coverage
   - Conflicting with other evidence

5. **MAINTAIN SCIENTIFIC RIGOR**:
   - Use standard gene/protein nomenclature
   - Distinguish between correlation and causation
   - Avoid overinterpretation of weak associations

Please be scientifically rigorous, citing standard biological knowledge where appropriate, and avoid hallucinations.")

    return(base_prompt)

  }

}
