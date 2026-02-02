get_ppi_context_text <- function(genes, organism, ppi_limit) {
  if (length(genes) == 0) return(NULL)

  if(!is.null(ppi_limit)){input_for_ppi <- head(genes, ppi_limit)}else{input_for_ppi <- genes}

  # Try to determine taxID
  current_taxID <- "auto"
  if (!is.null(organism) && !is.list(organism)) {
    current_taxID <- getTaxID(organism)
  }

  ppi_res <- tryCatch({
    g <- getPPI(input_for_ppi, taxID = current_taxID, output = "igraph", network_type = "functional")

    if (!is.null(g)) {
      el <- igraph::as_data_frame(g, what="edges")
      if (nrow(el) > 0) {
        if ("score" %in% names(el)) el <- el[order(el$score, decreasing = TRUE), ]
        if (!is.null(ppi_limit)){el_subset <- head(el, ppi_limit)}else{el_subset <- el}
        edges_text <- apply(el_subset, 1, function(row) {
          score_info <- ""
          if ("score" %in% names(row)) score_info <- paste0(" (Score: ", row["score"], ")")
          paste0(row["from"], " -- ", row["to"], score_info)
        })
        paste(edges_text, collapse = "\n")
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, error = function(e) NULL)

  return(ppi_res)
}


networkParamsParser <- function(
    identifiers,
    species,
    required_score = NULL,
    network_type = "functional",
    add_nodes = 1,
    show_query_node_labels = 0,
    caller_identity = NULL
) {
  # Format the identifiers
  identifiers <- paste(identifiers, collapse = "\n")

  # Check parameters
  if (missing(species)) {
    stop("Please provide an NCBI taxon identifier for the species.")
  }

  # Create parameters list
  params <- list(
    identifiers = identifiers,
    species = species,
    required_score = required_score,
    network_type = network_type,
    add_nodes = add_nodes,
    show_query_node_labels = show_query_node_labels,
    caller_identity = caller_identity
  )

  # Remove NULL elements from the list
  filtered_params <- Filter(Negate(is.null), params)
  return(filtered_params)
}

getPPI <- function(
    x,
    ID=1,
    taxID = "auto",
    required_score = NULL,
    network_type = "functional",
    add_nodes = 0,
    show_query_node_labels = 0,
    output = 'igraph') {

  output <- match.arg(output, c("igraph", "data.frame"))
  if (taxID == "auto") {
    if (!inherits(x, 'enrichResult') || is.null(x@organism) || length(x@organism) == 0) {
      stop("Unable to determine taxonomy ID. You can use the `getTaxID()` to get the taxonomy ID from scientific species name.\n")
    }
    taxID <- getTaxID(x@organism)
  }
  if (inherits(x, 'enrichResult')) {
    if (is(ID, 'numeric')) {
      id <- x$ID[ID]
    } else {
      id <- ID
    }
    genes <- unlist(geneInCategory(x)[id])
  } else {
    genes <- x
  }

  networkParams <- networkParamsParser(
    identifiers = genes,
    species = taxID,
    required_score = required_score,
    network_type = network_type,
    add_nodes = add_nodes,
    show_query_node_labels = show_query_node_labels
  )

  # Set stringDB base URL
  address <- "https://string-db.org"

  # Validate the address
  #
  httr::stop_for_status(httr::GET(address))

  # read data from stringDB api
  response <- httr::GET(paste(address, "/api/tsv/network", sep = ""), query = networkParams)
  res <- yulab.utils::yread_tsv(response$url, params = list(header = TRUE))

  if (output == "data.frame") {
    return(res)
  }

  node <- unique(c(res$preferredName_A, res$preferredName_B))

  igraph::graph_from_data_frame(
    d = unique(res[,c(3,4,6)]),
    vertices = node,
    directed=FALSE
  )
}
