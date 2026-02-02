contact_with_LLM <- function(prompt,
                             api_key,
                             base_url,
                             model) {

  Sys.setenv(ARK_API_KEY = api_key)

  tryCatch({
    os <- reticulate::import("os")
    openai <- reticulate::import("openai")

    client <- openai$OpenAI(
      base_url = base_url,
      api_key = os$environ$get("ARK_API_KEY")
    )

    messages <- list(
      list(role = "system", content = "You are an expert biologist and bioinformatics researcher."),
      list(role = "user", content = prompt)
    )

    completion <- client$chat$completions$create(
      model = model,
      messages = messages
    )

    cat(completion$choices[[1]]$message$content)
    return(completion$choices[[1]]$message$content)

  }, error = function(e) {
    message("Error in API call: ", e$message)
    return(NULL)
  })
}
