response <- function(text, status = 200L) {
  structure(list(status_code = status, content = charToRaw(text),
                 headers = list("content-type" = "application/json"), url = "https://test.invalid"),
            class = "response")
}
