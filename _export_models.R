# Export models 
library(jsonlite)

models <- readRDS("workflow/models/models.rds")

extract_model <- function(m) {
  list(
    mem  = as.list(coef(m$mem)),
    time = as.list(coef(m$time))
  )
}

models_json <- lapply(models, extract_model)

write_json(
  models_json,
  path = "workflow/models/models.json",
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("Created workflow/models/models.json\n")