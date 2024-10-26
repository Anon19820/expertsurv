library(shiny)

required_packages <- c("expertsurv")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  devtools::install_github("Anon19820/expertsurv", ref  = "non-compile") 
}

expertsurv::elicit_surv(compile_mods = expertsurv::compiled_models_saved)

