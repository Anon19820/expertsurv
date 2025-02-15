rmarkdown::render("README.Rmd", output_format = "html_document", output_file = "inst/app/www/README.html")
rmarkdown::render("README.Rmd", output_format = "html_document", output_file = "README.html")
rmarkdown::render(paste0(getwd(),"/vignettes/ShinyExpertsurv-Vignette.Rmd"), output_format = "html_document", output_file = "C:/Users/phili/OneDrive/PhD/R_packages_2023/expertsurv/inst/app/www/ShinyExpertsurv-Vignette.html")

source_file <- "vignettes/ShinyExpertsurv-Vignette.html"
dest_file <- "inst/app/www/ShinyExpertsurv-Vignette.html"

# Ensure the destination directory exists
if (!dir.exists(dirname(dest_file))) {
  dir.create(dirname(dest_file), recursive = TRUE)
}

# Copy the file
file.copy(source_file, dest_file, overwrite = TRUE)

rmarkdown::render(paste0(getwd(),"/vignettes/ShinyExpertsurv-Vignette.Rmd"), output_format = "html_document")


