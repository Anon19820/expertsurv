
options(rmarkdown.html_vignette.check_title = FALSE)
# build_vignette.R
library(rmarkdown)
library(fs)

# Path to your README.Rmd
readme_path <- "README.Rmd"

# Path to your vignettes directory
vignette_dir <- "vignettes"

# Ensure the vignettes directory exists
if (!dir_exists(vignette_dir)) {
  dir_create(vignette_dir)
}

# Path to the new vignette
vignette_path <- file.path(vignette_dir, "Expertsurv-Vignette.Rmd")

if(FALSE){
  

# Copy README.Rmd to the vignettes directory
file_copy(readme_path, vignette_path, overwrite = TRUE)

# Read the contents of README.Rmd
readme_content <- readLines(vignette_path)

# Replace the specified text
start_end <- grep("^---$", readme_content)
replacement <- c(
  "---",
  "title: \"<img src=\\" `r system.file('figures/hexsticker.png', package = 'expertsurv')`"\\align='left' height='150'/> Introduction to Expertsurv\"",
  "output: rmarkdown::html_vignette",
  "bibliography: \"`r system.file('REFERENCES.bib', package = 'expertsurv')`\"",
  "vignette: >",
  "  %\\VignetteIndexEntry{Introduction to Expertsurv}",
  "  %\\VignetteEngine{knitr::rmarkdown}",
  "  %\\VignetteEncoding{UTF-8}",
  "---"
)


#knitr::include_graphics(system.file("image/Vignette_Example_1_DIC.png", package = "expertsurv"))

if (length(start_end) == 2) {
  readme_content <- c(replacement, readme_content[(start_end[2]+1):length(readme_content)])
}

# Create the replacement string

#replacement_string <- "# <img src=\"C:/Users/phili/OneDrive/PhD/R_packages_2023/expertsurv/inst/figures/hexsticker.png\" align=\"left\" height=\"150\"/> expertsurv"

string_rep <- "# <img src=\"inst/figures/hexsticker.png\" align=\"left\" height=\"150\"/> expertsurv"
line_num <- grep(string_rep,
                   readme_content)-1
readme_content[line_num]  <- "\n<br clear='all'/>\n"

replacement_string <- "# Overview"

readme_content <- gsub(string_rep,
                        replacement_string,
                          readme_content)


grepl("knitr::include_graphics(paste0(\"inst/image/\",img_temp))",
      #  "knitr::include_graphics(system.file(paste0(\"image/\",img_temp), package = \"expertsurv\"))",
                       readme_content, fixed = TRUE)
readme_content <- gsub("knitr::include_graphics(paste0(\"inst/image/\",img_temp))",
                       "knitr::include_graphics(system.file(paste0(\"image/\",img_temp), package = \"expertsurv\"))", fixed= TRUE,
                       readme_content)




# Write the modified content back to the file
writeLines(readme_content, vignette_path)

}
# You need to manually add the following lines
# 
# title: "<img src=\"`r system.file('figures/hexsticker.png', package = 'expertsurv')`\" align='left' height='150'/> Introduction to Expertsurv"
# output: rmarkdown::html_vignette
# bibliography: "`r system.file('REFERENCES.bib', package = 'expertsurv')`"
# vignette: >
#   %\VignetteIndexEntry{Introduction to Expertsurv}
# %\VignetteEngine{knitr::rmarkdown}
# %\VignetteEncoding{UTF-8}


# Render the vignette
rmarkdown::render(vignette_path, output_format = "html_vignette")
