require(rmarkdown)
render(file.path(getwd(), "Volcanic50.Rmd"), output_format="html_document")
file.rename(file.path(getwd(), "Volcanic50.md"), file.path(getwd(), "README.md"))
