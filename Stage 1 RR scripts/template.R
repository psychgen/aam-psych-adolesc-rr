# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
#'
#'# Your script title (`your_script_filename.R`) 
#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to .... 
#'
#'<!-- The script contains additional formatting to allow direct rendering to a 
#'code report via rmarkdown (specifically using knitr::spin). You can ignore it 
#'all, and just use it as an Rscript, add light-touch RMarkdown formatting like
#'highlighting code chunks, or go all-in and write a formatted walk-through of 
#'your code. 
#'
#'Key pointers:
#' - Text following the " #' " sign will be normal text in the report (and will 
#'simply be ignored by R)
#' - R code chunks are initiated by the "# ---- " combination, can 
#'include regular (#) comments to be included with the code chunk, and can have
#'options specified as per RMarkdown standards. They end when " #' " text resumes.
#' - Check out the RMarkdown cheatsheet for standard formatting options 
#' 
#'  (NB this text will be ignored by both R and knitr, so you can leave
#'  for reference or delete as you prefer) --> 
#'
#'&nbsp;
#'&nbsp;
#'
#'## Example section header
#'
#'### Example subsection header
#'
#'  * Example bullet list
#'
#'
#'## This is the first R code chunk
# ----  echo=F
head(cars)
# This is a 'true' R comment, to be included in the chunk
# The chunk will end with the next " #' " line...
summary(cars)
#'







# ---- eval=F,include=F
#This is the code to render and knit the script - called automatically by 
#generate_code_walkthrough()
ezknitr::ezspin("scripts/template.R", 
                out_dir="scripts/reports",
                keep_rmd = T, 
                verbose=T)
