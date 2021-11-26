# xx_script_mapper.R

library(DiagrammeR)
library(rsvg)
library(DiagrammeRsvg)
library(tidyverse)


create_codemap_files <- function(){
  
  scripts<- 
    suppressWarnings(list.files("scripts/")[grep(".*[0-9].*",list.files("scripts/"))] %>% 
    str_remove(".R") %>% 
    as.data.frame() %>% 
    `colnames<-`("name") %>% 
    mutate(fullname = name) %>% 
    separate(name ,into= c("number", "first"), sep="_") %>% 
    mutate(rank = as.integer(str_sub(number, end=2)),
           label = paste0(number,"_",first,"..."),
           numberspare = number) %>% 
    separate(numberspare, into=c("main","sub")  ) %>% 
    mutate(sub=ifelse(is.na(sub),0,sub)) %>% 
    arrange(main,sub) %>% 
    rowid_to_column() %>% 
    group_by(rank) %>% 
    arrange(number) %>% 
    mutate(id = row_number(),
           to = case_when(rank==1&id==1~0,  
                          id>1~1,
                          TRUE~1),
           from = case_when(id == 1 ~ n(),
                            TRUE~as.integer(0))) %>% 
    ungroup() %>% 
    mutate(from =ifelse(rank == max(rank)&id==1 , from-1, from),
           arrow = case_when(rank==1&id==1~as.character(NA),
                             id>1 ~ "crow",
                             TRUE ~ "normal")))
  
  froms <- rep(scripts$fullname,scripts$from)
  tos <- rep(scripts$fullname,scripts$to)
  
  
  # Create a simple NDF
  nodes_1 <-
    create_node_df(
      n = nrow(scripts), 
      label = scripts$label,
      fontcolor = "black",
      rank = scripts$rank,
      style = "filled",
      fillcolor = "grey",
      color="grey",
      shape = "circle")
  
  edges_1 <-
    create_edge_df(
      from = tibble(fullname= froms) %>% left_join(scripts) %>% .$rowid,
      to = tibble(fullname= tos) %>% left_join(scripts) %>% .$rowid,
      color = "grey",
      arrowhead = scripts$arrow %>% na.omit())
  
  #silently add a white first level node causing the first, and last level node
  #caused by the last as this is the only way I have found to create whitespace
  
  #NB this doesn't really work if the network has different numbers of nodes
  #on the different levels, and can needlessly shift the node alignment in these
  #cases, so probably come up with an alternative or drop it
  
  nodes_1 <- nodes_1 %>% 
    bind_rows(tibble(id=c(nrow(nodes_1)+1,nrow(nodes_1)+2),
                     type=rep(NA,2),
                     label=rep(NA,2),
                     fontcolor=rep("white",2),
                     rank=c(1,max(scripts$rank)),
                     style=rep("filled",2),
                     fillcolor=rep("white",2),
                     color=rep("white",2),
                     shape=rep("plain",2)))
  
  edges_1 <- edges_1 %>% 
    bind_rows(tibble(id=c(nrow(nodes_1)-1,nrow(nodes_1)),
                     from=c(nrow(nodes_1),nrow(nodes_1)-2),
                     to=c(1,nrow(nodes_1)),
                     color=rep("white",2),
                     arrowhead=rep("normal",2)))
  
  # Render code overview png
  
  graph <-
    create_graph(
      nodes_df = nodes_1,
      edges_df = edges_1,
      attr_theme = "tb"
    ) %>%
    set_node_attrs(
      node_attr = "fontname",
      values = "Helvetica"
    ) %>% 
    set_node_attrs(
      node_attr = "fontsize",
      values = 8
    ) %>%
    set_edge_attrs(
      edge_attr = "arrowsize",
      values = 1) %>%
    set_edge_attrs(
      edge_attr= "dir",
      values = "forward"
    )
  dir.create(file.path("scripts", "reports"), showWarnings = FALSE)
  
  export_graph(graph,
               file_name = paste0("./scripts/reports/overview_codemap.png"),
               file_type = "png",
               height = 1000
  )
  
  for (primary in scripts$fullname){
    
    pr_lab <- scripts %>% 
      filter(fullname==primary) %>% 
      .$label
    
    nodes_temp <- nodes_1 %>% 
      mutate(fillcolor=ifelse(label==pr_lab, "firebrick2", fillcolor),
             color=ifelse(label==pr_lab, "firebrick2", color),
             label=ifelse(label==pr_lab, label, ""))
    
    graph <-
      create_graph(
        nodes_df = nodes_temp,
        edges_df = edges_1,
        attr_theme = "tb"
      ) %>%
      set_node_attrs(
        node_attr = "fontname",
        values = "Helvetica"
      ) %>% 
      set_node_attrs(
        node_attr = "fontsize",
        values = 8
      ) %>%
      set_edge_attrs(
        edge_attr = "arrowsize",
        values = 1) %>%
      set_edge_attrs(
        edge_attr= "dir",
        values = "forward"
      )
    dir.create(file.path("scripts", "reports"), showWarnings = FALSE)
    
    export_graph(graph,
                 file_name = paste0("./scripts/reports/",primary, "_codemap.png"),
                 file_type = "png",
                 height = 1000
    )
  }
}

## knit into single "gitbook"

generate_code_walkthrough <- function(name = "guess",
                                      overview_text="",
                                      format = "gitbook"){
  
  if(name=="guess"){
   name= rprojroot::find_rstudio_root_file() %>% 
      stringr::str_split("/") %>% 
      .[[1]] %>% 
      dplyr::last()
  }
  
  # we'll work from the basis of numeric scripts as above for consistency
  sections<- list.files("scripts/")[grep(".*[0-9].*",list.files("scripts/"))] %>% 
    str_remove_all(".R") 
  
  sections <- suppressWarnings(tibble(
    ref= sections,
    name=sections) %>% 
    separate(name ,into= c("number"), sep="_")  %>% 
    separate(number, into=c("main","sub")  ) %>% 
    mutate(sub=ifelse(is.na(sub),0,sub)) %>% 
    arrange(main,sub) %>% 
    .$ref)

  ## Add images to scripts ahead of spinning
  
  for (script in sections){
    fileconnR<- file(paste0("scripts/",script,".R"), open="r")
    tmp <-readLines(fileconnR, warn = F)
    close(fileconnR)
    
    if(any(grepl("knitr::opts", tmp)) ==F){
      tmp <- c(paste0("# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'"),
               tmp)
    }
    if(any(grepl("_codemap", tmp)) ==F){
      tmp <- c(tmp[1:min(grep("#'# ",tmp))],
               paste0("#'
# ---- echo=F, eval=T, out.width='33%', out.extra = 'style=\"float:right;\"'
knitr::include_graphics(\"",script,"_codemap.png\", error = FALSE)
#'
" ),
               tmp[min(grep("#'# ",tmp))+1:length(tmp)])
    }
    fileconnWT<- file(paste0("scripts/",script,".R"), open="wt")
    writeLines(tmp,fileconnWT)
    close(fileconnWT)
    rm(tmp)
  }
  
  ## Spin all scripts
  for(script in sections){
    ezknitr::ezspin(paste0("scripts/",script,".R"), 
                    out_dir="scripts/reports",
                    keep_rmd = T, 
                    verbose=T)
  }
  
  
  #First make the index file if it doesn't exist - this will be user customisable
  
  if(!file.exists("./scripts/reports/index.rmd"))
  {
    
    fileConn<-file("./scripts/reports/index.rmd")
    writeLines(paste0(c(
      "--- 
title: \"",name,"\"
date: \"`r Sys.Date()`\"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
description: \"This is a code walkthrough for the project", name,".\"
---

```{r include=FALSE}
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

```
# Overview of project
&nbsp;
",overview_text,"&nbsp;

```{r echo=F, eval=T, out.width='33%', fig.align='center'}
knitr::include_graphics(\"overview_codemap.png\", error = FALSE)
```
")), fileConn)
    
    close(fileConn)
  }
  
  if(!file.exists("./scripts/reports/_bookdown.yml"))
  {
    fileConn<-file("./scripts/reports/_bookdown.yml")
    writeLines(paste0(c(
  "rmd_files: [\"index.Rmd\", ",paste0('"',sections,'.Rmd"', collapse=","),"]" ), collapse=""), fileConn)

    close(fileConn)
  }
  
  if(format=="pdf"){
  bookdown::render_book("./scripts/reports/", 
                        output_format = "bookdown::pdf_book", 
                        preview=FALSE, 
                        output_dir = "./combined/" )
  }else if(format=="gitbook"){
    bookdown::render_book("./scripts/reports/", 
                          output_format = "bookdown::gitbook", 
                          preview=FALSE, 
                          output_dir = "./combined/" )
  }
}

