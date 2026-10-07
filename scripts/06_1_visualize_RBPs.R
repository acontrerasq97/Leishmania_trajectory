library(httr2)
library(jsonlite)
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)

protist_mart <- useEnsemblGenomes(
    biomart = "protists_mart", 
    dataset = "lmajor_eg_gene"
)

genes <- c("LmjF.04.1170", "LmjF.07.0340")


req <- request("https://ensembl.org") %>%
  req_headers(
    "Content-Type" = "application/json", "Accept" = "application/json") %>%
  req_body_json(list(ids = genes)) %>% 
  req_perform()

annotations <- as.data.frame(t(simplify2array(resp_body_json(res))))

