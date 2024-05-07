#!/user/bin/Rscript
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))  
if (!requireNamespace("optparse", quietly = TRUE))
  # install.packages("optparse",repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  install.packages("optparse")
library(optparse)
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = FALSE,
              help = "You input file "),
  make_option(c("-s", "--species"), type = "character", default = FALSE,
              help = "You species in [human,mouse] "), 
  make_option(c("-o", "--outputDir"), type = "character", default = FALSE,
              help = "Your outputDir ")
)
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# 1.library package
library(modelr)
library(magrittr)
library(tidyverse)
library(rstudioapi)
options(na.action = na.warn)
# vignette("dplyr" ,package = "dplyr")
# 2.glob function
# 1.add analysis function
#get the file path
fun_path_ <- dirname(getSourceEditorContext()$path)
source(paste0(fun_path_, "/fun/basis_fun.R"), encoding = "UTF-8")
source(paste0(fun_path_, "/fun/unique_fun.R"), encoding = "UTF-8")

# 2.global variables

# input_  <-  path_clean(opt$input)
# species <- choose_species(opt$species)
# output_ <- path_clean(opt$outputDir)

## run chenyiping's m6a DEA
input_  <-  path_clean("./1-output/20221027_newRNA/input/")
species <- choose_species("human")
output_ <- path_clean("./1-output/20221027_newRNA/all_newRNA_results3")

run_data <- data.frame(
    group=c("newRNA"), 
    input = input_,
    output = output_) %>% 
    mutate(
        ipath = paste0(input, "/newRNA_counts_matrix.tsv"),
        newRNA_DE = pmap(list(ipath, output, group), newRNA_DE_Rd))

run_data <- data.frame(
    group=c("allRNA"), 
    input = input_,
    output = output_) %>% 
    mutate(
        ipath = paste0(input, "/newRNA_counts_matrix.tsv"),
        newRNA_DE = pmap(list(ipath, output, group), newRNA_DE_Rd))

run_data <- data.frame(
    group=c("allRNA"), 
    input = input_,
    output = output_) %>% 
    mutate(
        ipath = paste0(input, "/newRNA_counts_matrix.tsv"),
        newRNA_DE = pmap(list(ipath, output, group), newRNA_DE_Rd))