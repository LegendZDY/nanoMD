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
input_  <-  path_clean("./1-output/2022_08_check_m6A/input/")
species <- choose_species("human")
output_ <- path_clean("./1-output/2022_08_check_m6A/all_m6a_check_results2")

run_data <- data.frame(
    group=c("m6a_check"), 
    input = input_,
    output = output_) %>% 
    mutate(
        CvsT_matrix = map2(input, output, ~ merages_data_run(.x, .y, "_m6A.bed")),
        CvsT_name = map(group, ~ paste0(output_, "/", .x, "CVST")),
        Q1_DE = map2(CvsT_matrix, CvsT_name, DE_analysis_run))
## run full_join results
all_data_matrix <- all_data %>%
    dplyr::select(id, ends_with("counts")) %>% 
    dplyr::rename(setNames(colnames(.), gsub("_counts", "", colnames(.)))) %>% 
    dplyr::rename(Geneid = id) %>% 
    dplyr::relocate(., starts_with("control"), .after = Geneid)

all_data_matrix[is.na(all_data_matrix)] <- 0
write.table(all_data_matrix, paste0(output_, "/all_count_matrix.tsv"), sep = "\t", row.names = F, quote = F)

all_sample %>%
    dplyr::select(Geneid, ends_with("counts")) %>% 
    dplyr::rename(setNames(colnames(.), gsub("_counts", "", colnames(.)))) %>% 
    dplyr::relocate(., starts_with("control"), .after = Geneid) %>%
    write.table(., paste0(output_path, "/all_count_matrix.tsv"), sep = "\t", row.names = F, quote = F)
run_data2 <- run_data %>% 
    mutate(Q1_DE = map2(CvsT_matrix, CvsT_name, DE_analysis_run)) 
#run merage_GPX4 data
merage_GPX4_run_data <- data.frame(
    group=c("m6a_check"), 
    input = input_,
    output = output_) %>% 
    mutate(
        CvsT_matrix = paste0(output_, "/all_count_matrix.tsv"),
        CvsT_name = map(group, ~ paste0(output_, "/", .x, "CVST")),
        Q1_DE = map2(CvsT_matrix, CvsT_name, DE_analysis_run))

bookdown::render_book(input = ".", output_format = NULL, ..., clean = TRUE,
  envir = parent.frame(), clean_envir = !interactive(),
  output_dir = NULL, new_session = NA, preview = FALSE,
  config_file = "_bookdown.yml")

bookdown::render_book("index.Rmd", "bookdown::pdf_book")