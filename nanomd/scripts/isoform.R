#!/opt/conda/bin/Rscript
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))  
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
library(optparse)
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = FALSE,
                help = "input file path"),
    make_option(c("-o", "--output"), type = "character", default = FALSE,
              help = "output file path"),
    make_option(c("-p", "--prefix"), type = "character", default = FALSE,
              help = "prefix of output file"),
    make_option(c("-s", "--species"), type = "character", default = FALSE,
              help = "species name (e.g. human, mouse, rat etc.)"),
    make_option(c("-t", "--type"), type = "character", default = FALSE,
              help = "file type (e.g. salmon, polyA, etc.)"),
    make_option(c("-n", "--split_num"), type = "numeric", default = NA,
              help = "split num is control the number of split for the plot"),
    make_option(c("-f", "--filter_num"), type = "numeric", default = 1,
              help = "filter num is control the number of filter for the plot"),
    make_option(c("-m", "--model"), type = "character", default = "normal",
              help = "model is control the model for the plot, you can choose in [normal|filter]")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input  <-  opt$input
output <- opt$output
prefix <- opt$prefix
species <- opt$species
type <- opt$type
split_num <- as.integer(opt$split_num)
filter_num <- as.integer(opt$filter_num)
model <- opt$model

# library packages
library(dplyr)
library(purrr)
library(ggplot2)
library(magrittr)
library(baseModule)
library(legendBaseModel)

species <- baseModule::choose_species(species)

matrix_init <- function(input_file, output_path, type) {
    # test
    # input_file <- "./python/nanoMD/nanomd/scripts/matrix_count.tsv"
    # input_file <- "./python/nanoMD/nanomd/scripts/matrix_polyA.tsv"
    # output_path <- "./python/nanoMD/nanomd/scripts"
    ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))
    tryCatch({
    # START
    file_name <- basename(input_file) %>% tools::file_path_sans_ext()
    cat(crayon::green(paste0(file_name, " Start matrix initialization\n")))
    # read data
    matrix_row <- readr::read_delim(input_file, col_names = T, delim = "\t", comment = "#", show_col_types = FALSE)  %>% 
        dplyr::rename(isoform = 1) %>% 
        dplyr::mutate(across(where(is.numeric), floor)) %>% 
        dplyr::rename(setNames(colnames(.), gsub("_quant$|_polyA$", "", colnames(.)))) %>% 
        tidyr::separate(isoform, c("transcript_id", "V2", "V3", "V4", "V5", "SYMBOL", "V7", "V8", "V9"), sep = "[|]") %>% 
        dplyr::select(-c("V2", "V3", "V4", "V5", "V7", "V8", "V9"))
    
    matrix_isoform <- matrix_row %>% 
        tidyr::unite(isoform, c("transcript_id", "SYMBOL"), sep = "|")

    write.csv(matrix_isoform, file = paste0(output_path, "/", file_name, "_isoform_standard.csv"), quote = FALSE, row.names =F)
    write.table(matrix_isoform, file = paste0(output_path, "/", file_name, "_isoform_standard.tsv"), sep = "\t", quote = FALSE, row.names =F)
    write.table(matrix_isoform, file = paste0(output_path, "/", file_name, "_isoform_standard.xls"), sep = "\t", quote = FALSE, row.names =F)
    if (type == "salmon") {
        matrix_gene <- matrix_row %>% 
            dplyr::select(-transcript_id) %>% 
            dplyr::group_by(SYMBOL) %>%
            dplyr::summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE)))
        write.csv(matrix_gene, file = paste0(output_path, "/", file_name, "_gene_standard.csv"), quote = FALSE, row.names =F)
        write.table(matrix_gene, file = paste0(output_path, "/", file_name, "_gene_standard.tsv"), sep = "\t", quote = FALSE, row.names =F)
        write.table(matrix_gene, file = paste0(output_path, "/", file_name, "_gene_standard.xls"), sep = "\t", quote = FALSE, row.names =F)
    }
    cat(crayon::green(paste0(file_name, " Matrix initialization done\n")))
    # END
    return("sucess")},
    error = function(e){print(e$message);message(return("failled"))})
}

auto_plot <- function(input_file, output_path, prefix, species, type = "salmon", split_num = NA, filter_num = 1, model = "normal") {
    # test
    # input_file <- "./python/nanoMD/nanomd/scripts/matrix_count.tsv"
    # input_file <- "./python/nanoMD/nanomd/scripts/matrix_polyA.tsv"
    # output_path <- "./python/nanoMD/nanomd/scripts/output"
    # prefix <- "test"
    # species <- "mouse"
    # type <- "salmon"
    # split_num <- NA
    # filter_num <- 1
    # model <- "normal"
    # 2. check output path
    options(warn = -1)
    ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))
    tryCatch({
    # START
    # init matrix
    file_name <- basename(input_file) %>% tools::file_path_sans_ext()
    file_path <- dirname(input_file)
    cat(crayon::green(paste0(file_name, " Start auto plot\n")))
    matrix_init(input_file, output_path, type)
    # read isoform matrix
    isoform_path  <- paste0(output_path, "/", file_name, "_isoform_standard.tsv")
    isoform_output <- paste0(output_path, "/", "isoform")
    ifelse(dir.exists(isoform_output), "Dir exist alreadly!", dir.create(isoform_output, recursive = TRUE))
    isoform_matrix <- readr::read_delim(isoform_path, col_names = T, delim = "\t", comment = "#", show_col_types = FALSE) %>% 
        dplyr::distinct(isoform, .keep_all = TRUE) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var="isoform")
    legendBaseModel::deseq2_Rd(isoform_matrix, isoform_output, prefix, split_num, filter_num, model)
    # read gene matrix
    if (type == "salmon") {
        gene_path  <- paste0(output_path, "/", file_name, "_gene_standard.tsv")
        gene_output <- paste0(output_path, "/", "gene")
        ifelse(dir.exists(gene_output), "Dir exist alreadly!", dir.create(gene_output, recursive = TRUE))
        gene_matrix <- readr::read_delim(gene_path, col_names = T, delim = "\t", comment = "#", show_col_types = FALSE) %>% 
            dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
            tibble::remove_rownames() %>%
            tibble::column_to_rownames(var="SYMBOL")
        legendBaseModel::deseq2_Rd(gene_matrix, gene_output, prefix, split_num, filter_num, model)
        deseq2_file <- paste0(gene_output, "/", prefix, "_deseq2.tsv")
        enrichment_output <- paste0(gene_output, "/", "enrichment")
        ifelse(dir.exists(enrichment_output), "Dir exist alreadly!", dir.create(enrichment_output, recursive = TRUE))
        legendBaseModel::auto_enrichment(deseq2_file, enrichment_output, species)
    }
    # END
    return("sucess")},
    error = function(e){print(e$message);message(return("failled"))})
}


# main run
auto_plot(input, output, prefix, species, type, split_num, filter_num, model)