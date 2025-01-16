#!/user/bin/Rscript
options(repos <- structure(c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
library(optparse)
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = FALSE,
              help = "You input file "),
  make_option(c("-o", "--output"), type = "character", default = FALSE,
              help = "Your outputDir "),
  make_option(c("-p", "--prefix"), type = "character", default = FALSE,
              help = "You species in [human,mouse] ")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 1.library package
library(magrittr, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(legendBaseModel, quietly = TRUE)
# vignette("dplyr" ,package = "dplyr")
# 2.glob function
source(paste0(getwdR(2), "/fun/unique_fun.R"), encoding = "UTF-8")
# 3.global variables
input_path  <-  clean_path(opt$input)
output_path <- clean_path(opt$output)
prefix <- opt$prefix

input_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/01_inputs/gpx4.xlsx"
output_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/02_outputs/gpx4_M6A"

input_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/01_inputs/mapk2.xlsx"
output_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/02_outputs/mapk2_M6A"

ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))

zero_data <- readxl::read_xlsx(input_path, sheet = 3)

raw_data <- zero_data %>%
  dplyr::select(-c("引物")) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(sample = paste(sample, seq_along(sample), sep = "-")) %>%
  dplyr::ungroup() %>%
  tibble::column_to_rownames("sample") %>%
  t() %>% 
  as.data.frame() %>%
  dplyr::mutate(x = seq_len(nrow(.)))
raw_data$x <- factor(raw_data$x, levels = raw_data$x)

data <- colnames(raw_data)

nc_data <- data[grepl("^NC-", data)]  # 以 "NC-" 开头的元素
si_data <- data[grepl("^si-", data)]  # 以 "si-" 开头的元素

combinations <- expand.grid(nc_data, si_data) %>% as.data.frame()

lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 20, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))


for (i in seq_len(nrow(combinations))) {
  nc_value <- combinations[i, 1] %>% as.character()
  si_value <- combinations[i, 2] %>% as.character()
  
  prefix <- paste0("M6A_", nc_value, "_", si_value)
  data_subset <- raw_data %>%
    dplyr::select(c("x", nc_value, si_value)) %>%
    tidyr::pivot_longer(cols = c(nc_value, si_value), 
               names_to = "group", 
               values_to = "value")

  p1 <- ggplot2::ggplot(data_subset, aes(x, value, group = group, color = group)) +
    geom_line() +
    labs(x="", y= "RT") 
  ggplot2::ggsave(paste0(output_path, "/", prefix, ".png"), 
                p1, width=16, height=8)
}


