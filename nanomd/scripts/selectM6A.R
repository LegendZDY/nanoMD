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
library(cowplot)
# install.packages("magick")
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
output_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/02_outputs/mapk2_M6A_483-67C"

ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))

zero_data <- readxl::read_xlsx(input_path, sheet = 3) %>% 
  dplyr::filter(Assay=="MAP2-449-60C")

zero_data <- readxl::read_xlsx(input_path, sheet = 3) %>% 
  dplyr::filter(Assay=="MAP2-483-67C")

# 2022
input_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/01_inputs/20220613.xlsx"
output_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/02_outputs/gpx4_2022"

ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))

zero_data <- readxl::read_xlsx(input_path, sheet = 3)

zero_data <- readxl::read_xlsx(input_path, sheet = 3) %>% 
  dplyr::filter(Assay=="GPX4-690") %>% 
  dplyr::rename(sample = Sample)


raw_data <- zero_data %>%
  dplyr::select(-c("Assay")) %>%
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
si_data <- data[grepl("^si", data)]  # 以 "si-" 开头的元素

combinations <- expand.grid(nc_data, si_data) %>% as.data.frame()

lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 40, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))


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

# paper
## GPX4-650
input_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/01_inputs/20220613.xlsx"

output_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/02_outputs/paper"

ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))

zero_data <- readxl::read_xlsx(input_path, sheet = 3) %>% 
  dplyr::filter(Assay=="GPX4-690") %>% 
  dplyr::rename(sample = Sample)


raw_data <- zero_data %>%
  dplyr::select(-c("Assay")) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(sample = paste(sample, seq_along(sample), sep = "-")) %>%
  dplyr::ungroup() %>%
  tibble::column_to_rownames("sample") %>%
  t() %>% 
  as.data.frame() %>%
  dplyr::mutate(x = seq_len(nrow(.)))
raw_data$x <- factor(raw_data$x, levels = raw_data$x)

nc_value <- "siNC-1-1"
si_value <- "siMettl3-1-1"

prefix <- "GPX4"


sub <- raw_data %>%
  dplyr::select(all_of(c("x", nc_value, si_value)))

data_subset <- sub %>%
  tidyr::pivot_longer(cols = c(nc_value, si_value), 
              names_to = "group", 
              values_to = "value") %>% 
  dplyr::mutate(group2 = case_when(
              grepl("^siNC", group) ~ "NC",
              TRUE ~ "siMettl3"))

bar_data <- as.data.frame(
  group <- c("NC", "siMettl3")) %>%
  dplyr::mutate(values = case_when(
    group == "NC" ~ sub[which(sub[, 2] > 0.1)[1], "x"],
    group == "siMettl3" ~ sub[which(sub[, 3] > 0.1)[1], "x"]
  ))

color_set <- c("siMettl3"="#3e7aaf", "NC"="#c37f54")

p <- ggplot2::ggplot(data_subset, aes(x, value, group = group2, color = group2)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = color_set) +  # 使用 scale_color_manual
  labs(x = "GPX4 650 site", y = "RT") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # 去掉背景网格线
    legend.title = element_blank(),  # 去掉图例标题
    legend.position = c(0.98, 0.1),  # 设置图例位置（右下角）
    legend.justification = c(1, 0),
    axis.title.x = element_text(size = 20, face = "bold"),  # 横坐标加大且加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵坐标加大且加粗
    axis.text.x = element_text(size = 12, face = "bold"),  # 横坐标刻度加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标刻度加粗
    legend.text = element_text(size = 12, face = "bold"),  # 图例文字加粗
    axis.line = element_line(size = 1),  # 坐标轴线条加粗
    panel.border = element_blank()
  ) 

bar <- ggplot(bar_data, aes(fill=group, y=values, x=group)) + 
    geom_bar(position = position_stack(), stat="identity", width = .7) +
    scale_fill_manual(values=color_set) +
    geom_text(aes(label = ifelse(group == "siMettl3", "***", ""),
                vjust = -0.5),
            size = 6,
            color = "black") +
    guides(fill = guide_legend(title = NULL)) +
    xlab("") +
    ylab(bquote(bold(qPCR) ~ C[T])) +
    theme_bw() +
    theme(
    panel.grid = element_blank(),  # 去掉背景网格线
    legend.title = element_blank(),  # 去掉图例标题
    axis.title.x = element_text(size = 20, face = "bold"),  # 横坐标加大且加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵坐标加大且加粗
    axis.text.x = element_text(size = 20, face = "bold"),  # 横坐标刻度加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标刻度加粗
    legend.text = element_text(size = 12, face = "bold"),  # 图例文字加粗
    axis.line = element_line(size = 1),  # 坐标轴线条加粗
    panel.border = element_blank()
    )

p2 <- p + draw_plot(bar, x = 0, y = 0.02, hjust = 0.25, vjust = 0, width = 50, height = 0.9, scale = 0.4)
p2
ggplot2::ggsave(paste0(output_path, "/", prefix, "-650.png"), 
              p2, width=12, height=8)
ggplot2::ggsave(paste0(output_path, "/", prefix, "-650.pdf"), 
              p2, width=12, height=8)
ggplot2::ggsave(paste0(output_path, "/", prefix, "-650.tiff"), 
              p2, width=12, height=8)

# p1 <- plot_grid(p, bar, labels = c("A", "B"), rel_widths = c(2, 1), label_size = 12)

# ggplot2::ggsave(paste0(output_path, "/", prefix, ".png"), 
#               p1, width=16, height=8)
input_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/01_inputs/mapk2.xlsx"

output_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/02_outputs/paper"

ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))

zero_data <- readxl::read_xlsx(input_path, sheet = 3)


raw_data <- zero_data %>%
  dplyr::select(-c("Assay")) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(sample = paste(sample, seq_along(sample), sep = "-")) %>%
  dplyr::ungroup() %>%
  tibble::column_to_rownames("sample") %>%
  t() %>% 
  as.data.frame() %>%
  dplyr::mutate(x = seq_len(nrow(.)))
raw_data$x <- factor(raw_data$x, levels = raw_data$x)

# MAPKAPK2
nc_value <- "NC-1-2"
si_value <- "si-mettl3-1-4"

prefix <- "MAPKAPK2"


sub <- raw_data %>%
  dplyr::select(all_of(c("x", nc_value, si_value)))

data_subset <- sub %>%
  tidyr::pivot_longer(cols = c(nc_value, si_value), 
              names_to = "group", 
              values_to = "value") %>% 
  dplyr::mutate(group2 = case_when(
              grepl("^NC", group) ~ "NC",
              TRUE ~ "siMettl3"))

bar_data <- as.data.frame(
  group <- c("NC", "siMettl3")) %>%
  dplyr::mutate(values = case_when(
    group == "NC" ~ sub[which(sub[, 2] > 0.2)[1], "x"],
    group == "siMettl3" ~ sub[which(sub[, 3] > 0.2)[1], "x"]
  ))

color_set <- c("siMettl3"="#3e7aaf", "NC"="#c37f54")

p <- ggplot2::ggplot(data_subset, aes(x, value, group = group2, color = group2)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = color_set) +  # 使用 scale_color_manual
  labs(x = "MAPKAPK2 483 site", y = "RT") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # 去掉背景网格线
    legend.title = element_blank(),  # 去掉图例标题
    legend.position = c(0.98, 0.1),  # 设置图例位置（右下角）
    legend.justification = c(1, 0),
    axis.title.x = element_text(size = 20, face = "bold"),  # 横坐标加大且加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵坐标加大且加粗
    axis.text.x = element_text(size = 12, face = "bold"),  # 横坐标刻度加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标刻度加粗
    legend.text = element_text(size = 12, face = "bold"),  # 图例文字加粗
    axis.line = element_line(size = 1),  # 坐标轴线条加粗
    panel.border = element_blank()
  ) 

bar <- ggplot(bar_data, aes(fill=group, y=values, x=group)) + 
    geom_bar(position = position_stack(), stat="identity", width = .7) +
    scale_fill_manual(values=color_set) +
    geom_text(aes(label = ifelse(group == "siMettl3", "***", ""),
                vjust = -0.5),
            size = 6,
            color = "black") +
    guides(fill = guide_legend(title = NULL)) +
    xlab("") +
    ylab(bquote(bold(qPCR) ~ C[T])) +
    theme_bw() +
    theme(
    panel.grid = element_blank(),  # 去掉背景网格线
    legend.title = element_blank(),  # 去掉图例标题
    axis.title.x = element_text(size = 20, face = "bold"),  # 横坐标加大且加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵坐标加大且加粗
    axis.text.x = element_text(size = 20, face = "bold"),  # 横坐标刻度加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标刻度加粗
    legend.text = element_text(size = 12, face = "bold"),  # 图例文字加粗
    axis.line = element_line(size = 1),  # 坐标轴线条加粗
    panel.border = element_blank()
    )

p2 <- p + draw_plot(bar, x = 0, y = 0.06, hjust = 0.27, vjust = 0, width = 50, height = 0.8, scale = 0.4)
p2
ggplot2::ggsave(paste0(output_path, "/", prefix, "-483.png"), 
              p2, width=12, height=8)
ggplot2::ggsave(paste0(output_path, "/", prefix, "-483.pdf"), 
              p2, width=12, height=8)
ggplot2::ggsave(paste0(output_path, "/", prefix, "-483.tiff"), 
              p2, width=12, height=8)

# MAPKAPK2-449
input_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/01_inputs/mapk2.xlsx"

output_path <- "/media/epibiotek/zdy/projects/n3_nanomd/outputs/00_OUT_papers/02_outputs/paper"

ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))

zero_data <- readxl::read_xlsx(input_path, sheet = 3)


raw_data <- zero_data %>%
  dplyr::select(-c("Assay")) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(sample = paste(sample, seq_along(sample), sep = "-")) %>%
  dplyr::ungroup() %>%
  tibble::column_to_rownames("sample") %>%
  t() %>% 
  as.data.frame() %>%
  dplyr::mutate(x = seq_len(nrow(.)))
raw_data$x <- factor(raw_data$x, levels = raw_data$x)

# MAPKAPK2
nc_value <- "NC-1-2"
si_value <- "si-mettl3-1-4"

prefix <- "MAPKAPK2"


sub <- raw_data %>%
  dplyr::select(all_of(c("x", nc_value, si_value)))

data_subset <- sub %>%
  tidyr::pivot_longer(cols = c(nc_value, si_value), 
              names_to = "group", 
              values_to = "value") %>% 
  dplyr::mutate(group2 = case_when(
              grepl("^NC", group) ~ "NC",
              TRUE ~ "siMettl3"))

bar_data <- as.data.frame(
  group <- c("NC", "siMettl3")) %>%
  dplyr::mutate(values = case_when(
    group == "NC" ~ sub[which(sub[, 2] > 0.2)[1], "x"],
    group == "siMettl3" ~ sub[which(sub[, 3] > 0.2)[1], "x"]
  ))

color_set <- c("siMettl3"="#3e7aaf", "NC"="#c37f54")

p <- ggplot2::ggplot(data_subset, aes(x, value, group = group2, color = group2)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = color_set) +  # 使用 scale_color_manual
  labs(x = "MAPKAPK2 483 site", y = "RT") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # 去掉背景网格线
    legend.title = element_blank(),  # 去掉图例标题
    legend.position = c(0.98, 0.1),  # 设置图例位置（右下角）
    legend.justification = c(1, 0),
    axis.title.x = element_text(size = 20, face = "bold"),  # 横坐标加大且加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵坐标加大且加粗
    axis.text.x = element_text(size = 12, face = "bold"),  # 横坐标刻度加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标刻度加粗
    legend.text = element_text(size = 12, face = "bold"),  # 图例文字加粗
    axis.line = element_line(size = 1),  # 坐标轴线条加粗
    panel.border = element_blank()
  ) 

bar <- ggplot(bar_data, aes(fill=group, y=values, x=group)) + 
    geom_bar(position = position_stack(), stat="identity", width = .7) +
    scale_fill_manual(values=color_set) +
    geom_text(aes(label = ifelse(group == "siMettl3", "***", ""),
                vjust = -0.5),
            size = 6,
            color = "black") +
    guides(fill = guide_legend(title = NULL)) +
    xlab("") +
    ylab(bquote(bold(qPCR) ~ C[T])) +
    theme_bw() +
    theme(
    panel.grid = element_blank(),  # 去掉背景网格线
    legend.title = element_blank(),  # 去掉图例标题
    axis.title.x = element_text(size = 20, face = "bold"),  # 横坐标加大且加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵坐标加大且加粗
    axis.text.x = element_text(size = 20, face = "bold"),  # 横坐标刻度加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标刻度加粗
    legend.text = element_text(size = 12, face = "bold"),  # 图例文字加粗
    axis.line = element_line(size = 1),  # 坐标轴线条加粗
    panel.border = element_blank()
    )

p2 <- p + draw_plot(bar, x = 0, y = 0.06, hjust = 0.27, vjust = 0, width = 50, height = 0.8, scale = 0.4)
p2
ggplot2::ggsave(paste0(output_path, "/", prefix, "-483.png"), 
              p2, width=12, height=8)
ggplot2::ggsave(paste0(output_path, "/", prefix, "-483.pdf"), 
              p2, width=12, height=8)
ggplot2::ggsave(paste0(output_path, "/", prefix, "-483.tiff"), 
              p2, width=12, height=8)


# MAPKAPK2-449
nc_value <- "NC-3-3"
si_value <- "si-mettl3-1-6"

prefix <- "MAPKAPK2"


sub <- raw_data %>%
  dplyr::select(all_of(c("x", nc_value, si_value)))

data_subset <- sub %>%
  tidyr::pivot_longer(cols = c(nc_value, si_value), 
              names_to = "group", 
              values_to = "value") %>% 
  dplyr::mutate(group2 = case_when(
              grepl("^NC", group) ~ "NC",
              TRUE ~ "siMettl3"))

bar_data <- as.data.frame(
  group <- c("NC", "siMettl3")) %>%
  dplyr::mutate(values = case_when(
    group == "NC" ~ sub[which(sub[, 2] > 0.2)[1], "x"],
    group == "siMettl3" ~ sub[which(sub[, 3] > 0.2)[1], "x"]
  ))

color_set <- c("siMettl3"="#3e7aaf", "NC"="#c37f54")

p <- ggplot2::ggplot(data_subset, aes(x, value, group = group2, color = group2)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = color_set) +  # 使用 scale_color_manual
  labs(x = "MAPKAPK2 449 site", y = "RT") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # 去掉背景网格线
    legend.title = element_blank(),  # 去掉图例标题
    legend.position = c(0.98, 0.1),  # 设置图例位置（右下角）
    legend.justification = c(1, 0),
    axis.title.x = element_text(size = 20, face = "bold"),  # 横坐标加大且加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵坐标加大且加粗
    axis.text.x = element_text(size = 12, face = "bold"),  # 横坐标刻度加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标刻度加粗
    legend.text = element_text(size = 12, face = "bold"),  # 图例文字加粗
    axis.line = element_line(size = 1),  # 坐标轴线条加粗
    panel.border = element_blank()
  ) 

bar <- ggplot(bar_data, aes(fill=group, y=values, x=group)) + 
    geom_bar(position = position_stack(), stat="identity", width = .7) +
    scale_fill_manual(values=color_set) +
    geom_text(aes(label = ifelse(group == "siMettl3", "***", ""),
                vjust = -0.5),
            size = 6,
            color = "black") +
    guides(fill = guide_legend(title = NULL)) +
    xlab("") +
    ylab(bquote(bold(qPCR) ~ C[T])) +
    theme_bw() +
    theme(
    panel.grid = element_blank(),  # 去掉背景网格线
    legend.title = element_blank(),  # 去掉图例标题
    axis.title.x = element_text(size = 20, face = "bold"),  # 横坐标加大且加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵坐标加大且加粗
    axis.text.x = element_text(size = 20, face = "bold"),  # 横坐标刻度加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标刻度加粗
    legend.text = element_text(size = 12, face = "bold"),  # 图例文字加粗
    axis.line = element_line(size = 1),  # 坐标轴线条加粗
    panel.border = element_blank()
    )

p2 <- p + draw_plot(bar, x = 0, y = 0.06, hjust = 0.27, vjust = 0, width = 50, height = 0.74, scale = 0.4)
p2
ggplot2::ggsave(paste0(output_path, "/", prefix, "-449.png"), 
              p2, width=12, height=8)
ggplot2::ggsave(paste0(output_path, "/", prefix, "-449.pdf"), 
              p2, width=12, height=8)
ggplot2::ggsave(paste0(output_path, "/", prefix, "-449.tiff"), 
              p2, width=12, height=8)


# p1 <- plot_grid(p, bar, labels = c("A", "B"), rel_widths = c(2, 1), label_size = 12)

# ggplot2::ggsave(paste0(output_path, "/", prefix, ".png"), 
#               p1, width=16, height=8)





