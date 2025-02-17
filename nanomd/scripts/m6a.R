filter_modification_mRNA <- function() {
    # 0.test
    input_path <- "./1-output/20250107_mm_2sample/00_inputs"
    output_path <- input_path
    prefix <- "ZHY_"
    # 1.library package
    library(clusterProfiler)
    # 2. read files
    options(warn = -1)
    ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))
    # 3. run
    ifelse(length(list.files(output_path, pattern = "WashU.txt")) > 0, return("sucess"),
    tryCatch({
    # START
    mRNA_list <- read.delim(paste0(input_path, "/protein.gene"), header = F) %>% 
        dplyr::mutate(type = "protein")
    # Treat
    mod_raw_m1a <- read.delim(paste0(input_path, "/Treat/m1A.plus.bed"), header = F) %>% 
        tidyr::separate(V5, into = c("genid", "SYMBOL", "cds", "site"), sep = "[|]", remove = F) %>% 
        left_join(mRNA_list, by = c("SYMBOL" = "V1")) %>% 
        dplyr::filter(type == "protein") %>%
        dplyr::select(c("V1","V2","V3","V4","V5","V6","V7","V8")) %>% 
        write.table(., paste0(input_path, "/Treat/m1A.plus.mrna.bed"), sep = "\t", row.names = F, col.names = F, quote = F)
    
    mod_raw_m6a <- read.delim(paste0(input_path, "/Treat/m6A.plus.bed"), header = F) %>% 
        tidyr::separate(V5, into = c("genid", "SYMBOL", "cds", "site"), sep = "[|]", remove = F) %>% 
        left_join(mRNA_list, by = c("SYMBOL" = "V1")) %>% 
        dplyr::filter(type == "protein") %>%
        dplyr::select(c("V1","V2","V3","V4","V5","V6","V7","V8")) %>% 
        write.table(., paste0(input_path, "/Treat/m6A.plus.mrna.bed"), sep = "\t", row.names = F, col.names = F, quote = F)

    mod_raw_m5c <- read.delim(paste0(input_path, "/Treat/m5C.plus.bed"), header = F) %>% 
        tidyr::separate(V5, into = c("genid", "SYMBOL", "cds", "site"), sep = "[|]", remove = F) %>% 
        left_join(mRNA_list, by = c("SYMBOL" = "V1")) %>% 
        dplyr::filter(type == "protein") %>%
        dplyr::select(c("V1","V2","V3","V4","V5","V6","V7","V8")) %>% 
        write.table(., paste0(input_path, "/Treat/m5C.plus.mrna.bed"), sep = "\t", row.names = F, col.names = F, quote = F)

    mod_raw_psi <- read.delim(paste0(input_path, "/Treat/psi.plus.bed"), header = F) %>% 
        tidyr::separate(V5, into = c("genid", "SYMBOL", "cds", "site"), sep = "[|]", remove = F) %>% 
        left_join(mRNA_list, by = c("SYMBOL" = "V1")) %>% 
        dplyr::filter(type == "protein") %>%
        dplyr::select(c("V1","V2","V3","V4","V5","V6","V7","V8")) %>% 
        write.table(., paste0(input_path, "/Treat/psi.plus.mrna.bed"), sep = "\t", row.names = F, col.names = F, quote = F)
    # WT
    mod_raw_m1a <- read.delim(paste0(input_path, "/WT/m1A.plus.bed"), header = F) %>% 
        tidyr::separate(V5, into = c("genid", "SYMBOL", "cds", "site"), sep = "[|]", remove = F) %>% 
        left_join(mRNA_list, by = c("SYMBOL" = "V1")) %>% 
        dplyr::filter(type == "protein") %>%
        dplyr::select(c("V1","V2","V3","V4","V5","V6","V7","V8")) %>% 
        write.table(., paste0(input_path, "/WT/m1A.plus.mrna.bed"), sep = "\t", row.names = F, col.names = F, quote = F)
    
    mod_raw_m6a <- read.delim(paste0(input_path, "/WT/m6A.plus.bed"), header = F) %>% 
        tidyr::separate(V5, into = c("genid", "SYMBOL", "cds", "site"), sep = "[|]", remove = F) %>% 
        left_join(mRNA_list, by = c("SYMBOL" = "V1")) %>% 
        dplyr::filter(type == "protein") %>%
        dplyr::select(c("V1","V2","V3","V4","V5","V6","V7","V8")) %>% 
        write.table(., paste0(input_path, "/WT/m6A.plus.mrna.bed"), sep = "\t", row.names = F, col.names = F, quote = F)

    mod_raw_m5c <- read.delim(paste0(input_path, "/WT/m5C.plus.bed"), header = F) %>% 
        tidyr::separate(V5, into = c("genid", "SYMBOL", "cds", "site"), sep = "[|]", remove = F) %>% 
        left_join(mRNA_list, by = c("SYMBOL" = "V1")) %>% 
        dplyr::filter(type == "protein") %>%
        dplyr::select(c("V1","V2","V3","V4","V5","V6","V7","V8")) %>% 
        write.table(., paste0(input_path, "/WT/m5C.plus.mrna.bed"), sep = "\t", row.names = F, col.names = F, quote = F)

    mod_raw_psi <- read.delim(paste0(input_path, "/WT/psi.plus.bed"), header = F) %>% 
        tidyr::separate(V5, into = c("genid", "SYMBOL", "cds", "site"), sep = "[|]", remove = F) %>% 
        left_join(mRNA_list, by = c("SYMBOL" = "V1")) %>% 
        dplyr::filter(type == "protein") %>%
        dplyr::select(c("V1","V2","V3","V4","V5","V6","V7","V8")) %>% 
        write.table(., paste0(input_path, "/WT/psi.plus.mrna.bed"), sep = "\t", row.names = F, col.names = F, quote = F)


    
    # END
    return("sucess")},
    error = function(e){print(e$message);message(return("failled"))}))
}

plot_metagene_Rd <- function(input_data, output_path, prefix) {
    ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))
    tryCatch({
    # START
    # cat(crayon::green(paste0(prefix, "Start ploting Metagene\n")))
    colour_fill <- c("black", "red")
    names(colour_fill) <- input_data$group  %>% unique()
    p1 <- ggplot(input_data, aes(x=value))+
        geom_line(aes(colour = group), stat = "density", adjust = 2) +
        scale_colour_manual(values=colour_fill) +
        geom_vline(xintercept = 1:2, col = "red", linetype="dashed") + 
        annotate("text", x = 0.5, y = -0.14, label = "5'UTR") +
        annotate("text", x = 1.5, y = -0.14, label = "CDS") +
        annotate("text", x = 2.5, y = -0.14, label = "3'UTR")  +
        annotate("rect", xmin = 0, xmax = 1, ymin = -0.08, ymax = -0.04, alpha = .99, colour = "black")+
        annotate("rect", xmin = 2, xmax = 3, ymin = -0.08, ymax = -0.04, alpha = .99, colour = "black")+
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.12, ymax = 0, alpha = .2, colour = "black") +
        xlab("psi metagene") +
        ylab("Frequency") +
        guides(colour = guide_legend(title = NULL)) +
        theme_bw() +
        theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
        theme(legend.background = element_blank()) 
    ggsave(paste0(output_path, "/", prefix, "_metagene.pdf"), p1, width = 8, height = 6)
    ggsave(paste0(output_path, "/", prefix, "_metagene.png"), p1, width = 8, height = 6)
    
    # cat(crayon::green(paste0(prefix, "Metagene analysis done\n")))
    # END
    return(p1)},
    error = function(e){print(e$message);message(return("fail"))})
}
restructure_coord <- function(input_path, group) {
    m6a.dist <- read.delim(input_path, header = T) %>%
        dplyr::filter(cds_size > utr3_size & utr3_size > utr5_size) %>% 
        dplyr::mutate(
            trx_len = utr5_size + cds_size + utr3_size)

    temp.df <- m6a.dist %>%
        dplyr::select(c("gene_name", "refseqID", "trx_len")) %>%
        dplyr::rename(setNames(colnames(.), c("gene_name", "gid", "trx_len"))) %$%
        .[order(.$gene_name, .$gid, -.$trx_len), ] %$%
        .[!duplicated(.$gene_name), ]

    m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid, ]
    utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
    utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
    # print(paste0("utr5.SF: ", utr5.SF, " utr3.SF: ", utr3.SF))
    utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
    cds.m6a.dist <- m6a.dist[m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
    utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2 & m6a.dist$rel_location < 3, ]
    utr5.m6a.dist$rel_location <- scales::rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0, 1))
    utr3.m6a.dist$rel_location <- scales::rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2, 3))
    m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)

    m6a_data <- data.frame("value"= m6a.metagene.coord) %>% 
        dplyr::mutate("group"= group)
    return(m6a_data)
}
auto_meta_plot <- function(input_path, output_path, prefix) {
    # test
    input_path <- "./1-output/20250107_mm_2sample/00_inputs/all"
    output_path <- "./1-output/20250107_mm_2sample/01_output/psi"
    prefix <- "ZHY_"
    # 2. check output path
    options(warn = -1)
    ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))
    # START
    sample_all <- legendBaseModel::make_sample_sheet(input_path, ".psi.dist.measures.txt") %>%
        dplyr::mutate(
            group = sample,
            dist_measures = purrr::map2(sample_path, group, restructure_coord)
        ) %>% 
        dplyr::mutate(
            output_plot = output_path,
            prefix_plot = paste0(prefix, group),
            plot_meta <- purrr::pmap(list(dist_measures, output_plot, prefix_plot), plot_metagene_Rd)
        )
     p1 <- legendBaseModel::rbind_merge(sample_all, "dist_measures") %>% 
        plot_metagene_Rd(., output_path, prefix)
    # END
    return(p1)
}
