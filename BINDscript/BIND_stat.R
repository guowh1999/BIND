#! /usr/bin/env Rscript
library(argparse)
library(dplyr)
library(grid)

parser <- ArgumentParser()

parser$add_argument("--data_mtx_dir", nargs=1, help="Data matrix file")
parser$add_argument("--grp_info_dir", nargs=1, help="Group infomation file")
parser$add_argument("--grp_name", nargs=1, default="Nogrpname", help="Name of group column")
parser$add_argument("--output_dir", nargs=1, help="Output dir")

args <- parser$parse_args()

data_mtx <- read.csv(args$data_mtx_dir, header = T, check.names = FALSE)

#data_mtx <- read.csv("E:/EV_ML/raw_data/Batch_PPI/RDATA/Webserver/ccle/smpdata.csv", header = T, check.names = FALSE)

data_mtx <- data_mtx[!duplicated(data_mtx[[1]]),]
data_mtx <- data_mtx[is.na(data_mtx[[1]]) == F,]
rownames(data_mtx) <- data_mtx[,1]
data_mtx <- data_mtx[,-1]

bind_stat <- function(df){
  exp_df <- df
  res <- list()
  df[!is.na(df)] <- 1
  df[is.na(df)] <- 0
  row_NA_count <- rowSums(df == 0)
  res[['binary_df']] <- df[order(row_NA_count, decreasing = TRUE), ]
  res[['row_na_stat']] <- as.data.frame(row_NA_count)
  res[['row_na_stat']]$row_NA_ratio <- res[['row_na_stat']]$row_NA_count / ncol(df)
  res[['row_na_stat']]$means_without_na <- rowMeans(exp_df, na.rm = T)
  res[['row_na_stat']]$na_ratio_cat <- ifelse(res[['row_na_stat']]$row_NA_ratio > 0.8, 1, ifelse(res[['row_na_stat']]$row_NA_ratio < 0.2, -1, 0))
  return(res)
}

data_mtx_stat <- bind_stat(data_mtx)

write.csv(data_mtx_stat[['binary_df']], file = paste0(args$output_dir, "/binary_mtx.csv"))
write.csv(data_mtx_stat[['row_na_stat']], file = paste0(args$output_dir, "/overall_stat.csv"))

library(pheatmap)
p <- pheatmap(data_mtx_stat[['binary_df']], color = c("#1c3f58","#A12D45"),cluster_rows = FALSE, cluster_cols = F,  show_colnames = T,
                show_rownames = FALSE, legend = F, fontsize_col = 3)


png(filename =  paste0(args$output_dir, "/NA_stat_heatmap.png"), width = 8, height = 6,units = 'in', res = 300)
#jpeg(paste0(args$output_dir, "/NA_stat_heatmap.jpg"), width = 8, height = 6,units = 'in', res = 300)
#pdf(paste0(args$output_dir, "/NA_stat_heatmap.pdf"), width = 8, height = 6)
print(p)
dev.off()                        


###group stat
info_mtx <- read.csv(args$grp_info_dir, header = T, check.names = F)
#info_mtx <- cellline_info; group_name <- "Tissue_type"
rownames(info_mtx) <- info_mtx[[1]]

group_name <- args$grp_name

if(group_name != "Nogrpname"){

bind_group_stat <- function(df, meta_df, group_name){
  #df <- nc_nci_df;meta_df <- organ_info_nci;group_name <- "tissue"
  groups <- table(meta_df[,colnames(meta_df) %in% group_name])
  groups <- names(groups[groups>=5])
  grp_na_rto <- list()
  grp_na_stat <- list()
  nogrp_na_rto <- list()
  nogrp_na_stat <- list()
  for(grp in groups){
    meta_df[[group_name]]
    smp <- meta_df %>% dplyr::filter(.[[group_name]] == grp) %>% rownames()
    grp_df <- df[,colnames(df) %in% smp]
    nogrp_df <- df[,colnames(df) %in% smp == F]
    
    grp_na_stat[[grp]] <- bind_stat(grp_df)
    grp_na_rto[[grp]] <- grp_na_stat[[grp]]$row_na_stat$row_NA_ratio
    
    nogrp_na_stat[[paste0("no_",grp)]] <- bind_stat(nogrp_df)
    nogrp_na_rto[[paste0("no_",grp)]] <- nogrp_na_stat[[paste0("no_",grp)]]$row_na_stat$row_NA_ratio
  }
  grp_comb <- do.call(cbind, grp_na_rto)
  nogrp_comb <- do.call(cbind, nogrp_na_rto)

  rownames(grp_comb) <- rownames(df)
  rownames(nogrp_comb) <- rownames(df)
  return(list(grp_comb, nogrp_comb))
}

grp_stat <- bind_group_stat(data_mtx, info_mtx, group_name)
write.csv(grp_stat[[1]], file = paste0(args$output_dir, "/group_stat.csv"))
write.csv(do.call(cbind, grp_stat), file = paste0(args$output_dir, "/img_grpvsnogrp_stat.csv"))

dist_matrix <- dist(t(grp_stat[[1]]), method = 'euclidean')
hclust_result <- hclust(dist_matrix)
library(ape)
my_tree <- as.phylo(hclust_result)
write.tree(phy=my_tree, file=paste0(args$output_dir,"/img_grp_naproportion_tree.newick"))


stat <- bind_stat(data_mtx)
dist_matrix <- dist(t(stat[['binary_df']]), method = 'binary')
hclust_result <- hclust(dist_matrix)
my_tree <- as.phylo(hclust_result)
write.tree(phy=my_tree, file=paste0(args$output_dir,"/img_binary_tree.newick"))

} else{
  print("No group information provided.")
}
# bind_pie <- function(df){
#   stat <- bind_stat(df)
#   stat_df <- stat[['row_na_stat']]
#   pie_df <- data.frame("no NA"=nrow(stat_df %>% dplyr::filter(.$row_NA_ratio == 0)),
#              "NA_proportion less than 0.2"=nrow(stat_df %>% dplyr::filter(.$row_NA_ratio > 0 & .$row_NA_ratio <= 0.2)),
#              "NA_proportion 0.2 to 0.8"=nrow(stat_df %>% dplyr::filter(.$row_NA_ratio <= 0.8 & .$row_NA_ratio >= 0.2)),
#              "NA_proportion more than 0.8"=nrow(stat_df %>% dplyr::filter(.$row_NA_ratio >= 0.8)))
#   return(pie_df)
# }
# pie_stat <- bind_pie(data_mtx)
# write.csv(pie_stat, file = paste0(args$output_dir, "/img_pie.csv"))
