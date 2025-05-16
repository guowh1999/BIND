library(argparse)
library(dplyr)

parser <- ArgumentParser()

parser$add_argument("--data_mtx_dir", nargs=1, help="Data matrix file")
parser$add_argument("--grp_info_dir", nargs=1, help="Group infomation file")
parser$add_argument("--grp_name", nargs=1, default="Nogrpname", help="Name of group column")
parser$add_argument("--output_dir", nargs=1, help="Output dir")
parser$add_argument("--grp", nargs=1, default="Nogrp",help="The group you are intreasted")
parser$add_argument("--k", nargs=1, type='integer', help="The number of neighbour")
parser$add_argument("--th", nargs=1, type='double', help="The threshold of classification")

###loading data

args <- parser$parse_args()
data_mtx <- read.csv(args$data_mtx_dir, header = T, check.names = FALSE)
info_mtx <- read.csv(args$grp_info_dir, header = T, check.names = F)
group_name <- args$grp_name
target_group <- args$grp
k <- args$k
th <- args$th

#data_mtx <- read.csv("E:/EV_ML/raw_data/Batch_PPI/RDATA/Webserver/ccle/smpdata1.csv", check.names = F)
#info_mtx <- read.csv("E:/EV_ML/raw_data/Batch_PPI/RDATA/Webserver/ccle/smpinfo1.csv", check.names = F)
# # #data_mtx <- nc_nci_df
# # #info_mtx <- nci_info
 
#group_name <- "Tissue_type"
#target_group <- "Haematopoietic and Lymphoid"
#k <- 5;th <- 3

###preprocess of data matrix
data_mtx <- data_mtx[!duplicated(data_mtx[[1]]),]
data_mtx <- data_mtx[is.na(data_mtx[[1]]) == F,]
rownames(data_mtx) <- data_mtx[,1]
data_mtx <- data_mtx[,-1]

rownames(info_mtx) <- info_mtx[[1]]

bind_stat_sub <- function(df){
  exp_df <- df
  res <- list()
  df[!is.na(df)] <- 1
  df[is.na(df)] <- 0
  row_NA_count <- rowSums(df == 0)
  res[['binary_df']] <- df #df[order(row_NA_count, decreasing = TRUE), ]
  res[['row_na_stat']] <- as.data.frame(row_NA_count)
  res[['row_na_stat']]$row_NA_ratio <- res[['row_na_stat']]$row_NA_count / ncol(df)
  res[['row_na_stat']]$means_without_na <- rowMeans(exp_df, na.rm = T)
  res[['row_na_stat']]$na_ratio_cat <- ifelse(res[['row_na_stat']]$row_NA_ratio > 0.8, 1, ifelse(res[['row_na_stat']]$row_NA_ratio < 0.2, -1, 0))
  return(res)
}

proc_res <- bind_stat_sub(data_mtx)
proc_df <- proc_res[['binary_df']]


ft_binarydf <- proc_df[rowSums(proc_df == 0) / ncol(proc_df) < 0.9,]
ft_expdf <- data_mtx[rowSums(is.na(data_mtx)) / ncol(data_mtx) < 0.9,]
ft_expdf[is.na(ft_expdf)] <- 0

### distance calculation
dist_matrix <- dist(t(ft_expdf), method = 'euclidean') %>% as.matrix()

gaussian_kernel <- function(dist_matrix, sigma) {
  sim_matrix <- exp(-dist_matrix^2 / (2 * sigma^2))
  return(sim_matrix)
}

library(progress)
na_classify <- function(na_df, exp_df, dist_mtx, class_cut, k){
  sub_af <- exp_df
  na_pt <- which(na_df == 0, arr.ind = TRUE)
  pb <- progress_bar$new(total = nrow(na_pt))
  imp_vc = c()
  for(i in 1:nrow(na_pt)){
    na_pt[i,'col']
    sample = colnames(na_df)[na_pt[i,'col']]
    protein = rownames(na_df)[na_pt[i,'row']]
    sum_total=0
    sum_sub=0
    
    sigma <- dist_mtx[sample,names(sort(dist_mtx[sample,])[k+1])]
    sim_mtx <- gaussian_kernel(dist_matrix, sigma) %>% as.matrix()
    
    for (smp_loop in names(sort(dist_mtx[sample,])[1:k+1])){
      sum_total = sum_total + sim_mtx[smp_loop,sample] * exp_df[protein, smp_loop]
      sum_sub = sum_sub + sim_mtx[smp_loop,sample]
    }
    imp_vl = sum_total / sum_sub
    imp_vc <- c(imp_vc, imp_vl)
    
    if(abs(imp_vl) < class_cut){
      sub_af[protein, sample] <- "biological"
    }
    else{
      sub_af[protein, sample] <- "technical"
    }
    pb$tick()
  }
  na_pt = na_pt %>% as.data.frame()
  na_pt$imp_vl <- imp_vc
  res = list(natp_df = as.data.frame(sub_af),
             imp_vl = na_pt)
  return(res)
}

cl_test <- na_classify(ft_binarydf, ft_expdf, dist_matrix, th, k)

cl_resdf <- cl_test[['natp_df']]
proteins <- rownames(cl_resdf)
write.csv(cl_resdf, file = paste0(args$output_dir, "/aftercalssification_label.csv"))

cl_resdf[cl_resdf != 'biological' & cl_resdf != 'technical'] <- 1
cl_resdf[cl_resdf == 'biological'] <- 0
cl_resdf[cl_resdf == 'technical'] <- 0.5
cl_resdf <- lapply(as.data.frame(cl_resdf), as.numeric) %>% as.data.frame()
rownames(cl_resdf) <- proteins

heat_df <- cl_resdf[order(rowSums(cl_resdf != 1), decreasing = T),]

write.csv(heat_df, file = paste0(args$output_dir, "/aftercalssification_mtx.csv"))
library(pheatmap)
p <- pheatmap(heat_df, color = c("#1c3f58","#fbf70c","#A12D45"), cluster_rows = FALSE, cluster_cols = FALSE,  show_colnames = T,
              show_rownames = FALSE,use_raster = F, legend = F, fontsize_col = 3)
png(filename =  paste0(args$output_dir, "/NA_classification_heatmap.png"), width = 8, height = 6,units = 'in', res = 300)
print(p)
dev.off()
unlink("Rplots.pdf")

###find markers
# info_mtx$target_grp <- ifelse(info_mtx[[group_name]] == target_group, target_group, "Other")
# write.csv(info_mtx, file = paste0(args$output_dir, "/groupinfo_labeltargetgrp.csv"))

if(group_name != "Nogrpname"){
target_df <- cl_resdf[,rownames(info_mtx[info_mtx[[group_name]] == target_group,])]
nontarget_df <- cl_resdf[,rownames(info_mtx[info_mtx[[group_name]] != target_group,])]
tg_allna_rto <- rowSums(target_df == 0) / ncol(target_df)
nontg_allna_rto <- rowSums(nontarget_df == 0) / ncol(nontarget_df)

#tg_tech_perc <- rowSums(target_df == 0.5) / rowSums(target_df != 1)
#nontg_tech_perc <- rowSums(nontarget_df == 0.5) / rowSums(nontarget_df != 1)

tg_nontg_df <- data.frame(tg = tg_allna_rto,
                          nontg = nontg_allna_rto,
                          #tg_techperc = tg_tech_perc,
                          #nontg_techperc = nontg_tech_perc,
                          proteins = rownames(target_df))

tg_nontg_df$diff <- tg_nontg_df$nontg - tg_nontg_df$tg

tg_nontg_df <- tg_nontg_df %>% dplyr::filter(.$tg < 0.7 & .$nontg > 0.3)
tg_nontg_df <- tg_nontg_df[order(tg_nontg_df$diff, decreasing = T),][1:min(50, nrow(tg_nontg_df)),]

#tg_nontg_df$techdiff <- tg_nontg_df$tg_techperc - tg_nontg_df$nontg_techperc

#tg_mks <- intersect(rownames(tg_nontg_df[order(tg_nontg_df$techdiff, decreasing = T),][1:200,]), rownames(tg_nontg_df[order(tg_nontg_df$diff, decreasing = T),][1:200,]))
#tg_mks <- rownames(tg_nontg_df[order(tg_nontg_df$diff, decreasing = T),][1:floor(nrow(tg_nontg_df)*0.2),])
tg_mks <- rownames(tg_nontg_df)

#tg_nontg_df <- tg_nontg_df[tg_mks,]
#tg_nontg_df <- tg_nontg_df %>% dplyr::filter(.$techdiff>0)

#tg_mks <- rownames(tg_nontg_df[order(tg_nontg_df$techdiff, decreasing = T),][1:min(50, nrow(tg_nontg_df)),])

mks_df <- cbind(target_df, nontarget_df)
mks_df <- mks_df[tg_mks,]

#group <- "Breast"
#info_mtx$groupsub <- ifelse(info_mtx$Tissue.of.origin == group, "Breast", "Other")
#annotation_col <- data.frame(Group = info_mtx$groupsub)
#rownames(annotation_col) <- rownames(info_mtx)

write.csv(mks_df, file = paste0(args$output_dir, "/aftercalssification_groupmarkers.csv"))
} else{
  print("No group information provided.")
}