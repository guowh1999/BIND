library(argparse)
library(dplyr)
library(data.table)
library(biomaRt)
library(tidyr)
library(doFuture)
library(progressr)

parser <- ArgumentParser()

parser$add_argument("--data_mtx_dir", nargs=1, help="Data matrix file")
parser$add_argument("--grp_info_dir", nargs=1, help="Group infomation file")
parser$add_argument("--grp_name", nargs=1, default="Nogrpname", help="Name of group column")
parser$add_argument("--output_dir", nargs=1, help="Output dir")
parser$add_argument("--grp", nargs=1, default="Nogrp",help="The group you are intreasted")
parser$add_argument("--reward", nargs=1,type='double', default=0.01,help="Reward")
parser$add_argument("--penalty", nargs=1,type='double', default=0.01,help="Penalty")
parser$add_argument("--NAtype_mtx_dir", nargs=1,help="Data matrix after NA classification")
parser$add_argument("--protein_list_dir", nargs=1,help="Input protein list dir")
parser$add_argument("--ppidb_dir", nargs=1,help="PPI database dir")

###loading data

args <- parser$parse_args()
data_mtx <- read.csv(args$data_mtx_dir, header = T, check.names = FALSE)
info_mtx <- read.csv(args$grp_info_dir, header = T, check.names = F)
group_name <- args$grp_name
group <- args$grp
reward <- args$reward
penalty <- args$penalty
type_mtx <- read.csv(args$NAtype_mtx_dir, header = T, check.names = F)
tg_pros <- read.csv(args$protein_list_dir, header = F, check.names = F)

pro9606phy <- read.table(args$ppidb_dir, header = T)

### preprocess data
data_mtx <- data_mtx[!duplicated(data_mtx[[1]]),]
data_mtx <- data_mtx[is.na(data_mtx[[1]]) == F,]
rownames(data_mtx) <- data_mtx[,1]
data_mtx <- data_mtx[,-1]

rownames(type_mtx) <- type_mtx[,1]
type_mtx <- type_mtx[,-1]

rownames(info_mtx) <- info_mtx[[1]]


tg_pros <- tg_pros[[1]]
tg_pros <- intersect(tg_pros, rownames(data_mtx))

ft_expdf <- data_mtx[rownames(type_mtx),]

###creating input data for grp/nogrp

if(args$grp == 'Nogrp'){
  ft_exp_tg <- ft_expdf[tg_pros,]
  ft_bio_tg <- type_mtx[tg_pros,]
} else{
  #filter of less than 0.5 NA in group and nongroup
  grp_exp_df <- ft_expdf[,rownames(info_mtx[info_mtx[[group_name]]==group,])]
  grp_exp_df <- grp_exp_df[rowSums(is.na(grp_exp_df)) / ncol(grp_exp_df) < 0.5,]
  
  if(ncol(grp_exp_df)<10) {print("sample number in group less than 10!")}
  
  nogrp_exp_df <- ft_expdf[,rownames(info_mtx[info_mtx[[group_name]]!=group,])]
  nogrp_exp_df <- nogrp_exp_df[rowSums(is.na(nogrp_exp_df)) / ncol(nogrp_exp_df) < 0.5,]
  
  if(ncol(nogrp_exp_df)<10) {print("sample number in other sample except target group less than 10!")}
  
  #intersect of grp and nogrp <0.5NA proteins
  inter_grp <- intersect(rownames(grp_exp_df), rownames(nogrp_exp_df))
  inter_tg <- intersect(inter_grp, tg_pros)
  
  #final data exp matrix
  grp_exp_df <- grp_exp_df[inter_tg,]
  nogrp_exp_df <- nogrp_exp_df[inter_tg,]
  
  grp_bio_df <- type_mtx[inter_tg,colnames(grp_exp_df)]
  nogrp_bio_df <- type_mtx[inter_tg, colnames(nogrp_exp_df)]
}

### get local PPI database

#pro9606phy <- pro9606phy %>% dplyr::filter(.$combined_score > 900)
pro9606phy$protein1 <- pro9606phy$protein1 %>% gsub("9606.","",.)
pro9606phy$protein2 <- pro9606phy$protein2 %>% gsub("9606.","",.)

ensembl <- useMart("ensembl")
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
result <- getBM(attributes = c("ensembl_peptide_id", "hgnc_symbol"),
                filters = "ensembl_peptide_id",
                values = unique(c(pro9606phy$protein1,pro9606phy$protein2)),
                mart = dataset)
pro1_list <- result[match(pro9606phy$protein1, result[,"ensembl_peptide_id"]),]
pro2_list <- result[match(pro9606phy$protein2, result[,"ensembl_peptide_id"]),]

pro9606phy$node1 <- pro1_list$hgnc_symbol
pro9606phy$node2 <- pro2_list$hgnc_symbol
pro9606phy <- pro9606phy %>% dplyr::filter(is.na(.$node1) == F & is.na(.$node2) == F)

### define function
cor_with_NA_type <- function(my_mat, my_mat_type, x, y, reward = NULL, punish = NULL) {
  if (is.null(reward)) { reward = 0.01 }
  if (is.null(punish)) { punish = 0.01 }
  
  cor_input <- my_mat[c(x, y),] %>% t()
  cor_input_type <- my_mat_type[c(x, y),] #%>% t()
  
  num_rows_all_na <- sum(apply(cor_input, 1, function(z) all(is.na(z))))
  num_rows_some_na <- sum(apply(cor_input, 1, function(z) any(is.na(z)) && !all(is.na(z))))
  
  
  tt <- length(which(cor_input_type[1,] == 'technical' & cor_input_type[2,] == 'technical'))
  bb <- length(which(cor_input_type[1,] == 'biological' & cor_input_type[2,] == 'biological'))
  #tb <- length(which(cor_input_type[1,] == 'biological' & cor_input_type[2,] == 'technical')) + length(which(cor_input_type[1,] == 'technical' & cor_input_type[2,] == 'biological'))
  vt <- length(which(((cor_input_type[1,] %in% c('biological','technical') == F) & cor_input_type[2,] == 'technical') | ((cor_input_type[2,] %in% c('biological','technical') == F) & cor_input_type[1,] == 'technical')))
  #vb <- length(which(((cor_input_type[1,] %in% c('biological','technical') == F) & cor_input_type[2,] == 'biological') | ((cor_input_type[2,] %in% c('biological','technical') == F) & cor_input_type[1,] == 'biological')))
  vv <- length(which((cor_input_type[1,] %in% c('biological','technical') == F) & (cor_input_type[2,] %in% c('biological','technical') == F )))
  
  tp_re <- tt+bb+vt
  
  cor_input_without_na <- cor_input[complete.cases(cor_input),]
  
  xinv <- try({
    cor_res <- cor.test(x = cor_input_without_na[, 1], y = cor_input_without_na[, 2], method = "spearman")
    weighted_rho = cor_res$estimate[[1]] * (1 + num_rows_all_na * reward) * (1 - num_rows_some_na * punish) * (1 + tp_re * reward)
    if(weighted_rho >= 0){weighted_rho = min(weighted_rho,1)}
    else {weighted_rho = max(weighted_rho,-1)}
    cor_df <- data.frame(
      protein_1 = x,
      protein_2 = y,
      raw_rho = cor_res$estimate[[1]],
      weighted_rho =weighted_rho,
      raw_p = cor_res$p.value
    )
  }, silent = TRUE)
  
  if ("try-error" %in% class(xinv)) {
    cor_df <- data.frame(
      protein_1 = x,
      protein_2 = y,
      raw_rho = NA,
      weighted_rho = NA,
      raw_p = NA
    )
  }
  return(cor_df)
}
options(future.globals.maxSize= 7*1024*1024^2)

verify_spearman_type <- function(PPI, protein_mtx, na_type_mtx, cores_use, reward, punish){
  Res <- list()
  Res[["Input_data"]] <- list("PPI" = PPI, "Protein_matrix" = protein_mtx)

  message("Preprocessing protein matrix inputted...")
  shared_protn <- intersect(unique(c(PPI[['node1']],PPI[['node2']])),rownames(protein_mtx))
  protein_mtx <- Res[["Preprocessed_data"]][["Protein_matrix"]] <- protein_mtx[shared_protn,]
  PPI <- Res[["Preprocessed_data"]][["PPI"]] <- PPI %>% subset(PPI[['node1']] %in% shared_protn & PPI[['node2']] %in% shared_protn)
  
  message("Runing weighted Spearman in parallel...")
  registerDoFuture()
  plan(multisession,workers = cores_use)
  with_progress({
    p <- progressor(along = shared_protn) 
    all_cor_res <- foreach(p_x = shared_protn,
                           .packages = c("purrr","data.table","tidyverse"),
                           .combine = rbind) %dopar% {
                             p()
                             x_cor_res <- purrr::map(shared_protn, function(p_y){
                               if(p_x == p_y){
                                 temp <- data.frame(
                                   protein_1 = NA,
                                   protein_2 = NA,
                                   raw_rho = NA,
                                   weighted_rho = NA,
                                   raw_p = NA)
                               }else{
                                 temp <- cor_with_NA_type(my_mat = protein_mtx, my_mat_type = na_type_mtx, x = p_x, y = p_y, reward = reward, punish = punish)
                               }
                               return(temp)
                             }) %>% dplyr::bind_rows()
                             
                           }
    
  })
  
  message("Curating results...")
  PPI[['Pair']] <- paste(PPI[['node1']],PPI[['node2']], sep = " ")
  Res[['Results']][['RAW_Cor_Res']] <- all_cor_res 
  
  Res[['Results']][['Final_Cor_Res']] <- all_cor_res %>%
    # dplyr::filter(if_any(everything(), ~ !is.na(.))) #del all NA row only
    tidyr::drop_na() %>% # del any row with NA 
    {.[['Pair']] <- paste(.[[1]],.[[2]], sep = " ");
    .[['IN_PPI']] <- ifelse(.[['Pair']] %in% PPI[['Pair']],"IN","OUT");
    .[["PPI_score"]] <- lapply(.[['Pair']], function(pr){ifelse(pr %in% PPI[['Pair']], 
                                                                return(subset(PPI, PPI[['Pair']] == pr)[['combined_score']]),return(NA))});.}
  message("Done.")
  return(Res)
}

rm_pair <- function(df){
node <- unique(c(df$protein_1,df$protein_2))
node_comb <- vector()
k=0
for(i in 1:length(node)){
  for(j in 1:length(node)){
    if(i < j){
      k = k + 1
      node_comb[k] <- paste0(node[i], " ", node[j])
    }}}
return(df[df[['Pair']] %in% node_comb,])
}

### conduct ppi analysis
if(args$grp == 'Nogrp'){
ppi_type_res <- verify_spearman_type(PPI = pro9606phy, 
                                     protein_mtx = ft_exp_tg, 
                                     na_type_mtx = ft_bio_tg, 
                                     cores_use = 8, 
                                     reward = reward, 
                                     punish = penalty)
ppi_resdf <- ppi_type_res[['Results']][['Final_Cor_Res']]
ppi_resdf$PPI_score <- ppi_resdf$PPI_score %>% unlist()

ppi_resdf <- rm_pair(ppi_resdf)
ppi_resdf$weighted_raw_change <- abs(ppi_resdf$weighted_rho) - abs(ppi_resdf$raw_rho)
mean_para <- mean(reward, penalty)
ppi_resdf_ft <- ppi_resdf %>% dplyr::filter(.$weighted_rho > 0.5)
ppi_resdf_ft <- ppi_resdf_ft[order(ppi_resdf_ft$weighted_raw_change, decreasing = T),][0:min(200,nrow(ppi_resdf_ft)),]

write.csv(ppi_resdf, file = paste0(args$output_dir, "/BINDppi_withoutgroup.csv"))
write.csv(ppi_resdf_ft, file = paste0(args$output_dir, "/BINDppi_withoutgroup_filtered.csv"))
} else{
  ppi_type_grp <- verify_spearman_type(PPI = pro9606phy, 
                                       protein_mtx = grp_exp_df, 
                                       na_type_mtx = grp_bio_df, 
                                       cores_use = 8, 
                                       reward = reward, 
                                       punish = penalty)
  ppi_resdf_grp <- ppi_type_grp[['Results']][['Final_Cor_Res']]
  ppi_resdf_grp$PPI_score <- ppi_resdf_grp$PPI_score %>% unlist()
  ppi_resdf_grp <- rm_pair(ppi_resdf_grp)
  
  ppi_type_nogrp <- verify_spearman_type(PPI = pro9606phy, 
                                        protein_mtx = nogrp_exp_df, 
                                        na_type_mtx = nogrp_bio_df, 
                                        cores_use = 8, 
                                        reward = reward, 
                                        punish = penalty)
  ppi_resdf_nogrp <- ppi_type_nogrp[['Results']][['Final_Cor_Res']]
  ppi_resdf_nogrp$PPI_score <- ppi_resdf_nogrp$PPI_score %>% unlist()
  ppi_resdf_nogrp <- rm_pair(ppi_resdf_nogrp)
  
  grp_integ_df <- data.frame(protein_1 = ppi_resdf_grp$protein_1,
                             protein_2 = ppi_resdf_grp$protein_2,
                             raw_rho_group = ppi_resdf_grp$raw_rho,
                             raw_rho_nogroup = ppi_resdf_nogrp$raw_rho,
                             weighted_rho_group = ppi_resdf_grp$weighted_rho,
                             weighted_rho_nogroup = ppi_resdf_nogrp$weighted_rho,
                             Pair = ppi_resdf_grp$Pair,
                             IN_PPI = ppi_resdf_grp$IN_PPI,
                             PPI_score = ppi_resdf_grp$PPI_score)
  grp_integ_df$weighted_rho_diff <- grp_integ_df$weighted_rho_group - grp_integ_df$weighted_rho_nogroup
  grp_integ_df$raw_rho_diff <- grp_integ_df$raw_rho_group - grp_integ_df$raw_rho_nogroup
  #grp_integ_df$group_weighted_raw_change <- grp_integ_df$weighted_rho_group - grp_integ_df$raw_rho_group
  #grp_integ_df$nogroup_weighted_raw_change <- grp_integ_df$weighted_rho_nogroup - grp_integ_df$raw_rho_nogroup
  
  mean_para <- mean(reward, penalty)
  grp_integ_ft <- grp_integ_df %>% dplyr::filter(abs(.$weighted_rho_nogroup) < 0.5 & abs(.$weighted_rho_group) > 0.5 & abs(.$weighted_rho_diff) - abs(.$raw_rho_diff) > mean_para)
  
  grp_integ_ft <- grp_integ_ft[order(grp_integ_ft$weighted_rho_diff, decreasing = T),][0:min(200,nrow(grp_integ_ft)),]
  
  write.csv(grp_integ_df, file = paste0(args$output_dir, "/BINDppi_groupdiff.csv"))
  write.csv(grp_integ_ft, file = paste0(args$output_dir, "/BINDppi_groupdiff_filtered.csv"))
}

