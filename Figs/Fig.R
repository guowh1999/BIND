setwd("E:/EV_ML/raw_data/Batch_PPI/RDATA/Figcode")

libraries <- c("dplyr", "ggplot2", "tidyr", "stringr")
lapply(libraries, library, character.only = TRUE)

#data reading and preprocessing-------------------------------------------------
# ccle 949 cell line dataset
cellline_df <- read.csv("./dataset/ccle949/cell_2022.csv", header = T)
cellline_df$Project_Identifier <- as.data.frame(str_split_fixed(cellline_df$Project_Identifier, ';', 2))[[2]]
rownames(cellline_df) <- cellline_df[,1]
cellline_df <- cellline_df[,-1]
cellline_df <- t(cellline_df)
cellline_df <- cellline_df %>% as.data.frame()

cellline_info <- read.csv("./dataset/ccle949/celllines_info.csv")
rownames(cellline_info) <- cellline_info[,1]

# nci 60 cell line data
nc_nci_df <- read.csv("./dataset/nci60/NCI60.csv", header = T)
gn_process <- nc_nci_df$Gene_symbol
nc_nci_df$Gene_symbol <- sapply(strsplit(gn_process, ";"), function(x) x[1])
nc_nci_df <- nc_nci_df[!duplicated(nc_nci_df$Gene_symbol),]
nc_nci_df <- nc_nci_df[is.na(nc_nci_df$Gene_symbol) == F,]
rownames(nc_nci_df) <- nc_nci_df[,1]
nc_nci_df <- nc_nci_df[,-1]

nci_info <- read.csv("./dataset/nci60/NCI60_info.csv")
organ_info_nci <- data.frame(nci_info[,2])
rownames(organ_info_nci) <- nci_info[,1]

# crc 65 cell line data

nc_crc_df <- read.csv("./dataset/crc65/CRC65.csv", header = T)
gn_process <- nc_crc_df$Gene_symbol
nc_crc_df$Gene_symbol <- sapply(strsplit(gn_process, ";"), function(x) x[1])
nc_crc_df <- nc_crc_df[!duplicated(nc_crc_df$Gene_symbol),]
nc_crc_df <- nc_crc_df[is.na(nc_crc_df$Gene_symbol) == F,]
rownames(nc_crc_df) <- nc_crc_df[,1]
nc_crc_df <- nc_crc_df[,-1]

# define NA transform function
process_na <- function(df){
  res <- list()
  df[!is.na(df)] <- 1
  df[is.na(df)] <- 0
  row_ones_count <- rowSums(df == 0)
  res[['sorted_df']] <- df[order(row_ones_count, decreasing = TRUE), ]
  res[['row_na_stat']] <- as.data.frame(row_ones_count)
  res[['row_na_stat']]$row_ones_ratio <- res[['row_na_stat']]$row_ones_count / ncol(df)
  return(res)
}

process_na_nosort <- function(df){
  res <- list()
  df[!is.na(df)] <- 1
  df[is.na(df)] <- 0
  row_ones_count <- rowSums(df == 0)
  res[['sorted_df']] <- df
  res[['row_na_stat']] <- as.data.frame(row_ones_count)
  res[['row_na_stat']]$row_ones_ratio <- res[['row_na_stat']]$row_ones_count / ncol(df)
  return(res)
}

#Figure 2---------------------------------------------------------------------
##Figure 2A
library(pheatmap)
cl_res <- process_na(cellline_df)
cl_res_df <- cl_res[['sorted_df']]
pheatmap(cl_res_df, color = c("#1c3f58","#A12D45"), cluster_rows = FALSE, cluster_cols = F,  show_colnames = FALSE,
         show_rownames = FALSE) #ccle 949

na_res_nci <- process_na(nc_nci_df)
pheatmap(na_res_nci[['sorted_df']], color = c("#1c3f58","#A12D45"), cluster_rows = FALSE, cluster_cols = FALSE,  show_colnames = FALSE,
         show_rownames = FALSE) #nci 60

na_res_crc <- process_na(nc_crc_df)
pheatmap(na_res_crc[['sorted_df']], color = c("#1c3f58","#A12D45"), cluster_rows = FALSE, cluster_cols = F,  show_colnames = FALSE,
              show_rownames = FALSE) # crc65

##Figure2B
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
pieplot <- function(df) {
  cl_pie <- data.frame(
    "no_na"      = nrow(df %>% dplyr::filter(row_ones_ratio == 0)),
    "na_0p2"     = nrow(df %>% dplyr::filter(row_ones_ratio > 0 & row_ones_ratio <= 0.2)),
    "na_0p2_0p8" = nrow(df %>% dplyr::filter(row_ones_ratio > 0.2 & row_ones_ratio <= 0.8)),
    "na_0p8"     = nrow(df %>% dplyr::filter(row_ones_ratio > 0.8))
  ) %>% 
    t() %>% 
    as.data.frame() %>%
    dplyr::rename(count = V1)
  
  cl_pie <- cl_pie %>%
    dplyr::mutate(
      cat = rownames(cl_pie),
      percentage = paste0(round(count / nrow(df) * 100, 1), "%"),
      label = paste0(count, " (", percentage, ")")
    )
  
  cl_pie <- cl_pie %>%
    dplyr::arrange(desc(cat)) %>%
    dplyr::mutate(
      y_position = cumsum(count) - 0.5 * count
    )
  
  pieplot <- ggplot(cl_pie, aes(x = "", y = count, fill = cat)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    geom_text(
      aes(y = y_position, label = label), 
      color = "black",
      size = 3, 
      fontface = "bold" 
    ) +
    scale_fill_manual(values = c("#5b7c99", "#88d3db", "#93a8ac", "#3e517a")) +
    blank_theme + 
    theme(
      axis.text.x = element_blank(),
      axis.title = element_blank()
    ) +
    labs(fill = "Categories") 
  
  return(pieplot)
}

pieplot(cl_res[['row_na_stat']]) #ccle 949
pieplot(na_res_nci[['row_na_stat']]) #nci60
pieplot(na_res_crc[['row_na_stat']]) #crc65

##Figure 2C

ggplot(cl_res[['row_na_stat']] %>% dplyr::filter(.$row_ones_ratio>0), aes(x = row_ones_ratio))+
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), fill = '#1c3f58', alpha = 0.8)+
  theme_bw()+theme(panel.grid=element_blank()) #ccle 949

ggplot(na_res_nci[['row_na_stat']]%>% dplyr::filter(.$row_ones_ratio >0), aes(x = row_ones_ratio))+
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),fill = '#1c3f58', alpha = 0.8)+
  theme_bw()+theme(panel.grid=element_blank()) #nci 60

ggplot(na_res_crc[['row_na_stat']]%>% dplyr::filter(.$row_ones_ratio >0), aes(x = row_ones_ratio))+
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),boundary = min(na_res_crc[['row_na_stat']]$row_ones_ratio), fill = '#1c3f58', alpha = 0.8, bins = 30)+
  theme_bw()+theme(panel.grid=element_blank()) #crc 65
 
##Figure 2D
cl_res_stat <- cl_res[['row_na_stat']]
cl_res_stat$means_without_na <- rowMeans(cellline_df, na.rm = T)

fit <- lm(means_without_na ~ row_ones_ratio, data = cl_res_stat)
slope <- coef(fit)["row_ones_ratio"]
intercept <- coef(fit)["(Intercept)"]
r <- round(-sqrt(summary(fit)$r.squared), 2)


ggplot(cl_res_df,aes(x=row_ones_ratio, y=means_without_na))+
  geom_point(color='#A12D45', alpha=0.7)+
  geom_smooth(method = "lm", se=TRUE, color="#1C3F58")+
  annotate("text", x = Inf, y = -Inf, label = paste("r=", r,";","pvalue<2.2e-16"), hjust=1, vjust=-18, size=8)+
  theme_bw()+theme(panel.grid=element_blank())

##Figure 2E
library(ggvenn)
no_na_cl <- cl_res[['row_na_stat']] %>% dplyr::filter(.$row_ones_ratio ==0)
no_na_nci <- na_res_nci[['row_na_stat']] %>% dplyr::filter(.$row_ones_ratio == 0)
no_na_crc <- na_res_crc[['row_na_stat']] %>% dplyr::filter(.$row_ones_ratio == 0)

ggvenn(data = list(cellline_949 = rownames(no_na_cl), nci60 = rownames(no_na_nci), crc65 = rownames(no_na_crc)), 
       fill_color = c("#A12D45","#1C3F58","#8E736C"), stroke_color = "white", fill_alpha = 0.7)

##Figure 2F
nona_all <- intersect(intersect(rownames(no_na_cl),rownames(no_na_nci)),rownames(no_na_crc))
library(clusterProfiler)
library(org.Hs.eg.db)
library(DO.db)

go_plot <- function(genes){
diff_go <- enrichGO(gene = genes,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = 'ALL',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 1)
group_go_table <- as.data.frame(diff_go)

barplot(diff_go, x = "GeneRatio",
        fill = "pvalue",
        showCategory =5, 
        split="ONTOLOGY") + 
  scale_fill_gradient(low = "#8e736c",high ='#1C3F58')+
  facet_grid(ONTOLOGY~., scale='free') + theme_bw()+theme(panel.grid=element_blank())
}

go_plot(nona_all)

##Figure 2G
almost_na <- cl_res[['row_na_stat']] %>% dplyr::filter(.$row_ones_ratio >= 0.8)
almost_na_nci <- na_res_nci[['row_na_stat']] %>% dplyr::filter(.$row_ones_ratio >= 0.8)
almost_na_crc <- na_res_crc[['row_na_stat']] %>% dplyr::filter(.$row_ones_ratio >= 0.8)
ggvenn(data = list(cellline_949 = rownames(almost_na), nci60 = rownames(almost_na_nci), crc65 = rownames(almost_na_crc)), 
       fill_color = c("#A12D45","#1C3F58","#8E736C"), stroke_color = "white", fill_alpha = 0.7)

##Figure 2H
ccle_almost <- setdiff(rownames(almost_na),union(rownames(almost_na_crc), rownames(almost_na_nci)))
go_plot(ccle_almost)

# Figure 3----------------------------------------------------------------------
##Figure 3A
cl_res_df <- cl_res[['sorted_df']]
subtypes <- c("Lung","Breast","Haematopoietic and Lymphoid","Large Intestine","Stomach")
subcancer_process <- function(df, can){row_ones_count = rowSums(df == 0); stat = as.data.frame(row_ones_count); stat$row_ones_ratio = stat$row_ones_count / ncol(df); stat$cancer_from = can; return(stat)}
subtype_df <- list()
subtype_comb <- list()
for (c in subtypes){
subtype_df[[c]] <- subcancer_process(cl_res_df[,cellline_info[cellline_info$Tissue_type == c,]$Cell_line], c)
subtype_comb[[c]] <- subtype_df[[c]]$row_ones_ratio
}

subtype_comb <- do.call(cbind, subtype_comb)
rownames(subtype_comb) <- rownames(cellline_df)
corr_mtx <- cor(subtype_comb[rowSums(subtype_comb == 0) != ncol(subtype_comb),])
pheatmap(corr_mtx,cellwidth = 20, cellheight = 20,
         breaks = seq(0.8, 1, length.out = 100),  
         color = colorRampPalette(c("#1C3F58", "white", "#A12D45"))(100))

##Figure 3B
library(patchwork)
ggplot(rbind(subtype_df[['Lung']], subtype_df[['Haematopoietic and Lymphoid']]), aes(x = row_ones_ratio, color = cancer_from))+
  geom_density(alpha=0.5, linewidth=1.5)+
  scale_color_manual(name = 'cancer_from',
                     guide = 'legend',
                     values = c("Lung" = "#fdb813",
                                "Haematopoietic and Lymphoid" = "#1C3F58"))+theme_bw() +theme(panel.grid=element_blank())+
  ggplot(rbind(subtype_df[['Lung']], subtype_df[['Haematopoietic and Lymphoid']]), aes(x = row_ones_ratio, fill = cancer_from))+
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), alpha=0.9)+
  scale_fill_manual(name = 'cancer_from',
                    guide = 'legend',
                    values = c("Lung" = "#fdb813",
                               "Haematopoietic and Lymphoid" = "#1C3F58"))+theme_bw() +theme(panel.grid=element_blank()) +plot_layout(nrow = 2)

##Figure 3C
### left panel(jaccard distance)

sub_cellline <- cl_res_df[,rownames(cellline_info[cellline_info$Tissue_type %in% c("Lung","Haematopoietic and Lymphoid"),])]
dist_matrix <- dist(t(sub_cellline), method = 'binary')
hclust_result <- hclust(dist_matrix)

library(ggtree);library(ggtreeExtra);library(ggnewscale)
ggtree(hclust_result,layout = "fan")+ #%<+% cellline_info +
  #geom_tippoint(aes(col=Cancer_type), size=1)
  geom_fruit(
    data=cellline_info,
    geom=geom_tile,  
    mapping=aes(y=Cell_line, fill=Tissue_type),  
    color="black",
    width=0.08,
    #width=22,
    stat="identity",
    offset = 0.15
  )+scale_fill_manual(values=c( "#3e517a","#fdb813"))+
  new_scale_fill()+
  geom_fruit(
    data=cellline_info,
    geom=geom_tile,  
    mapping=aes(y=Cell_line, fill=Cancer_type),  
    color="black",
    width=0.08,
    #width=22,
    stat="identity",
    offset = 0.15,
  )

###right panel(euclidean distance)
sub_cellline <- cellline_df[,rownames(cellline_info[cellline_info$Tissue_type %in% c("Lung","Haematopoietic and Lymphoid"),])]
dist_matrix <- dist(t(sub_cellline), method = 'euclidean')
hclust_result <- hclust(dist_matrix)
ggtree(hclust_result,layout = "fan")+ #%<+% cellline_info +
  #geom_tippoint(aes(col=Cancer_type), size=1)
  geom_fruit(
    data=cellline_info,
    geom=geom_tile,  
    mapping=aes(y=Cell_line, fill=Tissue_type),  
    color="black",
    #width=0.08,
    width=22,
    stat="identity",
    offset = 0.15
  )+scale_fill_manual(values=c( "#3e517a","#fdb813"))+
  new_scale_fill()+
  geom_fruit(
    data=cellline_info,
    geom=geom_tile,  
    mapping=aes(y=Cell_line, fill=Cancer_type),  
    color="black",
    #width=0.08,
    width=22,
    stat="identity",
    offset = 0.15,
  )

##Figure 3D
tg <- 'Lung'

non_target_df <- cl_res[['sorted_df']][,cellline_info[cellline_info$Tissue_type != tg,]$Cell_line]
non_target_stat <- subcancer_process(non_target_df, paste0('non_',tg))
target_stat <- subcancer_process(cl_res[['sorted_df']][,cellline_info[cellline_info$Tissue_type == tg,]$Cell_line], tg)
comp_target <- data.frame('ratio_target' = target_stat$row_ones_ratio,
                          'ratio_non_target' = non_target_stat$row_ones_ratio)
rownames(comp_target) <- rownames(target_stat)

fit <- lm(ratio_non_target ~ ratio_target, data = comp_target)
summary(fit)$r.squared
slope <- coef(fit)["ratio_target"]
intercept <- coef(fit)["(Intercept)"]
cut_itt <- 0.3

x_min <- (cut_itt-intercept)/slope;x_max <- 1;y_min <- 0;y_max <- slope * x_max + intercept - cut_itt
target_na_diff <- rownames(comp_target[comp_target$ratio_target >= x_min & comp_target$ratio_target <= x_max & comp_target$ratio_non_target >= y_min & comp_target$ratio_non_target <= y_max & ((comp_target$ratio_non_target - y_min) / (comp_target$ratio_target - x_min) < slope), ] )

lung_na_diff <- target_na_diff

x_min <- 0;x_max <- ((1-cut_itt-intercept)/slope);y_min <- slope * x_min + intercept + cut_itt;y_max <- 1
target_no_na_diff <- rownames(comp_target[comp_target$ratio_target >= x_min & comp_target$ratio_target <= x_max & comp_target$ratio_non_target >= y_min & comp_target$ratio_non_target <= y_max & ((comp_target$ratio_non_target - y_min) / (comp_target$ratio_target - x_min) > slope), ] )

lung_no_na_diff <- target_no_na_diff

ggplot()+
  geom_point(data = comp_target, aes(x=ratio_target, y=ratio_non_target),fill="#f1f1f1", alpha = 0.6,
             color="grey", size = 2.5, shape = 21)+
  geom_abline(intercept = intercept, slope = slope, color = 'black', linewidth = 1, linetype = "dashed")+
  geom_abline(intercept = intercept + cut_itt, slope = slope, color = '#1C3F58', linewidth = 1, linetype = "dashed")+geom_abline(intercept = intercept - cut_itt, slope = slope,color = '#A12D45', size = 1, linetype = "dashed")+
  geom_point(data = comp_target[target_na_diff,], aes(x=ratio_target, y=ratio_non_target,color='tg_unique'),alpha = 0.6,size = 2.5)+
  geom_point(data = comp_target[target_no_na_diff,], aes(x=ratio_target, y=ratio_non_target,color='nontg_unique'),alpha = 0.6, size = 2.5)+
  scale_color_manual(name = "proteins",
                     guide = "legend",
                     values = c("tg_unique" = "#A12D45",
                                "nontg_unique" = "#1C3F58")) +
  theme_bw()+theme(panel.grid=element_blank())

##Figure 3E
tg <- 'Haematopoietic and Lymphoid'

non_target_df <- cl_res[['sorted_df']][,cellline_info[cellline_info$Tissue_type != tg,]$Cell_line]
non_target_stat <- subcancer_process(non_target_df, paste0('non_',tg))
target_stat <- subcancer_process(cl_res[['sorted_df']][,cellline_info[cellline_info$Tissue_type == tg,]$Cell_line], tg)
comp_target <- data.frame('ratio_target' = target_stat$row_ones_ratio,
                          'ratio_non_target' = non_target_stat$row_ones_ratio)
rownames(comp_target) <- rownames(target_stat)

fit <- lm(ratio_non_target ~ ratio_target, data = comp_target)
summary(fit)$r.squared
slope <- coef(fit)["ratio_target"]
intercept <- coef(fit)["(Intercept)"]
cut_itt <- 0.3

x_min <- (cut_itt-intercept)/slope;x_max <- 1;y_min <- 0;y_max <- slope * x_max + intercept - cut_itt
target_na_diff <- rownames(comp_target[comp_target$ratio_target >= x_min & comp_target$ratio_target <= x_max & comp_target$ratio_non_target >= y_min & comp_target$ratio_non_target <= y_max & ((comp_target$ratio_non_target - y_min) / (comp_target$ratio_target - x_min) < slope), ] )

haema_na_diff <- target_na_diff

x_min <- 0;x_max <- ((1-cut_itt-intercept)/slope);y_min <- slope * x_min + intercept + cut_itt;y_max <- 1
target_no_na_diff <- rownames(comp_target[comp_target$ratio_target >= x_min & comp_target$ratio_target <= x_max & comp_target$ratio_non_target >= y_min & comp_target$ratio_non_target <= y_max & ((comp_target$ratio_non_target - y_min) / (comp_target$ratio_target - x_min) > slope), ] )

haema_no_na_diff <- target_no_na_diff

ggplot()+
  geom_point(data = comp_target, aes(x=ratio_target, y=ratio_non_target),fill="#f1f1f1", alpha = 0.6,
             color="grey", size = 2.5, shape = 21)+
  geom_abline(intercept = intercept, slope = slope, color = 'black', linewidth = 1, linetype = "dashed")+
  geom_abline(intercept = intercept + cut_itt, slope = slope, color = '#1C3F58', linewidth = 1, linetype = "dashed")+geom_abline(intercept = intercept - cut_itt, slope = slope,color = '#A12D45', size = 1, linetype = "dashed")+
  geom_point(data = comp_target[target_na_diff,], aes(x=ratio_target, y=ratio_non_target,color='tg_unique'),alpha = 0.6,size = 2.5)+
  geom_point(data = comp_target[target_no_na_diff,], aes(x=ratio_target, y=ratio_non_target,color='nontg_unique'),alpha = 0.6, size = 2.5)+
  scale_color_manual(name = "proteins",
                     guide = "legend",
                     values = c("tg_unique" = "#A12D45",
                                "nontg_unique" = "#1C3F58")) +
  theme_bw()+theme(panel.grid=element_blank())

##Figure 3F
go_plot(haema_na_diff)

##Figure 3G
go_plot(haema_no_na_diff)

#Figure 4-----------------------------------------------------------------------
##Figure 4B
load("./dataset/ksecreenres_new_1.rdata")
k_df <- do.call(rbind, k_lst)
k_df$kgrp <- as.factor(k_df$kgrp)
ggplot(data = k_df %>% dplyr::filter(.$kgrp %in% c(3,5,10,100)), aes(x = cut, y = acc, color = kgrp)) + 
  geom_point(size = 1) + 
  geom_line(linewidth = 1)+theme_bw()+theme(panel.grid=element_blank())

##Figure 4C
ggplot() + 
  geom_line(data = k_df%>% dplyr::filter(.$kgrp == 10), aes(x = fpr, y = tpr, color="kgrp"),size=1.5, color="#a12d45")+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",color = "grey",size = 0.8)+
  coord_cartesian(xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), expand = F)+
  theme_bw()+theme(panel.grid=element_blank())





##Figure 4D

### NA classify function
gaussian_kernel <- function(dist_matrix, sigma) {
  sim_matrix <- exp(-dist_matrix^2 / (2 * sigma^2))
  return(sim_matrix)
}

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
    sim_mtx <- gaussian_kernel(dist_mtx, sigma) %>% as.matrix()
    
    for (smp_loop in names(sort(dist_mtx[sample,])[1:(k+1)])){
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

## simulation NA type dataframe
library(progress)
load("./dataset/cclesim_natype_mtx_250426.rdata")

twotp_na_df <- na_type_df
twotp_na_df[(twotp_na_df != 'biological') & (twotp_na_df != 'technical')] = 1
twotp_na_df[twotp_na_df == 'biological'] = 0
twotp_na_df[twotp_na_df == 'technical'] = 0
twotp_na_df <- matrix(as.numeric(twotp_na_df), nrow=nrow(twotp_na_df))
rownames(twotp_na_df) <- rownames(na_type_df)
colnames(twotp_na_df) <- colnames(na_type_df)

twotp_zero <- na_type_df
twotp_zero[twotp_zero == 'biological'] = 0
twotp_zero[twotp_zero == 'technical'] = 0
twotp_zero <- matrix(as.numeric(twotp_zero), nrow=nrow(twotp_zero))
rownames(twotp_zero) <- rownames(na_type_df)
colnames(twotp_zero) <- colnames(na_type_df)

dist_matrix <- dist(t(twotp_na_df), method = 'euclidean') %>% as.matrix()

test_lst <- na_classify(twotp_na_df, twotp_zero, dist_matrix, 3, 10)

test_af <- test_lst[['natp_df']]

test_rows <- rownames(test_af)
test_cols <- colnames(test_af)
test_af[test_af != 'biological' & test_af != 'technical'] <- 1
test_af[test_af == 'biological'] <- 0
test_af[test_af == 'technical'] <- 1
test_af <- lapply(as.data.frame(test_af), as.numeric) %>% as.data.frame()
rownames(test_af) <- test_rows
colnames(test_af) <- test_cols

dist_matrix <- dist(t(twotp_na_df), method = 'binary')
dist_matrix2 <- dist(t(test_af), method = 'binary')

hclust_result <- hclust(dist_matrix)
hclust_result2 <- hclust(dist_matrix2)

ggtree(hclust_result)+#,layout = "fan")+ #%<+% cellline_info +
  #geom_tippoint(aes(col=Cancer_type), size=1)
  geom_fruit(
    data=cellline_info,
    geom=geom_tile,  
    mapping=aes(y=Cell_line, fill=Tissue_type),  
    color="black",
    width=0.13,
    stat="identity",
    offset = -2.5
  )+scale_fill_manual(values=c("#3e517a","#fdb813"))+
  theme(legend.position = "none")+
  ggtree(hclust_result2)+#,layout = "fan")+ #%<+% cellline_info +
  #geom_tippoint(aes(col=Cancer_type), size=1)
  geom_fruit(
    data=cellline_info,
    geom=geom_tile,  
    mapping=aes(y=Cell_line, fill=Tissue_type),  
    color="black",
    width=0.13,
    stat="identity",
    offset = -1.5
  )+scale_fill_manual(values=c( "#3e517a","#fdb813"))  + scale_x_reverse()


##Figure 4E
load("./dataset/ccle949_full_classNA_dismtxeculid_k10_th2p5.rdata")

cl_resdf <- cl_natp[['natp_df']]
cl_resdf[cl_resdf != 'biological' & cl_resdf != 'technical'] <- 1
cl_resdf[cl_resdf == 'biological'] <- 0
cl_resdf[cl_resdf == 'technical'] <- 0.5
pros <- rownames(cl_resdf)
smps <- colnames(cl_resdf)
cl_resdf <- lapply(as.data.frame(cl_resdf), as.numeric) %>% as.data.frame()
rownames(cl_resdf) <- pros
colnames(cl_resdf) <- smps

na_num <- sum(cl_resdf < 1)
tech_num_all <- sum(cl_resdf == 0.5)
bio_num_all <- sum(cl_resdf == 0)

natp_stat <- data.frame("na_rt" = (rowSums(cl_resdf < 1) / ncol(cl_resdf)))
natp_stat$nagroup <- cut(natp_stat$na_rt, 
                         breaks = seq(0, 1, by = 0.1), 
                         include.lowest = TRUE,
                         right = FALSE)
natp_stat$nagroup <- as.character(natp_stat$nagroup) %>% as.factor()

natp_stack <- data.frame("group" = levels(natp_stat$nagroup))

rownames(natp_stack) <- natp_stack$group
natp_stack$tech <- NA
natp_stack$bio <- NA
natp_stack$pronum <- NA
tech_num <- 0
bio_num <- 0
pro_num <- 0

for (grp in rownames(natp_stack)){
  row = rownames(natp_stat[natp_stat$nagroup == grp,])
  
  tech_num = tech_num + sum(cl_resdf[row,] == 0.5)
  bio_num = bio_num + sum(cl_resdf[row,] == 0)
  
  #tech_num = sum(cl_resdf[row,] == 0.5)
  #bio_num = sum(cl_resdf[row,] == 0)
  
  natp_stack[grp,"tech"] <- tech_num
  natp_stack[grp,"bio"] <- bio_num
  
  pro_num = pro_num + nrow(natp_stat[natp_stat$nagroup == grp,])
  #pro_num = nrow(natp_stat[natp_stat$nagroup == grp,])
  natp_stack[grp,"pronum"] <- pro_num
  
}

natp_stack$techfrac <- natp_stack$tech / tech_num_all

natp_stack$xmean <- c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95)

ggplot(natp_stack[,c("xmean",'techfrac')], aes(x=xmean, y=techfrac)) +
  geom_col(aes(x = xmean, y = techfrac), fill="#fdb813") +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, color = "#a12d45") +
  labs(x = "groups", y = "Values", title = "Tech NA Stacked Barplot") +
  #scale_y_log10()+
  #scale_fill_manual(values = c("tech" = "#fdb813", "bio" = "#1c3f58")) +
  theme_minimal()

##Figure 4F

##cellline df use uniport entry name, convert it to gene symbol if necessary
cellline_pro <- read.csv("./dataset/cellline_prosmap.csv", header = T)
cellline_pro$From <- sapply(strsplit(cellline_pro$From, "_HUMAN"), function(x) x[1])
cellline_pro <- cellline_pro[!duplicated(cellline_pro[,1]), ] 
rownames(cellline_pro) <- cellline_pro$From
map_pros <- function(x, tp){
  if(tp == "uni"){return(cellline_pro[x,]$To)}
  if(tp == "gn"){return(rownames(cellline_pro %>% dplyr::filter(.$To %in% x)))}
}
cellline_info$Haema <- ifelse(cellline_info$Tissue_type == 'Haematopoietic and Lymphoid','Haematopoietic and Lymphoid',"Other")

displays <- map_pros(c("ITGA3","ITGAV","ITGA2","LAMB1","LAMC1"), "gn")


haema_df <- cl_resdf[displays,rownames(cellline_info %>% dplyr::filter(.$Haema == 'Haematopoietic and Lymphoid'))]
nohaema_df <- cl_resdf[displays,rownames(cellline_info %>% dplyr::filter(.$Haema == 'Other'))]

display_df <- cbind(haema_df, nohaema_df)
  
my_color <- c("#1C3F58","#fbf70c","#A12D45")
cellline_info$Haema <- ifelse(cellline_info$Tissue_type == 'Haematopoietic and Lymphoid','Haematopoietic and Lymphoid',"Other")
annotation_col <- cellline_info %>% dplyr::filter(.$Cell_line %in% colnames(unique_nona_haema_df))
smp_anno <- rownames(annotation_col)
annotation_col <- data.frame(Haema = annotation_col[,'Haema'])
rownames(annotation_col) <- smp_anno

annotation_colors <- list(
  Haema = c("Haematopoietic and Lymphoid" = "#7b68ee", Other = "#2e8b57")
)

pheatmap(display_df, color = my_color, cluster_rows = F, cluster_cols = F,  show_colnames = F,
         show_rownames = T, annotation_col = annotation_col, annotation_colors = annotation_colors)

##Figure 4G
displays <- map_pros(c("INPP5D","PIK3R1","PTPN6","PTPRC","HCLS1","RASAL3"), "gn")
haema_df <- cl_resdf[displays,rownames(cellline_info %>% dplyr::filter(.$Haema == 'Haematopoietic and Lymphoid'))]
nohaema_df <- cl_resdf[displays,rownames(cellline_info %>% dplyr::filter(.$Haema == 'Other'))]

display_df <- cbind(haema_df, nohaema_df)

pheatmap(display_df, color = my_color, cluster_rows = F, cluster_cols = F,  show_colnames = F,
         show_rownames = T, annotation_col = annotation_col, annotation_colors = annotation_colors)

##Figure 5B
load("./dataset/cellline949_haema_resdf_corr.rdata")
ggplot(haema_type_resdf, aes(x=IN_PPI, y=weighted_rho, fill=IN_PPI))+
  geom_violin(width=1)+
  geom_boxplot(fill='white', width=0.1)+
  coord_cartesian(ylim = c(-0.5, 1))+
  scale_fill_manual(values = c("IN" = "#A12D45", "OUT" = "#1C3F58"))+
  theme_bw()+theme(panel.grid=element_blank())

##Figure 5C
haema_ft <- haema_type_resdf %>% dplyr::filter(.$IN_PPI == 'IN')
ggplot()+
  geom_point(data=haema_ft %>% dplyr::filter(( abs(.$weighted_rho) > abs(.$raw_rho)) & (sign(weighted_rho) == sign(raw_rho))), aes(x=raw_rho, y=weighted_rho, color='large_same'), alpha = 0.5)+
  geom_point(data=haema_ft %>% dplyr::filter(( abs(.$weighted_rho) < abs(.$raw_rho)) & (sign(weighted_rho) == sign(raw_rho))), aes(x=raw_rho, y=weighted_rho, color='small_same'), alpha = 0.5)+  #geom_point(data=haema_ft %>% dplyr::filter((.$weighted_rho < .$raw_rho) & (.$weighted_rho > 0  & .$raw_rho > 0)), aes(x=raw_rho, y=weighted_rho, color='small_pos'), alpha = 0.5)+
  geom_abline(color = 'black', linewidth = 1, linetype = "dashed")+
  scale_color_manual(name = "type",
                     guide = "legend",
                     values = c("large_same" = "#A12D45",
                                "small_same" = "#1C3F58"))+
  coord_cartesian(xlim=c(-1,1),ylim = c(-1, 1))+
  theme_bw()+theme(panel.grid=element_blank())

##Figure 5D-E
library(ggcorrplot)
cplx_sub <- c("EXOC1","EXOC2","EXOC3","EXOC4","EXOC7","EXOC5","EXOC6","EXOC8")
cplx_grid <- expand.grid(cplx_sub,cplx_sub)
colnames(cplx_grid) <- c("protein_1","protein_2")
cplx_grid$Pair <- paste0(cplx_grid$protein_1, " ", cplx_grid$protein_2)

cplx_grid <- merge(cplx_grid, haema_type_resdf, by = 'Pair')
cplx_grid[,c("protein_1.x","protein_2.x","raw_rho")]
cplx_mtx <- spread(cplx_grid[,c("protein_1.x","protein_2.x","raw_rho")] , key = "protein_1.x",
                   value = "raw_rho")
rownames(cplx_mtx) <- cplx_mtx[[1]]
cplx_mtx <- cplx_mtx[,-1]

cplx_mtx_wei <- spread(cplx_grid[,c("protein_1.x","protein_2.x","weighted_rho")] , key = "protein_1.x",
                       value = "weighted_rho")
rownames(cplx_mtx_wei) <- cplx_mtx_wei[[1]]
cplx_mtx_wei <- cplx_mtx_wei[,-1]

for(i in 1:nrow(cplx_mtx)){
  for(j in 1:nrow(cplx_mtx)){
    if(j < i){cplx_mtx_wei[i,j] <- cplx_mtx[i,j]}
  }
}
ggcorrplot(cplx_mtx_wei, colors = c("#1C3F58","white" ,"#A12D45"),
           legend.title = "complex_weighted_rho",ggtheme = ggplot2::theme_bw())

##Figure 6A
ft_ev <- read.csv("./dataset/ev/filtered_ev_exp.csv", row.names = 1)
evinfo_ft <- read.csv("./dataset/ev/evinfo_ft.csv", row.names = 1)

ev_raw_na <- process_na(ft_ev)
pheatmap(ev_raw_na[['sorted_df']], color = c("#1c3f58","#A12D45"), cluster_rows = FALSE, cluster_cols = FALSE,  show_colnames = FALSE, show_rownames = FALSE, use_raster = F)

##Figure 6B
ggplot(ev_raw_na[['row_na_stat']], aes(x = row_ones_ratio))+
  geom_histogram(fill = '#1C3F58', alpha = 0.7)+theme_bw()+theme(panel.grid=element_blank())

##Figure 6C
load("./dataset/ev_full_classNA_dismtxeculid_k5_th5.rdata")
cl_resdf <- cl_natp[['natp_df']]
cl_resdf[cl_resdf != 'biological' & cl_resdf != 'technical'] <- 1
cl_resdf[cl_resdf == 'biological'] <- 0
cl_resdf[cl_resdf == 'technical'] <- 0.5
pros <- rownames(cl_resdf)
smps <- colnames(cl_resdf)
cl_resdf <- lapply(as.data.frame(cl_resdf), as.numeric) %>% as.data.frame()
rownames(cl_resdf) <- pros
colnames(cl_resdf) <- smps

cl_resdf <- cl_resdf[rowSums(cl_resdf < 1) < ncol(cl_resdf)*0.9,]
cl_resdf <- cl_resdf[order(rowSums(cl_resdf < 1), decreasing = TRUE), ]

pheatmap(cl_resdf, color = c("#1c3f58","#fbf70c","#A12D45"), cluster_rows = FALSE, cluster_cols = FALSE,  show_colnames = FALSE,
         show_rownames = FALSE,use_raster = F)

##Figure 6D

dist_matrix <- dist(t(cl_resdf), method = 'binary')
hclust_result <- hclust(dist_matrix)

ggtree(hclust_result,layout = "fan")+ 
  geom_fruit(
    data=evinfo_ft,
    geom=geom_tile,  
    mapping=aes(y=file, fill=Label),  
    color="black",
    width=0.2,
    stat="identity",
    offset = 0.3
  )+scale_fill_manual(values=c("#fdb813","#3e517a"))

##Figure 6E
lung_newmks_df <- read.csv("./dataset/lung_newfind_markers.csv", row.names = 1)

my_color <- colorRampPalette(c("#1C3F58","#fdb813","#A12D45"))(100)
annotation_col <- data.frame(Group = ev_info_subc$subcat)
rownames(annotation_col) <- ev_info_subc$file
annotation_colors <- list(
  Group = c(lung_c = "#6a0573", lung_cn = "#8e44ad", lung_fn = "#b285d3")
)
pheatmap(lung_newmks_df, annotation_col = annotation_col, annotation_colors = annotation_colors, color = my_color, cluster_rows = FALSE, cluster_cols = F,  show_colnames = F,
              show_rownames = T,use_raster = F)

##Figure 6F
inhouse_evdf <- read.csv("./dataset/ev/exosome_evosepLFQ_20210528.csv", header = T)
process_gn_column <- function(df, column_name) {
  df %>%
    filter(str_detect(!!sym(column_name), "GN=")) %>%
    mutate(!!column_name := str_extract(!!sym(column_name), "(?<=GN=)[^ ]+"))
}

inhouse_evdf <- process_gn_column(inhouse_evdf,"Description")
rownames(inhouse_evdf) <- inhouse_evdf[[1]]
inhouse_evdf <- inhouse_evdf[,-1]
inhouse_evdf <- inhouse_evdf[rowSums(is.na(inhouse_evdf)) != ncol(inhouse_evdf),]

ih_info <- data.frame("Sample" = colnames(inhouse_evdf),
                      "Type" = c(rep("benigh", 10),rep("cancer",10)))
rownames(ih_info) <- ih_info$Sample

protein_violin <- function(tg){
  inhouse_tg <- inhouse_evdf[tg,] %>% t() %>% as.data.frame() %>% log1p()
  inhouse_tg$group <- ih_info$Type
  
  ggplot(data=inhouse_tg, aes(x=group, y=get(tg)))+
    geom_violin(aes(fill=group),trim = F)+
    geom_boxplot(aes(), width=0.2, fill="white")+
    scale_fill_manual(values = c("benigh"="#1c3f58",
                                 "cancer"="#a12d45")) +
    theme_bw()+theme(panel.grid=element_blank())
}



protein_violin('RPL8')

protein_violin('ISG15')
