#########################################################################################
# PREPARING FOR A SHINY APP 
# 
# Color,   epidem/notepidem,    final names,     
# 
# Beatriz Urda García 2020
#########################################################################################


library(stringr)
library(SummarizedExperiment)   # NEW
# library("backbone")
library(disparityfilter)   # there is a BUG in the package # https://stackoverflow.com/questions/43387316/igraph-and-disparityfilter-package-issues-with-characters-and-large-numbers 
library(igraph)
library(stringr)
source("../analysis_library.R")


setwd("~/Desktop/ANALYSIS/shiny_network_app/")

#################
# ADD BACKBONE TO THE NETWORKS
dsn<-read.csv2("final_pairwise_union_spearman_distance_sDEGs_network.txt",stringsAsFactors = F,sep="\t",header=T)
ssn<-read.csv2("final_metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt",stringsAsFactors = F,sep="\t",header=T)

# Open the backbones
metric_dsn = read.csv2("additional_data/backbones/dsn_metric.txt",stringsAsFactors = F,sep="\t",header=T)
ultrametric_dsn = read.csv2("additional_data/backbones/dsn_ultrametric.txt",stringsAsFactors = F,sep="\t",header=T)

metric_ssn = read.csv2("additional_data/backbones/ssn_metric.txt",stringsAsFactors = F,sep="\t",header=T)
ultrametric_ssn = read.csv2("additional_data/backbones/ssn_ultrametric.txt",stringsAsFactors = F,sep="\t",header=T)

# Select the interactions that belong to the backbone
metric_dsn = metric_dsn[metric_dsn$is_metric == 'True', ]; colnames(metric_dsn)[1:2] = c("Dis1", "Dis2")
ultrametric_dsn = ultrametric_dsn[ultrametric_dsn$is_ultrametric == 'True', ]; colnames(ultrametric_dsn)[1:2] = c("Dis1", "Dis2")
metric_ssn = metric_ssn[metric_ssn$is_metric == 'True', ]; colnames(metric_ssn)[1:2] = c("Dis1", "Dis2")
ultrametric_ssn = ultrametric_ssn[ultrametric_ssn$is_ultrametric == 'True', ]; colnames(ultrametric_ssn)[1:2] = c("Dis1", "Dis2")

metric_dsn = sort_icd_interactions(metric_dsn); metric_dsn$pair = paste(metric_dsn$Dis1, metric_dsn$Dis2, sep="_")
ultrametric_dsn = sort_icd_interactions(ultrametric_dsn); ultrametric_dsn$pair = paste(ultrametric_dsn$Dis1, ultrametric_dsn$Dis2, sep="_")
metric_ssn = sort_icd_interactions(metric_ssn); metric_ssn$pair = paste(metric_ssn$Dis1, metric_ssn$Dis2, sep="_")
ultrametric_ssn = sort_icd_interactions(ultrametric_ssn); ultrametric_ssn$pair = paste(ultrametric_ssn$Dis1, ultrametric_ssn$Dis2, sep="_")

metric_dsn$pair = paste(metric_dsn$Dis1, metric_dsn$Dis2, sep="_")
ultrametric_dsn$pair = paste(ultrametric_dsn$Dis1, ultrametric_dsn$Dis2, sep="_")
metric_ssn$pair = paste(metric_ssn$Dis1, metric_ssn$Dis2, sep="_")
ultrametric_ssn$pair = paste(ultrametric_ssn$Dis1, ultrametric_ssn$Dis2, sep="_")

# Add backbone to DSN
dsn$is_metric = rep(FALSE, nrow(dsn))
dsn$is_ultrametric = rep(FALSE, nrow(dsn))
for(k in 1:nrow(dsn)){
  # k = 1
  cint = paste(sort(c(dsn$Dis1[k], dsn$Dis2[k])), collapse ="_")
  if(cint %in% metric_dsn$pair){
    dsn$is_metric[k] = TRUE
  }
  if(cint %in% ultrametric_dsn$pair){
    dsn$is_ultrametric[k] = TRUE
  }
}

# Add backbone to SSN
ssn$is_metric = rep(FALSE, nrow(ssn))
ssn$is_ultrametric = rep(FALSE, nrow(ssn))
for(k in 1:nrow(ssn)){
  # k = 1
  cint = paste(sort(c(ssn$Dis1[k], ssn$Dis2[k])), collapse ="_")
  if(cint %in% metric_ssn$pair){
    ssn$is_metric[k] = TRUE
  }
  if(cint %in% ultrametric_ssn$pair){
    ssn$is_ultrametric[k] = TRUE
  }
}

# Checking that it is okay

# DSN
# Metric: 632
# Ultrametric: 88
table(dsn$is_metric)
table(dsn$is_ultrametric)

# SSN
# Metric: 4056
# Ultrametric: 331
table(ssn$is_metric)
table(ssn$is_ultrametric)

# Saving the new final networks
dsn_filename = "pairwise_union_spearman_distance_sDEGs_network.txt"
ssn_filename = "metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt"
write.table(dsn, file=paste0("final_",dsn_filename), quote=FALSE,sep="\t", row.names = FALSE)
write.table(ssn, file=paste0("final_",ssn_filename), quote=FALSE,sep="\t", row.names = FALSE)

#####################
# ADD SHARED GENES AND PATHWAYS TO THE NETWORKS
dsn<-read.csv2("final_pairwise_union_spearman_distance_sDEGs_network.txt",stringsAsFactors = F,sep="\t",header=T)
ssn<-read.csv2("final_metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt",stringsAsFactors = F,sep="\t",header=T)
extended_node_names = read.csv2("additional_data/extended_node_names.txt",stringsAsFactors = F,sep="\t",header=T)

nodes_dsn = union(dsn$Dis1, dsn$Dis2); length(nodes_dsn) # 45
nodes_ssn = union(ssn$Dis1, ssn$Dis2); length(nodes_ssn) # 161

genelist = read.csv2("additional_data/gene_list_df.txt",stringsAsFactors = F,sep="\t",header=T)

add_genes_pathways <- function(dsn){
  # dsn = ssn  # to test
  all_pathways = c()
  dsn$genes_up = rep(NA, nrow(dsn)); dsn$genes_down = rep(NA, nrow(dsn))
  dsn$pathways_up = rep(NA, nrow(dsn)); dsn$pathways_down = rep(NA, nrow(dsn))
  
  for(k in 1:nrow(dsn)){
    # k = 1; k=16
    dis1 = str_trim(dsn$Dis1[k]); dis2 = str_trim(dsn$Dis2[k])
    
    # Get old names
    odis1 = extended_node_names[extended_node_names$final_names == dis1, ]$old_names
    odis2 = extended_node_names[extended_node_names$final_names == dis2, ]$old_names
    
    # Get FC table for each node
    genes1 = read.csv2(paste0("additional_data/Genes/",odis1,"_DEGs.txt"),stringsAsFactors = F,sep="\t",header=T)
    genes2 = read.csv2(paste0("additional_data/Genes/",odis2,"_DEGs.txt"),stringsAsFactors = F,sep="\t",header=T)
    # Select sDEGs
    genes1 = genes1[genes1$adj.p.value <= 0.05, ]; genes2 = genes2[genes2$adj.p.value <= 0.05, ]; nrow(genes1); nrow(genes2)
    genes1up=genes1[genes1$logFC > 0, ]$Gene.Symbol; genes1down=genes1[genes1$logFC < 0, ]$Gene.Symbol; 
    genes2up=genes2[genes2$logFC > 0, ]$Gene.Symbol; genes2down=genes2[genes2$logFC < 0, ]$Gene.Symbol;
    
    # Get Pathway table for each node
    paths1 = read.csv2(paste0("additional_data/FE/",odis1,"_pathways.txt"),stringsAsFactors = F,sep="\t",header=T)
    paths2 = read.csv2(paste0("additional_data/FE/",odis2,"_pathways.txt"),stringsAsFactors = F,sep="\t",header=T)
    all_pathways = union(all_pathways, union(paths1$NAME, paths2$NAME)); # print(length(all_pathways))
    
    # Select significantly altered pathways
    paths1 = paths1[paths1$FDR.q.val <= 0.05, ]; paths2 = paths2[paths2$FDR.q.val <= 0.05, ]; nrow(paths1); nrow(paths2)
    paths1up=paths1[paths1$NES > 0, ]$NAME; paths1down=paths1[paths1$NES < 0, ]$NAME; 
    paths2up=paths2[paths2$NES > 0, ]$NAME; paths2down=paths2[paths2$NES < 0, ]$NAME;
    
    if(dsn$Distance[k] >= 0){
      # Genes
      dsn$genes_up[k] = list(intersect(genes1up, genes2up)) # up_up
      dsn$genes_down[k] = list(intersect(genes1down, genes2down)) # down_down
      # Pathways
      dsn$pathways_up[k] = list(intersect(paths1up, paths2up)) # up_up
      dsn$pathways_down[k] = list(intersect(paths1down, paths2down)) # down_down
    }else{
      # Genes
      dsn$genes_up[k] = list(intersect(genes1up, genes2up)) # up_down
      dsn$genes_down[k] = list(intersect(genes1down, genes2up)) # down_up
      # Pathways
      dsn$pathways_up[k] = list(intersect(paths1up, paths2down)) # up_up
      dsn$pathways_down[k] = list(intersect(paths1down, paths2up)) # down_down
    }
  }
  
  return(list(dsn, all_pathways))

}

# Add genes and pathways to the DSN
outp_dsn = add_genes_pathways(dsn)
fdsn = outp_dsn[[1]]; all_pathways = outp_dsn[[2]]

# Add genes and pathways to the SSN
outp_ssn = add_genes_pathways(ssn)
fssn = outp_ssn[[1]]; all_pathways_ssn = outp_ssn[[2]]

# Add an id to the interactions
fdsn$id = 1:nrow(fdsn)
fssn$id = 1:nrow(fssn)

dsn_filename = "pairwise_union_spearman_distance_sDEGs_network.rds"
ssn_filename = "metap_dis_pairwise_union_spearman_distance_sDEGs_network.rds"
saveRDS(fdsn, file=paste0("genes_paths_final_",dsn_filename))
saveRDS(fssn, file=paste0("genes_paths_final_",ssn_filename))

dsn_filename = "pairwise_union_spearman_distance_sDEGs_network.txt"
ssn_filename = "metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt"
# Saving the networks with the id of the pairs
to_save_dsn = subset(fdsn, select = -c(genes_up, genes_down, pathways_up, pathways_down))
to_save_ssn = subset(fssn, select = -c(genes_up, genes_down, pathways_up, pathways_down))
write.table(to_save_dsn, file=paste0("final_",dsn_filename), quote=FALSE,sep="\t", row.names = FALSE)
write.table(to_save_ssn, file=paste0("final_",ssn_filename), quote=FALSE,sep="\t", row.names = FALSE)

all_pathways = union(all_pathways, all_pathways_ssn); length(all_pathways) # 4997
saveRDS(all_pathways, "additional_data/all_pathways.rds")

# CREATE GENES DICTIONARY
create_genes_pathways_dictionaries <- function(fdsn){
  # Creates genes and paths df for a given network (DSN or SSN)
  genesdf = c() # For each gene, the interactions that share the gene up or down
  for(gene in genelist$ensemble){
    # gene = genelist$ensemble[1]; gene = "ENSG00000163453" # to test
    cup = c()
    cdown = c()
    for(k in 1:nrow(fdsn)){
      if(gene %in% fdsn$genes_up[[k]]){
        cup = append(cup, fdsn$id[k])
      }
      if(gene %in% fdsn$genes_down[[k]]){
        cdown = append(cdown, fdsn$id[k])    # cup por cdown?????
      }
    }
    genesdf[[gene]]$up = cup
    genesdf[[gene]]$down = cdown
  }
  
  # CREATE PATHWAYS DICTIONARY
  pathwaysdf = c() # For each pathway, the interactions that share the pathway up or down
  for(path in all_pathways){
    # gene = genelist$ensemble[1]; gene = "ENSG00000163453" # to test
    cup = c()
    cdown = c()
    for(k in 1:nrow(fdsn)){
      if(path %in% fdsn$pathways_up[[k]]){
        cup = append(cup, fdsn$id[k])
      }
      if(path %in% fdsn$pathways_down[[k]]){
        cdown = append(cdown, fdsn$id[k])       # cup por cdown?????
        
        
      }
    }
    pathwaysdf[[path]]$up = cup
    pathwaysdf[[path]]$down = cdowns
  }
  
  return(list("genesdf" = genesdf, "pathwaysdf" = pathwaysdf))
  
}

dics_dsn = create_genes_pathways_dictionaries(fdsn)
dics_ssn = create_genes_pathways_dictionaries(fssn)


saveRDS(dics_dsn, "additional_data/DSN_genes_pathways_dics.rds")
saveRDS(dics_ssn, "additional_data/SSN_genes_pathways_dics.rds")


# EXAMPLES
dics_fdsn$genesdf[["ENSG00000104093"]]$down
dics_fdsn$genesdf[["ENSG00000104093"]]


# EXAMPLES
dics_fssn$pathwaysdf[["REACTOME_SMOOTH_MUSCLE_CONTRACTION"]]
dics_fssn$pathwaysdf[["REACTOME_SMOOTH_MUSCLE_CONTRACTION"]]$down

# Change pathway names to lowercase

# # Put it the other way around --> as data frame
# genesdf = data.frame(matrix(NA, nrow = 0, ncol=3))
# for(gene in genelist$ensemble){
#   gene = genelist$ensemble[1]
#   gene = "ENSG00000163453"
#   cup = c()
#   cdown = c()
#   for(k in 1:nrow(fdsn)){
#     if(gene %in% fdsn$genes_up[[k]]){
#       cup = append(cup, fdsn$id[k])
#     }
#     if(gene %in% fdsn$genes_down[[k]]){
#       cup = append(cup, fdsn$id[k])
#     }
#     nrow = data.frame("Gene"=gene, "Up"=list(cup), "Down"=list(cdown)); nrow
#     genesdf = rbind(genesdf, data.frame(gene, cup, cdown))
#   }
#   
# }



#####################
# GET COMPLETE LIST OF GENES
se <- readRDS("../GREIN_SE_good_quality/BipolarDisorder_grein_se.rds")
genes = rownames(assays(se)$counts)
length(genes)  # 28089

ensemble = elementMetadata(se)$ensemble; length(ensemble)
gene_symbol = elementMetadata(se)$gene_symbol; length(gene_symbol)
genedf = data.frame("ensemble" = ensemble, "gene_symbol" = gene_symbol)
write.table(genedf, file="additional_data/gene_list_df.txt", sep="\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

### SAVE A FILE WITH THE REACTOME PATHWAY CATEGORIES
pcats<-read.csv2("../Reactome/Reactome_parents.txt",stringsAsFactors = F,sep="\t",header=T)
reactome_cats <- unique(pcats$Parent); length(reactome_cats) # 27
write.table(reactome_cats, file="additional_data/reactome_pathway_categories.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
colnames(pcats) = c("pathway","category")
write.table(pcats, file="additional_data/pcats.txt", sep="\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

### PREPARE LIST WITH THE GENES AND PATHWAYS OF THE DISEASES AND META-PATIENTS
library(readr)

meta <- read.csv("new_disease_metadata_final_names.txt",header=T, sep="\t",stringsAsFactors = F); dim(meta)
ndis = length(meta$disease_name)

ssn <- read.csv("final_metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt",header=T, sep="\t",stringsAsFactors = F)
unique_nodes <- str_trim(unique(union(ssn$Dis1, ssn$Dis2)))
nunique_nodes <- length(unique_nodes); nunique_nodes # 161

setClass("disease", slots = list(genes_up="character", genes_down="character", 
                                 pathways_up="character", pathways_down="character"))
diseases_dic = vector(mode = "list")

extended_node_names = c()


for (k in 1:nunique_nodes){
  print(k)
  # k = 2 # To comment
  if(grepl("\\d+", unique_nodes[k])){     # META-PATIENTS
    is_metap=TRUE
    # Get metap number
    cnumber = parse_number(gsub("-","",unique_nodes[k]))
    # Get metap name (old one)
    cdis = str_trim(gsub(cnumber,"",unique_nodes[k]))
    cdis = meta[meta$final_disease_name == cdis, ]$disease_name
    
    #Get genes
    cdir = paste0("../Metapatients/with_entire_count_matrix/DEAnalysis/",cdis,"/")
    cdir2 = list.files(cdir,paste0(cdis,"_DEGs_",".*_",cnumber,"\\.txt"))
    cdf = read.csv(paste0(cdir,cdir2),header=T, sep="\t",stringsAsFactors = F); dim(cdf)
    outputpath_g = paste0("aditional_data/Genes/",cdis,"_",cnumber,"_DEGs.txt")
    
    csign <- dim(cdf[cdf$adj.P.Val < 0.05, ])[1]
    extended_node_names <- rbind(extended_node_names, c(unique_nodes[k], paste(cdis,cnumber,sep="_"), csign))
    # print(tail(extended_node_names,1))
    
    
  }else{   # DISEASES
    is_metap=FALSE
    cdis = meta[meta$final_disease_name == unique_nodes[k], ]$disease_name
    #Get genes
    cdir = paste0("../DL_RESULTS_GREIN/",cdis,"/")
    cdf = read.csv(paste0(cdir,cdis,"_DEGs.txt"),header=T, sep="\t",stringsAsFactors = F); dim(cdf)
    outputpath_g = paste0("aditional_data/Genes/",cdis,"_DEGs.txt")
    
    csign <- dim(cdf[cdf$adj.P.Val < 0.05, ])[1]
    extended_node_names <- rbind(extended_node_names, c(unique_nodes[k], cdis, csign))
    # print(tail(extended_node_names,1))
  } 
  
  # Format genes
  cdf_tosafe = cdf[, c("symbol","logFC","adj.P.Val")]; colnames(cdf_tosafe) <- c("Gene Symbol","logFC","adj p-value")
  write.table(cdf_tosafe, file = outputpath_g, sep="\t", row.names = FALSE, quote = FALSE)
  
  cdf <- cdf[cdf$adj.P.Val <=0.05, ]
  genes_up=cdf[cdf$logFC>0, ]$symbol
  genes_down=cdf[cdf$logFC<0, ]$symbol
  
  # Get pathways
  if(is_metap == TRUE){ # 
    cdir = paste0("../Metapatients/with_entire_count_matrix/FE_results/")
    if(cdis == "ThyroidCancer_Papillary"){
      cdisdir = list.files(cdir,paste0('ThyroidCancerPapillary',"_.*_",cnumber,".*"), include.dirs=TRUE)[1]; cdisdir
      outputpath_p = paste0("aditional_data/FE/",'ThyroidCancer_Papillary',"_",cnumber,"_pathways.txt")
    }else{
      cdisdir = list.files(cdir,paste0(cdis,"_.*_",cnumber,".*"), include.dirs=TRUE)[1]; cdisdir
      outputpath_p = paste0("aditional_data/FE/",cdis,"_",cnumber,"_pathways.txt")
    }
    dir_up = list.files(paste0(cdir,cdisdir),"gsea_report_for_na_pos_.*\\.xls"); dir_up
    dir_down = list.files(paste0(cdir,cdisdir),"gsea_report_for_na_neg_.*\\.xls"); dir_down
    
    df_up = read.csv(paste0(cdir,cdisdir,"/",dir_up),header=T, sep="\t",stringsAsFactors = F); head(df_up)
    df_down = read.csv(paste0(cdir,cdisdir,"/",dir_down),header=T, sep="\t",stringsAsFactors = F); head(df_down)
    
  }else{ #
    cdir = paste0("../FE_results/Pathways_DEGs_diseases/")
    cdisdir = list.files(cdir,paste0(cdis,"_.+"), include.dirs=TRUE)[1]; cdisdir
    dir_up = list.files(paste0(cdir,cdisdir),"gsea_report_for_na_pos_.*\\.xls"); dir_up
    dir_down = list.files(paste0(cdir,cdisdir),"gsea_report_for_na_neg_.*\\.xls"); dir_down
    
    df_up = read.csv(paste0(cdir,cdisdir,"/",dir_up),header=T, sep="\t",stringsAsFactors = F); head(df_up)
    df_down = read.csv(paste0(cdir,cdisdir,"/",dir_down),header=T, sep="\t",stringsAsFactors = F); head(df_down)
    
    outputpath_p = paste0("aditional_data/FE/",cdis,"_pathways.txt")
    
  }
  # print("Here")
  df_all = rbind(df_up, df_down)
  df_all = df_all[, c("NAME","NES","FDR.q.val")]
  write.table(df_all, file = outputpath_p, sep="\t", row.names = FALSE, quote = FALSE)
  
  df_up = df_up[df_up$FDR.q.val <= 0.05, ]$NAME
  df_down = df_down[df_down$FDR.q.val <= 0.05, ]$NAME
  
  
  obj <- new("disease", 
             genes_up=genes_up, genes_down=genes_down, 
             pathways_up=df_up, pathways_down=df_down)
  
  cfinal_name <- unique_nodes[k]
  diseases_dic[[cfinal_name]] <- obj
  # print(k)
  # print("done")
}

saveRDS(diseases_dic, file="disease_dic.rds")
extended_node_names <- as.data.frame(extended_node_names); colnames(extended_node_names) <- c("final_names", "old_names", "sDEGs")
write.table(extended_node_names, file="aditional_data/extended_node_names.txt", sep="\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



### PREPARE THE NETWORKS WITH ALL THE NEEDED INFORMATION

##### Select network based on icd9
sort_icd_interactions <- function(df){
  # Sort interactions such that icd of Dis1 < icd of Dis2
  
  for (k in 1:length(df$Dis1)){
    if (df$Dis1[k] > df$Dis2[k]){
      dis1 <- df$Dis1[k]
      df$Dis1[k] <- df$Dis2[k]
      df$Dis2[k] <- dis1
    }
  }
  return(df)
}


dsn_filename = "pairwise_union_spearman_distance_sDEGs_network.txt"
ssn_filename = "metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt"
dsn <- read.csv(dsn_filename,header=T, sep="\t",stringsAsFactors = F)
ssn <- read.csv(ssn_filename,header=T, sep="\t",stringsAsFactors = F)


meta <- read.csv("new_disease_metadata_final_names.txt",header=T, sep="\t",stringsAsFactors = F); dim(meta)
discats <- unique(meta$disease_cat)
dis_icd <- meta[, c("disease_name", "icd9", "disease_cat", "new_dis_cat_colors", "final_disease_name")]

# Add disease names to ssn
ssn$Corr_Dis1 <- gsub("_\\d+","",ssn$Dis1)
ssn$Corr_Dis2 <- gsub("_\\d+","",ssn$Dis2)

# COLOR + CORRECT NAMES
# DSN
dim(dsn)
dsn <- merge(dsn, dis_icd, by.x = "Dis1", by.y = "disease_name", all.x = TRUE, all.y = FALSE); dim(dsn)
colnames(dsn)[6:9] <- c("Dis1_icd9", "Dis1_discat", "Dis1_catcolor", "fDis1_name") 
dsn <- merge(dsn, dis_icd, by.x = "Dis2", by.y = "disease_name", all.x = TRUE, all.y = FALSE); dim(dsn)
colnames(dsn)[10:13] <- c("Dis2_icd9", "Dis2_discat", "Dis2_catcolor", "fDis2_name")

# SSN
dim(ssn)
ssn <- merge(ssn, dis_icd, by.x = "Corr_Dis1", by.y = "disease_name", all.x = TRUE, all.y = FALSE); dim(ssn)
colnames(ssn)[8:11] <- c("Dis1_icd9", "Dis1_discat", "Dis1_catcolor", "fDis1_name") 
ssn <- merge(ssn, dis_icd, by.x = "Corr_Dis2", by.y = "disease_name", all.x = TRUE, all.y = FALSE); dim(ssn)
colnames(ssn)[12:15] <- c("Dis2_icd9", "Dis2_discat", "Dis2_catcolor", "fDis2_name")

library(readr)
n_dis1 <- parse_number(ssn$Dis1); n_dis1[is.na(n_dis1)] <- ""
n_dis2 <- parse_number(ssn$Dis2); n_dis2[is.na(n_dis2)] <- ""

ssn$fDis1_name <- paste(ssn$fDis1_name, n_dis1, sep=" ")
ssn$fDis2_name <- paste(ssn$fDis2_name, n_dis2, sep=" ")

# EPIDEM / NOT EPIDEM
barabasi <- read.csv("aditional_data/PDN_3_digits.net",header=F, sep="\t",stringsAsFactors = F)
colnames(barabasi) <- c('Dis1','Dis2','Prev1','Prev2','Co-ocurrence','RR','RR_99perc_left_bound','RR_99perc_right_bound','phi_corr','p_value')
barabasi <- barabasi[,c('Dis1','Dis2','RR_99perc_left_bound')]
barabasi <- barabasi[which(barabasi$RR_99perc_left_bound > 1 ), ]
dim(barabasi); str(barabasi)
barabasi <- sort_icd_interactions(barabasi)
barabasi$interactions <- paste(barabasi$Dis1, barabasi$Dis2, sep="_")

174 %in% barabasi$Dis1; 174 %in% barabasi$Dis2;
common_icds <- intersect(union(barabasi$Dis1, barabasi$Dis2), union(dsn$Dis1_icd9, dsn$Dis2_icd9)); length(common_icds)
b_barab <- barabasi[barabasi$Dis1 == 174 | barabasi$Dis2 == 174, ]
with_b_barab <- union(b_barab$Dis1, b_barab$Dis2)

b_dsn <- dsn[dsn$Dis1_icd9 == 174 | dsn$Dis2_icd9 == 174, ]
with_b_dsn <- union(b_dsn$Dis1_icd9, b_dsn$Dis2_icd9)
intersect(with_b_dsn, with_b_barab)

b_ssn <- ssn[ssn$Dis1_icd9 == 174 | ssn$Dis2_icd9 == 174, ]
with_b_ssn <- union(b_ssn$Dis1_icd9, b_ssn$Dis2_icd9)
intersect(with_b_ssn, with_b_barab)

# DSN
cunique_int <- c()
dsn$in_epidem <- FALSE
for(k in 1:length(dsn$Dis1)){
  cinteract <- sort(c(dsn$Dis1_icd9[k], dsn$Dis2_icd9[k]))
  cinteract <- paste(cinteract, collapse="_")
  if(cinteract %in% barabasi$interactions){
    dsn$in_epidem[k] <- TRUE
    cunique_int <- append(cunique_int, cinteract)
  }
}
length(unique(cunique_int)) # 235 interactions
dsn[dsn$Distance < 0, ]$in_epidem <- FALSE
dsn_in_epidem = dim(dsn[dsn$in_epidem == TRUE, ])[1]; dsn_in_epidem  # 197
dsn_size = dim(dsn[dsn$Distance>0, ])[1]; dsn_size # 417
(dsn_in_epidem/dsn_size)*100 # 47.24
# NO SALEN LOS MISMOS NÜMEROS PORQUE AQUÏ HAY INTERACCIONES DE ICD9 REPETIDAS

# SSN
# ---- all interactions: D-D, D-M, M-M
cunique_int <- c()
ssn$in_epidem <- FALSE
for(k in 1:length(ssn$Dis1)){
  cinteract <- sort(c(ssn$Dis1_icd9[k], ssn$Dis2_icd9[k]))
  cinteract <- paste(cinteract, collapse="_")
  if(cinteract %in% barabasi$interactions){
    ssn$in_epidem[k] <- TRUE
    cunique_int <- append(cunique_int, cinteract)
  }
}
length(unique(cunique_int)) # 309
ssn[ssn$Distance < 0, ]$in_epidem <- FALSE
ssn_in_epidem = dim(ssn[ssn$in_epidem == TRUE, ])[1]; ssn_in_epidem  # 2542
ssn_size = dim(ssn[ssn$Distance>0, ])[1]; ssn_size # 5869
(ssn_in_epidem/ssn_size)*100  # 43.31232

# ---- Only M-D & D-D interactions
ssn$in_epidem_dd_dm <- ssn$in_epidem

# Put False if Dis1 is M and Dis2 is M
ssn[(ssn$Corr_Dis1 != ssn$Dis1) & (ssn$Corr_Dis2 != ssn$Dis2), ]$in_epidem_dd_dm <- FALSE
# hey = ssn[(ssn$Corr_Dis1 != ssn$Dis1) & (ssn$Corr_Dis2 != ssn$Dis2), ]


### SAVE THE CORRECT NETWORKS
fdsn <- dsn[, c("fDis1_name","fDis2_name", "Dis1_icd9","Dis2_icd9", 
                "Dis1_discat","Dis2_discat","Dis1_catcolor","Dis2_catcolor",
                "Distance","pvalue","adj_pvalue",
                "in_epidem")]
colnames(fdsn)[1:2] <- c("Dis1", "Dis2")

fssn <- ssn[, c("fDis1_name","fDis2_name","Dis1_icd9","Dis2_icd9",
                "Dis1_discat","Dis2_discat","Dis1_catcolor","Dis2_catcolor",
                "Corr_Dis1","Corr_Dis2",
                "Distance","pvalue","adj_pvalue",
                "in_epidem","in_epidem_dd_dm")]
colnames(fssn)[1:2] <- c("Dis1", "Dis2")

# Add disease names to ssn
fssn$Corr_Dis1 <- str_trim(gsub(" \\d+","",fssn$Dis1))
fssn$Corr_Dis2 <- str_trim(gsub(" \\d+","",fssn$Dis2))
unique(fssn$Corr_Dis1)

write.table(fdsn, file=paste0("final_",dsn_filename), quote=FALSE,sep="\t", row.names = FALSE)
write.table(fssn, file=paste0("final_",ssn_filename), quote=FALSE,sep="\t", row.names = FALSE)

# Checking
ssn_filename = "metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt"
ssn <- read.csv(paste0("final_",ssn_filename),header=T, sep="\t",stringsAsFactors = F)



########## OLD ###########

# # DSN
# final_dsn <- read.csv("aditional_data/pairwise_union_spearman_distance_sDEGs_pos_allpos_B_TRUE_final_network.txt",header=T, sep="\t",stringsAsFactors = F)
# final_barabasi <- read.csv("aditional_data/pairwise_union_spearman_distance_sDEGs_pos_allpos_B_TRUE_final_barabasi.txt",header=T, sep="\t",stringsAsFactors = F)
# 
# final_dsng <- graph_from_data_frame(final_dsn, directed=FALSE)
# final_barabasig <- graph_from_data_frame(final_barabasi, directed=FALSE)
# intersect <- intersection(final_dsng,final_barabasig)
# gsize(intersect) # 157
# 
# intersectg <- as_data_frame(intersect)
# epidem_interacts_dd <- unique(paste(intersectg$from, intersectg$to, sep="_")); length(epidem_interacts_dd) # 157 = N Overlapping interactions
# 
# cunique_int <- c()
# dsn$in_epidem <- FALSE
# for(k in 1:length(dsn$Dis1)){
#   cinteract <- sort(c(dsn$Dis1_icd9[k], dsn$Dis2_icd9[k]))
#   cinteract <- paste(cinteract, collapse="_")
#   if(cinteract %in% epidem_interacts_dd){
#     dsn$in_epidem[k] <- TRUE
#     cunique_int <- append(cunique_int, cinteract)
#   }
# }
# length(unique(cunique_int)) # 106 interactions
# 
# # SSN
# #--- Considering all interactions: D-D, D-M, M-M
# all_barabasi <- read.csv("aditional_data/metap_dis_pairwise_union_spearman_distance_sDEGs_pos_all_TRUE_final_barabasi.txt",header=T, sep="\t",stringsAsFactors = F) 
# all_barabasi$interact <- paste(all_barabasi$Dis1, all_barabasi$Dis2, sep="_")
# epidem_interacts_ssn <- union(all_barabasi$interact, epidem_interacts_dd); length(epidem_interacts_ssn) # 522
# 
# cunique_int <- c()
# ssn$in_epidem <- FALSE
# for(k in 1:length(dsn$Dis1)){
#   cinteract <- sort(c(ssn$Dis1_icd9[k], ssn$Dis2_icd9[k]))
#   cinteract <- paste(cinteract, collapse="_")
#   if(cinteract %in% epidem_interacts_ssn){
#     ssn$in_epidem[k] <- TRUE
#     cunique_int <- append(cunique_int, cinteract)
#   }
# }
# length(unique(cunique_int)) # 77
# length(ssn[ssn$in_epidem == TRUE, ]$in_epidem) # 417
# 
# #--- Considering only: D-D and D-M interactions
# final_ssn <- read.csv("aditional_data/metap_dis_pairwise_union_spearman_distance_sDEGs_pos_metap_dis_TRUE_final_network.txt",header=T, sep="\t",stringsAsFactors = F)
# final_barab_ssn <- read.csv("aditional_data/metap_dis_pairwise_union_spearman_distance_sDEGs_pos_metap_dis_TRUE_final_barabasi.txt",header=T, sep="\t",stringsAsFactors = F)
# 
# final_ssng <- graph_from_data_frame(final_ssn, directed=FALSE)
# final_barab_ssng <- graph_from_data_frame(final_barab_ssn, directed=FALSE)
# intersect <- intersection(final_ssng,final_barab_ssng)
# gsize(intersect)
# gsize(final_ssng)
# 
# intersectg_ssn <- as_data_frame(final_ssng);  
# epidem_interacts_ssn <- paste(intersectg_ssn$from, intersectg_ssn$to, sep="_")
# epidem_interacts_all <- union(epidem_interacts, epidem_interacts_ssn); length(epidem_interacts_all) # 563
# 
# cunique_int <- c()
# ssn$in_epidem_dd_dm <- FALSE
# for(k in 1:length(dsn$Dis1)){
#   if((ssn$Corr_Dis1[k] == ssn$Dis1[k]) | (ssn$Corr_Dis2[k] == ssn$Dis2[k])){ # if D-D or D-M = si Dis1 o Dis2 son una enfermedad (no contienen número)
#     # Check if it is in epidemiology
#     cinteract <- sort(c(ssn$Dis1_icd9[k], ssn$Dis2_icd9[k]))
#     cinteract <- paste(cinteract, collapse="_")
#     if(cinteract %in% epidem_interacts_all){
#       ssn$in_epidem_dd_dm[k] <- TRUE
#       cunique_int <- append(cunique_int, cinteract)
#     }
#   }
# }
# length(unique(cunique_int)) # 67
# length(ssn[ssn$in_epidem_dd_dm == TRUE, ]$Corr_Dis2) # 200
# 
# ### SAVE THE CORRECT NETWORKS
# fdsn <- dsn[, c("Dis1","Dis2", "Dis1_icd9","Dis2_icd9", 
#                 "Dis1_discat","Dis2_discat","Dis1_catcolor","Dis2_catcolor",
#                 "Distance","pvalue","adj_pvalue",
#                 "in_epidem")]
# 
# fssn <- ssn[, c("fDis1_name","fDis2_name","Dis1_icd9","Dis2_icd9",
#                 "Dis1_discat","Dis2_discat","Dis1_catcolor","Dis2_catcolor",
#                 "Corr_Dis1","Corr_Dis2",
#                 "Distance","pvalue","adj_pvalue",
#                 "in_epidem","in_epidem_dd_dm")]
# colnames(fssn)[1:2] <- c("Dis1", "Dis2")
# 
# write.table(fdsn, file=paste0("final_",dsn_filename), quote=FALSE,sep="\t", row.names = FALSE)
# write.table(fssn, file=paste0("final_",ssn_filename), quote=FALSE,sep="\t", row.names = FALSE)


# Get network backbone  # LIBRARY WITH A BUG
# fdsn = read.csv("final_pairwise_union_spearman_distance_sDEGs_network.txt",header=T, sep="\t",stringsAsFactors = F)
# fssn = read.csv("final_metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt",header=T, sep="\t",stringsAsFactors = F)
# 
# fdsn = fdsn[fdsn$Distance > 0, c("Dis1", "Dis2", "Distance")]
# fssn = fssn[fssn$Distance > 0, c("Dis1", "Dis2", "Distance")]
# 
# colnames(fdsn)[3] = "weight"
# colnames(fssn)[3] = "weight"
# 
# fdsn$weight = fdsn$weight*100
# fssn$weight = fssn$weight*100
# 
# fdsn_graph = graph_from_data_frame(fdsn, directed = FALSE)
# fssn_graph = graph_from_data_frame(fssn, directed = FALSE)
# is_weighted(fdsn_graph) # FALSE; 
# is_weighted(fssn_graph) # FALSE
# 
# # epsilon=0.0000001
# # E(fdsn_graph)$weight = fdsn$Distance; plot(fdsn_graph)
# # E(fssn_graph)$weight = fssn$Distance; plot(fssn_graph)
# 
# # Backbone: edge represents the strenght of the interaction
# # E(fdsn_graph)$weight = sample(1:25, ecount(fdsn_graph), replace = TRUE); plot(fdsn_graph)
# bkbone_dsn = backbone(fdsn_graph); plot(bkbone_dsn)
# bkbone_ssn = backbone(fssn_graph); plot(bkbone_ssn)
# plot(backbone(fdsn_graph))
# 
# 
# g <- sample_pa(n = 250, m = 5, directed = FALSE)
# E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)
# plot(backbone(g))
# 
# g <- sample_pa(n = 250, m = 5, directed = FALSE)
# E(g)$weight <- rep(0.5,ecount(g))
# plot(backbone(g))
