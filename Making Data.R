
############### Set Repositories################################################
setRepositories(ind = 1:7)

setwd("/disk4/bilje/2022_CI_final")
############### Load library####################################################
library(dplyr)
library(data.table)
library(stringr)
library(plyr)
############### data 불러오기 (meta)################################################################################################################

metadata1_5138 <- fread("/disk4/bilje/2022_CI_final/MGYS00005138/MGYS00005138_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
srarun_5138 <- fread("/disk4/bilje/2022_CI_final/MGYS00005138/MGYS00005138_SraRunTable.txt", header = TRUE, stringsAsFactors = FALSE)

metadata1_5184 <- fread("/disk4/bilje/2022_CI_final/MGYS00005184/MGYS00005184_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
srarun_5184 <- fread("/disk4/bilje/2022_CI_final/MGYS00005184/MGYS00005184_SraRunTable.txt", header = TRUE, stringsAsFactors = FALSE)

metadata1_5194 <- fread("/disk4/bilje/2022_CI_final/MGYS00005194/MGYS00005194_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
srarun_5194 <- fread("/disk4/bilje/2022_CI_final/MGYS00005194/MGYS00005194_SraRunTable.txt", header = TRUE, stringsAsFactors = FALSE)

# 2  metadata와 SraRunTable 합치기 - match되지 않은 자료 제거
meta_5138 <- merge(metadata1_5138, srarun_5138, by.x='run_accession', by.y='Run', all=F)
meta_5184 <- merge(metadata1_5184, srarun_5184, by.x='run_accession', by.y='Run', all=F)
meta_5194 <- merge(metadata1_5194, srarun_5194, by.x='run_accession', by.y='Run', all=F)


# 5184, 5194에 host disease 추가
meta_5184$Host_disease <- rep("Healthy", nrow(meta_5184))
meta_5194$Host_disease <- rep("Chronic Rhinosinusitis", nrow(meta_5194))

############### meta cleansing################################################################################################################
############### meta 5138########################################################

# host_disease 가 ""인 행 찾기 & 버리기
index_notdisease <- which(meta_5138$Host_disease == "")
sample_notdisease <- meta_5138$run_accession[index_notdisease]

meta_5138_update <- meta_5138[-index_notdisease,]

# NA인 열 제거
idxNa <- c()
for (i in 1:ncol(meta_5138_update)){
  idxNa[i] <- sum(is.na(meta_5138_update[,..i]))
}

idxNa_col <- which(idxNa == nrow(meta_5138_update))
meta_5138_update <- meta_5138_update[,-..idxNa_col]

# posix인 열 찾기
col_class <- matrix(lapply(meta_5138_update, function(x) paste(class(x), collapse = ',')))
idx_posix <- which(str_detect(col_class,'POSIX'))
meta_notposix <- meta_5138_update[,-..idx_posix]

# ""인 열 제거
idx_empty <- c()
for (i in 1:ncol(meta_notposix)){
  idx_empty[i] <- sum("" == meta_notposix[,..i])
}

idx_empty_col <- which(idx_empty == nrow(meta_notposix))
idx_empty_colname <- colnames(meta_notposix[,..idx_empty_col])

meta_5138_update <- meta_5138_update[,-..idx_empty_colname]

# 동일 col : `sample_collection-date` & collection_date
meta_5138_update <- meta_5138_update[,-"collection_date"]

# 동일 col : `sample_host taxid` & `sample_host-tax-id`
meta_5138_update <- meta_5138_update[,-"sample_host-tax-id"]

# 동일 col : `sample_sample-name & `Sample Name`
meta_5138_update <- meta_5138_update[,-"Sample Name"]

############### meta 5184########################################################
# 동일 col : `sample_collection date` & `sample_collection-date`
meta_5184_update <- meta_5184[,-"sample_collection date"]

# 동일 col : `sample_environment-material` & `sample_environment (material)` & "environment_(material)"
meta_5184_update <- meta_5184_update[,-"sample_environment (material)"]
meta_5184_update <- meta_5184_update[,-"environment_(material)"]

# 동일 col : `sample_host sex` & `host_sex`
meta_5184_update <- meta_5184_update[,-"host_sex"]


# 동일 col : `sample_geographic location (country and/or sea,region)` & `geo_loc_name_country` & "geographic_location_(country_and/or_sea)"
meta_5184_update <- meta_5184_update[,-"geo_loc_name_country"]
meta_5184_update <- meta_5184_update[,-"geographic_location_(country_and/or_sea)"]

# 동일 col : `sample_environmental package` & `human_gut_environmental_package`
meta_5184_update <- meta_5184_update[,-"human_gut_environmental_package"]

# 동일 col : `sample_investigation type` & `Investigation_type` 
# col 이름 변경 : library_source
meta_5184_update <- meta_5184_update[,-"Investigation_type"]
names(meta_5184_update)[names(meta_5184_update)=="sample_investigation type"] <- "library_source"

# 동일 col : `sample_longitude` & `sample_geographic location (longitude)` & "geographic_location_(longitude)"
meta_5184_update <- meta_5184_update[,-"sample_geographic location (longitude)"]
meta_5184_update <- meta_5184_update[,-"geographic_location_(longitude)"]

# 동일 col : `sample_latitude` & `sample_geographic location (latitude)` & "geographic_location_(latitude)"
meta_5184_update <- meta_5184_update[,-"sample_geographic location (latitude)"]
meta_5184_update <- meta_5184_update[,-"geographic_location_(latitude)"]

# 동일 col : `analysis_instrument-platform` & `Platform`
meta_5184_update <- meta_5184_update[,-"analysis_instrument-platform"]

############### meta 5194########################################################

# 동일 col : `Sample_Name` & `SRA_accession`
meta_5194_update <- meta_5194[,-"Sample_Name"]

# 동일 col : `sample_species` & `sample_host scientific name`
meta_5194_update <- meta_5194_update[,-"sample_host scientific name"]

############### KEY 열 형성######################################################
meta_5138_update$KEY <- meta_5138_update$analysis_accession
meta_5184_update$KEY <- meta_5184_update$run_accession
meta_5194_update$KEY <- meta_5194_update$run_accession
############### meta merge#######################################################
# 3가지의 metadata: meta_5138_update, meta_5184_update, meta_5194_update
meta_first <- merge(meta_5184_update, meta_5194_update, by = intersect(names(meta_5184_update), names(meta_5194_update)), all=T)
name <- intersect(names(meta_first), names(meta_5138_update))

name_first <- meta_first[,..name]
name_5138 <- meta_5138_update[,..name]

# class 다른 col 이름 변경
names(meta_5138_update)[names(meta_5138_update)=="Library Name"] <- "Library_Name"
names(meta_5138_update)[names(meta_5138_update)=="host_Age"] <- "host age"

name <- intersect(names(meta_first), names(meta_5138_update))
meta_final <- merge(meta_first, meta_5138_update, by = name, all=T)

# csv 파일
#write.csv(meta_final, "/disk4/bilje/2022_CI_final/meta.csv")

############### data 불러오기 (OTU)###########################################################################################################
otu_5138 <- fread("/disk4/bilje/2022_CI_final/MGYS00005138/MGYS00005138_SSU_OTUtable.csv", header = TRUE, stringsAsFactors = FALSE)
otu_5184 <- fread("/disk4/bilje/2022_CI_final/MGYS00005184/MGYS00005184_OTU_Table.csv", header = TRUE, stringsAsFactors = FALSE)
otu_5194 <- fread("/disk4/bilje/2022_CI_final/MGYS00005194/MGYS00005194_OTU_Table.csv", header = TRUE, stringsAsFactors = FALSE)
############### OTU cleansing################################################################################################################
taxonomic <- c("sk", "k", "p", "c", "o", "f", "g")

# 1: "sk", "k", "s" 열 제거
otu_5138_1 <- otu_5138[,-c("sk", "k", "s")]
otu_5184_1 <- otu_5184[,-c("sk", "k", "s")]
otu_5194_1 <- otu_5194[,-c("sk", "k", "s")]
############### aggregate########################################################
## phylum

# p1: "c", "o", "f", "g" 열 제거
otu_5138_p1 <- otu_5138_1[,-c("c", "o", "f", "g")]
otu_5184_p1 <- otu_5184_1[,-c("c", "o", "f", "g")]
otu_5194_p1 <- otu_5194_1[,-c("c", "o", "f", "g")]

# p2: p열에 ""이 존재하는 행 제거
otu_5138_p2 <- otu_5138_p1[-which(otu_5138_p1$p == ""),]
otu_5184_p2 <- otu_5184_p1[-which(otu_5184_p1$p == ""),]
otu_5194_p2 <- otu_5194_p1[-which(otu_5194_p1$p == ""),]

# p3: 행 중복 합하기
otu_5138_p3 <- aggregate(otu_5138_p2[,-1], by=list(otu_5138_p2$p), FUN=sum)
otu_5184_p3 <- aggregate(otu_5184_p2[,-1], by=list(otu_5184_p2$p), FUN=sum)
otu_5194_p3 <- aggregate(otu_5194_p2[,-1], by=list(otu_5194_p2$p), FUN=sum)

# p4: p__ 행 제거
otu_5138_p4 <- otu_5138_p3[-1,]
otu_5184_p4 <- otu_5184_p3[-1,]
otu_5194_p4 <- otu_5194_p3[-1,]


# 최종 phylum OTU
otu_p_first <- merge(otu_5138_p4, otu_5184_p4, by='Group.1', all=T)
otu_p <- merge(otu_p_first, otu_5194_p4, by='Group.1', all=T)
otu_p[is.na(otu_p)] <- 0


## class

# c1: "p", "o", "f", "g" 열 제거
otu_5138_c1 <- otu_5138_1[,-c("p", "o", "f", "g")]
otu_5184_c1 <- otu_5184_1[,-c("p", "o", "f", "g")]
otu_5194_c1 <- otu_5194_1[,-c("p", "o", "f", "g")]

# c2: c열에 ""이 존재하는 행 제거
otu_5138_c2 <- otu_5138_c1[-which(otu_5138_c1$c == ""),]
otu_5184_c2 <- otu_5184_c1[-which(otu_5184_c1$c == ""),]
otu_5194_c2 <- otu_5194_c1[-which(otu_5194_c1$c == ""),]

# c3: 행 중복 합하기
otu_5138_c3 <- aggregate(otu_5138_c2[,-1], by=list(otu_5138_c2$c), FUN=sum)
otu_5184_c3 <- aggregate(otu_5184_c2[,-1], by=list(otu_5184_c2$c), FUN=sum)
otu_5194_c3 <- aggregate(otu_5194_c2[,-1], by=list(otu_5194_c2$c), FUN=sum)

# c4: c__ 행 제거
otu_5138_c4 <- otu_5138_c3[-1,]
otu_5184_c4 <- otu_5184_c3[-1,]
otu_5194_c4 <- otu_5194_c3[-1,]

## 최종 class OTU
otu_c_first <- merge(otu_5138_c4, otu_5184_c4, by='Group.1', all=T)
otu_c <- merge(otu_c_first, otu_5194_c4, by='Group.1', all=T)
otu_c[is.na(otu_c)] <- 0


## order

# o1: "p", "c", "f", "g" 열 제거
otu_5138_o1 <- otu_5138_1[,-c("p", "c", "f", "g")]
otu_5184_o1 <- otu_5184_1[,-c("p", "c", "f", "g")]
otu_5194_o1 <- otu_5194_1[,-c("p", "c", "f", "g")]

# o2: o열에 ""이 존재하는 행 제거
otu_5138_o2 <- otu_5138_o1[-which(otu_5138_o1$o == ""),]
otu_5184_o2 <- otu_5184_o1[-which(otu_5184_o1$o == ""),]
otu_5194_o2 <- otu_5194_o1[-which(otu_5194_o1$o == ""),]

# o3: 행 중복 합하기
otu_5138_o3 <- aggregate(otu_5138_o2[,-1], by=list(otu_5138_o2$o), FUN=sum)
otu_5184_o3 <- aggregate(otu_5184_o2[,-1], by=list(otu_5184_o2$o), FUN=sum)
otu_5194_o3 <- aggregate(otu_5194_o2[,-1], by=list(otu_5194_o2$o), FUN=sum)

# o4: o__ 행 제거
otu_5138_o4 <- otu_5138_o3[-1,]
otu_5184_o4 <- otu_5184_o3[-1,]
otu_5194_o4 <- otu_5194_o3[-1,]


# 최종 order OTU
otu_o_first <- merge(otu_5138_o4, otu_5184_o4, by='Group.1', all=T)
otu_o <- merge(otu_o_first, otu_5194_o4, by='Group.1', all=T)
otu_o[is.na(otu_o)] <- 0


## family

# f1: "p", "c", "o", "g" 열 제거
otu_5138_f1 <- otu_5138_1[,-c("p", "c", "o", "g")]
otu_5184_f1 <- otu_5184_1[,-c("p", "c", "o", "g")]
otu_5194_f1 <- otu_5194_1[,-c("p", "c", "o", "g")]

# f2: o열에 ""이 존재하는 행 제거
otu_5138_f2 <- otu_5138_f1[-which(otu_5138_f1$f == ""),]
otu_5184_f2 <- otu_5184_f1[-which(otu_5184_f1$f == ""),]
otu_5194_f2 <- otu_5194_f1[-which(otu_5194_f1$f == ""),]

# f3: 행 중복 합하기
otu_5138_f3 <- aggregate(otu_5138_f2[,-1], by=list(otu_5138_f2$f), FUN=sum)
otu_5184_f3 <- aggregate(otu_5184_f2[,-1], by=list(otu_5184_f2$f), FUN=sum)
otu_5194_f3 <- aggregate(otu_5194_f2[,-1], by=list(otu_5194_f2$f), FUN=sum)

# f4: f__ 행 제거
otu_5138_f4 <- otu_5138_f3[-1,]
otu_5184_f4 <- otu_5184_f3[-1,]
otu_5194_f4 <- otu_5194_f3[-1,]

# 최종 family OTU
otu_f_first <- merge(otu_5138_f4, otu_5184_f4, by='Group.1', all=T)
otu_f <- merge(otu_f_first, otu_5194_f4, by='Group.1', all=T)
otu_f[is.na(otu_f)] <- 0


## genus

# g1: "p", "c", "o", "f" 열 제거
otu_5138_g1 <- otu_5138_1[,-c("p", "c", "o", "f")]
otu_5184_g1 <- otu_5184_1[,-c("p", "c", "o", "f")]
otu_5194_g1 <- otu_5194_1[,-c("p", "c", "o", "f")]

# g2: o열에 ""이 존재하는 행 제거
otu_5138_g2 <- otu_5138_g1[-which(otu_5138_g1$g == ""),]
otu_5184_g2 <- otu_5184_g1[-which(otu_5184_g1$g == ""),]
otu_5194_g2 <- otu_5194_g1[-which(otu_5194_g1$g == ""),]

# g3: 행 중복 합하기
otu_5138_g3 <- aggregate(otu_5138_g2[,-1], by=list(otu_5138_g2$g), FUN=sum)
otu_5184_g3 <- aggregate(otu_5184_g2[,-1], by=list(otu_5184_g2$g), FUN=sum)
otu_5194_g3 <- aggregate(otu_5194_g2[,-1], by=list(otu_5194_g2$g), FUN=sum)

# g4: g__ 행 제거
otu_5138_g4 <- otu_5138_g3[-1,]
otu_5184_g4 <- otu_5184_g3[-1,]
otu_5194_g4 <- otu_5194_g3[-1,]

# 최종 genus OTU
otu_g_first <- merge(otu_5138_g4, otu_5184_g4, by='Group.1', all=T)
otu_g <- merge(otu_g_first, otu_5194_g4, by='Group.1', all=T)
otu_g[is.na(otu_g)] <- 0


############### 전치 & 최종 OTU##################################################
otu_p_t <- t(otu_p)
otu_c_t <- t(otu_c)
otu_o_t <- t(otu_o)
otu_f_t <- t(otu_f)
otu_g_t <- t(otu_g)

colnames(otu_p_t) <-  otu_p_t["Group.1",]
otu_p_final <- otu_p_t[-1,]
colnames(otu_c_t) <-  otu_c_t["Group.1",]
otu_c_final <- otu_c_t[-1,]
colnames(otu_o_t) <-  otu_o_t["Group.1",]
otu_o_final <- otu_o_t[-1,]
colnames(otu_f_t) <-  otu_f_t["Group.1",]
otu_f_final <- otu_f_t[-1,]
colnames(otu_g_t) <-  otu_g_t["Group.1",]
otu_g_final <- otu_g_t[-1,]


otu_p_final <- data.frame(otu_p_final) %>% mutate_all(as.numeric)
otu_c_final <- data.frame(otu_c_final) %>% mutate_all(as.numeric)
otu_o_final <- data.frame(otu_o_final) %>% mutate_all(as.numeric)
otu_f_final <- data.frame(otu_f_final) %>% mutate_all(as.numeric)
otu_g_final <- data.frame(otu_g_final) %>% mutate_all(as.numeric)

dim(otu_p_final)
dim(otu_c_final)
dim(otu_o_final)
dim(otu_f_final)
dim(otu_g_final)

#write.csv(otu_p_final, "/disk4/bilje/2022_CI_final/OTU/OTU_phylum.csv")
#write.csv(otu_c_final, "/disk4/bilje/2022_CI_final/OTU/OTU_class.csv")
#write.csv(otu_o_final, "/disk4/bilje/2022_CI_final/OTU/OTU_order.csv")
#write.csv(otu_f_final, "/disk4/bilje/2022_CI_final/OTU/OTU_family.csv")
#write.csv(otu_g_final, "/disk4/bilje/2022_CI_final/OTU/OTU_genus.csv")

############### meta & OTU merge###############################################################################################################
meta_KEY <- c("KEY", "Host_disease")
extract_meta_KEY <- meta_final[,..meta_KEY]

metaOTU_p <- merge(otu_p_final, extract_meta_KEY, by.x = 0, by.y = "KEY", all = F)
metaOTU_c <- merge(otu_c_final, extract_meta_KEY, by.x = 0, by.y = "KEY", all = F)
metaOTU_o <- merge(otu_o_final, extract_meta_KEY, by.x = 0, by.y = "KEY", all = F)
metaOTU_f <- merge(otu_f_final, extract_meta_KEY, by.x = 0, by.y = "KEY", all = F)
metaOTU_g <- merge(otu_g_final, extract_meta_KEY, by.x = 0, by.y = "KEY", all = F)

rownames(metaOTU_p) <- metaOTU_p[,1]
metaOTU_p <- metaOTU_p[,-1]
rownames(metaOTU_c) <- metaOTU_c[,1]
metaOTU_c <- metaOTU_c[,-1]
rownames(metaOTU_o) <- metaOTU_o[,1]
metaOTU_o <- metaOTU_o[,-1]
rownames(metaOTU_f) <- metaOTU_f[,1]
metaOTU_f <- metaOTU_f[,-1]
rownames(metaOTU_g) <- metaOTU_g[,1]
metaOTU_g <- metaOTU_g[,-1]
############### 차원축소######################################################################################################################
library(Rtsne)
library(pheatmap)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(caret)
library(tidyverse)
library(tsne)
library(Rtsne)
library(scales)
library(factoextra)
library(FactoMineR)
############### normalization (TMM)#############################################
#질병정보 제거한 metaOTU data 형성
metaOTU_p_sub <- subset(metaOTU_p, select=-Host_disease)
metaOTU_c_sub <- subset(metaOTU_c, select=-Host_disease)
metaOTU_o_sub <- subset(metaOTU_o, select=-Host_disease)
metaOTU_f_sub <- subset(metaOTU_f, select=-Host_disease)
metaOTU_g_sub <- subset(metaOTU_g, select=-Host_disease)

# TSS
# normalization
# nor_minmax = function(x){
#   result = (x - min(x)) / (max(x) - min(x))
#   return(result)
# }
# metaOTU_p_nor <- metaOTU_p_sub
# metaOTU_p_nor$library_size <- rowSums(metaOTU_p_sub)
# for (i in c(1:nrow(metaOTU_p_nor))){
#   metaOTU_p_nor[i,] <- nor_minmax(metaOTU_p_nor[i,])
# }
# metaOTU_p_nor <- metaOTU_p_nor[,-91]

# TMM
library(edgeR)

OTU_p_count <- DGEList(counts=t(metaOTU_p_sub))
OTU_p_NormFac <- calcNormFactors(OTU_p_count, method= "TMM")
logCPM_p <- cpm(OTU_p_NormFac, log = T)
OTU_p_norm <- as.data.frame(t(logCPM_p))
OTU_p_norm$Host_disease <- metaOTU_p$Host_disease

OTU_c_count <- DGEList(counts=t(metaOTU_c_sub))
OTU_c_NormFac <- calcNormFactors(OTU_c_count, method= "TMM")
logCPM_c <- cpm(OTU_c_NormFac, log = T)
OTU_c_norm <- as.data.frame(t(logCPM_c))
OTU_c_norm$Host_disease <- metaOTU_c$Host_disease

OTU_o_count <- DGEList(counts=t(metaOTU_o_sub))
OTU_o_NormFac <- calcNormFactors(OTU_o_count, method= "TMM")
logCPM_o <- cpm(OTU_o_NormFac, log = T)
OTU_o_norm <- as.data.frame(t(logCPM_o))
OTU_o_norm$Host_disease <- metaOTU_o$Host_disease

OTU_f_count <- DGEList(counts=t(metaOTU_f_sub))
OTU_f_NormFac <- calcNormFactors(OTU_f_count, method= "TMM")
logCPM_f <- cpm(OTU_f_NormFac, log = T)
OTU_f_norm <- as.data.frame(t(logCPM_f))
OTU_f_norm$Host_disease <- metaOTU_f$Host_disease

OTU_g_count <- DGEList(counts=t(metaOTU_g_sub))
OTU_g_NormFac <- calcNormFactors(OTU_g_count, method= "TMM")
logCPM_g <- cpm(OTU_g_NormFac, log = T)
OTU_g_norm <- as.data.frame(t(logCPM_g))
OTU_g_norm$Host_disease <- metaOTU_g$Host_disease

############### PCA #########################################
# phylum 분산 0인거 제거 & 주성분 분석 & PCA 
metrix_p_be <- nearZeroVar(metaOTU_p, saveMetrics=TRUE)
idx_p_be <- which(metrix_p_be$zeroVar == TRUE)
OTU_p_be <- metaOTU_p[,-idx_p_be]

pca_p_be <- prcomp(subset(OTU_p_be, select=-Host_disease))
summary(pca_p_be)

plot(pca_p_be, type='lines', main="Phylum")
biplot(pca_p_be, main="Phylum")

pca_p2_be <- PCA(subset(OTU_p_be, select=-Host_disease), graph = T)
OTU_p_be$Host_disease <- as.factor(OTU_p_be$Host_disease)
pca1_be <- fviz_pca_ind(pca_p2_be, habillage = OTU_p_be$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Phylum")  +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# class 분산 0인거 제거 & 주성분 분석 & PCA 
metrix_c_be <- nearZeroVar(metaOTU_c, saveMetrics=TRUE)
idx_c_be <- which(metrix_c_be$zeroVar == TRUE)
OTU_c_be <- metaOTU_c[,-idx_c_be]

pca_c_be <- prcomp(subset(OTU_c_be, select=-Host_disease))
summary(pca_c_be)

plot(pca_c_be, type='lines', main="Class")
biplot(pca_c_be, main="Class")

pca_c2_be <- PCA(subset(OTU_c_be, select=-Host_disease), graph = T)
OTU_c_be$Host_disease <- as.factor(OTU_c_be$Host_disease)
pca2_be <- fviz_pca_ind(pca_c2_be, habillage = OTU_c_be$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Class")  +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# order 분산 0인거 제거 & 주성분 분석 & PCA 
metrix_o_be <- nearZeroVar(metaOTU_o, saveMetrics=TRUE)
idx_o_be <- which(metrix_o_be$zeroVar == TRUE)
OTU_o_be <- metaOTU_o[,-idx_o_be]

pca_o_be <- prcomp(subset(OTU_o_be, select=-Host_disease))
summary(pca_o_be)

plot(pca_o_be, type='lines', main="Order")
biplot(pca_o_be, main="Order")

pca_o2_be <- PCA(subset(OTU_o_be, select=-Host_disease), graph = T)
OTU_o_be$Host_disease <- as.factor(OTU_o_be$Host_disease)
pca3_be <- fviz_pca_ind(pca_o2_be, habillage = OTU_o_be$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Order")  +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# family 분산 0인거 제거 & 주성분 분석 & PCA 
metrix_f_be <- nearZeroVar(metaOTU_f, saveMetrics=TRUE)
idx_f_be <- which(metrix_f_be$zeroVar == TRUE)
OTU_f_be <- metaOTU_f[,-idx_f_be]

pca_f_be <- prcomp(subset(OTU_f_be, select=-Host_disease))
summary(pca_f_be)

plot(pca_f_be, type='lines', main="Family")
biplot(pca_f_be, main="Family")

pca_f2_be <- PCA(subset(OTU_f_be, select=-Host_disease), graph = T)
OTU_f_be$Host_disease <- as.factor(OTU_f_be$Host_disease)
pca4_be <- fviz_pca_ind(pca_f2_be, habillage = OTU_f_be$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Family")  +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# genus 분산 0인거 제거 & 주성분 분석 & PCA 
metrix_g_be <- nearZeroVar(metaOTU_g, saveMetrics=TRUE)
idx_g_be <- which(metrix_g_be$zeroVar == TRUE)
OTU_g_be <- metaOTU_g[,-idx_g_be]

pca_g_be <- prcomp(subset(OTU_g_be, select=-Host_disease))
summary(pca_g_be)

plot(pca_g_be, type='lines', main="Genus")
biplot(pca_g_be, main="Genus")

pca_g2_be <- PCA(subset(OTU_g_be, select=-Host_disease), graph = T)
OTU_g_be$Host_disease <- as.factor(OTU_g_be$Host_disease)
pca5_be <- fviz_pca_ind(pca_g2_be, habillage = OTU_g_be$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Genus")  +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

library(patchwork)
pca1_be + pca2_be + pca3_be + pca4_be + pca5_be + 
  theme(legend.position = c(1.6,0.5),
        legend.background = element_rect(color = "black", fill = "gainsboro"), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(2.6, 'cm'),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16),
        legend.key.size = unit(10, 'cm')) + 
  guides(disease = guide_legend(title.position = "top"), color = guide_legend(override.aes = list(size = 5)))

############### PCA after normalization#########################################
# phylum 분산 0인거 제거 & 주성분 분석 & PCA 
metrix_p <- nearZeroVar(OTU_p_norm, saveMetrics=TRUE)
idx_p <- which(metrix_p$zeroVar == TRUE)
OTU_p_norm_sub <- OTU_p_norm[,-idx_p]

pca_p <- prcomp(subset(OTU_p_norm_sub, select=-Host_disease))
summary(pca_p)

plot(pca_p, type='lines', main="Phylum")
biplot(pca_p, main="Phylum")

pca_p2 <- PCA(subset(OTU_p_norm_sub, select=-Host_disease), graph = T)
OTU_p_norm_sub$Host_disease <- as.factor(OTU_p_norm_sub$Host_disease)
pca1 <- fviz_pca_ind(pca_p2, habillage = OTU_p_norm_sub$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Phylum")  +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# class 분산 0인거 제거 & 주성분 분석 & PCA 
metrix_c <- nearZeroVar(OTU_c_norm, saveMetrics=TRUE)
idx_c <- which(metrix_c$zeroVar == TRUE)
OTU_c_norm_sub <- OTU_c_norm[,-idx_c]

pca_c <- prcomp(subset(OTU_c_norm_sub, select=-Host_disease))
summary(pca_c)

plot(pca_c, type='lines', main="Class")
biplot(pca_c, main="Class")

pca_c2 <- PCA(subset(OTU_c_norm_sub, select=-Host_disease), graph = T)
OTU_c_norm_sub$Host_disease <- as.factor(OTU_c_norm_sub$Host_disease)
pca2 <- fviz_pca_ind(pca_c2, habillage = OTU_c_norm_sub$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Class")  +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# order 분산 0인거 제거 & 주성분 분석 & PCA 
metrix_o <- nearZeroVar(OTU_o_norm, saveMetrics=TRUE)
idx_o <- which(metrix_o$zeroVar == TRUE)
OTU_o_norm_sub <- OTU_o_norm[,-idx_o]

pca_o <- prcomp(subset(OTU_o_norm_sub, select=-Host_disease))
summary(pca_o)

plot(pca_o, type='lines', main="Order")
biplot(pca_o, main="Order")

pca_o2 <- PCA(subset(OTU_o_norm_sub, select=-Host_disease))
OTU_o_norm_sub$Host_disease <- as.factor(OTU_o_norm_sub$Host_disease)
pca3 <- fviz_pca_ind(pca_o2, habillage = OTU_o_norm_sub$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Order") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# faimly 분산 0인거 제거 & 주성분 분석 & PCA
metrix_f <- nearZeroVar(OTU_f_norm, saveMetrics=TRUE)
idx_f <- which(metrix_f$zeroVar == TRUE)
OTU_f_norm_sub <- OTU_f_norm[,-idx_f]

pca_f <- prcomp(subset(OTU_f_norm_sub, select=-Host_disease))
summary(pca_f)

plot(pca_f, type='lines', main="Family")
biplot(pca_f, main="Family")

pca_f2 <- PCA(subset(OTU_f_norm_sub, select=-Host_disease), graph = T)
OTU_f_norm_sub$Host_disease <- as.factor(OTU_f_norm_sub$Host_disease)
pca4 <- fviz_pca_ind(pca_f2, habillage = OTU_f_norm_sub$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Family") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# genus 분산 0인거 제거 & 주성분 분석 & PCA
metrix_g <- nearZeroVar(OTU_g_norm, saveMetrics=TRUE)
idx_g <- which(metrix_g$zeroVar == TRUE)
OTU_g_norm_sub <- OTU_g_norm[,-idx_g]

pca_g <- prcomp(subset(OTU_g_norm_sub, select=-Host_disease))
summary(pca_g)

plot(pca_g, type='lines', main="Genus")
biplot(pca_g, main="Genus")

pca_g2 <- PCA(subset(OTU_g_norm_sub, select=-Host_disease), graph = T)
OTU_g_norm_sub$Host_disease <- as.factor(OTU_g_norm_sub$Host_disease)
pca5 <- fviz_pca_ind(pca_g2, habillage = OTU_g_norm_sub$Host_disease, label="none", geom.ind = c("point"), geom = c("point"), title="PCA - Genus") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

library(patchwork)
pca1 + pca2 + pca3 + pca4 + pca5 + 
  theme(legend.position = c(1.6,0.5),
        legend.background = element_rect(color = "black", fill = "gainsboro"), 
        legend.key.size = unit(3, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(2.6, 'cm'),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16)) + 
  guides(disease = guide_legend(title.position = "top"), color = guide_legend(override.aes = list(size = 5)))

############### tsne############################################################
#phylum
#OTU_p_be는 분산제거
p_tsne_be <- Rtsne(subset(OTU_p_be, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

p_tsne_df_be <- data.frame(tsne_x = p_tsne_be$Y[,1], 
                           tsne_y = p_tsne_be$Y[,2], 
                           Host_disease = OTU_p_be$Host_disease)

plot_pt_be <- p_tsne_df_be %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Phylum") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

# metaOTU_p_nor$Host_disease <- metaOTU_p$Host_disease
# p_tsne <- Rtsne(metaOTU_p_nor[,-91], PCA = TRUE, dims = 2, max_iter = 500, perplexity = 40)
# 
# p_tsne_df <- data.frame(tsne_x = p_tsne$Y[,1],
#                         tsne_y = p_tsne$Y[,2],
#                         Host_disease = metaOTU_p_nor$Host_disease)
# 
# p_tsne_df %>%
#   ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) +
#   geom_point(size = 1) +
#   labs(title="tsne phylum")


#class
c_tsne_be <- Rtsne(subset(OTU_c_be, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

c_tsne_df_be <- data.frame(tsne_x = c_tsne_be$Y[,1], 
                          tsne_y = c_tsne_be$Y[,2], 
                          Host_disease = OTU_c_be$Host_disease)

plot_ct_be <- c_tsne_df_be %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Class") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

#order
o_tsne_be <- Rtsne(subset(OTU_o_be, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

o_tsne_df_be <- data.frame(tsne_x = o_tsne_be$Y[,1], 
                           tsne_y = o_tsne_be$Y[,2], 
                           Host_disease = OTU_o_be$Host_disease)

plot_ot_be <- o_tsne_df_be %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Order") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

#family
f_tsne_be <- Rtsne(subset(OTU_f_be, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

f_tsne_df_be <- data.frame(tsne_x = f_tsne_be$Y[,1], 
                           tsne_y = f_tsne_be$Y[,2], 
                           Host_disease = OTU_f_be$Host_disease)

plot_ft_be <- f_tsne_df_be %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Family") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

#genus
g_tsne_be <- Rtsne(subset(OTU_g_be, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

g_tsne_df_be <- data.frame(tsne_x = g_tsne_be$Y[,1], 
                           tsne_y = g_tsne_be$Y[,2], 
                           Host_disease = OTU_g_be$Host_disease)

plot_gt_be <- g_tsne_df_be %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  labs(title="tSNE - Genus") +
  theme_bw() +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

plot_pt_be + plot_ct_be + plot_ot_be + plot_ft_be + plot_gt_be +
  theme(legend.position = c(1.6,0.5),
        legend.background = element_rect(color = "black", fill = "gainsboro"), 
        legend.key.size = unit(3, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(2.6, 'cm'),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16)) +
        guides(disease = guide_legend(title.position = "top"), color = guide_legend(override.aes = list(size = 7)))

############### tsne after normalization########################################
# phylum
# OTU_p_norm_sub는 분산 제거
p_tsne <- Rtsne(subset(OTU_p_norm_sub, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

p_tsne_df <- data.frame(tsne_x = p_tsne$Y[,1], 
                        tsne_y = p_tsne$Y[,2], 
                        Host_disease = OTU_p_norm_sub$Host_disease)

plot_p <- p_tsne_df %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Phylum") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

#class
c_tsne <- Rtsne(subset(OTU_c_norm_sub, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

c_tsne_df <- data.frame(tsne_x = c_tsne$Y[,1], 
                         tsne_y = c_tsne$Y[,2], 
                         Host_disease = OTU_c_norm_sub$Host_disease)

plot_c <- c_tsne_df %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Class") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

#order
o_tsne <- Rtsne(subset(OTU_o_norm_sub, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

o_tsne_df <- data.frame(tsne_x = o_tsne$Y[,1], 
                         tsne_y = o_tsne$Y[,2], 
                         Host_disease = OTU_o_norm_sub$Host_disease)

plot_o <- o_tsne_df %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Order") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

#family
f_tsne <- Rtsne(subset(OTU_f_norm_sub, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

f_tsne_df <- data.frame(tsne_x = f_tsne$Y[,1], 
                        tsne_y = f_tsne$Y[,2], 
                        Host_disease = OTU_f_norm_sub$Host_disease)

plot_f <- f_tsne_df %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Family") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

#genus
g_tsne <- Rtsne(subset(OTU_g_norm_sub, select=-Host_disease), PCA = TRUE, dims = 2, max_iter = 500)

g_tsne_df <- data.frame(tsne_x = g_tsne$Y[,1], 
                         tsne_y = g_tsne$Y[,2], 
                         Host_disease = OTU_g_norm_sub$Host_disease)

plot_g <- g_tsne_df %>% 
  ggplot(aes(x = tsne_x, y = tsne_y, color = Host_disease)) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(title="tSNE - Genus") +
  theme(legend.position = "none", plot.title=element_text(size=20), axis.text.x = element_text(size =10), axis.text.y = element_text(size = 10))

plot_p + plot_c + plot_o + plot_f + plot_g + 
  theme(legend.position = c(1.6,0.5),
        legend.background = element_rect(color = "black", fill = "gainsboro"), 
        legend.key.size = unit(3, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(2.6, 'cm'),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16)) +
  guides(disease = guide_legend(title.position = "top"), color = guide_legend(override.aes = list(size = 7)))

############### modeling######################################################################################################################
#normalization 이후
############### class 패키지 - knn##############################################
#tsne이용
library(class)

## phylum
p_tsne_sample <- sample(2, nrow(p_tsne_df), replace=T, prob=c(0.8,0.2))
train_pt <- p_tsne_df[p_tsne_sample==1,]
test_pt <- p_tsne_df[p_tsne_sample==2,]

# train & test 확인
plot(formula = tsne_y ~ tsne_x, data = train_pt, pch = 20, col = c("red", "green", "blue")[train_pt$Host_disease], main = "Classification")
points(formula = tsne_y ~ tsne_x, data = train_pt, pch = 20, cex = 1, col = "purple")
points(formula = tsne_y ~ tsne_x, data = test_pt, pch = 20, cex = 1, col = "pink")

# k = 3 일 때
set.seed(1234)
knn_pt <- knn(train = train_pt[,-3],
             test = test_pt[,-3],
             cl = train_pt[,3],
             k = 3)

# train_pt 산점도 그리기
plot(formula = tsne_y ~ tsne_x,
     data = train_pt,
     col = alpha(c("red", "green", "blue"), 0.7)[train_pt$Host_disease],
     main = "KNN (k = 3) - Phylum", pch = 20)

# knn test 결과 표시하기
points(formula = tsne_y ~ tsne_x,
       data = test_pt,
       pch = 20,
       cex = 1,
       col = alpha(c("red", "green", "blue"), 0.7)[knn_pt])


# 분류 정확도 계산하기
sum(knn_pt == test_pt[,3]) / length(test_pt[,3])

## class
c_tsne_sample <- sample(2, nrow(c_tsne_df), replace=T, prob=c(0.8,0.2))
train_ct <- c_tsne_df[c_tsne_sample==1,]
test_ct <- c_tsne_df[c_tsne_sample==2,]

# train & test 확인
plot(formula = tsne_y ~ tsne_x, data = train_ct, pch = 20, col = c("red", "green", "blue")[train_ct$Host_disease], main = "Classification")
points(formula = tsne_y ~ tsne_x, data = train_ct, pch = 20, cex = 1, col = "purple")
points(formula = tsne_y ~ tsne_x, data = test_ct, pch = 20, cex = 1, col = "pink")

# k = 3 일 때
set.seed(1234)
knn_ct <- knn(train = train_ct[,-3],
              test = test_ct[,-3],
              cl = train_ct[,3],
              k = 3)

# train_ct 산점도 그리기
plot(formula = tsne_y ~ tsne_x,
     data = train_ct,
     col = alpha(c("red", "green", "blue"), 0.7)[train_ct$Host_disease],
     main = "KNN (k = 3) - Class", pch = 20)

# knn test 결과 표시하기
points(formula = tsne_y ~ tsne_x,
       data = test_ct,
       pch = 20,
       cex = 1,
       col = alpha(c("red", "green", "blue"), 0.7)[knn_ct])

# 분류 정확도 계산하기
sum(knn_ct == test_ct[,3]) / length(test_ct[,3])

## order
o_tsne_sample <- sample(2, nrow(o_tsne_df), replace=T, prob=c(0.8,0.2))
train_ot <- o_tsne_df[o_tsne_sample==1,]
test_ot <- o_tsne_df[o_tsne_sample==2,]

# train & test 확인
plot(formula = tsne_y ~ tsne_x, data = train_ot, pch = 20, col = c("red", "green", "blue")[train_ot$Host_disease], main = "Classification")
points(formula = tsne_y ~ tsne_x, data = train_ot, pch = 20, cex = 1, col = "purple")
points(formula = tsne_y ~ tsne_x, data = test_ot, pch = 20, cex = 1, col = "pink")

# k = 3 일 때
set.seed(1234)
knn_ot <- knn(train = train_ot[,-3],
              test = test_ot[,-3],
              cl = train_ot[,3],
              k = 3)

# train_ot 산점도 그리기
plot(formula = tsne_y ~ tsne_x,
     data = train_ot,
     col = alpha(c("red", "green", "blue"), 0.7)[train_ot$Host_disease],
     main = "KNN (k = 3) - Order", pch = 20)

# knn test 결과 표시하기
points(formula = tsne_y ~ tsne_x,
       data = test_ot,
       pch = 20,
       cex = 1,
       col = alpha(c("red", "green", "blue"), 0.7)[knn_ot])

# 분류 정확도 계산하기
sum(knn_ot == test_ot[,3]) / length(test_ot[,3])

## Family
f_tsne_sample <- sample(2, nrow(f_tsne_df), replace=T, prob=c(0.8,0.2))
train_ft <- f_tsne_df[f_tsne_sample==1,]
test_ft <- f_tsne_df[f_tsne_sample==2,]

# train & test 확인
plot(formula = tsne_y ~ tsne_x, data = train_ft, pch = 20, col = c("red", "green", "blue")[train_ft$Host_disease], main = "Classification")
points(formula = tsne_y ~ tsne_x, data = train_ft, pch = 20, cex = 1, col = "purple")
points(formula = tsne_y ~ tsne_x, data = test_ft, pch = 20, cex = 1, col = "pink")

# k = 3 일 때
set.seed(1234)
knn_ft <- knn(train = train_ft[,-3],
              test = test_ft[,-3],
              cl = train_ft[,3],
              k = 3)

# train_ft 산점도 그리기
plot(formula = tsne_y ~ tsne_x,
     data = train_ft,
     col = alpha(c("red", "green", "blue"), 0.7)[train_ft$Host_disease],
     main = "KNN (k = 3) - Family", pch = 20)

# knn test 결과 표시하기
points(formula = tsne_y ~ tsne_x,
       data = test_ft,
       pch = 20,
       cex = 1,
       col = alpha(c("red", "green", "blue"), 0.7)[knn_ft])

# 분류 정확도 계산하기
sum(knn_ft == test_ft[,3]) / length(test_ft[,3])

## Genus
g_tsne_sample <- sample(2, nrow(g_tsne_df), replace=T, prob=c(0.8,0.2))
train_gt <- g_tsne_df[g_tsne_sample==1,]
test_gt <- g_tsne_df[g_tsne_sample==2,]

# train & test 확인
plot(formula = tsne_y ~ tsne_x, data = train_gt, pch = 20, col = c("red", "green", "blue")[train_gt$Host_disease], main = "Classification")
points(formula = tsne_y ~ tsne_x, data = train_gt, pch = 20, cex = 1, col = "purple")
points(formula = tsne_y ~ tsne_x, data = test_gt, pch = 20, cex = 1, col = "pink")

# k = 3 일 때
set.seed(1234)
knn_gt <- knn(train = train_gt[,-3],
              test = test_gt[,-3],
              cl = train_gt[,3],
              k = 3)

# train_gt 산점도 그리기
plot(formula = tsne_y ~ tsne_x,
     data = train_gt,
     col = alpha(c("red", "green", "blue"), 0.7)[train_gt$Host_disease],
     main = "KNN (k = 3) - Genus", pch = 20)

# knn test 결과 표시하기
points(formula = tsne_y ~ tsne_x,
       data = test_gt,
       pch = 20,
       cex = 1,
       col = alpha(c("red", "green", "blue"), 0.7)[knn_gt])

# 분류 정확도 계산하기
sum(knn_gt == test_gt[,3]) / length(test_gt[,3])

############### caret 패키지####################################################
library(caret)
# train : train_np, train_nc, train_no, train_nf, train_ng
# test : test_np, test_nc, test_no, test_nf, test_ng
# OTU_p_norm_sub는 정규화 이후 분산 0인거 제거한것
sample_np <- sample(2, nrow(OTU_p_norm_sub), replace=T, prob=c(0.8,0.2))
train_np <- OTU_p_norm_sub[sample_np==1,]
test_np <- OTU_p_norm_sub[sample_np==2,]

sample_nc <- sample(2, nrow(OTU_c_norm_sub), replace=T, prob=c(0.8,0.2))
train_nc <- OTU_c_norm_sub[sample_nc==1,]
test_nc <- OTU_c_norm_sub[sample_nc==2,]

sample_no <- sample(2, nrow(OTU_o_norm_sub), replace=T, prob=c(0.8,0.2))
train_no <- OTU_o_norm_sub[sample_no==1,]
test_no <- OTU_o_norm_sub[sample_no==2,]

sample_nf <- sample(2, nrow(OTU_f_norm_sub), replace=T, prob=c(0.8,0.2))
train_nf <- OTU_f_norm_sub[sample_nf==1,]
test_nf <- OTU_f_norm_sub[sample_nf==2,]

sample_ng <- sample(2, nrow(OTU_g_norm_sub), replace=T, prob=c(0.8,0.2))
train_ng <- OTU_g_norm_sub[sample_ng==1,]
test_ng <- OTU_g_norm_sub[sample_ng==2,]

############### knn#########
#modeling
model_np_knn <- train(Host_disease~., data = train_np, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

model_nc_knn <- train(Host_disease~., data = train_nc, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

model_no_knn <- train(Host_disease~., data = train_no, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

model_nf_knn <- train(Host_disease~., data = train_nf, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

model_ng_knn <- train(Host_disease~., data = train_ng, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

# perdict
preds_np_knn <- predict(model_np_knn, newdata = test_np)
confusionMatrix(preds_np_knn, as.factor(test_np$Host_disease))

preds_nc_knn <- predict(model_nc_knn, newdata = test_nc)
confusionMatrix(preds_nc_knn, as.factor(test_nc$Host_disease))

preds_no_knn <- predict(model_no_knn, newdata = test_no)
confusionMatrix(preds_no_knn, as.factor(test_no$Host_disease))

preds_nf_knn <- predict(model_nf_knn, newdata = test_nf)
confusionMatrix(preds_nf_knn, as.factor(test_nf$Host_disease))

preds_ng_knn <- predict(model_gt_knn, newdata = test_ng)
confusionMatrix(preds_ng_knn, as.factor(test_ng$Host_disease))

# Accuracy
plot(model_np_knn, xlab = "k", main = "Accuracy by KNN - Phylum") 
plot(model_nc_knn, xlab = "k", main = "Accuracy by KNN - Class")
plot(model_no_knn, xlab = "k", main = "Accuracy by KNN - Order")
plot(model_nf_knn, xlab = "k", main = "Accuracy by KNN - Family")
plot(model_ng_knn, xlab = "k", main = "Accuracy by KNN - Genus")

model_np_knn_ac <- data.frame(model_np_knn$resample[,1])
model_nc_knn_ac <- data.frame(model_nc_knn$resample[,1])
model_no_knn_ac <- data.frame(model_no_knn$resample[,1])
model_nf_knn_ac <- data.frame(model_nf_knn$resample[,1])
model_ng_knn_ac <- data.frame(model_ng_knn$resample[,1])

boxplot(model_np_knn_ac)
boxplot(model_nc_knn_ac)
boxplot(model_no_knn_ac)
boxplot(model_nf_knn_ac)
boxplot(model_ng_knn_ac)

############### svmRadial#########
set.seed(132245)
model_np_svm <-  train(Host_disease ~., data = train_np, method = 'svmRadial', preProcess = c("center","scale"),
                  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")

model_nc_svm <-  train(Host_disease ~., data = train_nc, method = 'svmRadial', preProcess = c("center","scale"),
                  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")

model_no_svm <-  train(Host_disease ~., data = train_no, method = 'svmRadial', preProcess = c("center","scale"),
                  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")

model_nf_svm <-  train(Host_disease ~., data = train_nf, method = 'svmRadial', preProcess = c("center","scale"),
                  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")

model_ng_svm <-  train(Host_disease ~., data = train_ng, method = 'svmRadial', preProcess = c("center","scale"),
                  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")

preds_np_svm <- predict(model_np_svm, newdata = test_np)
preds_nc_svm <- predict(model_nc_svm, newdata = test_nc)
preds_no_svm <- predict(model_no_svm, newdata = test_no)
preds_nf_svm <- predict(model_nf_svm, newdata = test_nf)
preds_ng_svm <- predict(model_ng_svm, newdata = test_ng)

confusionMatrix(data = preds_np_svm, test_np$Host_disease)
confusionMatrix(data = preds_nc_svm, test_nc$Host_disease)
confusionMatrix(data = preds_no_svm, test_no$Host_disease)
confusionMatrix(data = preds_nf_svm, test_nf$Host_disease)
confusionMatrix(data = preds_ng_svm, test_ng$Host_disease)

# Accuracy
plot(model_np_svm, xlab = "C", main = "Accuracy by svmRadial - Phylum") 
plot(model_nc_svm, xlab = "C", main = "Accuracy by svmRadial - Class") 
plot(model_no_svm, xlab = "C", main = "Accuracy by svmRadial - Order") 
plot(model_nf_svm, xlab = "C", main = "Accuracy by svmRadial - Family") 
plot(model_ng_svm, xlab = "C", main = "Accuracy by svmRadial - Genus") 

model_np_svm_ac <- data.frame(model_np_svm$resample[,1])
model_nc_svm_ac <- data.frame(model_nc_svm$resample[,1])
model_no_svm_ac <- data.frame(model_no_svm$resample[,1])
model_nf_svm_ac <- data.frame(model_nf_svm$resample[,1])
model_ng_svm_ac <- data.frame(model_ng_svm$resample[,1])

boxplot(model_np_svm_ac)
boxplot(model_nc_svm_ac)
boxplot(model_no_svm_ac)
boxplot(model_nf_svm_ac)
boxplot(model_ng_svm_ac)

############### nb 나이브 베이즈#########
#modeling
model_np_nb <-  train(Host_disease ~., data = train_np, method = 'nb', preProcess = c("center","scale"),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")

#predict 
preds_np_nb <- predict(model_np_nb, newdata = test_np)

confusionMatrix(data = preds_np_nb, test_np$Host_disease)

# Accuracy
plot(model_np_nb, xlab = "usekernel", main = "Accuracy by nb - Phylum") 

model_np_nb_ac <- data.frame(model_np_nb$resample[,1])

boxplot(model_np_nb_ac)

############### rpart#########
#modeling
model_np_rpart <-  train(Host_disease ~., data = train_np, method = 'rpart', preProcess = c("center","scale"),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")
#predict 
preds_np_rpart <- predict(model_np_rpart, newdata = test_np)

confusionMatrix(data = preds_np_rpart, test_np$Host_disease)

# Accuracy
plot(model_np_rpart, xlab = "cp", main = "Accuracy by rpart - Phylum") 

model_np_rpart_ac <- data.frame(model_np_rpart$resample[,1])

boxplot(model_np_rpart_ac)

############### lda Linear Discriminant Analysis#########
set.seed(7)
# modeling
model_np_lda <-  train(Host_disease ~., data = train_np, method = 'lda', preProcess = c("center","scale"),
                         trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")
# predict 
preds_np_lda <- predict(model_np_lda, newdata = test_np)

confusionMatrix(data = preds_np_lda, test_np$Host_disease)

# Accuracy
model_np_lda_ac <- data.frame(model_np_lda$resample[,1])

boxplot(model_np_lda_ac)


# GLMNET
############### rpart#########
set.seed(7)
fit.glmnet <- train(diabetes~., data=dataset, method="glmnet", metric=metric, preProc=c("center", "scale"), trControl=control)

# CART
############### rpart#########
set.seed(7)
fit.cart <- train(diabetes~., data=dataset, method="rpart", metric=metric, trControl=control)
# C5.0
############### rpart#########
set.seed(7)
fit.c50 <- train(diabetes~., data=dataset, method="C5.0", metric=metric, trControl=control)
# Bagged CART
############### rpart#########
set.seed(7)
fit.treebag <- train(diabetes~., data=dataset, method="treebag", metric=metric, trControl=control)
# 랜덤포레스트
############### rpart#########
set.seed(7)
fit.rf <- train(diabetes~., data=dataset, method="rf", metric=metric, trControl=control)
# Gradient Boosting
############### rpart#########
set.seed(7)
fit.gbm <- train(diabetes~., data=dataset, method="gbm", metric=metric, trControl=control, verbose=FALSE)

# 에러
############### glm 로지스틱 회귀분석#########
# set.seed(7)
# # modeling
# model_np_glm <-  train(Host_disease ~., data = train_np, method = 'glm', preProcess = c("center","scale"),
#                        trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3), metric = "Accuracy")
# # predict 
# preds_np_glm <- predict(model_np_glm, newdata = test_np)
# 
# confusionMatrix(data = preds_np_glm, test_np$Host_disease)
# 
# # Accuracy
# plot(model_np_glm, xlab = "cp", main = "Accuracy by glm - Phylum") 
# 
# model_np_glm_ac <- data.frame(model_np_glm$resample[,1])
# 
# boxplot(model_np_glm_ac)
# 






model_p_knn <- data.frame(rep("knn", times = 30), model_np_knn_ac)
colnames(model_p_knn) <- c( "Method", "Accuracy") 

model_p_svm <- data.frame(rep("svmRadial", times = 30), model_np_svm_ac)
colnames(model_p_svm) <- c( "Method", "Accuracy") 

rbind(model_p_knn, model_p_svm, )




# 참고자료
# https://partrita.github.io/posts/ggpubr/
# https://topepo.github.io/caret/available-models.html
# https://blogs.sas.com/content/saskorea/2017/08/22/%EC%B5%9C%EC%A0%81%EC%9D%98-%EB%A8%B8%EC%8B%A0%EB%9F%AC%EB%8B%9D-%EC%95%8C%EA%B3%A0%EB%A6%AC%EC%A6%98%EC%9D%84-%EA%B3%A0%EB%A5%B4%EA%B8%B0-%EC%9C%84%ED%95%9C-%EC%B9%98%ED%8A%B8/  https://rooney-song.tistory.com/80
# https://blog.daum.net/buillee/908
# https://m.blog.naver.com/PostView.naver?isHttpsRedirect=true&blogId=tjdudwo93&logNo=220980156432














# rattle 패키지의 fancyRpartPlot() 함수를 이용하여 그래프의 모양을 좋게 함
require(rattle)
fancyRpartPlot(model$finalModel)


#normalization 이전
############### caret 패키지####################################################
library(caret)
# 분산 0(nzv)인거 제거
metrix_zero_p <- nearZeroVar(metaOTU_p, saveMetrics=TRUE)
idx_pp <- which(metrix_zero_p$nzv == TRUE)
metaOTU_p_delete <- metaOTU_p[,-idx_pp]

metrix_zero_c <- nearZeroVar(metaOTU_c, saveMetrics=TRUE)
idx_cc <- which(metrix_zero_c$nzv == TRUE)
metaOTU_c_delete <- metaOTU_c[,-idx_cc]

metrix_zero_o <- nearZeroVar(metaOTU_o, saveMetrics=TRUE)
idx_oo <- which(metrix_zero_o$nzv == TRUE)
metaOTU_o_delete <- metaOTU_o[,-idx_oo]

metrix_zero_f <- nearZeroVar(metaOTU_f, saveMetrics=TRUE)
idx_ff <- which(metrix_zero_f$nzv == TRUE)
metaOTU_f_delete <- metaOTU_f[,-idx_ff]

metrix_zero_g <- nearZeroVar(metaOTU_g, saveMetrics=TRUE)
idx_gg <- which(metrix_zero_g$nzv == TRUE)
metaOTU_g_delete <- metaOTU_g[,-idx_gg]

# train : train_pp, train_cc, train_oo, train_ff, train_gg
# test : test_pp, test_cc, test_oo, test_ff, test_gg
sample_p <- sample(2, nrow(metaOTU_p_delete), replace=T, prob=c(0.8,0.2))
train_pp <- metaOTU_p_delete[sample_p==1,]
test_pp <- metaOTU_p_delete[sample_p==2,]

sample_c <- sample(2, nrow(metaOTU_c_delete), replace=T, prob=c(0.8,0.2))
train_cc <- metaOTU_c_delete[sample_c==1,]
test_cc <- metaOTU_c_delete[sample_c==2,]

sample_o <- sample(2, nrow(metaOTU_o_delete), replace=T, prob=c(0.8,0.2))
train_oo <- metaOTU_o_delete[sample_o==1,]
test_oo <- metaOTU_o_delete[sample_o==2,]

sample_f <- sample(2, nrow(metaOTU_f_delete), replace=T, prob=c(0.8,0.2))
train_ff <- metaOTU_f_delete[sample_f==1,]
test_ff <- metaOTU_f_delete[sample_f==2,]

sample_g <- sample(2, nrow(metaOTU_g_delete), replace=T, prob=c(0.8,0.2))
train_gg <- metaOTU_g_delete[sample_g==1,]
test_gg <- metaOTU_g_delete[sample_g==2,]

dim(train_pp)
dim(train_cc)
dim(train_oo)
dim(train_ff)
dim(train_gg)
dim(test_pp)
dim(test_cc)
dim(test_oo)
dim(test_ff)
dim(test_gg)
############### knn#########
#tsne
model_pt_knn <- train(Host_disease~., data = train_p, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

model_ct_knn <- train(Host_disease~., data = train_c, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

model_ot_knn <- train(Host_disease~., data = train_o, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

model_ft_knn <- train(Host_disease~., data = train_f, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))

model_gt_knn <- train(Host_disease~., data = train_g, method = "knn", tuneGrid = expand.grid(k=1:30), 
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                      metric = "Accuracy",  preProcess = c("center","scale"))


preds_pt_knn <- predict(model_pt_knn, newdata = test_p)
confusionMatrix(preds_pt_knn, as.factor(test_p$Host_disease))

preds_ct_knn <- predict(model_ct_knn, newdata = test_c)
confusionMatrix(preds_ct_knn, as.factor(test_c$Host_disease))

preds_ot_knn <- predict(model_ot_knn, newdata = test_o)
confusionMatrix(preds_ot_knn, as.factor(test_o$Host_disease))

preds_ft_knn <- predict(model_ft_knn, newdata = test_f)
confusionMatrix(preds_ft_knn, as.factor(test_f$Host_disease))

preds_gt_knn <- predict(model_gt_knn, newdata = test_g)
confusionMatrix(preds_gt_knn, as.factor(test_g$Host_disease))

plot(model_pt_knn, xlab = "k", main = "Accuracy by KNN - phylum") 
plot(model_ct_knn, xlab = "k", main = "Accuracy by KNN - class")
plot(model_ot_knn, xlab = "k", main = "Accuracy by KNN - order")
plot(model_ft_knn, xlab = "k", main = "Accuracy by KNN - family")
plot(model_gt_knn, xlab = "k", main = "Accuracy by KNN - genus")

# 
model_p_knn <- train(Host_disease~., data = train_pp, method = "knn", tuneGrid = expand.grid(k=1:30), 
                     trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                     metric = "Accuracy",  preProcess = c("center","scale"))

model_c_knn <- train(Host_disease~., data = train_cc, method = "knn", tuneGrid = expand.grid(k=1:30), 
                     trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                     metric = "Accuracy",  preProcess = c("center","scale"))

model_o_knn <- train(Host_disease~., data = train_oo, method = "knn", tuneGrid = expand.grid(k=1:30), 
                     trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                     metric = "Accuracy",  preProcess = c("center","scale"))

model_f_knn <- train(Host_disease~., data = train_ff, method = "knn", tuneGrid = expand.grid(k=1:30), 
                     trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                     metric = "Accuracy",  preProcess = c("center","scale"))

model_g_knn <- train(Host_disease~., data = train_gg, method = "knn", tuneGrid = expand.grid(k=1:30), 
                     trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3),  
                     metric = "Accuracy",  preProcess = c("center","scale"))


preds_p_knn <- predict(model_p_knn, newdata = test_pp)
confusionMatrix(preds_p_knn, as.factor(test_pp$Host_disease))

preds_c_knn <- predict(model_c_knn, newdata = test_cc)
confusionMatrix(preds_c_knn, as.factor(test_cc$Host_disease))

preds_o_knn <- predict(model_o_knn, newdata = test_oo)
confusionMatrix(preds_o_knn, as.factor(test_oo$Host_disease))

preds_f_knn <- predict(model_f_knn, newdata = test_ff)
confusionMatrix(preds_f_knn, as.factor(test_ff$Host_disease))

preds_g_knn <- predict(model_g_knn, newdata = test_gg)
confusionMatrix(preds_g_knn, as.factor(test_gg$Host_disease))

plot(model_p_knn, xlab = "k", main = "Accuracy by KNN - phylum") 
plot(model_c_knn, xlab = "k", main = "Accuracy by KNN - class")
plot(model_o_knn, xlab = "k", main = "Accuracy by KNN - order")
plot(model_f_knn, xlab = "k", main = "Accuracy by KNN - family")
plot(model_g_knn, xlab = "k", main = "Accuracy by KNN - genus")
############### svmRadial#########
set.seed(132245)
model_p_svm <-  train(Host_disease ~., 
                      data = train_p, 
                      method = 'svmRadial', 
                      preProcess = c("center","scale"),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3))

model_c_svm <-  train(Host_disease ~., 
                      data = train_c, 
                      method = 'svmRadial', 
                      preProcess = c("center","scale"),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3))

model_o_svm <-  train(Host_disease ~., 
                      data = train_o, 
                      method = 'svmRadial', 
                      preProcess = c("center","scale"),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3))

model_f_svm <-  train(Host_disease ~., 
                      data = train_f, 
                      method = 'svmRadial', 
                      preProcess = c("center","scale"),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3))

model_g_svm <-  train(Host_disease ~., 
                      data = train_g, 
                      method = 'svmRadial', 
                      preProcess = c("center","scale"),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3))

preds_p_svm <- predict(model_p_svm, newdata = test_p)
preds_c_svm <- predict(model_c_svm, newdata = test_c)
preds_o_svm <- predict(model_o_svm, newdata = test_o)
preds_f_svm <- predict(model_f_svm, newdata = test_f)
preds_g_svm <- predict(model_g_svm, newdata = test_g)

confusionMatrix(data = preds_p_svm, test_p$Host_disease)
confusionMatrix(data = preds_c_svm, test_c$Host_disease)
confusionMatrix(data = preds_o_svm, test_o$Host_disease)
confusionMatrix(data = preds_f_svm, test_f$Host_disease)
confusionMatrix(data = preds_g_svm, test_g$Host_disease)

############### nb#########
model_p_nb <- train(train_p, train_p$Host_disease, method="nb", 
                    preProcess = c("center","scale"),
                    trControl=trainControl(method = "repeatedcv", number = 10, repeats = 3))

############### rpart#########
model_p_rpart <- train(Host_disease ~., 
                       data=train_p,
                       method="rpart",
                       preProcess = c("center","scale"),
                       trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3))


# rattle 패키지의 fancyRpartPlot() 함수를 이용하여 그래프의 모양을 좋게 함
require(rattle)
fancyRpartPlot(model$finalModel)


sample_c <- sample(2, nrow(metaOTU_c_delete), replace=T, prob=c(0.8,0.2))

train_cc <- metaOTU_c_delete[sample_c==1,]
test_cc <- metaOTU_c_delete[sample_c==2,]

# Linear Discriminant Analysis
set.seed(7)
fit.lda <- train(Host_disease~., 
                 data=train_cc, 
                 method="lda", 
                 metric="Accuracy", 
                 preProcess=c("center", "scale"), 
                 trControl=trainControl(method = "repeatedcv", number = 10, repeats = 3))
# 로지스틱 회귀분석
set.seed(seed)
fit.glm <- train(diabetes~., data=dataset, method="glm", metric=metric, trControl=control)
# GLMNET
set.seed(seed)
fit.glmnet <- train(diabetes~., data=dataset, method="glmnet", metric=metric, preProc=c("center", "scale"), trControl=control)
# SVM Radial
set.seed(seed)
fit.svmRadial <- train(diabetes~., data=dataset, method="svmRadial", metric=metric, preProc=c("center", "scale"), trControl=control, fit=FALSE)
# 나이브 베이즈
set.seed(seed)
fit.nb <- train(diabetes~., data=dataset, method="nb", metric=metric, trControl=control)
# CART
set.seed(seed)
fit.cart <- train(diabetes~., data=dataset, method="rpart", metric=metric, trControl=control)
# C5.0
set.seed(seed)
fit.c50 <- train(diabetes~., data=dataset, method="C5.0", metric=metric, trControl=control)
# Bagged CART
set.seed(seed)
fit.treebag <- train(diabetes~., data=dataset, method="treebag", metric=metric, trControl=control)
# 랜덤포레스트
set.seed(seed)
fit.rf <- train(diabetes~., data=dataset, method="rf", metric=metric, trControl=control)
# Gradient Boosting
set.seed(seed)
fit.gbm <- train(diabetes~., data=dataset, method="gbm", metric=metric, trControl=control, verbose=FALSE)

##############################################################################################################################################
#normalization 이후



