# Load libraries
library(raster)
library(sp)
library(sf)
library(terra)
library(caret)
library(FactoMineR)
library(randomForest)
library(parallel)
library(doParallel)

# Load data
setwd("C://Documents and Settings/gvanh/Documents/UFAL/mod/")
base <- raster("base.tif")
base_ind <- raster("base_ind_n.tif")
base_amp <- raster("base_amp_n.tif")
base_dec <- raster("base_dec_n.tif")
base_asp <- raster("base_asp_n.tif")
base_tpi <- raster("base_tpi_n.tif")
base_tri <- raster("base_tri_n.tif")
base_wet <- raster("base_wet_n.tif")
base_vis <- raster("base_vis_n.tif")
base_dist <- raster("base_dist_n.tif")
base_solo <- raster("base_solo.tif")
base_unit <- raster("base_unit.tif")
base_df <- stack(base_dec, base_tpi, base_wet, base_vis, base_dist, base_solo, base_unit)
#base_df <- stack(base_ind, base_amp, base_wet, base_vis)
alagoas <- read_sf("AL_UF_2022.shp")
#mde <- mask(x=base,mask=alagoas)

# Stack Correlation Check
#check <- layerStats(base_df, 'pearson', na.rm=T)
#corr_matrix=check$'pearson correlation coefficient'

# Parameters
size_rd <- 1000
t <- 100
n <- 100
level <- 0.95
config <- c(size_rd,t,n, level)
no_cores <- detectCores(logical = TRUE)
seq_id_all<- seq_along(1:config[2])
df_var <- c("dec","tpi","wet","vis", "dist", "solo", "unit")
#df_var <- c("ind", "amp", "wet", "vis")

# Define target
reg <- read_sf("AL_Municipios_2022.shp")
#random_pt <- sample.int(nrow(reg),1)
target <- crop(base_df,extent(reg[100,]))
target_tmp <- rasterToPoints(target)
target_pt <- SpatialPoints(target_tmp[,1:2])
target_df <- target_tmp[,3:9]
target_df <- as.data.frame(target_df)
target_df[is.na(target_df)] <- 0
colnames(target_df) <- df_var

# Load CNSA
#for (j in 1:2) {
cnsa <- read_sf("sitios545.shp")
opt <- levels(as.factor(cnsa$Cat2))
cnsa <- subset(cnsa, cnsa$Cat2==opt[7])

train <- subset(cnsa, cnsa$NM_MUN != "Traipu")
test <- subset(cnsa, cnsa$NM_MUN == "Traipu")


#smp_size <- floor(0.75 * nrow(cnsa))
#train_ind <- sample(seq_len(nrow(cnsa)), size = smp_size)
#train <- cnsa[train_ind,]
#test <- cnsa[-train_ind,]
train_df <- terra::extract(base_df,train)
test_df <- terra::extract(base_df,test)
train_df <- as.data.frame(train_df)
test_df <- as.data.frame(test_df)
colnames(train_df) <- df_var
colnames(test_df) <- df_var

# Preparation
randompoint <- st_sample(reg[100,],t)
randombuf <- st_buffer(randompoint,config[1])
res_df <- list()
for (i in 1:length(randombuf)) {
  mdex_alt <- crop(base_df, as(randombuf[i], "Spatial"))
  tmp_pt <- rasterToPoints(mdex_alt)
  pt <- SpatialPoints(tmp_pt[,1:2])
  res_tmp <- tmp_pt[,3:9]
  colnames(res_tmp) <- df_var
  res_tmp[is.na(res_tmp)] <- 0
  res_tmp <- as.data.frame(res_tmp)
  res_tmp$solo <- as.factor(as.character(res_tmp$solo))
  res_tmp$unit <- as.factor(as.character(res_tmp$unit))
  res_df[[i]] <- res_tmp
}

# Functions
mdia_hcpc <- function(x) {
  library(FactoMineR)
  respca <- PCA(x,graph=F)
  reshcpc <- HCPC(respca, graph=F)
  res_clust <- reshcpc$data.clust$clust
  return(res_clust)
}


mdia_famd <- function(x) {
  library(FactoMineR)
  resfamd <- FAMD(x, graph=F)
  reshcpc <- HCPC(resfamd, graph=F)
  res_clust <- reshcpc$data.clust$clust
  return(res_clust)
}

mdia_nbclust <- function(x) {
  library(NbClust)
  # Classification
  res_clust <- NbClust(data=x, distance="euclidean",method="ward.D2", min.nc=2, max.nc=6)
  res_clust <- as.factor(res_clust$Best.partition)
  return(res_clust)
}

mdia_fuzzy <- function(x) {
  library(fclust)
  res_fkm <- FKM(X=x, stand=1, RS=10)
  res_clust <- as.factor(res_fkm$clus[,1])
  return(res_clust)
}

model_general <- function(x) {
  library(randomForest)
  library(caret)
  tuneGrid <- expand.grid(.mtry = c(1: 10))
  trControl <- trainControl(method = "cv", number = 10, search = "grid")
  model <- train(x[,1:7],x[,8],method="rf", metric="Accuracy", tuneGrid=tuneGrid, trControl=trControl, importance=T,nodesize=10,ntree=config[3])
  return(model)
}

predict_sites <- function(x) {
  library(randomForest)
  library(caret)
  tmp <- vector()
  train_df <- train_df[(train_df$solo %in% x$trainingData$solo),]
  train_df$solo <- as.factor(as.character(train_df$solo))
  levels(train_df$solo) <- levels(x$trainingData$solo)
  train_df <- train_df[(train_df$unit %in% x$trainingData$unit),]
  train_df$unit <- as.factor(as.character(train_df$unit))
  levels(train_df$unit) <- levels(x$trainingData$unit)
  if(nrow(train_df) > 0) {
    train_predict <- predict(x$finalModel, train_df,progress="text", predict.all=TRUE)
    for(i in 1:nrow(train_df)) {
      tmp[i] <- ifelse(max(as.numeric(table(train_predict$individual[i,])))>=config[3]*config[4],1,NA)
    }
    prob <- as.data.frame(table(as.numeric(train_predict$aggregate)*tmp))
    prob$def <- prob$Freq/nrow(train_df)*100
    return(prob) 
  }
}

predict_zones <- function(x) {
  library(randomForest)
  library(caret)
  tmp <- vector()
  target_df$nrow <- seq.int(nrow(target_df))
  target_df <- target_df[(target_df$solo %in% x$trainingData$solo),]
  target_df$solo <- as.factor(as.character(target_df$solo))
  levels(target_df$solo) <- levels(x$trainingData$solo)
  target_df <- target_df[(target_df$unit %in% x$trainingData$unit),]
  target_df$unit <- as.factor(as.character(target_df$unit))
  levels(target_df$unit) <- levels(x$trainingData$unit)
  res2 <- target_df
  target_df <- target_df[,-9]
  if(nrow(target_df) > 0) {
    target_predict <- predict(x$finalModel, target_df, progress="text", predict.all=TRUE)
    for (i in 1:nrow(target_predict$individual)) {
      tmp[i] <- ifelse(max(as.numeric(table(target_predict$individual[i,])))>=config[3]*config[4],1,NA)
    }
    res <- vector()
    for (i in 1:nrow(target_predict$individual)) {
      res[i] <- as.numeric(target_predict$aggregate[i]) * tmp[i]
    }
    return(list(res, res2))
  }
}



# Parallel Functions
cl <- makeCluster(4)
registerDoParallel(cl)
start_time <- Sys.time()
clusterExport(cl,list('mdia_hcpc','config'))
results<- c(parLapply(cl,res_df,fun=mdia_hcpc))
stopCluster(cl)
end_time <- Sys.time()
end_time - start_time
for (i in 1:t) {
  res_df[[i]]$Cluster <- results[[i]]
}

rm(results,mdex_alt,target_pt,target_tmp,base,base_amp,base_asp,base_dec,base_ind,base_tpi,base_tri,base_vis,base_wet)

cl <- makeCluster(4)
registerDoParallel(cl)
start_time <- Sys.time()
clusterExport(cl,list('model_general','config'))
results_pred<- c(parLapply(cl,res_df, fun=model_general))
stopCluster(cl)
end_time <- Sys.time()
end_time - start_time

cl <- makeCluster(4)
registerDoParallel(cl)
start_time <- Sys.time()
clusterExport(cl,list('predict_sites','config', 'train_df'))
results_sites<- c(parLapply(cl,results_pred, fun=predict_sites))
stopCluster(cl)
end_time <- Sys.time()
end_time - start_time

cl <- makeCluster(4)
registerDoParallel(cl)
start_time <- Sys.time()
clusterExport(cl,list('predict_zones','config', 'target_df'))
results_zones<- c(parLapply(cl,results_pred, fun=predict_zones))
stopCluster(cl)
end_time <- Sys.time()
end_time - start_time

# End
for (i in 1:100) {
  for (k in 1:max(results_zones[[i]])) {
    if(!is.na(match(k, results_sites[[i]]$Var1))) {
      results_zones[[i]][results_zones[[i]] == k] <- results_sites[[i]]$def[results_sites[[i]]$Var1 == k]}
    else {
      results_zones[[i]][results_zones[[i]] == k] <- NA
    }
  }
}

# Run full cycle
cl <- makeCluster(7)
registerDoParallel(cl)
start_time <- Sys.time()
clusterExport(cl,list('mdia_famd','config'))
results<- c(parLapply(cl,res_df,fun=mdia_famd))
for (i in 1:t) {
  res_df[[i]]$Cluster <- results[[i]]
}
clusterExport(cl,list('model_general','config'))
results_pred<- c(parLapply(cl,res_df, fun=model_general))
clusterExport(cl,list('predict_sites','config', 'train_df'))
results_sites<- c(parLapply(cl,results_pred, fun=predict_sites))
clusterExport(cl,list('predict_zones','config', 'target_df'))
results_zones<- c(parLapply(cl,results_pred, fun=predict_zones))
stopCluster(cl)
for (i in 1:100) {
  for (k in 1:max(results_zones[[i]][[1]], na.rm=T)) {
    if(!is.na(match(k, results_sites[[i]]$Var1))) {
      results_zones[[i]][[1]][results_zones[[i]][[1]] == k] <- results_sites[[i]]$def[results_sites[[i]]$Var1 == k]}
    else {
      results_zones[[i]][[1]][results_zones[[i]][[1]] == k] <- 0
    }
  }
}
final_df <- vector()
target_df$nrow <- seq.int(nrow(target_df))
for(i in 1:100) {
  tmp <- results_zones[[i]]
  output1 <- tmp[[1]]
  output2 <- tmp[[2]]$nrow
  output <- cbind(output1, output2)
  tmp2 <- output[,1][match(target_df$nrow,output[,2])]
  final_df <- cbind(final_df,tmp2)
}
#final_df[is.na(final_df)] <- 55
r <- raster()
ncol(r) <- ncol(target)
nrow(r) <- nrow(target)
crs(r) <- crs(target)
extent(r) <- extent(target)
values(r) <- rowMeans(as.data.frame(final_df))
end_time <- Sys.time()
end_time - start_time

#ROC-AUC
library(pROC)
res_arq <- as.data.frame(extract(r, test))
res_rd <- as.data.frame(extract(r, as(st_sample(reg[100,], 14), "Spatial")))
res_arq$response <- 1
res_rd$response <- 0
res_arq$Type <- "Arq"
res_rd$Type <- "rd"
colnames(res_arq) <- c("predictor", "response", "type")
colnames(res_rd) <- c("predictor", "response", "type")
res_fin <- rbind(res_arq, res_rd)
res_auc <- roc(res_fin$response, res_fin$predictor, plot=T, auc=T)