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
library(rgeos)

# Load data
setwd("c://Users/gvanh/Documents/UFAL/mod/")
base <- raster("base.tif")
base_ind <- raster("base_ind_n.tif")
#base_amp <- raster("base_amp_n.tif")
base_dec <- raster("base_dec_n.tif")
base_asp <- raster("base_asp_n.tif")
base_tpi <- raster("base_tpi_n.tif")
#base_tri <- raster("base_tri_n.tif")
base_wet <- raster("base_wet_n.tif")
base_vis <- raster("base_vis_n.tif")
base_dist <- raster("base_dist_n.tif")
#base_df <- stack(base_ind, base_amp, base_dec, base_tpi, base_tri, base_wet, base_vis)
base_df <- stack(base, base_ind, base_dec, base_tpi, base_wet, base_vis, base_dist)
alagoas <- read_sf("AL_UF_2022.shp")
mde <- mask(x=base,mask=alagoas)

# Parameters
size_rd <- 1000
t <- 100
config <- c(size_rd,t)
no_cores <- detectCores(logical = TRUE)
seq_id_all<- seq_along(1:config[2])
df_var <- c("alt", "ind", "dec", "tpi", "wet", "vis", "dist")

# Make State wide
#target_tmp <- rasterToPoints(base_df)
#target_pt <- SpatialPoints(target_tmp[,1:2])
#target_df <- target_tmp[,3:9]
#target_df[is.na(target_df)] <- 0
#colnames(target_df) <- df_var

# Define random target
mun <- read_sf("area2.shp")
#random_pt <- sample.int(nrow(mun),1)
target <- crop(base_df,extent(mun))
target_tmp <- rasterToPoints(target)
target_pt <- SpatialPoints(target_tmp[,1:2])
target_df <- target_tmp[,3:9]
target_df <- target_tmp[,3:9]
target_df[is.na(target_df)] <- 0
colnames(target_df) <- df_var

# Load CNSA
#for (j in 1:2) {
  cnsa <- read_sf("sitios2.shp")
  opt <- levels(as.factor(cnsa$Class_cron))
  cnsa <- subset(cnsa, cnsa$Class_cron==opt[2])
  smp_size <- floor(0.75 * nrow(cnsa))
  train_ind <- sample(seq_len(nrow(cnsa)), size = smp_size)
  cnsa_df <- terra::extract(base_df,cnsa[train_ind, ])
  colnames(cnsa_df) <- df_var

  # Prep
  randompoint <- st_sample(alagoas,100)
  randombuf <- st_buffer(randompoint,config[1])
  res_df <- list()
  for (i in 1:length(randombuf)) {
    mdex_alt <- crop(base_df, as(randombuf[i], "Spatial"))
    tmp_pt <- rasterToPoints(mdex_alt)
    pt <- SpatialPoints(tmp_pt[,1:2])
    #res_tmp <- tmp_pt[,3:9]
    res_tmp <- tmp_pt[,3:9]
    colnames(res_tmp) <- df_var
    res_tmp[is.na(res_tmp)] <- 0
    res_tmp <- as.data.frame(res_tmp)
    res_df[[i]] <- res_tmp
  }
  
  # Prep Alt
  randombuf <- list()
  for (i in 1:100) {
    randombuf[i] <- gBuffer(sampleRandom(mde, size=1, sp=T), width=0.01)  
  }
  res_df <- list()
  for (i in 1:length(randombuf)) {
    mdex_alt <- crop(base_df, randombuf[[i]])
    tmp_pt <- rasterToPoints(mdex_alt)
    pt <- SpatialPoints(tmp_pt[,1:2])
    res_tmp <- tmp_pt[,3:9]
    colnames(res_tmp) <- df_var
    res_tmp[is.na(res_tmp)] <- 0
    res_tmp <- as.data.frame(res_tmp)
    res_df[[i]] <- res_tmp
  }
  
  # Function
  calculate_mdia <- function(x) {
    library(FactoMineR)
    library(randomForest)
    library(caret)
    library(NbClust)
    library(fclust)
    # Classification
    #res_clustA <- NbClust(data=x, distance="euclidean",method="ward.D2", min.nc=2, max.nc=6)
    #res_clustA <- as.factor(res_clustA$Best.partition)
    respca <- PCA(x,graph=F)
    reshcpc <- HCPC(respca, graph=F)
    res_clustB <- reshcpc$data.clust$clust
    #res_fkm <- FKM(X=x, stand=1, RS=10)
    #res_clustB <- as.factor(res_fkm$clus[,1])
    # Training
    tuneGrid <- expand.grid(.mtry = c(1: 10))
    trControl <- trainControl(method = "cv", number = 10, search = "grid")
    model <- train(x,res_clustB,method="rf", metric="Accuracy", tuneGrid=tuneGrid, trControl=trControl, importance=T,nodesize=10,ntree=100)
    # Prediction
    #final <- predict(model$finalModel, cnsa_df[,1:7],progress="text")
    final <- predict(model$finalModel, cnsa_df[,1:7],progress="text", predict.all=TRUE)
    tmp <- vector()
    for(i in 1:nrow(cnsa_df)) {
      tmp[i] <- ifelse(max(as.numeric(table(final$individual[i,])))>=95,1,0)
    }
    prob <- as.data.frame(table(as.numeric(final$aggregate) * tmp))
    prob$def <- prob$Freq/nrow(cnsa_df)*100
    #final2 <- predict(model$finalModel, target_df[,1:7], progress="text")
    final2 <- predict(model$finalModel, target_df[,1:7], progress="text", predict.all=TRUE)
    # Get probabilities
    tmp <- vector()
    for(i in 1:nrow(target_df)) {
      tmp[i] <- ifelse(max(as.numeric(table(final2$individual[i,])))>=95,1,0)
    }
    res <- as.numeric(final2$aggregate) * tmp
    #res_cnsa <- cbind(cnsa_df,final)
    #prob <- as.data.frame(table(final))
    #for (i in 1:max(as.integer(final))) {
     # prob$weight[i] <- mean(dist(subset(res_cnsa,final==i)))
      #prob$weightinv[i] <- 1/mean(dist(subset(res_cnsa,final==i)))
    #}
    #prob$def <- prob$Freq*prob$weightinv
    #prob$def[!is.finite(prob$def)] <- 0
    # Insert probabilities to area
    #res <- as.numeric(final2)
    for (k in 1:max(res)) {
      if ((sum(res[res == k]) != 0) && (sum(prob$def[prob$Var1 == k]) != 0)) {res[res == k] <- prob$def[prob$Var1 == k]}
    }
    #res[is.na(res)] <- 0
    #res <- ifelse(res==max(res),1,0)
    return(res)
    # End for-loop
  }
  
  cl <- makeCluster(no_cores-1)
  registerDoParallel(cl)
  start_time <- Sys.time()
  clusterExport(cl,list('calculate_mdia','config','target_df','alagoas','base_df','cnsa_df','df_var'))
  results<- c(parLapply(cl,res_df,fun=calculate_mdia))
  stopCluster(cl)
  end_time <- Sys.time()
  end_time - start_time

  
  # Make final
  r <- raster()
  ncol(r) <- ncol(target)
  nrow(r) <- nrow(target)
  crs(r) <- crs(base_df)
  extent(r) <- extent(target)
  values(r) <- rowMeans(as.data.frame(results))
  df_minmax <- minmax(rast(r))
  df_norm <- (r - df_minmax[1,]) / (df_minmax[2,] - df_minmax[1,])
  #df_norm <- mask(x=df_norm,mask=mun[random_pt,2])
  plot(df_norm, zlim=c(0,1), main=mun[random_pt,2]$NM_MUN)
