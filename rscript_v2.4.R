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
setwd("c://Users/gvanh/Documents/UFAL/mod/")
base <- raster("base.tif")
base_ind <- raster("base_ind_n.tif")
base_amp <- raster("base_amp_n.tif")
base_dec <- raster("base_dec_n.tif")
base_asp <- raster("base_asp_n.tif")
base_tpi <- raster("base_tpi_n.tif")
base_tri <- raster("base_tri_n.tif")
base_wet <- raster("base_wet_n.tif")
base_vis <- raster("base_vis_n.tif")
base_df <- stack(base_ind, base_amp, base_dec, base_tpi, base_tri, base_wet, base_vis)
alagoas <- read_sf("AL_UF_2022.shp")
mde <- mask(x=base,mask=alagoas)

# Parameters
size_rd <- 1000
t <- 100
config <- c(size_rd,t)
no_cores <- detectCores(logical = TRUE)
cl <- makeCluster(no_cores-1)
registerDoParallel(cl)
seq_id_all<- seq_along(1:config[2])
df_var <- c("ind", "amp", "dec", "tpi", "tri", "wet", "vis")

# Make State wide
#target_tmp <- rasterToPoints(base_df)
#target_pt <- SpatialPoints(target_tmp[,1:2])
#target_df <- target_tmp[,3:9]
#target_df[is.na(target_df)] <- 0
#colnames(target_df) <- df_var

# Define random target
mun <- read_sf("AL_Municipios_2022.shp")
random_pt <- sample.int(nrow(mun),1)
target <- crop(base_df,extent(mun[random_pt,1]))
target_tmp <- rasterToPoints(target)
target_pt <- SpatialPoints(target_tmp[,1:2])
target_df <- target_tmp[,3:9]
target_df[is.na(target_df)] <- 0
colnames(target_df) <- df_var

# Load CNSA
#for (j in 1:2) {
  cnsa <- read_sf("sitios2.shp")
  opt <- levels(as.factor(cnsa$Class_cron))
  cnsa <- subset(cnsa, cnsa$Class_cron==opt[4])
  smp_size <- floor(0.75 * nrow(cnsa))
  train_ind <- sample(seq_len(nrow(cnsa)), size = smp_size)
  cnsa_df <- terra::extract(base_df,cnsa[train_ind, ])
  colnames(cnsa_df) <- df_var

  
  # Function
  calculate_mdia <- function(...) {
    library(raster)
    library(FactoMineR)
    library(sf)
    library(sp)
    library(randomForest)
    library(caret)
    library(terra)
    library(NbClust)
    library(fclust)
    # Get random coordinates
    i <- (...)
    randompoint <- st_sample(alagoas,1)
    randombuf <- st_buffer(randompoint,config[1])
    # Crop raster
    mdex_alt <- crop(base_df, as(randombuf,"Spatial"))
    # Extraction
    tmp_pt <- rasterToPoints(mdex_alt)
    pt <- SpatialPoints(tmp_pt[,1:2])
    res_df <- tmp_pt[,3:9]
    colnames(res_df) <- df_var
    res_df[is.na(res_df)] <- 0
    res_df <- as.data.frame(res_df)
    # Classification
    #res_clustA <- NbClust(data=res_df, distance="euclidean",method="ward.D2", min.nc=2, max.nc=6)
    #res_clustA <- as.factor(res_clustA$Best.partition)
    #respca <- PCA(res_df,graph=F)
    #reshcpc <- HCPC(respca, graph=F)
    #res_clustB <- reshcpc$data.clust$clust
    res_fkm <- FKM(X=res_df, stand=1, RS=10)
    res_clustC <- as.factor(res_fkm$clus[,1])
    # Training
    tuneGrid <- expand.grid(.mtry = c(1: 10))
    trControl <- trainControl(method = "cv", number = 10, search = "grid")
    model <- train(res_df,res_clustC,method="rf", metric="Accuracy", tuneGrid=tuneGrid, trControl=trControl, importance=T,nodesize=10,ntree=100)
    # Prediction
    final <- predict(model$finalModel, cnsa_df[,1:7],progress="text")
    final2 <- predict(model$finalModel, target_df[,1:7], progress="text")
    # Get probabilities
    prob <- as.data.frame(table(final))
    prob$def <- prob$Freq/nrow(cnsa_df)*100
    # Insert probabilities to area
    res <- as.numeric(final2)
    for (k in 1:max(res)) {
      res[res == k] <- prob$def[k]
    }
    return(res)
    # End for-loop
  }
  
  # Calculate
  start_time <- Sys.time()
  clusterExport(cl,list('calculate_mdia','config','target_df','alagoas','base_df','cnsa_df','df_var'))
  results<- c(parLapply(cl,seq_id_all,fun=my.function))
  end_time <- Sys.time()
  end_time - start_time
  
  # Make final
  r <- raster()
  ncol(r) <- ncol(target)
  nrow(r) <- nrow(target)
  crs(r) <- crs(target)
  extent(r) <- extent(target)
  values(r) <- rowMeans(as.data.frame(results))
  df_minmax <- minmax(rast(r))
  df_norm <- (r - df_minmax[1,]) / (df_minmax[2,] - df_minmax[1,])
  df_norm <- mask(x=df_norm,mask=mun[random_pt,2])
  plot(df_norm, main=mun[random_pt,2]$NM_MUN)
