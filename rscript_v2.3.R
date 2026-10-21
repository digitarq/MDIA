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
start_time <- Sys.time()
setwd("c://Users/Kleython/Documents/greg/modelo/")
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

# Prepare cores
no_cores <- detectCores(logical = TRUE)
cl <- makeCluster(no_cores-1)
registerDoParallel(cl)
seq_id_all<- seq_along(1:config[2])

# Make State wide
target_tmp <- rasterToPoints(base_df)
target_pt <- SpatialPoints(target_tmp[,1:2])
target_df <- target_tmp[,3:9]
#target_df <- terra::extract(base_df,target_pt)
target_df[is.na(target_df)] <- 0
colnames(target_df) <- c("ind", "amp", "dec", "tpi", "tri", "wet", "vis")

# Define random target
#mun <- read_sf("AL_Municipios_2022.shp")
#random_pt <- sample.int(nrow(mun),1)
#target <- crop(base_df,extent(mun[random_pt,1]))
#target_pt <- SpatialPoints(rasterToPoints(target))
#target_df <- as.data.frame(extract(target,target_pt))
#colnames(target_df) <- c("ind", "amp", "dec", "tpi", "tri", "wet", "vis")

# Load CNSA
#for (j in 1:2) {
  cnsa <- read_sf("sitios2.shp")
  opt <- c("PH", "H")
  cnsa <- subset(cnsa, cnsa$Class_cron==opt[1])
  smp_size <- floor(0.75 * nrow(cnsa))
  train_ind <- sample(seq_len(nrow(cnsa)), size = smp_size)
  cnsa_df <- terra::extract(base_df,cnsa[train_ind, ])
  colnames(cnsa_df) <- c("ind", "amp", "dec", "tpi", "tri", "wet", "vis")
  
  # Parameters
  size_rd <- 1000
  t <- 10
  nb_cluster <- 30
  config <- c(size_rd,t,nb_cluster)
  
  # Function
  my.function <- function(...) {
    library(raster)
    library(FactoMineR)
    library(sf)
    library(sp)
    library(randomForest)
    library(caret)
    library(terra)
    # Get random coordinates
    i <- (...)
    randompoint <- st_sample(alagoas,1)
    randombuf <- st_buffer(randompoint,config[1])
    # Crop raster
    mdex_alt <- crop(base_df, as(randombuf,"Spatial"))
    # Extraction
    pt <- SpatialPoints(rasterToPoints(mdex_alt))
    res_df <- terra::extract(mdex_alt,pt)
    colnames(res_df) <- c("ind", "amp", "dec", "tpi", "tri", "wet", "vis")
    res_df[is.na(res_df)] <- 0
    res_df <- as.data.frame(res_df)
    # Classification
    respca <- PCA(res_df,graph=F)
    reshcpc <- HCPC(respca, graph=F)
    res_clust <- reshcpc$data.clust$clust
    # Training
    tuneGrid <- expand.grid(.mtry = c(1: 10))
    trControl <- trainControl(method = "cv", number = 10, search = "grid")
    model <- train(res_df,res_clust,method="rf", metric="Accuracy", tuneGrid=tuneGrid, trControl=trControl, importance=T,nodesize=10,ntree=100)
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
  clusterExport(cl,list('my.function','config','target_df','alagoas','base_df','cnsa_df'))
  results<- c(parLapply(cl,seq_id_all,fun=my.function))
  end_time <- Sys.time()
  end_time - start_time
  
  # Make final
  res_r <- stack()
  r <- raster()
  ncol(r) <- ncol(mde)
  nrow(r) <- nrow(mde)
  crs(r) <- crs(mde)
  extent(r) <- extent(mde)
  for (i in 1:config[2]){
  values(r) <- results[[i]]
  res_r <- stack(res_r,r)
  }
  df_mean <- calc(res_r,mean)
  df_minmax <- minmax(rast(df_mean))
  df_norm <- (df_mean - df_minmax[1,]) / (df_minmax[2,] - df_minmax[1,])
