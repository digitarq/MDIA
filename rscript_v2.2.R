# Load libraries
library(raster)
library(sp)
library(sf)
library(caret)
library(FactoMineR)
library(rgeos)
library(randomForest)
library(viridis)
library(spatialEco)
library(terra)

# Load data
setwd("c://Users/Kleython/Documents/greg/modelo/")
base <- raster("base.tif")
base_dec <- raster("base_dec.tif")
base_asp <- raster("base_asp.tif")
base_tpi <- raster("base_tpi.tif")
base_tri <- raster("base_tri.tif")
base_fete <- raster("base_fete20.tif")
base_df <- stack(base, base_dec, base_asp, base_tpi, base_tri, base_fete)
alagoas <- read_sf("AL_UF_2022.shp")
mde <- mask(x=base,mask=alagoas)

# Define random target
mun <- read_sf("AL_Municipios_2022.shp")
random_pt <- sample.int(nrow(mun),1)
target <- crop(base_df,extent(mun[random_pt,1]))
target_pt <- SpatialPoints(rasterToPoints(target))
target_df <- as.data.frame(extract(target,target_pt))
colnames(target_df) <- c("mde", "dec", "asp", "tpi", "tri", "fete")

# Load CNSA
for (j in 1:2) {
  cnsa <- read_sf("sitios2.shp")
  #opt <- c("PH", "H")
  #cnsa <- subset(cnsa, cnsa$Class_cron==opt[j])
  cnsa <- cnsa[,3]
  cnsa_df <- extract(base_df,cnsa)
  colnames(cnsa_df) <- c("mde", "dec", "asp", "tpi", "tri", "fete")
  
  # Parameters
  size_rd <- 2000
  t <- 10
  nb_cluster <- 5
  
  # Create stack
  res <- stack()
  
  # Start for-loop
  pb <- txtProgressBar(min = 0, max = t, initial = 0, style=3) 
  for (i in 1:t) {
    
    # Get random coordinates
    randompoint <- st_sample(alagoas,1)
    randombuf <- st_buffer(randompoint,size_rd)
    
    # Crop raster
    mdex_alt <- crop(base, as(randombuf,"Spatial"))
    
    #if(ncol(mdex_alt) <= size_rd) {
      
      # Extraction
      pt <- SpatialPoints(rasterToPoints(mdex_alt))
      res_df <- extract(base_df,pt)
      colnames(res_df) <- c("mde", "dec", "asp", "tpi", "tri", "fete")
      res_df[is.na(res_df)] <- 0
      res_df <- as.data.frame(res_df)
      
      # Check-in
      if(sum(res_df$mde) > 0) {
        
        # Classification
        respca <- PCA(res_df,graph=F)
        reshcpc <- HCPC(respca, nb.clust=nb_cluster, graph=F)
        res_clust <- reshcpc$data.clust$clust
        res_df <- cbind(res_df,res_clust)
        rm(res_clust)
        
        # Training
        x <- res_df[,-7]
        y <- as.factor(res_df[,7])
        tuneGrid <- expand.grid(.mtry = c(1: 10))
        trControl <- trainControl(method = "cv", number = 10, search = "grid")
        model <- train(x,y,method="rf", metric="Accuracy", tuneGrid=tuneGrid, trControl=trControl, importance=T,nodesize=10,ntree=100)
        
        # Prediction
        final <- predict(model$finalModel, cnsa_df[,1:6],progress="text")
        final2 <- predict(model$finalModel, target_df[,1:6], progress="text")
        
        # Get probabilities
        prob <- as.data.frame(table(final))
        prob$def <- prob$Freq/nrow(cnsa_df)*100
        
        # Insert probabilities to area
        xxx <- as.numeric(final2)
        for (k in 1:nb_cluster) {
          xxx[xxx == k] <- prob$def[k]
        }
        
        # Rasterize and add to Stack
        r <- raster()
        ncol(r) <- ncol(target)
        nrow(r) <- nrow(target)
        crs(r) <- crs(target)
        extent(r) <- extent(target)
        r <- rasterize(target_pt@coords[,1:2],r,xxx)
        names(r) <- paste('run',i)
        res <- stack(res,r)
        
        # Print progress
        setTxtProgressBar(pb,i)
      }
 #} 
    # End for-loop
  }
  res2 <- mask(x=res,mask=mun[random_pt,1])
  png(file=paste("res/",mun$NM_MUN[random_pt]," - ",opt[j],sprintf("%02d", i),".png",sep=""), width = 1200, height = 1000)
  par(mfrow=c(1,2), oma=c(0,0,2,0))
  plot(calc(res2,median), zlim=c(0,100),col=viridis(1e3), main="Mediana", sub=opt[j])
  plot((calc(res2,median)/maxValue(calc(res2,median)))*100, zlim=c(0,100),col=viridis(1e3), main="Mediana normalizada", sub=opt[j])
  title(mun$NM_MUN[random_pt], line=0, side=3, outer=TRUE, cex=2)
  dev.off()
  close(pb)
}