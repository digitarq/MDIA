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
library(tidyr)
library(scales)
library(pROC)

# Load data
setwd("C://Documents and Settings/User/Documents/UFPI/Pesquisa/Alagoas/2025/v09/")
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
alagoas <- read_sf("AL_UF_2022.shp")

# Parameters
size_rd <- 1000
t <- 10
n <- 100
level <- 0.75
config <- c(size_rd,t,n,level)
nb_mun <- 68
no_cores <- detectCores(logical = TRUE)
seq_id_all<- seq_along(1:config[2])
df_var <- c("dec","tpi","wet","vis", "dist", "solo", "unit")

# Functions
mdia_hcpc <- function(x) {
  library(FactoMineR)
  respca <- PCA(x,graph=F)
  reshcpc <- HCPC(respca, graph=F)
  res_clust <- 
  return(res_hcpc)
}

mdia_famd <- function(x) {
  library(FactoMineR)
  resfamd <- FAMD(x, graph=F)
  reshcpc <- HCPC(resfamd, graph=F)
  return(reshcpc)
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
    return(list(prob, train_predict, train_df))
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
    return(list(res, res2, target_predict))
  }
}

create_raster <- function(x) {
  x <- raster(ncol=ncol(target), nrow=nrow(target), crs=crs(target), ext=extent(target))
}

# Define target
reg <- read_sf("AL_Municipios_2022.shp")
t_mun <- reg[nb_mun,]
target <- crop(base_df,extent(t_mun))
target_tmp <- rasterToPoints(target)
target_pt <- SpatialPoints(target_tmp[,1:2])
target_df <- target_tmp[,3:9]
target_df <- as.data.frame(target_df)
target_df[is.na(target_df)] <- 0
colnames(target_df) <- df_var
  
  # Preparation
  randompoint <- st_sample(t_mun,t)
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
  
  cl <- makeCluster(no_cores - 1)
  registerDoParallel(cl)
  clusterExport(cl,list('mdia_famd','config'))
  results<- c(parLapply(cl,res_df,fun=mdia_famd))
  for (i in 1:t) {
    res_df[[i]]$Cluster <- results[[i]]$data.clust$clust
  }
  clusterExport(cl,list('model_general','config'))
  results_pred<- c(parLapply(cl,res_df, fun=model_general))
  clusterExport(cl,list('predict_zones','config', 'target_df'))
  init_zones<- c(parLapply(cl,results_pred, fun=predict_zones))
  
  infodata <- list()
  list_auc <- vector()
  final_df <- list()
  r_fin <- stack()
  oob_error <- list()
  table_sites <- list()
  table_zones <- list()
  cnsa <- read_sf("df_cnsa.shp")
  opt <- levels(as.factor(cnsa$L2))
  opt <- opt[opt != "NA"]
  for (k in 1:length(opt)) {
  # Load CNSA
  infodata[[k]] <- vector()
  train <- subset(cnsa, cnsa$L2==opt[k])
  train$dup <- duplicated(train[c("Chave", "L2")])
  train <- subset(train, train$dup == FALSE)
  infodata[[k]][1] <- nrow(train)
  train <- subset(train, train$NM_MUN != t_mun$NM_MUN)
  infodata[[k]][2] <- nrow(train)
  test <- subset(cnsa, cnsa$L1==opt[k])
  test$dup <- duplicated(test[c("Chave", "L2")])
  test <- subset(test, test$dup == FALSE)
  test <- subset(test, test$NM_MUN == t_mun$NM_MUN)
  infodata[[k]][3] <- nrow(test)
  train_df <- as.data.frame(terra::extract(base_df,train))
  colnames(train_df) <- df_var
  if(nrow(test)>0){
  test_df <- terra::extract(base_df,test)
  test_df <- as.data.frame(test_df)
  colnames(test_df) <- df_var
  }
  
  # Run full cycle
  clusterExport(cl,list('predict_sites','config', 'train_df'))
  results_sites<- c(parLapply(cl,results_pred, fun=predict_sites))
  results_zones <- init_zones
  for (i in 1:t) {
    if(!is.null(results_zones[[i]]) & !is.na(max(results_zones[[i]][[1]], na.rm=T)) & length(results_zones[[i]][[1]]) != sum(is.na(results_zones[[i]][[1]]))) {
      for (j in 1:max(results_zones[[i]][[1]], na.rm=T)) {
        if(!is.na(match(j, results_sites[[i]][[1]]$Var1))) {
          results_zones[[i]][[1]][results_zones[[i]][[1]] == j] <- results_sites[[i]][[1]]$def[results_sites[[i]][[1]]$Var1 == j]}
        else {
          results_zones[[i]][[1]][results_zones[[i]][[1]] == j] <- 0
        }
      }
      }
    else {
      results_zones[[i]][[1]] <- rep(0, nrow(target_df))
      results_zones[[i]][[2]] <- seq_along(1:nrow(target_df))      
    }
  }
  final_df[[k]] <- vector()
  target_df$nrow <- seq.int(nrow(target_df))
  for(i in 1:t) {
    tmp <- results_zones[[i]]
    output1 <- tmp[[1]]
    if(length(tmp[[2]]) == 8) {output2 <- tmp[[2]]$nrow} else {output2 <- tmp[[2]]}
    output <- cbind(output1, output2)
    tmp2 <- output[,1][match(target_df$nrow,output[,2])]
    final_df[[k]] <- cbind(final_df[[k]],tmp2)
  }
  r_stack <- stack()
  r_res <- create_raster(x)
  r_sig <- create_raster(x)
  r_temp <- create_raster(x)
  values(r_res) <- apply(final_df[[k]], 1, FUN=mean, na.rm=T)
  sig_df <- rowSums(as.data.frame(!is.na(final_df[[k]])))
  values(r_sig) <- sig_df
  r_stack <- stack(r_res, r_sig)
  names(r_stack) <- c("Resultados", "Controle")
  df <- as.data.frame(values(r_stack))
  df$rank1 <- scale(df$Resultados)
  df$rank1 <- rescale(df$rank1, to=c(-1,1), from=range(df$rank1, na.rm=T))
  df$rank2 <- with(df, match(df$Controle, sort(unique(df$Controle, na.rm=T))))
  df$rank2 <- rescale(df$rank2, to=c(0,1), from=c(0,config[2]), na.rm=T)
  df$rank3 <- df$rank1*df$rank2
  values(r_temp) <- df$rank3
  r_fin <- stack(r_fin, r_temp)
  if(nrow(test) > 0 && !is.null(test)) {
    res_arq <- as.data.frame(terra::extract(r_temp, as(test, "Spatial")))
    res_rd <- as.data.frame(terra::extract(r_temp, as(st_sample(t_mun, nrow(test)), "Spatial")))
    res_arq$response <- 1
    res_rd$response <- 0
    res_arq$Type <- "Arq"
    res_rd$Type <- "rd"
    colnames(res_arq) <- c("predictor", "response", "type")
    colnames(res_rd) <- c("predictor", "response", "type")
    res_fin <- rbind(res_arq, res_rd)
    res_fin$predictor[is.na(res_fin$predictor)] <- 0
    res_auc <- roc(res_fin$response, res_fin$predictor, plot=F, auc=T)
    list_auc[k] <- res_auc$auc
  } else {
    list_auc[k] <- 0
  }
  tmp3 <- vector()
  for (i in 1:t) {
    tmp3[i] <- 1-(sum(diag(results_pred[[i]]$finalModel$confusion[,-ncol(results_pred[[i]]$finalModel$confusion)])/sum(results_pred[[i]]$finalModel$confusion[,-ncol(results_pred[[i]]$finalModel$confusion)])))
  }
  oob_error[[k]] <- as.data.frame(tmp3)
  tmp4 <- vector()
  for (i in 1:t) {
    if (is.null(results_sites[[i]][[1]])) {
      tmp4[i] <- NA
    }
    else {
      tmp4[i] <- sum(results_sites[[i]][[1]]$Freq)
    }
  }
  table_sites[[k]] <- tmp4
  tmp5 <- vector()
  for (j in 1:t) {
    if (is.null(results_zones[[j]][[1]])) {
      tmp5[j] <- NA
    }
    else {
      tmp5[j] <- round((length(na.omit(results_zones[[j]][[1]]))/ncell(target))*100)
    }
  }
  table_zones[[k]] <- tmp5
  }
  stopCluster(cl)
  names(r_fin) <- opt
  r_res <- r_fin
  r_def <- mean(r_res, na.rm=T)
  infodata <- t(as.data.frame(infodata))
  infodata <- as.data.frame(infodata)
  rownames(infodata) <- opt
  colnames(infodata) <- c("Total", "Train", "Test")  
  infodata$Proporção <- round((infodata$Total/sum(infodata$Total))*100)
  for (i in 1:6){
    r_res[i] <- r_res[i] * infodata[[4]][[i]]
  }
  coords_sample <- as(randompoint, "Spatial")@coords
  coords_sample <- cbind(seq(1:t), coords_sample)
  colnames(coords_sample) <- c("Pt", "X", "Y")
  localidade <- reg$NM_MUN[nb_mun]
  
save(infodata, t_mun, randompoint, list_auc, table_sites, table_zones, mdex_alt, target, config, target_df, target_tmp, train_df, coords_sample, oob_error, r_res, r_def, file="results.Rdata")
