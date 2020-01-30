#This code is for randomly selecting n-cells and running findMarkers in order to figure out
#whether a given cluster is real. 

packrat::init("Seuratv3", restart = TRUE, infer.dependencies = FALSE)

library(Seurat)
library(doParallel)

load("MP_IPC_ForPresentation.RData")

#First, I will permute the data. 
#second, I need to take a random subset of cells and make it identity 1 and the rest Identity 2
#Then I can run FindAllMarkers on this subset.
#I will permute the data 100 times, and for each permutation, do 100 sampling. 
#And I want this code in parallel. 
#Then I want to do it for the relative size of each cluster. 

#Set up Parallel processes
cl <- makeCluster(10)
registerDoParallel(cl)
#export data matrix

data <- as.matrix(MP_IPC3_2@assays$RNA@data)
scale.data <- as.matrix(MP_IPC3_2@assays$RNA@scale.data)
counts <- as.matrix(MP_IPC3_2@assays$RNA@counts)

object1 <- MP_IPC3_2

ThreshMarkers <- foreach(j=1:1000, .combine='rbind', .packages=c('foreach','Seurat')) %dopar%{
  
  #Duplicate each of the matrices for permutation. This is just being extra cautious. the parallelize foreach loop 
  #should create an isolated environment whereby changes will not affect the global variable. 
  
  #For each for loop, I'm generating a column resampling once. 
  #I then apply that to both counts and data, but not scaleData
  #For scale data (which is smaller) -> I find the gene that's being permuted
  #and apply the resampling for the same gene, if it's there in scale.data
  data1 <- data
  scale.data1 <- scale.data
  counts1 <- counts
  
  for(k in 1:dim(data)[1]){
    resample <-sample(length(data[k,]))
    data1[k,] <- data[k,resample]
    if(rownames(data)[k] %in% row.names(scale.data)){
      scale.data1[rownames(data)[k],] <- scale.data[rownames(data)[k],resample]
    }
  
  }
  object1@assays$RNA@counts = counts1
  object1@assays$RNA@data = data1
  object1@assays$RNA@scale.data = scale.data1
  object1 <- FindVariableFeatures(object1, selection.method = "vst", nfeatures = 2000)
  object1 <- RunPCA(object1, features = VariableFeatures(object = object1))
  object1 <- FindNeighbors(object1, dims = 1:10)
  
  foreach(m=seq(0.5, 4, by = 0.5), .combine='rbind', .packages=c('Seurat')) %do% {
    object1 <- FindClusters(object1, resolution = m, n.iter = 30)
    test <- FindAllMarkers(object1, logfc.threshold = 0.25, only.pos = TRUE)
    CellNum <- table(object1@active.ident)
    Resolution <-rep(m,dim(test)[1])
    Iteration <-rep(j, dim(test)[1])
    NumberCells <- as.character(test$cluster)
    NumberCells <- unlist(lapply(1:dim(test)[1], function (x) NumberCells[x] <- CellNum[test$cluster[x]]))
    Combined <- paste(Iteration, test$cluster, Resolution)
    cbind(test, Resolution, Iteration, Combined, NumberCells)
  }
  
}

write.csv(ThreshMarkers, file = "PermuteCells_ClusterResolution_ThreshMarkers.csv")

quit()
