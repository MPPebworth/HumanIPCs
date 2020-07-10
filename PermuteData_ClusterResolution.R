#This R Code is meant to take a single cell dataset, scrambled it 1000 times, and generate
# a list of all clusters markers, along with metadata. 
#This scrambled cluster marker list can then be used as a statistical threshold 
#bioinformatic clusters and cluster markers in the original dataset. 

#This code was run on Seurat version three.  

packrat::init("Seuratv3", restart = TRUE, infer.dependencies = FALSE)

library(Seurat)
library(doParallel)

#Seurat version 3 object with all IPCs, named IPC3
load("IPCs.RData")

#Step 1. I will permute the data for each gene across all cells, keeping the gene distribution the same, 
#but destroying the gene expression pattern of each cell. 
#Step 2: Generate clusters at different cluster resolutions to see how cluster size influences markers generated.
#Step 3: Run FindAllMarkers to find positive gene markers for these random clusters.

#Steps 1-4 will be repeated 1000x, in parallel across 30 processors. 

#Set up Parallel processes (30)
cl <- makeCluster(30, outfile = "")
registerDoParallel(cl)

#Extra data matrices from Seurat object for scrambling. 
data <- as.matrix(IPC3@assays$RNA@data)
scale.data <- as.matrix(IPC3@assays$RNA@scale.data)
counts <- as.matrix(IPC3@assays$RNA@counts)

#Create an object for manipulation in the foreach loop.
object1 <- IPC3

ThreshMarkers <- foreach(j=1:1000, .combine='rbind', .errorhandling = "remove",.packages=c('foreach','Seurat')) %dopar%{
  
  #Duplicate each of the matrices for permutation. This is just being extra cautious. The parallelize foreach loop 
  #should create an isolated environment whereby changes will not affect the global variable. 
  
  data1 <- data
  scale.data1 <- scale.data
  counts1 <- counts

  #For each for loop below, I'm generating a column resampling once. 
  #I then apply that to both counts and data, but not scaleData
  #For scale data (which is smaller) -> I find the gene that's being permuted
  #and apply the resampling for the same gene, if it's there in scale.data
  
  for(k in 1:dim(data)[1]){
    resample <-sample(length(data[k,]))
    data1[k,] <- data[k,resample]
    if(rownames(data)[k] %in% row.names(scale.data)){
      scale.data1[rownames(data)[k],] <- scale.data[rownames(data)[k],resample]
    }
  
  }

  warning("Now overwriting with scrambled data")
  object1@assays$RNA@counts = counts1
  object1@assays$RNA@data = data1
  object1@assays$RNA@scale.data = scale.data1

  warning("Now finding variable features for this scrambling")
  object1 <- FindVariableFeatures(object1, selection.method = "vst", nfeatures = 2000)
  object1 <- RunPCA(object1, features = VariableFeatures(object = object1))
  object1 <- FindNeighbors(object1, dims = 1:10)
  
  warning("Now finding clusters and their markers over a wide variety of resolutions.")
  foreach(m=c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4), .combine='rbind', .errorhandling ="remove", .packages=c('Seurat')) %do% {
    warning(paste("Finding Clusters & Markers for Resolution ", m, sep = ""))    
    object1 <- FindClusters(object1, resolution = m, n.iter = 30)
    test <- FindAllMarkers(object1, logfc.threshold = 0.25, only.pos = TRUE)
    CellNum <- table(object1@active.ident)
    Resolution <-rep(m,dim(test)[1])
    Iteration <-rep(j, dim(test)[1])
    NumberCells <- as.character(test$cluster)
    warning("Calculating number of cells")
    NumberCells <- unlist(lapply(1:dim(test)[1], function (x) NumberCells[x] <- CellNum[test$cluster[x]]))
    Combined <- paste(Iteration, test$cluster, Resolution)
    cbind(test, Resolution, Iteration, Combined, NumberCells)
   
  }
  
}

#Export
write.csv(ThreshMarkers, file = "PermuteCells_ClusterResolution_ThreshMarkers.csv")

quit()
