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

table(MP_IPC3_2@active.ident)

#Set up Parallel processes
cl <- makeCluster(10)
registerDoParallel(cl)
#export data matrix
data1 <- as.matrix(MP_IPC3_2@assays$RNA@data)

cellnumber = seq(from = 10, to = 510, by = 25)
object1 <- MP_IPC3_2
object1[["old.ident"]] <- Idents(object = object1)
NumCells <-dim(data1)[2]

start_time <- Sys.time()
ThreshMarkers<- foreach(m = 1:length(cellnumber),.combine='rbind', .packages=c('foreach','Seurat','doParallel')) %do%{
  
  foreach(j=1:100, .combine='rbind', .packages=c('foreach','Seurat','doParallel')) %dopar%{
    
    #Permute data
    shuf <- data1
    shuf<- t(sapply(1:dim(shuf)[1], function (row) shuf[row,] <<- sample(shuf[row,],replace = FALSE)))
    rownames(shuf) = rownames(data1)
    colnames(shuf) = colnames(data1)
    object1@assays$RNA@data = shuf
    
    foreach(i=1:100, .combine = 'rbind', .packages=c('Seurat')) %do% {
      idents1 <- as.character(Idents(object1))
      idents1[sample(1:NumCells, size= cellnumber[m], replace = FALSE)] <- rep("Identity1",cellnumber[m])
      names(idents1) = names(Idents(object1))
      Idents(object1) <- factor(idents1)
      test <- FindMarkers(object1, ident.1 = "Identity1", logfc.threshold = 0, 
                          min.pct = 0, only.pos = TRUE)
      Iteration <- rep(i,dim(test)[1])
      Permute <- rep(j,dim(test)[1])
      CellNum <- rep(cellnumber[m],dim(test)[1])
      test <- cbind(test,Iteration,Permute,CellNum)
      Idents(object = object1) <- object1[["old.ident"]]
      test
    }
    
  }
}
Sys.time()- start_time
stopCluster(cl)

write.csv(ThreshMarkers, file ="Permute_ThreshMarkers.csv")
