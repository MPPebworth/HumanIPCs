##Comparing the null distribution

##
setwd("OrganoidReportCard")
packrat::init("SeuratV3", restart = TRUE, infer.dependencies = FALSE)
packrat::status()
packrat::on()
.libPaths()
library(tidyverse)

NullDist <- read.csv("PermuteCells_ClusterResolution_ThreshMarkers.csv")

# #####Looking at Marker genes per cluster across different settin --------
#####Looking at Marker genes per cluster across different settings
table(NullDist$cluster)
#Look at the average number of markers in each cluster across all settings
hist(table(NullDist$Combined))
#What's the average number of marker genes per cluster at each resolution?
select(NullDist,cluster,Resolution,Iteration) %>% 
  group_by(Resolution, Iteration) %>%
  summarise(MarkPerCluster = n()) %>%
  summarize(meanMarkers = mean(MarkPerCluster),
            RangeMarkers = sd(MarkPerCluster))



#####Looking at the number of cells per cluster at different settings
#Look at the nmber of cells in each cluster. Looks like a poisson distribution of some kind.
hist(NullDist$NumberCells)

#Let's graph the number of cells in each cluster for each resolution value? 
plot(NullDist[c(9,12)])
#What's the average number of cells/cluster across iterations at each resolution value?
select(NullDist,Resolution,Iteration,cluster, NumberCells) %>%
  group_by(Resolution, Iteration) %>%
  distinct() %>%
  summarize(Itermean = mean(NumberCells, na.rm = TRUE)) %>%
  group_by(Resolution) %>%
  summarize(MeanCellsClusters = mean(Itermean, na.rm = TRUE),
            RangeCellsCluster = sd(Itermean))

select(NullDist,Resolution,Iteration,cluster,NumberCells) %>%
  group_by(Resolution, Iteration) %>%
  distinct() %>% 
  ggplot(aes(x= Resolution, y = NumberCells)) +
  geom_point() + ylab("Cells per cluster")

#Here's the best way to illustrate the effect of resolution of cluster size.
select(NullDist,Resolution,Iteration,cluster,NumberCells) %>%
  group_by(Resolution, Iteration) %>%
  distinct() %>% 
  summarize(MeanCells = mean(NumberCells, na.rm = TRUE), RangeCellsCluster = sd(NumberCells)) %>%
  ggplot(aes(x= Resolution, y = MeanCells)) +
  geom_point() + ylab("Average Cluster Size [# of cells]") + 
  geom_errorbar(aes(ymin = MeanCells - RangeCellsCluster, ymax = MeanCells + RangeCellsCluster)) +
  ggtitle("Effect of Resolution on Average Cluster Size")

#####Looking at the p-values (after correction) across different values
hist(NullDist$p_val_adj)
hist(-log(NullDist$p_val_adj))

sum(NullDist$p_val_adj < 0.05)/length(NullDist$p_val_adj)
#9.5% of markers across all resolutions were significant. 

#Let's look at significance by resolution
select(NullDist,Resolution,cluster,p_val_adj) %>%
  group_by(Resolution) %>%
  summarize(Sig = sum(p_val_adj < 0.05)/sum(p_val_adj >= 0),
            MinP = min(p_val_adj), MaxP = max(p_val_adj))
#lower resolution values produced a higher percentage of false positive p-values. 

#How does the adjuste p-value related to the number of cells detected in each cluster?
select(NullDist, Combined, p_val_adj, NumberCells) %>%
  group_by(Combined, NumberCells) %>%
  summarize(Sig = sum(p_val_adj < 0.05)/sum(p_val_adj >= 0),
            MinP = min(p_val_adj), MaxP = max(p_val_adj))
#Plot the following:
#-- NumberCells by the % of significant P-val-adj in each cluster
#-- NumberCells by the minimum p-value generated
#-- Number of cells y the maximum p-value generated

ggplot(select(NullDist, p_val_adj, NumberCells), aes(x = p_val_adj, y = NumberCells, )) + stat_density_2d()

#It may be more helpful to look at adjusted p-value as compared to number of cells per cluster
#My computer can't handle the numbers here. So I'm going to find the minimum adjusted p-value for each cluster
#I'm also going to find the mean

select(NullDist,Resolution, Iteration, cluster,NumberCells, p_val_adj) %>%
  group_by(Resolution, Iteration, NumberCells) %>%
  summarize(minP = min(p_val_adj), maxP = max(p_val_adj))


NullWide <- NullDist[c(9,10,12)] %>% spread(Resolution, NumberCells)



# Comparing to Original Dataset -------------------------------------------

load("MP_IPC_ForPresentation.RData")
load('MP_IPC3.RData')
MP_IPC3_n <- MP_IPC3_2
table(MP_IPC3@active.ident)
hist(table(MP_IPC3_n@meta.data$Old.Cluster))
Idents(MP_IPC3_n) <- MP_IPC3_n@meta.data$Old.Cluster
table(Idents(MP_IPC3_n))[1]
OrigMarkers <-FindAllMarkers(MP_IPC3_n, only.pos = TRUE)
write.csv(OrigMarkers, "MP_IPC3_Markers_Presentation.csv")
OrigMarkers<- read.csv("MP_IPC3_Markers_Presentation.csv")


#Plot Original compared to Null Distribution
OrigMarkersPlot <- OrigMarkers %>% mutate(NumberCells = table(MP_IPC3_2@active.ident)[cluster]) 
OrigMarkersPlot1 <- OrigMarkers1 %>%filter(p_val_adj < 0.05)
NullDistPlot <- select(NullDist, NumberCells, p_val_adj)%>%filter(p_val_adj < 0.05)
#Cluster size by marker significance
ggplot(NullDistPlot, aes(x = -log(p_val_adj), y = NumberCells), colour = 'black') + 
  geom_point() +
  geom_point(data=OrigMarkersPlot1, colour = 'red') +
  ylab("Cluster Size") + xlab("-Log of Corrected P-values")

#Number of Marks by Cluster Size
NullDistPlot1 <- NullDist%>% group_by(Combined, NumberCells)%>% filter(p_val_adj < 0.05) %>% 
  summarise(ClusterNum = n())
OrigMarkersPlot2 <- OrigMarkersPlot %>% group_by(cluster,NumberCells)%>%
  filter(p_val_adj < 0.05) %>% summarise(ClusterNum = n())
ggplot(NullDistPlot1, aes(x = NumberCells, y= ClusterNum ), colour = 'black') + 
  geom_point() +
  geom_point(data=OrigMarkersPlot2, colour = 'red', label = OrigMarkersPlot2$cluster) + 
  geom_text(data=OrigMarkersPlot2,label=OrigMarkersPlot2$cluster) +
  xlab("Cluster Size") + ylab("Number of Markers") + theme_classic()


#Identify which Resolution will generate the most clusters of a comparable size to the original cluster.
MatchingSize <- select(NullDist,Resolution, Iteration, cluster, Combined, NumberCells) %>% 
  group_by(Resolution, Iteration) %>% filter(abs(NumberCells-table(Idents(MP_IPC3_n))[1])<= 50) %>%
  distinct() %>% group_by(Resolution) %>%
  summarize(n())

#Filter out all comparably sized clusters at the above resolution.
ClusterMatch <- select(NullDist, Resolution, Iteration, cluster, Combined, gene, NumberCells, p_val_adj) %>%
  group_by(Iteration, Resolution)%>% 
  filter(Resolution == MatchingSize$Resolution[which.max(MatchingSize$`n()`)]) %>% 
  filter(abs(NumberCells-table(Idents(MP_IPC3_n))[1])<= 50)

ggplot(select(ClusterMatch, NumberCells, p_val_adj), aes(x = p_val_adj, y = NumberCells)) +geom_point()+
  ylab("Cluster Size") + xlab("BH-corrected P-value")

ggplot(select(ClusterMatch, NumberCells, p_val_adj), aes(x = -log(p_val_adj), y = NumberCells)) + geom_point() +
  ylab("Cluster Size") + xlab("- Log of BH-corrected P-value")

#Let's look at the overlap of my FindAllMarkers geneset with the random geneset
sum(OrigMarkers$gene %in% ClusterMatch$gene)/length(OrigMarkers$gene) #19% overlap
Overlap <- OrigMarkers[OrigMarkers$gene %in% ClusterMatch$gene,] %>% filter(p_val_adj < 0.05) %>% distinct()
dim(Overlap) #440 marker genes overlap

#Find the most significant p-value for each iteration
topNullHits <- ClusterMatch %>% group_by(Iteration) %>% arrange(desc(-p_val_adj)) %>% slice(1) %>% arrange(desc(-p_val_adj))
topNullHits
range(topNullHits$p_val_adj) #Range is e-31 to e-11 for null
hist(-log(topNullHits$p_val_adj))
#Follows a very tight curve, wih almost the vast majority of interations producing a top marker around 7.7e-14
ggplot(topNullHits, aes(x = -log(p_val_adj))) + geom_histogram() + xlab("-Log of BH-corrected P-value") + 
  ggtitle("Distribution of Most Significant P-Value in Each Iteration")

#Compare to original marker lists. 
OrigMarkers %>% group_by(cluster) %>% arrange(desc(-p_val_adj)) %>% slice(1) %>% arrange(desc(-p_val_adj)) #Many clusters are below this threshold.
range(OrigMarkers$p_val_adj) #Range is up to e-202

##Hurray! p-values seems to converge. Let's take the mean to find a threshold.
mean(topNullHits$p_val_adj) #7.7e-14

#Test how many clusters have at least 50% of their significant markers above this threshold. 
test <- OrigMarkers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% mutate(Above = (p_val_adj <=  mean(topNullHits$p_val_adj))) %>%
  summarize(Percent = sum(Above)/n())
test %>% filter(Percent > 0.5)
hist(test$Percent)
# 2 initial clusters passed the 50% threshold.

#How about a statistically more sound proposition?
#Let's go back to the original ClusterMatch object, and filter for signifiance BH-correct p-values.
#Let's assume that we're measuring noise. So then we are sampling noise, which should hae an approximate mean.
#I want to compare that mean p-value of noise to the p-values of my other samples. 
#I'll then find the 95th confidence interval via bootstrap o fthe mean of significant p-values.
#I'll do a bias accelerated, bias corrected bootstrap to generate 95% confidence intervals
dim(ClusterMatch) #414,439
Thresh <- ClusterMatch %>% ungroup() %>% filter(p_val_adj < 0.05) 
  mutate(SRank = percent_rank(p_val_adj)) %>% 
  arrange(desc(-SRank)) %>%
  mutate(Near = abs(SRank-0.05)) %>%
  arrange(desc(-Near)) %>% slice(1)
  
#Thresh is about 0.0000816. (8.16e-05) -- but just for this cluster
test <- OrigMarkers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% 
  mutate(Above = (p_val_adj <=  mean(Thresh$p_val_adj))) %>%
  summarize(Percent = sum(Above)/n())
test
hist(test$Percent)
test %>% filter(Percent > 0.9)



# Find p-value threshold for cluster based on matched cluster size ----------------

#Let's create a chart of cluster size vs. p-value threshold by values of + or - 10 cells (bin by steps sizes of 20)
range(NullDist$NumberCells) #22 to 2269. Let's start the iteration at a cluster size of 33.
range(table(IPCs@active.ident))
#Let's figure out what bin size to use. 
#First, let's interate through a given bin size, figure out if any clusters show up there. 
n=50
for(i in seq(n, 2258, by = 25)){

  clustersleft <- select(NullDist,Resolution, Iteration, cluster, Combined, NumberCells) %>% 
    group_by(Resolution, Iteration) %>% filter(abs(NumberCells-i)<= n) %>% group_by(Resolution) %>%
    summarize(Num=n()) %>% arrange(desc(Num)) %>% slice(1)
  print(paste(i, clustersleft$Num,sep =" "))
  flush.console()
  if(dim(clustersleft)[1] == 0){ break}
}
#Good up to 1475 cells. That's pretty good.
SizeMatchSummary <- NULL
ClusterMatchSummary <- NULL
ThreshMatchSummary <- NULL

for(i in seq(n, 1475, by = n)){
  
  #Identify which Resolution will generate the most clusters of a comparable size to the original cluster.
  MatchingSize <- select(NullDist,Resolution, Iteration, cluster, Combined, NumberCells) %>% 
    group_by(Resolution, Iteration) %>% filter(abs(NumberCells-i)<= n ) %>%
    distinct() %>% group_by(Resolution) %>%
    summarize(ClustNum = n()) %>% mutate(ClusterSize = i)
  SizeMatchSummary <- rbind(SizeMatchSummary, MatchingSize)
  
  #Filter out all comparably sized clusters at the above resolution.
  ClusterMatch <- select(NullDist, Resolution, Iteration, cluster, Combined, gene, NumberCells, p_val_adj) %>%
    group_by(Iteration, Resolution)%>% 
    filter(Resolution == MatchingSize$Resolution[which.max(MatchingSize$ClustNum)]) %>% 
    filter(abs(NumberCells-i)<= n )
  #See how many markers are showing up in the matching cluster size
  ClusterMatch1 <- select(ClusterMatch,Resolution, Combined) %>% ungroup() %>% 
    group_by(Combined) %>% summarize(MarkNum = n()) %>% 
    mutate(ClusterSize = i)
  ClusterMatchSummary <- rbind(ClusterMatchSummary, ClusterMatch1)
  
  #Find Threshold for that cluster size
  #Let's take the distribution of p-values (adjusted) and calculated the accelerated-bias, bias-corrected bootstrap
  #
  print(dim(ClusterMatch))
  Thresh <- ClusterMatch  %>% filter(p_val_adj <= 0.05) %>% ungroup() %>% select(p_val_adj) %>% 
    mutate(logP = -log(p_val_adj)) 
  Thresh1 <- cbind(exp(-quantile(Thresh$logP, probs = 0.95)), i)
  colnames(Thresh1) <- c('Thresh','ClusterSize')
  ThreshMatchSummary <- rbind(ThreshMatchSummary,Thresh1)
}

ThreshMatchSummary <- as.data.frame(ThreshMatchSummary)
ClusterMatchSummary <- as.data.frame(ClusterMatchSummary)
SizehMatchSummary <- as.data.frame(SizeMatchSummary)

save(ThreshMatchSummary, SizeMatchSummary, ClusterMatchSummary, file = "NullDistribution_ThresholdvClusterSize.RData")
ggplot(as.data.frame(ThreshMatchSummary), aes(y = Thresh, x = ClusterSize)) + geom_point() + scale_y_continuous(trans='log2')

#beyond cluster size of 600, the algorthm  falls apart. Probably because it can't. 
ggplot(as.data.frame(ClusterMatchSummary), aes(y = MarkNum, x = ClusterSize)) + geom_point()
ggplot(as.data.frame(ClusterMatchSummary), aes(y = MarkNum, x = ClusterSize)) + geom_point()

ggplot(as.data.frame(SizeMatchSummary), aes(y = Resolution, x = ClusterSize)) + geom_point()
ggplot(as.data.frame(SizeMatchSummary), aes(y = Resolution, x = ClustNum)) + geom_point()
#Let's plot the

#Let's also check for every significant marker, what the highest detected p-value is in our scrambled dataset.
# 1) Find all genes from that cluster in that dataset
# 3) What's the 95% percentile of the p-values detected?
# 3) Add those values to an additional column that was made, but only by match that column. 
#Double check that identities and OrigMarker cluster names match
table(MP_IPC3_2@active.ident)
table(OrigMarkers$cluster)
OrigMarkers$cluster = as.character(OrigMarkers$cluster)
OrigMarkers$cluster[OrigMarkers$cluster=="Generic IPCs"] = "Archetypic IPCs"
OrigMarkers$cluster[OrigMarkers$cluster=="Neuronal1"] = "N1"
OrigMarkers$cluster[OrigMarkers$cluster=="Neuronal2"] = "N2"
OrigMarkers$cluster[OrigMarkers$cluster=="Neuronal3"] = "N3"
sum(names(table(MP_IPC3_2@active.ident)) %in% names(table(OrigMarkers$cluster)))
length(names(table(MP_IPC3_2@active.ident)))

Gene Threshold <- c(1:length(OrigMarkers$gene))
OrigMarkers2 <- cbind(OrigMarkers, Gene Threshold)
SizeMatch1 <- SizeMatchSummary %>% group_by(ClusterSize) %>% filter(max(ClustNum) == ClustNum)

for(i in 1:dim(table(MP_IPC3_2@active.ident))){
  #Identify cluster
  clusterN <- names(table(MP_IPC3_2@active.ident))[i]
  #Pull out cluster markres for that cluster
  NOrigMarkers<-OrigMarkers %>% filter(clusterN == cluster)
  clusterSize <- SizeMatch1 %>% ungroup() %>% slice(which.min(abs(table(MP_IPC3_2@active.ident)[i]- ClustNum)))
  
  clusterN1 <- NullDist %>% filter(Resolution == clusterSize$Resolution) %>% 
    group_by(gene) %>% slice(which.min(p_val_adj))
  
  clusterN2 <-  clusterN1 %>% filter(gene %in% NOrigMarkers$gene) %>% select(gene, p_val_adj)
  clusterN3 <- OrigMarkers %>% filter(cluster == clusterN) %>% full_join(clusterN2, by = 'gene',keep =TRUE)
  #Test to make sure order isn't messed up
  sum(as.character(NOrigMarkers$gene) == clusterN3$gene)
  dim(NOrigMarkers)
  sum( OrigMarkers2[OrigMarkers$cluster == clusterN,]$gene == clusterN3$gene)
  
  #Now set the size matched p-values into the OrigMarker file
  OrigMarkers2[OrigMarkers$cluster == clusterN,] = clusterN3
}

OrigMarkers2 %>% filter(p_val_adj < Gene Threshold) %>% dim()
dim(OrigMarkers) #462
sum(is.na(OrigMarkers2$Gene Threshold)) #1914
#462+1914
#2376/2436 markers are above the null distribution p-values
#Which genes aren't?
OrigMarkers2 %>% filter(p_val_adj >= Gene Threshold) %>% filter(!is.na(maxPvalue))
View(OrigMarkers2 %>% filter(p_val_adj < maxPvalue))
View(OrigMarkers2)

write.csv(OrigMarkers2, file = "ClusterSignificance_MarkerList.csv") #Written over after generating threshold for each

# Implement Recursive Algorithm to Merge Clusters -------------------------

MarkerGeneSpace <- select(OrigMarkers, p_val_adj, gene) %>% filter(p_val_adj < 0.05) %>% select(gene) %>% distinct()
MP_IPC3_n <- FindVariableFeatures(MP_IPC3_n) 
top10 <- head(VariableFeatures(MP_IPC3_n), 10)
plot1 <- VariableFeaturePlot(MP_IPC3_n)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
#4000 because there are about 3000 marker genes total, 
#I didn't use marker genes as the space for generating correlations because some of the marker genes might be explained
#by noise. Rather, I took the variable genes, because those are unbiased.

#Create all the variables needed for the While loop, without changing global variables

#testIPCs <- MP_IPC3_n
#Idents(testIPCs) <- IPCs@meta.data$New.Cluster1
#IPCs <- testIPCs
Idents(MP_IPC)
IPCs <- MP_IPC3_2
#Idents(IPCs) <- MP_IPC3_2@meta.data$Old.Cluster
scaleIPCs <- CreateSeuratObject(IPCs@assays$RNA@scale.data, assay = "RNA", meta.data = IPCs@meta.data)
scaleIPCs <- FindVariableFeatures(scaleIPCs, nfeatures = 3000)
Idents(scaleIPCs) <- IPCs@active.ident
newMarkers <- FindAllMarkers(scaleIPCs,features = scaleIPCs@assays$RNA@var.features)
newMarkers1 <- filter(newMarkers, p_val_adj < 0.05)

AvgExpr <- AverageExpression(IPCs,use.scale = TRUE, features = unique(newMarkers1$gene))
AvgExpr <- as.matrix(AvgExpr$RNA)
AvgExprCorr <- cor(AvgExpr, method = "pearson")
AvgExprCorr <- as.matrix(AvgExprCorr)
plot(AvgExprCorr)
heatmap(AvgExprCorr)
#Markers <- OrigMarkers
#Markers <- OrigMarkers2
#MarkersSaved <- FindAllMarkers(IPCs,only.pos = TRUE)
#Markers <- MarkersSaved
#save(Markers, file = "Markers_forPresentation.RData")
#load("Markers_forPresentation.RData")

#Test how many clusters are real
#We want more than 90% of markers to be real
#50% generated pretty small clusters and only some combined.
PercentP = 0.75

CustomThresh = table(IPCs@active.ident)
for(i in 1:length(table(IPCs@active.ident))){
  ClusterToThresh<- ThreshMatchSummary %>% 
    mutate(CloseSize = abs(ClusterSize-table(IPCs@active.ident)[i])) %>%
    slice(which.min(CloseSize))
  CustomThresh[i]=ClusterToThresh$Thresh
  
}

AvgpVal<- Markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% 
  mutate(Thresh = CustomThresh[cluster]) %>%
  mutate(Above = (p_val_adj <=  Thresh)) %>%
  summarize(Percent = sum(Above)/n()) %>% mutate(Pass = Percent > PercentP)
AvgpVal[order(as.character(AvgpVal$cluster)),]
table(newIPCs@active.ident)[order(names(table(IPCs@active.ident)))]
AvgpVal <- cbind(AvgpVal[order(as.character(AvgpVal$cluster)),], 
                 table(IPCs@active.ident)[order(names(table(IPCs@active.ident)))])
AvgpVal

if (any(AvgpVal$Freq > 400)){AvgpVal[AvgpVal$Freq > 400,3] = TRUE}

ggplot(AvgpVal, aes(x = Percent)) + geom_histogram(binwidth = 0.1) +
  ggtitle("Percentage of Markers Above Threshold by Cluster") + ylab('Number of Clusters') +
  xlab('Percent of Gene Markers that Pass Null Threshold')

#What if it was the number of marker genes that passed threshold that mattered? How many did?
PassedMarkers <- Markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% 
  mutate(Thresh = CustomThresh[cluster]) %>%
  filter(p_val_adj <=  Thresh) 
PassedMarkers
PassedMarkers %>%  group_by(cluster) %>% summarize(TotalMarkers = n())

AvgExpr1 <- as.matrix(AvgExpr) 
AvgExprCorr1 <- AvgExprCorr
heatmap(AvgExprCorr1)

diag(AvgExprCorr1) = 0
rownames(AvgExprCorr1)=names(table(IPCs@active.ident))
colnames(AvgExprCorr1)=names(table(IPCs@active.ident))

heatmap(AvgExprCorr1)

CorrFrame <- data.frame(row=rownames(AvgExprCorr1)[row(AvgExprCorr1)[upper.tri(AvgExprCorr1)]], 
           col=colnames(AvgExprCorr1)[col(AvgExprCorr1)[upper.tri(AvgExprCorr1)]], 
           corr=AvgExprCorr1[upper.tri(AvgExprCorr1)]) %>% arrange(desc(corr))
head(CorrFrame)



 #Here's the cut off valve. Stop combining when clusters get too big for my tests. (500 cells)

while(any(!AvgpVal$Pass)){
  
  #Find names of top two correlated clusters. If there's a tie, take the first. 
  bestCorr = rownames(AvgExprCorr1)[arrayInd(which.max(AvgExprCorr1), dim(AvgExprCorr1))]
  #To implement: test with both of those with high correlation are past the threshold found imperically. 
  #if so, then stop!
  BestClusters <- select(AvgpVal,cluster, Pass) %>% filter(cluster %in% bestCorr) %>% slice(c(1,2))
  
  #If both top clusters pass the threshold of a real cluster, then don't merge them
  #Instead, set their correlation to zero, and look for the next best cluster. 
  #This algorithm also chooses the first max correlation, so ties aren't an issue. 
  #Put limit on cluster size because I don't have values. 
  while(all(BestClusters$Pass) | any(table(IPCs@active.ident)[BestClusters$cluster] > 600)){
    AvgExprCorr1[arrayInd(which.max(AvgExprCorr1), dim(AvgExprCorr1))] = 0
    bestCorr = rownames(AvgExprCorr1)[arrayInd(which.max(AvgExprCorr1), dim(AvgExprCorr1))]
    BestClusters <- select(AvgpVal,cluster, Pass) %>% filter(cluster %in% bestCorr) %>% slice(c(1,2))
    
  }
  #Double check that at least one of the highest correlated clusters passes the significance threshold. 
  if(any(!BestClusters$Pass)){
    
    newIdent  <- paste(bestCorr[1], "_",bestCorr[2], sep ="")
    oldIdent1 <- bestCorr[1]
    oldIdent2 <- bestCorr[2]
    
    cells.use <- WhichCells(object = IPCs, idents = c(oldIdent1, oldIdent2))
    IPCs      <- SetIdent(IPCs,  cells =  cells.use, value = newIdent)

    newMarker <- FindMarkers(IPCs, ident.1 = newIdent, only.pos = TRUE)
    newMarker <- cbind(newMarker, rep(newIdent, dim(newMarker)[1]), row.names(newMarker))
    colnames(newMarker)[c(6,7)]= c('cluster','gene')
    Markers <- droplevels(rbind(Markers[!Markers$cluster %in% bestCorr,], newMarker))
    
    CustomThresh =  table(IPCs@active.ident)
    for(i in 1:length(table(IPCs@active.ident))){
      ClusterToThresh <- 
        ThreshMatchSummary %>% 
        mutate(CloseSize = abs(ClusterSize-table(IPCs@active.ident)[i])) %>%
        slice(which.min(CloseSize))
      CustomThresh[i]=ClusterToThresh$Thresh 
      
    }
    
    AvgpVal<- Markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% 
      mutate(Thresh = CustomThresh[cluster]) %>%
      mutate(Above = (p_val_adj <=  Thresh)) %>%
      summarize(Percent = sum(Above)/n()) %>% mutate(Pass = Percent > PercentP)
    AvgpVal <- cbind(AvgpVal[order(as.character(AvgpVal$cluster)),], 
                     table(IPCs@active.ident)[order(names(table(IPCs@active.ident)))])
    AvgpVal[AvgpVal$Freq > 500,3] = TRUE
    
    newCluster<- subset(IPCs, idents= newIdent)
    newAvg  <- AverageExpression(newCluster, features = unique(newMarkers1$gene), use.scale = TRUE)
    newAvg <- newAvg$RNA
    AvgExpr1 <- cbind(AvgExpr1[,!colnames(AvgExpr1) %in% c(oldIdent1, oldIdent2)], newAvg)
  
  }
  AvgExprCorr1 <- cor(AvgExpr1, method = c("pearson"))
  diag(AvgExprCorr1) = 0
  rownames(AvgExprCorr1)=names(table(IPCs@active.ident))
  colnames(AvgExprCorr1)=names(table(IPCs@active.ident))
  
}

TSNEPlot(IPCs)


# Guided reclustering of New Cluster1 ---------------------------------------

newIPCs <- MP_IPC3_2
table(newIPCs@meta.data$New.Cluster1)
Idents(newIPCs) <- newIPCs@meta.data$New.Cluster1

newMarkersSaved <- FindAllMarkers(newIPCs,only.pos = TRUE)
save(newMarkersSaved, file = "Markers_newCluster1.RData")

#Generate custom threshold
PercentP = 0.50
CustomThresh = table(newIPCs@active.ident)
for(i in 1:length(table(newIPCs@active.ident))){
  ClusterToThresh<- ThreshMatchSummary %>% 
    mutate(CloseSize = abs(ClusterSize-table(IPCs@active.ident)[i])) %>%
    slice(which.min(CloseSize))
  CustomThresh[i]=ClusterToThresh$Thresh
  
}
CustomThresh


#See if more than 65% of significant markers pass the custom threshold. 
AvgpVal<- newMarkersSaved %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% 
  mutate(Thresh = CustomThresh[cluster]) %>%
  mutate(Above = (p_val_adj <=  Thresh)) %>%
  summarize(Percent = sum(Above)/n()) %>% mutate(Pass = Percent > PercentP)
AvgpVal[order(as.character(AvgpVal$cluster)),]
table(IPCs@active.ident)[order(names(table(newIPCs@active.ident)))]
AvgpVal <- cbind(AvgpVal[order(as.character(AvgpVal$cluster)),], 
                 table(newIPCs@active.ident)[order(names(table(newIPCs@active.ident)))])
AvgpVal
if (any(AvgpVal$Freq > 500)){AvgpVal[AvgpVal$Freq > 500,3] = TRUE}

ggplot(AvgpVal, aes(x = Percent)) + geom_histogram(binwidth = 0.1) +
  ggtitle("Percentage of Markers Above Threshold by Cluster") + ylab('Number of Clusters') +
  xlab('Percent of Gene Markers that Pass Null Threshold')

AvgpVal %>% filter(!Pass)
#Clusters that don't pass
#astrocyte, div1, EGFR, lowquality1, EN-neuronal-DCX, IN-Neuronal-DCX, progenitor, neuronal2, novel
#None of the progenitor markers pass the threshold at all. 

AvgExpr1 <- AverageExpression(newIPCs,use.scale = TRUE)
AvgExpr1 <- as.matrix(AvgExpr1$RNA)
AvgExprCorr1 <- cor(AvgExpr1, method = "pearson")
plot(AvgExprCorr)
heatmap(AvgExprCorr)
diag(AvgExprCorr1)=0
heatmap(AvgExprCorr1)

CorrFrame <- data.frame(row=rownames(AvgExprCorr1)[row(AvgExprCorr1)[upper.tri(AvgExprCorr1)]], 
                        col=colnames(AvgExprCorr1)[col(AvgExprCorr1)[upper.tri(AvgExprCorr1)]], 
                        corr=AvgExprCorr1[upper.tri(AvgExprCorr1)]) %>% arrange(desc(corr))
failedCluster <-AvgpVal %>% filter(!Pass) %>% select(cluster) %>% unlist() %>% as.character()

CorrFrame1 <- CorrFrame[(CorrFrame$row %in% failedCluster | CorrFrame$col %in% failedCluster) & CorrFrame$corr > 0,]
CorrFrame1[c(1:15),]
#Recombine:
#loquality1 and outlier
#EGFR, astrocyte, Rgtoastrocyte, and novel
#IN-neuronal-DCX, 


#Let's combine these and try again
newIPCs1 <- RenameIdents(newIPCs, 'EGFR' = 'RG-Like', 'novel' = 'RG-Like', 'astrocyte' = 'RG-Like')
newIPCs1 <- RenameIdents(newIPCs1, "outlier" = 'Ribosomal','lowquality1' = 'Ribosomal')
table(newIPCs1@active.ident)

PercentP = 0.65
CustomThresh = table(newIPCs1@active.ident)
for(i in 1:length(table(newIPCs1@active.ident))){
  ClusterToThresh<- ThreshMatchSummary %>% 
    mutate(CloseSize = abs(ClusterSize-table(IPCs@active.ident)[i])) %>%
    slice(which.min(CloseSize))
  CustomThresh[i]=ClusterToThresh$Thresh
  
}

#See if more than 65% of significant markers pass the custom threshold. 
AvgpVal<- newMarkersSaved %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% 
  mutate(Thresh = CustomThresh[cluster]) %>%
  mutate(Above = (p_val_adj <=  Thresh)) %>%
  summarize(Percent = sum(Above)/n()) %>% mutate(Pass = Percent > PercentP)

AvgpVal <- cbind(AvgpVal[order(as.character(AvgpVal$cluster)),], 
                 table(newIPCs@active.ident)[order(names(table(newIPCs@active.ident)))])
AvgpVal

AvgExpr1 <- AverageExpression(newIPCs,use.scale = TRUE)
AvgExpr1 <- as.matrix(AvgExpr1$RNA)
AvgExprCorr1 <- cor(AvgExpr1, method = "pearson")
plot(AvgExprCorr)
heatmap(AvgExprCorr)
diag(AvgExprCorr1)=0
heatmap(AvgExprCorr1)

CorrFrame <- data.frame(row=rownames(AvgExprCorr1)[row(AvgExprCorr1)[upper.tri(AvgExprCorr1)]], 
                        col=colnames(AvgExprCorr1)[col(AvgExprCorr1)[upper.tri(AvgExprCorr1)]], 
                        corr=AvgExprCorr1[upper.tri(AvgExprCorr1)]) %>% arrange(desc(corr))

AvgpVal


# Filter Markers by Custom Threshold --------------------------------------

newIPCs1<- MP_IPC3
n=50
CustomThresh = table(newIPCs1@active.ident)

for(i in 1:length(CustomThresh)){
  
  print(names(CustomThresh)[i])
  
  #Identify which Resolution will generate the most clusters of a comparable size to the original cluster.
  MatchingSize <- select(NullDist,Resolution, Iteration, cluster, Combined, NumberCells) %>% 
    group_by(Resolution, Iteration) %>% filter(abs(NumberCells-table(newIPCs1@active.ident)[i])<= n ) %>%
    distinct() %>% group_by(Resolution) %>%
    summarize(ClustNum = n()) %>% arrange(desc(ClustNum)) %>% slice(1)
  
  #Filter for all comparably sized clusters at the above resolution.
  #Flter by p-value of 0.05 to look for false positives. We're testing to see whether we can get 
  #our markers to be able the 0.05 p-value
  ClusterMatch <- select(NullDist, Resolution, Iteration, cluster, Combined, gene, NumberCells, p_val,p_val_adj) %>%
    group_by(Iteration, Resolution)%>%
    filter(Resolution == MatchingSize$Resolution) %>% 
    filter(abs(NumberCells-table(newIPCs1@active.ident)[i])<= n ) %>% filter(p_val_adj < 1)
  
  #Caculate custom thresh
    CustomThresh[i]=exp(-quantile(-log(ClusterMatch$p_val_adj), probs = 0.95))
  
}
CustomThresh_adj = p.adjust(CustomThresh, method = "bonferroni")

#save(CustomThresh_adj, file = 'ClusterThreshold_NewCluster1.Rdata')

AvgpVal <- Markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% 
  mutate(Thresh = CustomThresh_adj[cluster]) %>%
  filter(p_val_adj <=  Thresh)

AvgpVal
View(AvgpVal)

OrigMarkers2 <- OrigMarkers2 %>% mutate(Threshold = CustomThresh_adj[cluster])


write.csv(OrigMarkers2 %>% select(-MinPvalue), file = "ClusterSignificance_MarkerList.csv") 

# Test significance of cluster by distribution of unfiltered p-val ------

IPCs <- MP_IPC3

Markers <- OrigMarkers2
#Markers <- FindAllMarkers(IPCs, only.pos = TRUE)

OrigClusterSize <- table(IPCs@active.ident)

n=50

ClusterSignificance <- NULL

for(i in 1:length(OrigClusterSize)){
  
  print(names(OrigClusterSize)[i])
  
  #Identify which Resolution will generate the most clusters of a comparable size to the original cluster.
  MatchingSize <- select(NullDist,Resolution, Iteration, cluster, Combined, NumberCells) %>% 
    group_by(Resolution, Iteration) %>% filter(abs(NumberCells-OrigClusterSize[i])<= n ) %>%
    distinct() %>% group_by(Resolution) %>%
    summarize(ClustNum = n()) %>% arrange(desc(ClustNum)) %>% slice(1)
  print(MatchingSize$Resolution)
  
  #Filter out all comparably sized clusters at the above resolution.
  ClusterMatch <- select(NullDist, Resolution, Iteration, cluster, Combined, gene, NumberCells, p_val_adj) %>%
    group_by(Iteration, Resolution)%>% 
    filter(Resolution == MatchingSize$Resolution) %>% 
    filter(abs(NumberCells-OrigClusterSize[i])<= n)
  print(dim(ClusterMatch))
 
  print(dim(ClusterMatch$cluster)[1])
  #Run a Mann Whitney U test on individual clusters to test if they are significant. 
  
  OrigClusterMarker1 <- Markers %>% filter(cluster == names(OrigClusterSize)[i])
  
  ClusterSignificance[[i]] <-wilcox.test(ClusterMatch$p_val_adj,OrigClusterMarker1$p_val_adj, 
              paired = FALSE, correct = FALSE)$p.value

}

names(ClusterSignificance) = names(OrigClusterSize)
write.csv(p.adjust(ClusterSignificance, method = "bonferroni"), file = "ClusterSignificance_Test.csv")


# Test significance of cluster by distribution of unfiltered p-val and recluster! ------

IPcs <- MP_IPC3_2
Idents(IPCs) <- IPCs@meta.data$Old.Cluster
Markers <- FindAllMarkers(IPCs, only.pos = TRUE)
write.csv(Markers, file ="OriginalClusterMarkers.csv")

Corr <- as.matrix(read.table('MP_IPC_corr.txt'))
heatmap(Corr)

OrigClusterSize <- table(Markers$cluster)
n=100
ClusterSignificance <- NULL

for(i in 1:length(OrigClusterSize)){
  
  print(names(OrigClusterSize)[i])
  
  #Identify which Resolution will generate the most clusters of a comparable size to the original cluster.
  MatchingSize <- select(NullDist,Resolution, Iteration, cluster, Combined, NumberCells) %>% 
    group_by(Resolution, Iteration) %>% filter(abs(NumberCells-OrigClusterSize[i])<= n ) %>%
    distinct() %>% group_by(Resolution) %>%
    summarize(ClustNum = n()) %>% arrange(desc(ClustNum)) %>% slice(1)
  #print(MatchingSize$Resolution)
  
  #Filter out all comparably sized clusters at the above resolution.
  ClusterMatch <- select(NullDist, Resolution, Iteration, cluster, Combined, gene, NumberCells, p_val_adj) %>%
    group_by(Iteration, Resolution)%>% 
    filter(Resolution == MatchingSize$Resolution) %>% 
    filter(abs(NumberCells-OrigClusterSize[i])<= n )
  #print(dim(ClusterMatch))
  
  #print(dim(ClusterMatch$cluster)[1])
  #Run a Mann Whitney U test on individual clusters to test if they are significant. 
  
  OrigClusterMarker1 <- Markers %>% filter(cluster == names(OrigClusterSize)[i])
  
  ClusterSignificance[[i]] <-wilcox.test(ClusterMatch$p_val_adj,OrigClusterMarker1$p_val_adj, 
   paired = FALSE, correct = FALSE)$p.value
  
}

names(ClusterSignificance) = names(OrigClusterSize)
SigClustAdj <- p.adjust(ClusterSignificance, method = "bonferroni")

inSigClust <- names(SigClustAdj[SigClustAdj >= 0.05])
view(inSigClust)

Corr1 <- Corr
diag(Corr1) = 0


IPCs1 <- IPCs
for(i in 1:length(inSigClust)){
  ToMix[i] <- rownames(Corr1)[which.max(Corr1[rownames(Corr1) == inSigClust[i],])]
  newId <- paste(ToMix, inSigClust[i], sep = "_")
  print(newId)
}

#cells.use <- WhichCells(object = IPCs, idents = c(ToMix,inSigClust[i]))
#IPCs1      <- SetIdent(IPCs1,  cells =  cells.use, value = newIdent)



newCluster<- subset(IPCs, idents= newIdent)
newAvg  <- AverageExpression(newCluster, features = unique(newMarkers1$gene), use.scale = TRUE)
newAvg <- newAvg$RNA
AvgExpr1 <- cbind(AvgExpr1[,!colnames(AvgExpr1) %in% c(oldIdent1, oldIdent2)], newAvg)

AvgExprCorr1 <- cor(AvgExpr1, method = c("pearson"))
diag(AvgExprCorr1) = 0
rownames(AvgExprCorr1)=names(table(IPCs@active.ident))
colnames(AvgExprCorr1)=names(table(IPCs@active.ident))

