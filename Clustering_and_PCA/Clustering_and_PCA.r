setwd(".")
set.seed(1122)
options("digits"=3)
library(factoextra)
library(ggplot2)
library(dplyr)
library(cluster)
text2csv <- function(in.file,out.file) {
  mammals.df <- read.csv(file= in.file,
                         header = FALSE,
                         sep = "",
                         dec = ".",
                         comment.char = "#"
                         
  )
  
  write.csv(mammals.df, file = out.file)
  mammals.df <- read.csv(file= out.file,
                         header = TRUE,
                         sep = ",",
                         dec = ".",
                         comment.char = "#",
                         skip = 4
  )
  write.csv(mammals.df[,2:10], file = out.file,row.names = F)
  mammals.df <- read.csv(out.file,row.names = 1)
  
  return( mammals.df)
}
mammals.df <- text2csv(in.file = "file19.txt", out.file = "mammals.csv")
summary(mammals.df)
mammals.scaled <- scale(mammals.df)
rm(mammals.df)
fviz_nbclust(mammals.scaled, kmeans, method="wss") # Elbow method minimizes total
fviz_nbclust(mammals.scaled, kmeans, method="silhouette") # Silhouette method
set.seed(1122)
kmu <- kmeans(mammals.scaled, centers=7,nstart=40)
#png(filename="mammalCluster.png", units = "px", width=2000, height=2000)
fviz_cluster(kmu, data=mammals.scaled)
#dev.off()
obsxclus <-as.data.frame(table(kmu$cluster))
names(obsxclus) <-c("Cluster","Observations")
print(obsxclus)
rm(obsxclus)
print(paste("Total SSE :", kmu$tot.withinss))
SSExclus <-c()
for (x in 1:max(kmu$cluster)){
  SSExclus <-rbind(SSExclus,c(x,kmu$withinss[x]))
}
colnames(SSExclus) <- c("Cluster", "SSE")
print(SSExclus)
rm(SSExclus,x)
for (x in 1:max(kmu$cluster)){
  cat(paste("Group ",x, ":\n",paste(names(which(kmu$cluster==x)),collapse = ", ")) ,"\n") 
}
rm(x)
set.seed(1122)
mammals.sample <-sample_n(read.csv("mammals.csv",row.names = 1),size = 35) 
hclust.single <- eclust(mammals.sample, FUNcluster = "hclust",hc_method="single")
hclust.complete <- eclust(mammals.sample, FUNcluster = "hclust",hc_method="complete")
hclust.average <- eclust(mammals.sample, FUNcluster = "hclust",hc_method="average")
fviz_dend(hclust.average,show_labels=T ,palette="jco",cex = .6, hang=-.5,horiz = T,main = "Cluster Dendogram Using Average Linkage")
fviz_dend(hclust.complete,show_labels=T ,palette="jco",cex = .6, hang=-.5,horiz = T, main = "Cluster Dendogram Using Complete Linkage")
fviz_dend(hclust.single,show_labels=T ,palette="jco",cex = .6, hang=-.5,horiz = T,main = "Cluster Dendogram Using Single Linkage")
print2sing <- function(hc) {
  for (x in 1:34){
    if(hc$merge[x,1] <0 & hc$merge[x,2] <0){
      cat(paste("\t{",rownames(mammals.sample[ - hc$merge[x,1],] ),
                ", ",rownames(mammals.sample[ - hc$merge[x,2],] ),"}\n"))
    }
  }
}
cat("Two singleton clusters formed using Single Linkage:\n")
print2sing(hclust.single)

cat("\nTwo singleton clusters formed using Complete Linkage\n")
print2sing(hclust.complete)

cat("\nTwo singleton clusters formed using Average Linkage\n")
print2sing(hclust.average)
plot(hclust.single,hang = -1,main = "Cutting Single Linkage Dendogram at height 2")
groups <- cutree(hclust.single, h=2) 
rect.hclust(hclust.single, h=2, border="red") 
hclust.single <- eclust(mammals.sample, FUNcluster = "hclust",hc_method="single",k=5)
hclust.complete <- eclust(mammals.sample, FUNcluster = "hclust",hc_method="complete",k=5)
hclust.average <- eclust(mammals.sample, FUNcluster = "hclust",hc_method="average",k=5)
fviz_dend(hclust.average,show_labels=T ,palette="jco",cex = .6, hang=-.5,horiz = T,main = "Cluster Dendogram Using Average Linkage")
fviz_dend(hclust.complete,show_labels=T ,palette="jco",cex = .6, hang=-.5,horiz = T, main = "Cluster Dendogram Using Complete Linkage")
fviz_dend(hclust.single,show_labels=T ,palette="jco",cex = .6, hang=-.5,horiz = T,main = "Cluster Dendogram Using Single Linkage")
pDunnSil <-function(clust){
  distance <-dist(scale(mammals.sample))
  tt<-fpc::cluster.stats(distance,clust$cluster)
  cat(paste("Dunn :",tt$dunn," ; Silhouette :",tt$avg.silwidth ,"\n"))
}

print("Printing Dunn and Silhouette for Single linkage")
pDunnSil(hclust.single)
print("Printing Dunn and Silhouette for Complete linkage")
pDunnSil(hclust.complete)
print("Printing Dunn and Silhouette for Average linkage")
pDunnSil(hclust.average)
rm(hclust.single,hclust.complete,hclust.average)
set.seed(1122)
htru.df <-read.csv(file = "HTRU_2-small.csv")
htru.pca <- prcomp(scale(htru.df[,1:8]))
htru.pca$sdev
summary(htru.df)
cat(paste("The first two components explain:\n",
          format(100*sum(htru.pca$sdev[1:2])/sum(htru.pca$sdev),digits=4),
          "% of the variance"
)
)
first2 <-htru.pca$x[,1:2]
plot(first2[,1],
     first2[,2],
     col=rainbow(n=1,v = c( htru.df$class==1,htru.df$class==0)),
     xlab = paste("PC1 (",
                  format(100*sum(htru.pca$sdev[1])/sum(htru.pca$sdev),digits=4),
                  "%)"),
     ylab = paste("PC2 (",
                  format(100*sum(htru.pca$sdev[2])/sum(htru.pca$sdev),digits=4),
                  "%)"),
     main="Plotting using PCA with first two components"
)
khtwu<- kmeans(scale(htru.df[,1:8]), centers=2,nstart=25)
fviz_cluster(khtwu, 
             data=scale(htru.df),
             main="Plotting using Kmeans Using all components")
cat(paste0("In the Left  cluster (2), the percentage of observations is: ",
           format(100*sum(khtwu$cluster==2)/sum(khtwu$cluster==khtwu$cluster),digits=5)," %\n",
           "In the Right cluster (1), the percentage of observations is: ",
           format(100*sum(khtwu$cluster==1)/sum(khtwu$cluster==khtwu$cluster),digits=5)," %\n"))
cat(paste0("The 'TRUE'  class label percentage of observations is:  ",
           format(100*sum(htru.df$class==1)/sum(htru.df$class==htru.df$class),digits=5),
           " %\n",
           "The 'FALSE' class label percentage of observations is: ",
           format(100*sum(htru.df$class==0)/sum(htru.df$class==htru.df$class),digits=5),
           " %\n"
)
)
cat(paste0("From the larger cluster, ",
           format(100*sum(htru.df[which(khtwu$cluster==1),9]==0)/sum(khtwu$cluster==1),
                  digits=4),
           " % belong to the majority class (0), while only ",
           format(100*sum(htru.df[which(khtwu$cluster==1),9]==1)/sum(khtwu$cluster==1),
                  digits=4),
           " % belong to the minority class (1)"
)
)
tn <- 100*sum(htru.df[which(khtwu$cluster==1),9]==0)
fn <- 100*sum(htru.df[which(khtwu$cluster==1),9]==1)
tp <- 100*sum(htru.df[which(khtwu$cluster==2),9]==1)
fp <- 100*sum(htru.df[which(khtwu$cluster==2),9]==0)
cat(paste0(
  "As a predictive model, it presents the following results:\n\taccuracy\t: ",
  format(100*(tn+tp)/(tp+tn+fp+fn),digits=4),
  "%\n\tError rate\t: ",
  format(100*(fn+fp)/(tp+tn+fp+fn),digits=4),
  "%\n\tPresicion\t: ",
  format(100*(tn)/(tn+fp),digits=4),
  "%\n\tSpecificity\t: ",
  format(100*(tp)/(tp+fp),digits=4),
  "%\n\t\tThis results are given that the positive class was the minority class (1)",
  "\nAlso, the balance accuracy (",
  format(100*((tp)/(tp+fp)+(tn)/(tn+fp))/2,digits=4)        ,
  "%) is pretty good for such an imbalanced sample: "
)
)
rm(tn,fn,tp,fp)
cat(paste0("The total Variance explained by the clustering is : ",
           format(100*(khtwu$tot.withinss/khtwu$totss),digits=4),
           " %")
)
    #cat(paste("Average Silhouette width of both of the cluster is:",format(clusStats$avg.silwidth,digits=4) ,"\n"))
silho<-silhouette(khtwu$cluster, dist(scale(htru.df[,1:8])))
cat(paste0("The average Silhouette width of both of the clusters is: ",
           format(sum(silho[,3])/sum(silho[,1]==silho[,1]),digits=4)
)
)
for (x in 1:2){
  cat(paste0("The average Silhouette width of Cluster ",x," is: ",
             format(sum(silho[which(silho[,1]==x),3])/sum(silho[,1]==x),digits=4)
             ,"\n"
  )
  )
}
set.seed(1122)
khtwu<- kmeans(scale(htru.pca$x[,1:2]), centers=2,nstart=25)
fviz_cluster(khtwu, 
             data=scale(htru.df),
             main="Plotting using Kmeans Using PCA 2 Principal components")
silho<-silhouette(khtwu$cluster, dist(scale(htru.pca$x[,1:8])))
cat(paste0("The average Silhouette width of both of the clusters is: ",
           format(sum(silho[,3])/sum(silho[,1]==silho[,1]),digits=4)
)
)
for (x in 1:2){
  cat(paste0("The average Silhouette width of cluster ",x," is: ",
             format(sum(silho[which(silho[,1]==x),3])/sum(silho[,1]==x),digits=4)
             ,"\n"
  )
  )
}
rm(list=ls())






