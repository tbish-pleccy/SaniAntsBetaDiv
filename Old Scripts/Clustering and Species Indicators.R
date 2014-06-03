###############
library(FD)
library(betapart)
library(vegan)
library(boot)
library(ecodist)
library(sp)
library(ape)
library(MuMIn)
library(cluster)




communityJ08 <- comm[comm$year == 2008 & comm$seasn == "J", ]
#rownames(communityJ08) <- paste(communityJ08$Site, communityJ08$rep, sep = "")

communityJ08 <- subset(communityJ08, select = - c(s.id, year, seasn, Site, rep))
communityJ08[communityJ08 > 0] <- 1
communityJ08 <- communityJ08[rowSums(communityJ08) > 4, colSums(communityJ08) > 0]
communityJ082 <- comm[comm$year == 2008 & comm$seasn == "J", ]
communityJ082 <- subset(communityJ082, select = - c(s.id, year, seasn, Site, rep))
communityJ082 <- communityJ082[rownames(communityJ082) %in% rownames(communityJ08),
                               colnames(communityJ082) %in% colnames(communityJ08)]



str(comm.dissim)
str(fun.dissim)
str(comm.dissim.bsim)
str(comm.dissim.bnes)
str(ele.dissim)
str(spa.dissim)


td <- comm.dissim[["2008Jbsim"]]
spe <- communityJ08
spe2 <- communityJ082
?hclust
td.sites <- rownames(as.matrix(td))
alts.s.id <- site.info[site.info$s.id %in% td.sites, c("altitude", "s.id")]
td.clust <- hclust(td, method = "single")
plot(td.clust)
plot(td.clust, labels = alts.s.id$altitude)

td.clust.sing <- hclust(td, method = "single")
td.clust.aver <- hclust(td, method = "average")
td.clust.comp <- hclust(td, method = "complete")
td.clust.ward <- hclust(td, method = "ward")

td.clust.sing.co <- cophenetic(td.clust.sing)
td.clust.aver.co <- cophenetic(td.clust.aver)
td.clust.comp.co <- cophenetic(td.clust.comp)
td.clust.ward.co <- cophenetic(td.clust.ward)

cor(td, td.clust.sing.co)
cor(td, td.clust.aver.co)
cor(td, td.clust.comp.co)
cor(td, td.clust.ward.co)
plot(td.clust.aver, labels = alts.s.id$altitude)

summary(td.clust.aver)




# Optimal number of clusters according to silhouette widths
# (Rousseeuw quality index)
# *********************************************************

# Plot average silhouette widths (using Ward clustering) for all partitions 
# except for the trivial partition in a single group (k=1)
# First, create an empty vector in which the asw values will be written
library(cluster)
asw <- numeric(nrow(alts.s.id ))
for (k in 2:(nrow(alts.s.id )-1)) {
  sil <- silhouette(cutree(td.clust.aver, k=k), td)
  asw[k] <- summary(sil)$avg.width
}
k.best <- which.max(asw)
# The plot is produced by function plot.silhouette {cluster}
#windows(title="Silhouettes - Ward - k = 2 to n-1")
plot(1:nrow(alts.s.id), asw, type="h", 
     main="Silhouette-optimal number of clusters, Ward", 
     xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")



# Optimal number of clusters according to Mantel statistic (Pearson)
# ******************************************************************

# Function to compute a binary distance matrix from groups
grpdist <- function(X)
{
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}

# Run based on the Ward clustering
kt <- data.frame(k=1:nrow(alts.s.id), r=0)
for (i in 2:(nrow(alts.s.id)-1)) {
  gr <- cutree(td.clust.aver, i)
  distgr <- grpdist(gr)
  mt <- cor(td, distgr, method="pearson")
  kt[i,2] <- mt
}
kt
k.best <- which.max(kt$r)
# The plot is produced by function plot.silhouette {cluster}
#windows(title="Optimal number of clusters - Mantel")
plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - Ward", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)


# Silhouette plot of the final partition
# **************************************

# Choose the number of clusters
k <- 5
# Silhouette plot
cutg <- cutree(td.clust.aver, k=k)
sil <- silhouette(cutg, td)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(spe)[attr(silo,"iOrd")]
#windows(title="Silhouette plot - Average - K=4")
plot(silo, main="Silhouette plot - Chord - Average", 
     cex.names=0.8, col=cutg+1, nmax.lab=100)



# Final dendrogram with the selected groups
# *****************************************
library(gclus)
spe.chwo <- reorder.hclust(td.clust.aver, td)
summary(spe.chwo)


# Plot reordered dendrogram with group labels
#windows(title="Final dendrogram",8,6)
plot(spe.chwo, hang=-1, xlab="5 groups", sub="", 
     ylab="Height", 
     main="Chord - Average (reordered)", 
     labels=alts.s.id$altitude
     #labels=cutree(spe.chwo, k=k)
)
rect.hclust(spe.chwo, k=k)

# Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
source("hcoplot.R")         # hcoplot.R must be in the working directory
hcoplot(spe.ch.ward, spe.ch, k=4)


cutree(spe.chwo, k = 5)


# Heat map
# ********

# Heat map of the distance matrix ordered with the dendrogram
dend <- as.dendrogram(spe.chwo)
#windows(title="Heatmap - sites")
heatmap(as.matrix(td), Rowv=dend, symm=TRUE, margin=c(3,3))

# Ordered community table
# Species are ordered by their weighted averages on site scores
or <- vegemite(spe, spe.chwo, scale = "log")

# Heat map of the doubly ordered community table, with dendrogram
library(RColorBrewer)
#windows(title="Heatmap - species")
heatmap(t(spe[rev(or$species)]), Rowv=NA, Colv=dend,
        col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4), 
        ylab="Species (weighted averages of sites)", xlab="Sites")




library(indicspecies)
ind.group1 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 1, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7,
                         func = "IndVal.g", nboot = 999, alpha = 0.05)
ind.group2 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 2, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7,
                         func = "IndVal.g", nboot = 999, alpha = 0.05)
ind.group3 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 3, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7,
                         func = "IndVal.g", nboot = 999, alpha = 0.05)
ind.group4 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 4, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7,
                         func = "IndVal.g", nboot = 999, alpha = 0.05)
ind.group5 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 5, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7,
                         func = "IndVal.g", nboot = 999, alpha = 0.05)

print(ind.group1, confint = TRUE)
print(ind.group2, confint = TRUE)
print(ind.group3, confint = TRUE)
print(ind.group4, confint = TRUE)
print(ind.group5, confint = TRUE)

