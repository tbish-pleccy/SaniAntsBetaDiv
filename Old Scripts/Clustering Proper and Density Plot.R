library(cluster)
library(gclus)


i <- 2008
j <- "J"

clustering.methods <- NULL
group.mem <- NULL
dendrogram.list <- list()
altitude.list <- list()
for (i in unique(comm$year)){
  for (j in unique(comm$seasn)){
    comm.temp <- comm[comm$year == i & comm$seasn == j, ]
    comm.temp <- subset(comm.temp, select = - c(s.id, year, seasn, Site, rep))
    comm.temp[comm.temp > 0] <- 1
    comm.temp <- comm.temp[rowSums(comm.temp) > 4, colSums(comm.temp) > 0]
        
    td <- comm.dissim[[paste(i, j, "bsim", sep = "")]]
    spe <- comm.temp
    td.sites <- rownames(as.matrix(td))
    alts.s.id <- site.info[site.info$s.id %in% td.sites, c("altitude", "s.id")]
    
    td.clust.sing <- hclust(td, method = "single")
    td.clust.aver <- hclust(td, method = "average")
    td.clust.mcqu <- hclust(td, method = "mcquitty")
    td.clust.cent <- hclust(td, method = "centroid")
    td.clust.medi <- hclust(td, method = "median")
    td.clust.comp <- hclust(td, method = "complete")
    td.clust.ward <- hclust(td, method = "ward")
    
    td.clust.sing.co <- cophenetic(td.clust.sing)
    td.clust.aver.co <- cophenetic(td.clust.aver)
    td.clust.mcqu.co <- cophenetic(td.clust.mcqu)
    td.clust.cent.co <- cophenetic(td.clust.cent)
    td.clust.medi.co <- cophenetic(td.clust.medi)
    td.clust.comp.co <- cophenetic(td.clust.comp)
    td.clust.ward.co <- cophenetic(td.clust.ward)
    
    si <- cor(td, td.clust.sing.co)
    av <- cor(td, td.clust.aver.co)
    mc <- cor(td, td.clust.mcqu.co)
    ce <- cor(td, td.clust.cent.co)
    me <- cor(td, td.clust.medi.co)
    co <- cor(td, td.clust.comp.co)
    wa <- cor(td, td.clust.ward.co)
    
    meth <- data.frame(si, av, mc, ce, me, co, wa)
    
    best.meth <- names(which.max(meth))
    coph.cor <- max(meth)
    
    
    asw <- numeric(nrow(alts.s.id ))
    for (m in 2:(nrow(alts.s.id )-1)) {
      sil <- silhouette(cutree(td.clust.aver, k = m), td)
      asw[m] <- summary(sil)$avg.width
    }
    k.best.sil <- which.max(asw)
    
    
    grpdist <- function(X)
    {
      require(cluster)
      gr <- as.data.frame(as.factor(X))
      distgr <- daisy(gr, "gower")
      distgr
    }
    
    kt <- data.frame(k = 1:nrow(alts.s.id), r = 0)
    for (n in 2:(nrow(alts.s.id) - 1)) {
      gr <- cutree(td.clust.aver, n)
      distgr <- grpdist(gr)
      mt <- cor(td, distgr, method="pearson")
      kt[n, 2] <- mt
    }
    k.best.man <- which.max(kt$r)
    
    spe.chwo <- reorder.hclust(td.clust.aver, td)
    dendrogram.list[[paste(i, j, sep = "")]] <- spe.chwo
    altitude.list[[paste(i, j, sep = "")]] <- alts.s.id$altitude
    
    
    group.mem.temp <- data.frame(cutree(spe.chwo, k = 5))
    names(group.mem.temp) <- "group.num"
    group.mem <- rbind(group.mem, group.mem.temp)
    
    clustering.temp <- data.frame(i, j, best.meth, coph.cor, k.best.sil, 
                                  k.best.man, 
                                  si, av, mc, ce, me, co, wa)
    clustering.methods <- rbind(clustering.methods, clustering.temp)
    
    
    
    
    

  }
}

clustering.methods
group.mem$s.id <- rownames(group.mem)
group.mem <- merge(group.mem, site.info, by = "s.id")
group.mem[group.mem$year == 2009 & group.mem$seasn == "S", ]

Mmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mean(clustering.methods$k.best.man)
Mmode(clustering.methods$k.best.man)
mean(clustering.methods$k.best.sil)
Mmode(clustering.methods$k.best.sil)

tapply(clustering.methods$k.best.man, clustering.methods$j, mean)
tapply(clustering.methods$k.best.sil, clustering.methods$j, mean)
tapply(clustering.methods$k.best.man, clustering.methods$j, Mmode)
tapply(clustering.methods$k.best.sil, clustering.methods$j, Mmode)

?mode


library(lattice)
xyplot(group.num ~ altitude | as.factor(year) + seasn, data = group.mem)

tapply(group.mem$group.num, group.mem$altitude, Mmode)
tapply(group.mem$altitude, group.mem$group.num, Mmode)


group.mem[group.mem$group.num == 1, ]

td.clust.aver.co

plot(dist(alts.s.id), td.clust.aver.co)


par(mfrow = c(5, 3))
for (i in unique(comm$year)){
  for (j in unique(comm$seasn)){
    dend <- dendrogram.list[[paste(i, j, sep = "")]]
    dend.alt <- altitude.list[[paste(i, j, sep = "")]]
    plot(dend, sub="", 
         ylab = "Height", 
         main = paste(i, j), hang = -0.5, 
         labels = dend.alt
    )
    rect.hclust(dend, k = 5, 
                #border = "red",
                which = c(1, 2, 3, 4, 5), 
                border = c("red", "blue", "pink", "darkgreen", "orange"), 
                cluster = cutree(dend, k = 5))
  }
}



group.mem[group.mem$group.num == 1, c("year", "seasn", "group.num", "altitude")]

group.mem[group.mem$year == 2010 & 
            group.mem$seasn == "S", c("year", "seasn", "group.num", "altitude")]



dend.d <- as.dendrogram(dend)
# vector of colors labelColors = c('red', 'blue', 'pink', 'darkgrey',
# 'purple')
labelColors <-  c("red", "blue", "pink", "darkgreen", "orange")
# cut dendrogram in 4 clusters
clusMember <- cutree(dend, 5)
# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    a$label <- dend.alt
    
  }
  n
}
# using dendrapply
clusDendro <- dendrapply(dend.d, colLab)
# make plot
plot(clusDendro, main = "Cool Dendrogram")



head(group.mem)
par(mfrow = c(1, 5))
for (i in (unique(group.mem$group.num))){
  gp <- group.mem[group.mem$group.num == i, ]
  plot(density(gp$altitude), main = paste("Group", i, sep = " "), 
       xlim = c(800, 3100))
}


library(sm)
library(RColorBrewer)
brewer.pal(5, "RdBu")
clus.col <- colorRampPalette(c("blue", "red"))
par(mfrow = c(1, 1))
sm.density.compare(group.mem$altitude, as.factor(group.mem$group.num), 
                   col = clus.col(5), 
                   xlab = "Elevation (m a.s.l.)", lwd = 1, lty = 1)
title(main = "Density plot of cluster membership with altitude")
#colfill<-c(2:(2+length(levels(as.factor(group.mem$group.num))))) 
par(usr = c(1,10,1,10))
#text(9, 9, "hi")
#legend(x = 9, y = 9, levels(as.factor(group.mem$group.num)), fill=clus.col(5))
legend(x = 9, y = 10, levels(as.factor(group.mem$group.num)), lty = seq(1, 5, 1),
       lwd = 1, col = clus.col(5), title = "Group", bty = "n")



for (i in (unique(group.mem$group.num))){
  gp <- group.mem[group.mem$group.num == i, ]
  hist(gp$altitude, main = paste("Group", i, sep = " "), 
       xlim = c(800, 3100))
}
hist(gp$altitude)

?density

