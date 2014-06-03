
par(mfrow = c(2, 4))
for (i in 1:8){
  require(ecodist)
  c.t <- comm.dissim.time[[paste(i, "bsim", sep = "")]]
  f.t <- fun.dissim.time[[paste(i, "bsim", sep = "")]]
  
  select <- site.info[site.info$Site == i & site.info$seasn == "J", "s.id"]
  
  tf <- labels(c.t) %in% select
  comm.dist <- as.dist(full(c.t)[tf, tf])
  fun.dist <- as.dist(full(f.t)[tf, tf])
  
  
  
  plot(comm.dist ~ fun.dist, ylim = c(0, 1), xlim = c(0, 1), main = i)
  abline(0, 1)
  
  
}




par(mfrow = c(4, 4))
for (seasn in c("J", "S")){
  for(year in unique(site.info$year)){
    require(ecodist)
    ele.dissim[[paste(year, seasn, sep = "")]]
    
    comm.dist <- comm.dissim[[paste(year, seasn, "bsim", sep = "")]]
    fun.dist <- fun.dissim[[paste(year, seasn, "bsim", sep = "")]]
    
    seasn.col <- ifelse(seasn == "J", "red", "skyblue")
    seasn.name <- ifelse(seasn == "J", "Wet Season", "Dry Season")
    plot(fun.dist ~ comm.dist, ylim = c(0, 1), xlim = c(0, 1),
         pch = 21, col = "white", bg = seasn.col,
         main = paste(year, sep = " "), 
         xlab = "Taxonomic Turnover", ylab = "Functional Turnover")
    abline(0, 1)
    
    
    
    
  }  
}






(which(ele.dissim[[paste(year, seasn, sep = "")]] == 300))
t.d <- (ele.dissim[[paste(year, seasn, sep = "")]])
t.df <- data.frame(as.matrix(t.d))
t.df[(t.df) == 300, (t.df) == 300]

x <- t.d
neighbour.comp <- function(x){
  row.c <- NULL
  for (i in 2:length(labels(x))){
    row.ct <- c(labels(x)[i:length(labels(x))])
    row.c <- c(row.c, row.ct)
  }
  t.d.v <- cbind(rep(labels(x)[1:length(labels(x)) - 1], 
                     times = as.vector(seq(length(labels(x)) - 1, 1, -1))),
                 row.c, as.vector(x))
  #t.d.v <- t.d.v[which(t.d.v[,3] == "300"),]
  t.d.v <- data.frame(t.d.v)
  
  names(t.d.v) <- c("col", "row", "diff")
  t.d.v$col <- as.numeric(levels(t.d.v$col))[t.d.v$col]
  t.d.v$row <- as.numeric(levels(t.d.v$row))[t.d.v$row]
  t.d.v$diff <- as.numeric(levels(t.d.v$diff)[t.d.v$diff])  
  col.temp <- site.info[site.info$s.id %in% t.d.v$col, c("s.id", "altitude")]
  colnames(col.temp) <- c("col", "col.elev")
  row.temp <- site.info[site.info$s.id %in% t.d.v$row, c("s.id", "altitude")]
  colnames(row.temp) <- c("row", "row.elev")
  t.d.v <- merge(merge(t.d.v, col.temp), row.temp)
  t.d.v
  
}

t.d.v <- neighbour.comp(t.d)

unique(c(as.vector(t.d.v$col), as.vector(t.d.v$row)))

rep(labels(t.d)[2:length(labels(t.d))], 
    times = as.vector(seq(length(labels(t.d)) - 1, 1, -1)))


comm.dist <- comm.dissim[[paste(2012, "J", "bsim", sep = "")]]
fun.dist <- fun.dissim[[paste(2012, "J", "bsim", sep = "")]]

test2 <- MRM(fun.dist ~ comm.dist, nperm = 1000)
test2
summary(lm(fun.dist ~ comm.dist))

plot(fun.dist ~ comm.dist)













par(mfrow = c(4, 4))
for (seasn in c("J", "S")){
  for(year in unique(site.info$year)){
    require(ecodist)
    ele.dist <- ele.dissim[[paste(year, seasn, sep = "")]]
    
    comm.dist <- comm.dissim[[paste(year, seasn, "bsim", sep = "")]]
    fun.dist <- fun.dissim[[paste(year, seasn, "bsim", sep = "")]]
    
    seasn.col <- ifelse(seasn == "J", "red", "skyblue")
    seasn.name <- ifelse(seasn == "J", "Wet Season", "Dry Season")
#     plot(fun.dist ~ comm.dist, ylim = c(0, 1), xlim = c(0, 1),
#          pch = 21, col = seasn.col, bg = seasn.col,
#          main = paste(year, sep = " "), 
#          xlab = "Taxonomic Turnover", ylab = "Functional Turnover")
#     abline(0, 1)
    
    c.neigh <- neighbour.comp(comm.dist)
    names(c.neigh)[3] <- "c.diff"
    f.neigh <- neighbour.comp(fun.dist)
    names(f.neigh)[3] <- "f.diff"
    e.neigh <- neighbour.comp(ele.dist)
    names(e.neigh)[3] <- "e.diff"
    neigh <- merge(merge(c.neigh, f.neigh), e.neigh)
    #neigh <- neigh[neigh$e.diff == 300, ]
    neigh$ratio <- neigh$c.diff / neigh$f.diff
    #points(f.diff ~ c.diff, data = neigh)
    neigh$comp <- as.factor(paste(neigh$row.elev, neigh$col.elev, sep = "-"))
    neigh$comp <- relevel(neigh$comp, "900-1200")
    #plot(log(ratio) ~ comp, data = neigh, main = paste(year, seasn, sep = " "))
plot(log(ratio) ~ e.diff, data = neigh, main = paste(year, seasn, sep = " "))

  }  
}


neigh$comp <- as.factor(paste(neigh$row.elev, neigh$col.elev, sep = "-"))
neigh$comp <- relevel(neigh$comp, "900-1200")
plot(log(ratio) ~ comp, data = neigh)

neigh$comp <- relevel(neigh$comp,
  c(sort(levels(neigh$comp))[7], sort(levels(neigh$comp))[1:6]))

neigh2 <- neigh[neigh$e.diff == 300, ]
neigh2$comp <- droplevels(neigh2$comp)


par(mfrow = c(1, 2), mar=c(5, 5, 5, 2))
plot(log(ratio) ~ comp, data = neigh2, xlab = "Elevational Comparison", 
     ylab = "log(Taxonomic Turnover / Functional Turnover)", pch = 21,
     bg = "black", xaxt = "n", tcl=0.2, mgp=c(2.5,0.5,0))
axis(1, labels = FALSE)
labels <- levels(neigh2$comp)
text(x =  seq_along(labels), y = par("usr")[3] - 0.2, srt = 25, adj = 1,
     labels = labels, xpd = TRUE, cex = 0.5)

plot(log(ratio) ~ e.diff, data = neigh, xlab = "Elevational Distance", 
     ylab = "log(Taxonomic Turnover / Functional Turnover)", pch = 21, 
     bg = "black")





par(mfrow=c(1,1),mar = c(5, 4, 4, 2) + 0.1)
NBN<-read.table("C:\\Users\\Tom Bishop\\Documents\\Research Projects\\Butterfly Project\\Data and Scripts for Manuscript\\Data\\NBN data.csv",header=T,sep=",")
NBN<-NBN[-13,]
NBN$rank<-rank((NBN$percent.greater))
NBN<-NBN[order(NBN$rank,decreasing=TRUE),]
NBN$fraction<-paste(NBN$greater.6500.records,"/",NBN$Total,sep="")
par(fg="black",mar = c(8, 6, 4, 2) + 0.1)
x<-barplot(height=NBN$percent.greater, names.arg=NBN$Taxon,xaxt="n",ylim=c(0,70),ylab="Percentage of species with greater\n than 6,500 records",bty="l",
           tcl=0.2,mgp=c(2.5,0.5,0))
text(x=x,y=-1.25,srt=45,adj=1,labels=NBN$Taxon,xpd=TRUE,cex=1)
text(x=x,y=NBN$percent.greater+3,labels=NBN$fraction,cex=0.8,srt=45,adj=c(0.2,0.1))
box(bty="l")
