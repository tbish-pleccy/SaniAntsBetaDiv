# Observed spatial taxonomic beta diversity

library(betapart)
tax.beta.obs <- lapply(community.list, beta.pair)
tax.beta.obs[["2006J"]][["beta.sne"]]

# Observed elevational difference

ele.dist.list <- list()
for(i in unique(site.info.reps$year)){
  for (j in unique(site.info.reps$seasn)){
    ele <- data.frame(site.info.reps[site.info.reps$year == i & 
                                       site.info.reps$seasn == j, "altitude"])
    rownames(ele) <- site.info.reps[site.info.reps$year == i & 
                                      site.info.reps$seasn == j, "s.id"]
    ele.dist.list[[paste(i, j, sep = "")]] <- dist(ele)
  }
}

# Observed spatial difference

spa.dist.list <- list()
for(i in unique(site.info.reps$year)){
  for (j in unique(site.info.reps$seasn)){
    spa <- as.matrix(site.info.reps[site.info.reps$year == i & 
                                       site.info.reps$seasn == j, c("x", "y")])
    spa.dist <- spDists(spa, longlat = TRUE)
    rownames(spa.dist) <- rownames(spa)
    colnames(spa.dist) <- rownames(spa)
    spa.dist.list[[paste(i, j, sep = "")]] <- as.dist(spa.dist)
  }
}


coord <- as.matrix(subset(site.info.temp, select = c(x, y)))
spa.dist <- as.dist(spDists(coord, longlat = TRUE))


plot(tax.beta.obs[["2008J"]][["beta.sim"]] ~ spa.dist.list[["2008J"]])
plot(ele.dist.list[["2008J"]] ~ spa.dist.list[["2008J"]])




