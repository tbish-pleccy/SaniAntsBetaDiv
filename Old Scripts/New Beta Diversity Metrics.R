# Create functional beta diversity function using Blonder's method
hv.list <- list()
for (i in unique(rownames(community))){
  c <- pcoa.axes[rownames(pcoa.axes) %in% names(which(as.matrix(community)[paste(i), ] > 0)), ]
  hv.list[i] <- hypervolume(c, reps = 1000, bandwidth = 0.3, warn = FALSE, 
                            verbose = F)
  print(paste(i, "out of", length(unique(rownames(community))), sep = " "))
}
x <- c(1, 20)
bsor.calc <- function (x){
  hv1 <- hv.list[[x[1]]]
  hv2 <- hv.list[[x[2]]]
  hv_set <- hypervolume_set(hv1, hv2, reduction_factor = 1, 
                            check_memory = FALSE, verbose = FALSE)
  overlaps <- data.frame(t(sapply(hv_set@HVList, function(y) {y@Volume})))
  bsor <- (overlaps$Unique_1 + overlaps$Unique_2) / 
    ((2 * overlaps$Intersection) + (overlaps$Unique_1 + overlaps$Unique_2))
  bsim <- min(overlaps$Unique_1, overlaps$Unique_2) / 
    (overlaps$Intersection + min(overlaps$Unique_1, overlaps$Unique_2))
  bsne <- ((max(overlaps$Unique_1, overlaps$Unique_2) - 
             min(overlaps$Unique_1, overlaps$Unique_2)) / 
    ((2 * overlaps$Intersection) + min(overlaps$Unique_1, overlaps$Unique_2) + 
       max(overlaps$Unique_1, overlaps$Unique_2))) * 
    (overlaps$Intersection / (overlaps$Intersection + 
                                min(overlaps$Unique_1, overlaps$Unique_2)))  
  
  t(c(bsor, bsim, bsne))
}

comp.list.t <- comp.list

cbind(comp.list[1, ], matrix(data = c(bsor, bsim, bsne),nrow = 1, ncol = 3))
cbind(c(1, 2), c(bsor, bsim, bsne))

# mat.creation <- function(x){
#   test.m[x[, 1], x[, 2]] <- as.numeric(x[, 3])
# }


fBd <- function(x){ # x must be community data
  blank.m <- matrix(nrow = nrow(test.comm.mod), ncol = nrow(test.comm.mod))
  rownames(blank.m) <- rownames(test.comm.mod)
  colnames(blank.m) <- rownames(test.comm.mod)  
  comp.list <- NULL
  comp.list <- (cbind(rep(colnames(blank.m), each = length(colnames(blank.m))),
                     rep(colnames(blank.m), length.out = length(blank.m))))
  comp.list <- cbind(comp.list, t(apply(comp.list, MARGIN = 1, FUN = bsor.calc)))
  
  fbsor <- blank.m
  for(i in 1:nrow(comp.list)){
    fbsor[comp.list[i, 1], comp.list[i, 2]] <- as.numeric(comp.list[i, 3])
  }
  fbsim <- blank.m
  for(i in 1:nrow(comp.list)){
    fbsim[comp.list[i, 1], comp.list[i, 2]] <- as.numeric(comp.list[i, 4])
  }
  fbsne <- blank.m
  for(i in 1:nrow(comp.list)){
    fbsne[comp.list[i, 1], comp.list[i, 2]] <- as.numeric(comp.list[i, 5])
  }
  
  fBd.results <- list()
  fBd.results[["fbsor"]] <- fbsor
  fBd.results[["fbsim"]] <- fbsim
  fBd.results[["fbsne"]] <- fbsne
  fBd.results
}


plot(fbsor ~ fbsim)
plot(fbsor ~ fbsne)
plot(fbsim ~ fbsne)

commBhv <- community
commBhv <- commBhv[rowSums(commBhv) > 1, ]
ts.list <- list()
for (year in unique(site.info$year)){
  for(seasn in unique(site.info$seasn)){
    sel <- site.info[site.info$year %in% year & 
                       site.info$seasn %in% seasn, "s.id"]
    c <- commBhv[rownames(commBhv) %in% sel, ]
    rowSums(c) > 1
    ts.list[[paste(year, seasn, sep = "")]] <- c
  }
}
str(ts.list)

trait.list <- list()
for (year in unique(site.info$year)){
  for(seasn in unique(site.info$seasn)){
    sel <- site.info[site.info$year %in% year & 
                       site.info$seasn %in% seasn, "s.id"]
    c <- commBhv[rownames(commBhv) %in% sel, ]
    c <- c[, colSums(c) > 0]
    

    trait.list[[paste(year, seasn, sep = "")]] <- 
      pcoa.axes[rownames(pcoa.axes) %in% colnames(c), ]
  }
}



library(betapart)
start <- Sys.time()
ts.list.beta <- lapply(ts.list, FUN = fBd)
end <- Sys.time()
end - start

str(ts.list.beta)

ts.list.beta[["2012S"]]$beta.sim

bsor.calc(c("1", "2"))

