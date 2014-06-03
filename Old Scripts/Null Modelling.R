require(ape)
require(vegan)
require(phytools)
require(geiger)
require(abind)
require(picante)
require(Rsundials)
require(nlme)
require(adephylo)
require(phylobase)
require(ecodist)
require(ade4)
require(bipartite)
require(geometry)
require(packfor)
require(GUniFrac)
require(SDMtools)
require(fBasics)
require(FD)
require(betapart)
?functional.beta.pair
###
# Null Modelling
###

null.model <- function(x, y){
  year <- (site.info[site.info$period == y, c("year", "seasn")])[1, 1]
  seasn <- ((site.info[site.info$period == y, 
                                         c("year", "seasn")])[1, 2])
  temp.comm.id <- site.info[site.info$year == year & site.info$seasn == seasn,
                            "s.id"] 
  temp.comm <- as.matrix(community[rownames(community) %in% temp.comm.id, ])
  temp.comm <- temp.comm[rowSums(temp.comm) > 4, ]
  temp.comm <- temp.comm[, colnames(temp.comm) %in% 
                           names(which(colSums(temp.comm) > 0))]
  rownames(x) <- sample(rownames(x), 
                        length(rownames(x)), replace = FALSE)
  x <- x[sort(rownames(x)), ]
  x <- as.matrix(x[rownames(x) %in% colnames(temp.comm), ])
  #x <- dist(x, method = "euclidean")
  #as.matrix(comdist(comm = temp.comm, dis = x))
  (functional.beta.pair(temp.comm, x))
  #beta.pair(temp.comm)
}

sim.null <- list()
sor.null <- list()
sne.null <- list()
for (period in 1:14){
  nulls.temp <- replicate(100, null.model(pcoa.axes, period))
  year <- (site.info[site.info$period == period, c("year", "seasn")])[1, 1]
  seasn <- (site.info[site.info$period == period, c("year", "seasn")])[1, 2]
  
  sim.array.temp <- simplify2array(lapply(nulls.temp["funct.beta.sim",], as.matrix))
  sor.array.temp <- simplify2array(lapply(nulls.temp["funct.beta.sor",], as.matrix))
  sne.array.temp <- simplify2array(lapply(nulls.temp["funct.beta.sne",], as.matrix))
  
  
  sim.null[[paste(year, seasn, sep = "")]] <- sim.array.temp
  sor.null[[paste(year, seasn, sep = "")]] <- sor.array.temp
  sne.null[[paste(year, seasn, sep = "")]] <- sne.array.temp
  
  
}

str(sim.null)


sim.obs <- list()
sor.obs <- list()
sne.obs <- list()
for (period in 1:14){
  year <- (site.info[site.info$period == period, c("year", "seasn")])[1, 1]
  seasn <- (site.info[site.info$period == period, c("year", "seasn")])[1, 2]
  sim.obs[[paste(year, seasn, sep = "")]] <- 
    as.matrix(fun.dissim[[paste(year, seasn, "bsim", sep = "")]])
  sor.obs[[paste(year, seasn, sep = "")]] <- 
    as.matrix(fun.dissim[[paste(year, seasn, "bsor", sep = "")]])
  sne.obs[[paste(year, seasn, sep = "")]] <- 
    as.matrix(fun.dissim[[paste(year, seasn, "bsne", sep = "")]])
}
sim.obs

observed.results <- list()
observed.results[["sor"]] <- sor.obs
observed.results[["sim"]] <- sim.obs
observed.results[["sne"]] <- sne.obs

save(observed.results, file = "C:/Users/Tom Bishop/Documents/Research Projects/PhD/2 Second Paper/Paper 2 Restart September/Data/observed.results.RData")

null.means <- apply(nulls, c(1:2), mean, na.rm = TRUE)
null.sd <- apply(nulls, c(1:2), sd, na.rm = TRUE)
ses <- (as.matrix(obs.func.dist) - null.means) / null.sd


test.ob <- site.info[site.info$year == 2006 & site.info$seasn == "J", ]
test.ob[order(as.numeric(test.ob$s.id)), ]

site.info[order(as.numeric(site.info$s.id)), ]
labels(fun.dissim[["2006Jbsor"]])
