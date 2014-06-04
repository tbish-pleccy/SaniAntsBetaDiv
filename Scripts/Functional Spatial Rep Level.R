###############################
# Functional beta diversity.  #
# Spatial.                    #
# Rep level.                  #
###############################


# Save objects for observed spatial rep level
# Communities
for (i in labels(community.list)){
  number <- which(labels(community.list) == i)
  condor.comm <- community.list[[i]]
  save(condor.comm, 
       file = paste("Condor Input/Spatial Rep Level/Comm", number - 1, 
                    ".RData", sep = "")) 
}
# Traits
for (i in labels(community.traits)){
  number <- which(labels(community.traits) == i)
  condor.traits <- community.traits[[i]]
  save(condor.traits, 
       file = paste("Condor Input/Spatial Rep Level/Traits", number - 1, 
                    ".RData", sep = "")) 
}

# Import completed observed spatial rep level
fun.beta.obs <- list()
empty.list <- list()
for (i in labels(community.list)){
  n <- i
  assign(n, empty.list)
  fun.beta.obs[[n]] <- get(n)
  
  for (j in c("funct.beta.sor", "funct.beta.sim", "funct.beta.sne")){
    assign(j, empty.list)
    fun.beta.obs[[n]][[j]] <- get(j)
  }  
}

for (i in 1:14){
  path <- paste("Condor Output/Spatial Rep Level/fun_beta_obs_rep", i - 1,
                ".RData", sep = "")
  possibleError <- tryCatch(load(path), 
                            error = function(e) e)
  if (!inherits(possibleError, "simpleError")){
    load(path)
    time.name <- labels(community.list)[which(labels(community.list) == i)]
    fun.beta.obs[[time.name]][["sor"]][["funct.beta.sor"]] <- as.matrix(beta.results[[3]])
    fun.beta.obs[[time.name]][["sim"]][["funct.beta.sim"]] <- as.matrix(beta.results[[1]])
    fun.beta.obs[[time.name]][["sne"]][["funct.beta.sne"]] <- as.matrix(beta.results[[2]])
  }
}