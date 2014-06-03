### Input data for Condor (convex hull method)
# Input community data

# Observed communities
ch.comm.list <- list()
for (year in unique(site.info$year)){
  for(seasn in unique(site.info$seasn)){
    sel <- site.info[site.info$year %in% year & 
                       site.info$seasn %in% seasn, "s.id"]
    c <- community[rownames(community) %in% sel, ]
    c <- c[rowSums(c) > 4, ]
    c <- c[, colSums(c) > 0]
    ch.comm.list[[paste(year, seasn, sep = "")]] <- c
  }
}

# Observed traits
ch.traits.list <- list()
for (i in names(ch.comm.list)){
  species <- names(ch.comm.list[[i]])
  traits.mod <- pcoa.axes[rownames(pcoa.axes) %in% species, ]
  ch.traits.list[[i]] <- traits.mod
}

# Randomised traits
library(abind)
ch.rand.traits.list <- list()
for (i in names(ch.comm.list)){
  randomisation <- NULL
  #randomisation <- list()
  for (j in 1:1000){
    traits.mod <- pcoa.axes
    rownames(traits.mod) <- sample(rownames(traits.mod), 
                                  length(rownames(traits.mod)), 
                                  replace = FALSE)
    traits.mod <- traits.mod[sort(rownames(traits.mod)), ]
    obs.species <- names(ch.comm.list[[i]])
    traits.mod <- traits.mod[rownames(traits.mod) %in% obs.species, ]  
    randomisation <- abind(randomisation, traits.mod, along = 3)
    
    #randomisation[[paste(i, "R", j, sep = "")]] <- traits.mod
  }
  ch.rand.traits.list[[i]] <- randomisation
}

# Save Observed Communities
save.frame <- data.frame(rep(names(ch.comm.list), each = 1000), seq(1, 14000, 1))
names(save.frame) <- c("comm.name", "index")
head(save.frame)

observed.communities.location <- 
  "C:/Users/Tom Bishop/Documents/Research Projects/PhD/2 Second Paper/Paper 2 Restart September/Data/Observed Communities Convex Hull"
for (i in 1:nrow(save.frame)){
  select <- save.frame[i, ]
  write.csv(as.matrix(ch.comm.list[[select$comm.name]]), row.names = TRUE, 
              paste(observed.communities.location, "/","Comm",  select$index - 1, ".csv", sep = ""))
}

save(ch.comm.list, file = "C:/Users/Tom Bishop/Documents/Research Projects/PhD/2 Second Paper/Paper 2 Restart September/Data/R Objects/ch.comm.list.RData")
save(ch.rand.traits.list, file = "C:/Users/Tom Bishop/Documents/Research Projects/PhD/2 Second Paper/Paper 2 Restart September/Data/R Objects/ch.rand.traits.list.RData")
save(ele.dissim, file = "C:/Users/Tom Bishop/Documents/Research Projects/PhD/2 Second Paper/Paper 2 Restart September/Data/R Objects/ele.dissim.RData")


# Save Randomised Traits
save.frame$rel.ind <- rep(seq(1, 1000, 1), 14)
random.traits.location <- 
  "C:/Users/Tom Bishop/Documents/Research Projects/PhD/2 Second Paper/Paper 2 Restart September/Data/Randomised Traits Convex Hull"
for (i in 1:nrow(save.frame)){
  select <- save.frame[i, ]
  randomisedtraits <- as.matrix(ch.rand.traits.list[[select$comm.name]][,,select$rel.ind])
  write.csv(randomisedtraits, file = 
           paste(random.traits.location, "/", "Traits", select$index - 1, 
                 ".csv", sep = ""))
}


test.save.object <- as.matrix(ch.rand.traits.list[[i]][,,j])
save(test.save.object, file =  
     "C:/Users/Tom Bishop/Documents/Research Projects/PhD/2 Second Paper/Paper 2 Restart September/Data/Condor Outputs/testhello.RData")
testhello <- load("C:/Users/Tom Bishop/Documents/Research Projects/PhD/2 Second Paper/Paper 2 Restart September/Data/Condor Outputs/testhello.RData")
testhello
######
# Condor R function

# R_script = functional_beta_condor.R
# input_files = community.csv
# indexed_input_files = randomisedtraits.csv
# indexed_output_files = betaresults.RData
# indexed_log = log
# total_jobs = 1000

assign("fandangle", seq(1, 100))
fandangle

