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


save(fun.dissim)