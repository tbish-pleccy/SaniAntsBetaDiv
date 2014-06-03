head(comm)
library(FD)

head(community)

func.div <- dbFD(tFD.temp, community)
func.div <- data.frame(func.div)
nrow(tFD.temp)
ncol(community)

pcoa.axes



## which communities did we actually use?
used.comms <- NULL
for (i in unique(comm$year)){
  for (j in unique(comm$seasn)){
    c.ab <- comm[comm$year == i & comm$seasn == j, ]  # c.ab is abundance data
    
    c.pa <- subset(c.ab, select = - c(s.id, Site, rep, year, seasn))
    rownames(c.pa) <- c.ab$s.id
    c.pa[c.pa > 0] <- 1  # from here, c.pa is presence/absence data
    c.pa <- c.pa[- which(rowSums(c.pa) < 5), ]  # use only comms with >5 species
    absent.species <- which(colSums(c.pa) == 0)
    c.pa <- c.pa[, - absent.species]  # remove empty species
    c.pa <- c.pa[as.vector(which(rowSums(c.pa) > 0)), ]  # remove empty communities
    
    used.comms <- c(used.comms, rownames(c.pa))
  }
}






comm.ascores <- NULL
for (i in used.comms){
  dat <- community[i, ]
  sp <- names(which(colSums(dat) > 0))
  if (length(sp) > 1){
    axis.scores <- data.frame(pcoa.axes[which(rownames(pcoa.axes) %in% sp), ])
    a1.max <- max(axis.scores[,1])  
    a1.min <- min(axis.scores[,1])
    a1.av <- mean(axis.scores[,1])
    
    a2.max <- max(axis.scores[,2])
    a2.min <- min(axis.scores[,2])
    a2.av <- mean(axis.scores[,2])
    
    a3.max <- max(axis.scores[,3])
    a3.min <- min(axis.scores[,3])
    a3.av <- mean(axis.scores[,3])
    
    a4.max <- max(axis.scores[,4])
    a4.min <- min(axis.scores[,4])
    a4.av <- mean(axis.scores[,4])
    
  }
  else {
    axis.scores <- data.frame(pcoa.axes[which(rownames(pcoa.axes) %in% sp), ])
    a1.max <- (axis.scores[1,])  
    a1.min <- (axis.scores[1,])
    a1.av <- mean(axis.scores[1,])
    
    a2.max <- (axis.scores[2,])
    a2.min <- (axis.scores[2,])
    a2.av <- mean(axis.scores[2,])
    
    a3.max <- (axis.scores[3,])
    a3.min <- (axis.scores[3,])
    a3.av <- mean(axis.scores[3,])
    
    a4.max <- (axis.scores[4,])
    a4.min <- (axis.scores[4,]) 
    a4.av <- mean(axis.scores[4,])
    
  }
  s.id <- as.numeric(i)
  temp <- data.frame(s.id, a1.max, a1.min, a1.av, a2.max, a2.min, a2.av, a3.max, a3.min, a3.av, a4.max, a4.min, a4.av)
  comm.ascores <- rbind(comm.ascores, temp)

}

comm.ascores$rangea1 <- abs(comm.ascores$a1.max - comm.ascores$a1.min)
comm.ascores$rangea2 <- abs(comm.ascores$a2.max - comm.ascores$a2.min)
comm.ascores$rangea3 <- abs(comm.ascores$a3.max - comm.ascores$a3.min)
comm.ascores$rangea4 <- abs(comm.ascores$a4.max - comm.ascores$a4.min)



comm.ascores$s.id

comm.ascores <- merge(comm.ascores, site.info, by = "s.id")
comm.ascores <- comm.ascores[order(comm.ascores$s.id), ]
rownames(comm.ascores) <- comm.ascores$s.id

library(lattice)
xyplot(rangea1 ~ altitude | year + seasn, data = comm.ascores)
xyplot(rangea2 ~ altitude | year + seasn, data = comm.ascores)
xyplot(rangea3 ~ altitude | year + seasn, data = comm.ascores)
xyplot(rangea4 ~ altitude | year + seasn, data = comm.ascores)
xyplot(a1.av ~ altitude | year + seasn, data = comm.ascores)
xyplot(a2.av ~ altitude | year + seasn, data = comm.ascores)
xyplot(a3.av ~ altitude | year + seasn, data = comm.ascores)
xyplot(a4.av ~ altitude | year + seasn, data = comm.ascores)





head(comm.ascores)
str(comm.ascores)
comm.ascores$altitude <- as.numeric(comm.ascores$altitude)
comm.ascores$alti2 <- comm.ascores$altitud^2
comm.ascores$fyear <- as.factor(comm.ascores$year)
comm.ascores$fsite <- as.factor(comm.ascores$Site)
comm.ascores$altiC <- scale(comm.ascores$altitude)
comm.ascores$alti2C <- scale(comm.ascores$alti2)
comm.ascores$seasnB <- rep(c(0, 1), length.out = nrow(comm.ascores))

##################################################
library(lme4)
a1ran.full <- lmer(rangea1 ~ altiC*alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
summary(a1ran.full)

# test 3 way interaction

a1ran.1.1 <- lmer(rangea1 ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
anova(a1ran.full, a1ran.1.1) # not significant. Remove 3 way and continue with a1ran.1.1

# test 2 way interactions
a1ran.2 <- lmer(rangea1 ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)

a1ran.2.1 <- lmer(rangea1 ~ altiC*alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
a1ran.2.2 <- lmer(rangea1 ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)
a1ran.2.3 <- lmer(rangea1 ~ altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)
anova(a1ran.2, a1ran.2.1) # alti2C*seasnB = 0.09825
anova(a1ran.2, a1ran.2.2) # altiC*seasnB = 0.1286 remove this interaction. 
anova(a1ran.2, a1ran.2.3) # altiC*alti2C = 0.0019** significant. 

# continue testing two way interactions. 
a1ran.3 <- lmer(rangea1 ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)

a1ran.3.1 <- lmer(rangea1 ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)
a1ran.3.2 <- lmer(rangea1 ~ altiC + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)
anova(a1ran.3, a1ran.3.1) # alti2*seasnB = 0.3145 remove this interaction. 
anova(a1ran.3, a1ran.3.2) # altiC*alti2C - 0.00214** significant

# test the last two interaction again. 
a1ran.4 <- lmer(rangea1 ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)

a1ran.4.1 <- lmer(rangea1 ~ altiC + alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
anova(a1ran.4, a1ran.4.1) # altiC*alti2C = 0.002215** significant. 
summary(a1ran.4)
# remove seasn

a1ran.4.2 <- lmer(rangea1 ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)
anova(a1ran.4, a1ran.4.2) # seasnB = 0.002717** significant. 

a1ran.final <- lmer(rangea1 ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                    control = list(maxIter = 5200), REML = FALSE)
summary(a1ran.final)

library(MuMIn)
r.squaredGLMM(x = a1ran.final)




comm.ascores$fit.a1ran <- fitted(a1ran.final)

par(mfrow = c(3, 3), mar = c(5, 5, 4, 2))
for (i in unique(comm.ascores$year)){
  raw <- comm.ascores[comm.ascores$seasn == "J" & comm.ascores$year == i, ]
  plot(rangea1 ~ altitude, data = raw, pch = 21, col = "white", bg = "orange", 
       bty = "l", 
       xlab = "Elevation (m a.s.l.)", ylab = "Range in A1", main = i, type = "p")
  
  raw2 <- comm.ascores[comm.ascores$seasn == "S" & comm.ascores$year == i, ]
  points(rangea1 ~ altitude, data = raw2, pch = 21, col = "white", bg = "lightblue")
  lines(fit.a1ran ~ altitude, data = raw, col = "red")
  lines(fit.a1ran ~ altitude, data = raw2, col = "blue")  
} 


# so the range in axis 1 shrinks with elevation at around 2400m.
# which complexes are we losing? The extreme positives? The extreme negatives?
# If so we would see a signifcant model with a direction. 
# Alternatively, we may be losing both positive and negative extremes. 
# In this case, an insignficant model would inidcate this. 

a1av.full <- lmer(a1.av ~ altiC*alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
summary(a1av.full)

# test three way
a1av.1.1 <- lmer(a1.av ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
anova(a1av.full, a1av.1.1) # marginally significant! but not convinced, lets remove it. 

# test 2 ways
a1av.2 <- lmer(a1.av ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a1av.2.1 <- lmer(a1.av ~ altiC*alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a1av.2.2 <- lmer(a1.av ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a1av.2.3 <- lmer(a1.av ~ altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

anova(a1av.2, a1av.2.1) # alti2C*seasnB = 0.1402 remove
anova(a1av.2, a1av.2.2) # altiC*seasnB = 0.1277
anova(a1av.2, a1av.2.3) # altiC*alti2C = 0.9666

# test last 2 two ways

a1av.3 <- lmer(a1.av ~ altiC*alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a1av.3.1 <- lmer(a1.av ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a1av.3.2 <- lmer(a1.av ~ alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a1av.3, a1av.3.1) # altiC*seasnB = 0.665 remove
anova(a1av.3, a1av.3.2) # altiC*alti2C = 0.09944

# test last 2 way
a1av.4 <- lmer(a1.av ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)

a1av.4.1 <- lmer(a1.av ~ altiC + alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a1av.4, a1av.4.1) # altiC*alti2C = 0.09955 remove

# test single effects
a1av.5 <- lmer(a1.av ~ altiC + alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a1av.5.1 <- lmer(a1.av ~ altiC + alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a1av.5.2 <- lmer(a1.av ~ altiC + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a1av.5.3 <- lmer(a1.av ~ alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
anova(a1av.5, a1av.5.1) # seasnB = 0.3805
anova(a1av.5, a1av.5.2) # alti2C = 0.995 remove
anova(a1av.5, a1av.5.3) # altiC = 0.5563

# test other singles
a1av.6 <- lmer(a1.av ~ altiC + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a1av.6.1 <- lmer(a1.av ~ altiC + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a1av.6.2 <- lmer(a1.av ~ seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a1av.6, a1av.6.1) # seasnB = 0.38 remvoe
anova(a1av.6, a1av.6.2) # altiC = 0.00000

# final test
a1av.7 <- lmer(a1.av ~ altiC + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a1av.7.1 <- lmer(a1.av ~ (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
anova(a1av.7, a1av.7.1) # altiC = 0.0000

# FINAL
a1av.final <- lmer(a1.av ~ altiC + (1 | fyear / seasn), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
summary(a1av.final)
library(MuMIn)
r.squaredGLMM(x = a1av.final)




comm.ascores$fit.a1av <- fitted(a1av.final)


par(mfrow = c(3, 3), mar = c(5, 5, 4, 2))
for (i in unique(comm.ascores$year)){
  raw <- comm.ascores[comm.ascores$seasn == "J" & comm.ascores$year == i, ]
  plot(a1.av ~ altitude, data = raw, pch = 21, col = "white", bg = "orange", 
       bty = "l", 
       xlab = "Elevation (m a.s.l.)", ylab = "Average A1 Score", main = i, type = "p", 
       ylim = c(-0.6, 0.6))
  
  raw2 <- comm.ascores[comm.ascores$seasn == "S" & comm.ascores$year == i, ]
  points(a1.av ~ altitude, data = raw2, pch = 21, col = "white", bg = "lightblue")
  lines(fit.a1av ~ altitude, data = raw, col = "red")
  lines(fit.a1av ~ altitude, data = raw2, col = "blue")  
} 


# Mild increase in average A1 scores with elevation. So an increase in those traits
# associated with surface foraging. But this is a weak effect with lots of variation. 


# So we see a shift from negative to positive scores for A1. This is weak however. 
# Indicating that there is a small shift from negative to positive A1 scores,
# but also a loss of the extreme values - giving us the reduced range as seen above. 

########################
####### Axis 2 #########
########################

a2ran.full <- lmer(rangea2 ~ altiC*alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
summary(a2ran.full)

# test three way
a2ran.1.1 <- lmer(rangea2 ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a2ran.full, a2ran.1.1) # take it out

# test two ways
a2ran.2 <- lmer(rangea2 ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

a2ran.2.1 <- lmer(rangea2 ~ altiC*alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
a2ran.2.2 <- lmer(rangea2 ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
a2ran.2.3 <- lmer(rangea2 ~ altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

anova(a2ran.2, a2ran.2.1) # alti2C*seasnB = 0.3119
anova(a2ran.2, a2ran.2.2) # altiC*seasnB = 0.3444 remove
anova(a2ran.2, a2ran.2.3) # altiC*alti2C = 0.0000

# test other 2 ways

a2ran.3 <- lmer(rangea2 ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

a2ran.3.1 <- lmer(rangea2 ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
a2ran.3.2 <- lmer(rangea2 ~ altiC + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

anova(a2ran.3, a2ran.3.1) #alit2C*seasnB = 0.6124 remove
anova(a2ran.3, a2ran.3.2) # altiC*alti2C = 0.00000

# test last 3 way again
a2ran.4 <- lmer(rangea2 ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

a2ran.4.1 <- lmer(rangea2 ~ altiC + alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a2ran.4, a2ran.4.1) # significant, keep altiC*alti2C

a2ran.4.2 <- lmer(rangea2 ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a2ran.4, a2ran.4.2) # seasnB = 0.09567

# keep going
a2ran.5 <- lmer(rangea2 ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

a2ran.5.1 <- lmer(rangea2 ~ altiC + alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a2ran.5, a2ran.5.1) # interaction is significant!


a2ran.final <- lmer(rangea2 ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

summary(a2ran.final)
r.squaredGLMM(x = a2ran.final)


comm.ascores$fit.a2ran <- fitted(a2ran.final)

par(mfrow = c(3, 3), mar = c(5, 5, 4, 2))
for (i in unique(comm.ascores$year)){
  raw <- comm.ascores[comm.ascores$seasn == "J" & comm.ascores$year == i, ]
  plot(rangea2 ~ altitude, data = raw, pch = 21, col = "white", bg = "orange", 
       bty = "l", 
       xlab = "Elevation (m a.s.l.)", ylab = "Range in A2", main = i, type = "p")
  
  raw2 <- comm.ascores[comm.ascores$seasn == "S" & comm.ascores$year == i, ]
  points(rangea2 ~ altitude, data = raw2, pch = 21, col = "white", bg = "lightblue")
  lines(fit.a2ran ~ altitude, data = raw, col = "red")
  lines(fit.a2ran ~ altitude, data = raw2, col = "blue")  
} 


# bell shaped curve of range in A2 with elevation.
# range dips around 1500m. Rise to a peak around 2100 - 2400m and then tails off. 


a2av.full <- lmer(a2.av ~ altiC*alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
summary(a2av.full)

# remove 3 way
a2av.1.1 <- lmer(a2.av ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
anova(a2av.full, a2av.1.1) # remove it. 


# test 2 ways
a2av.2 <- lmer(a2.av ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)

a2av.2.1 <- lmer(a2.av ~ altiC*alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a2av.2.2 <- lmer(a2.av ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a2av.2.3 <- lmer(a2.av ~ altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

anova(a2av.2, a2av.2.1) # alti2C*seasnB = 0.3826
anova(a2av.2, a2av.2.2) # altiC*seasnB = 0.4457 remove
anova(a2av.2, a2av.2.3) # altiC*alti2C = 0.00000


# test other 2 ways.
a2av.3 <- lmer(a2.av ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a2av.3.1 <- lmer(a2.av ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a2av.3.2 <- lmer(a2.av ~ altiC + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a2av.3, a2av.3.1) # alti2C*seasnB = 0.429 remove
anova(a2av.3, a2av.3.2) # altiC*seasnB = 00000

# test last two way
a2av.4 <- lmer(a2.av ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a2av.4.1 <- lmer(a2.av ~ altiC + alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a2av.4, a2av.4.1) # keep two way of altiC*alti2C


a2av.4.2 <- lmer(a2.av ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a2av.4, a2av.4.2) # seasnB = 0.71

# keep going

a2av.5 <- lmer(a2.av ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a2av.5.1 <- lmer(a2.av ~ altiC + alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a2av.5, a2av.5.1) # 0.000 keep interaction

a2av.final <- lmer(a2.av ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)


summary(a2av.final)
library(MuMIn)
r.squaredGLMM(x = a2av.final)




comm.ascores$fit.a2av <- fitted(a2av.final)
par(mfrow = c(3, 3), mar = c(5, 5, 4, 2))
for (i in unique(comm.ascores$year)){
  raw <- comm.ascores[comm.ascores$seasn == "J" & comm.ascores$year == i, ]
  plot(a2.av ~ altitude, data = raw, pch = 21, col = "white", bg = "orange", 
       bty = "l", 
       xlab = "Elevation (m a.s.l.)", ylab = "Average A2 Score", main = i, type = "p")
  
  raw2 <- comm.ascores[comm.ascores$seasn == "S" & comm.ascores$year == i, ]
  points(a2.av ~ altitude, data = raw2, pch = 21, col = "white", bg = "lightblue")
  lines(fit.a2av ~ altitude, data = raw, col = "red")
  lines(fit.a2av ~ altitude, data = raw2, col = "blue")  
} 

par(mfrow = c(1, 1))
plot(rangea2 ~ a2.av, data = comm.ascores)

par(mar  = c(5, 4, 4, 2))

# so increase in range of A2 at 2400m is due to a huge decrease in the average at
# this elevation. Average A2 score tends to be around 0 for most elevations,
# indicating that a broad range of strategies on this axis are present. 
# At 2400m, the average becomes very negative, indicating the introduction of 
# species with large bodies and small mandibles (large omnivores). Average A2
# scores then gradually increase back to around 0 with increasing elevation up
# to 3000m. 


########################
####### Axis 3 #########
########################

a3ran.full <- lmer(rangea3 ~ altiC*alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
summary(a3ran.full)

# test three way
a3ran.1.1 <- lmer(rangea3 ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a3ran.full, a3ran.1.1) # apparently it is significant.

# remove two ways
a3ran.2 <- lmer(rangea3 ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

a3ran.2.1 <- lmer(rangea3 ~ altiC*alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
a3ran.2.2 <- lmer(rangea3 ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
a3ran.2.3 <- lmer(rangea3 ~ altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a3ran.2, a3ran.2.1) # altiC2*seasnB = 1
anova(a3ran.2, a3ran.2.2) # altiC*seasnB = 1
anova(a3ran.2, a3ran.2.3) # altiC*alti2C = 0.0000

# test other two ways
a3ran.3 <- lmer(rangea3 ~ altiC*alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

a3ran.3.1 <- lmer(rangea3 ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
a3ran.3.2 <- lmer(rangea3 ~ alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a3ran.3, a3ran.3.1) # altiC*seasnB = 0.9
anova(a3ran.3, a3ran.3.2) # altiC*alti2C = 0.0000

# test last 2 way

a3ran.4 <- lmer(rangea3 ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)

a3ran.4.1 <- lmer(rangea3 ~ altiC + alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a3ran.4, a3ran.4.1) # altiC*alti2C = 0.0000

a3ran.4.2 <- lmer(rangea3 ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
anova(a3ran.4, a3ran.4.2) # seasnB = 0.01588

summary(a3ran.4)

a3ran.final <- lmer(rangea3 ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE, verbose = TRUE)
summary(a3ran.final)
library(MuMIn)
r.squaredGLMM(x = a3ran.final)


1comm.ascores$fit.a3ran <- fitted(a3ran.final)

par(mfrow = c(3, 3), mar = c(5, 5, 4, 2))
for (i in unique(comm.ascores$year)){
  raw <- comm.ascores[comm.ascores$seasn == "J" & comm.ascores$year == i, ]
  plot(rangea3 ~ altitude, data = raw, pch = 21, col = "white", bg = "orange", 
       bty = "l", 
       xlab = "Elevation (m a.s.l.)", ylab = "Range in A3", main = i, type = "p")
  
  raw2 <- comm.ascores[comm.ascores$seasn == "S" & comm.ascores$year == i, ]
  points(rangea3 ~ altitude, data = raw2, pch = 21, col = "white", bg = "lightblue")
  lines(fit.a3ran ~ altitude, data = raw, col = "red")
  lines(fit.a3ran ~ altitude, data = raw2, col = "blue")  
} 

# Range in A3 is greatest at 2400m-ish. Similar to A2. 


a3av.full <- lmer(a3.av ~ altiC*alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
summary(a3av.full)

# remove 3 way
a3av.1.1 <- lmer(a3.av ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
anova(a3av.full, a3av.1.1) # remove it. 

# test 2 ways
a3av.2 <- lmer(a3.av ~ altiC*alti2C + altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)


a3av.2.1 <- lmer(a3.av ~ altiC*alti2C + altiC*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a3av.2.2 <- lmer(a3.av ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a3av.2.3 <- lmer(a3.av ~ altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a3av.2, a3av.2.1) # alti2C*seasnB = 0.9
anova(a3av.2, a3av.2.2) # altiC*seasnB = 0.99 REMOVE
anova(a3av.2, a3av.2.3) # altiC*alti2C = 0.000


# test 2 ways
a3av.3 <- lmer(a3.av ~ altiC*alti2C + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a3av.3.1 <- lmer(a3.av ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a3av.3.2 <- lmer(a3.av ~ altiC + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a3av.3, a3av.3.1) # alti2C*seasnB = 0.7238 REMOVE
anova(a3av.3, a3av.3.2) # altiC*alti2C = 0.0000


# test last 2 way

a3av.4 <- lmer(a3.av ~ altiC*alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a3av.4.1 <- lmer(a3.av ~ altiC + alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a3av.4, a3av.4.1) # altiC*alti2C = 0.0000

a3av.4.2 <- lmer(a3.av ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a3av.4, a3av.4.2) # seasnB = 0.271 REMOVE

# test again....
a3av.5 <- lmer(a3.av ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)

a3av.5.1 <- lmer(a3.av ~ altiC + alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a3av.5, a3av.5.1) # altiC*alti2C = 0.000

a3av.final <- lmer(a3.av ~ altiC*alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
summary(a3av.final)

library(MuMIn)
r.squaredGLMM(x = a3av.final)

comm.ascores$fit.a3av <- fitted(a3av.final)
par(mfrow = c(3, 3), mar = c(5, 5, 4, 2))
for (i in unique(comm.ascores$year)){
  raw <- comm.ascores[comm.ascores$seasn == "J" & comm.ascores$year == i, ]
  plot(a3.av ~ altitude, data = raw, pch = 21, col = "white", bg = "orange", 
       bty = "l", 
       xlab = "Elevation (m a.s.l.)", ylab = "Average A3 Score", main = i, type = "p")
  
  raw2 <- comm.ascores[comm.ascores$seasn == "S" & comm.ascores$year == i, ]
  points(a3.av ~ altitude, data = raw2, pch = 21, col = "white", bg = "lightblue")
  lines(fit.a3av ~ altitude, data = raw, col = "red")
  lines(fit.a3av ~ altitude, data = raw2, col = "blue")  
} 

# Hm. Peaks at 2400m but is otherwise flat?

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


library(ape)
library(ade4)
scale.traits <- tFD.temp
scale.traits <- scale(scale.traits, center = TRUE, scale = TRUE)
rownames(scale.traits) <- seq(1:nrow(scale.traits))
rownames(scale.traits) <- rownames(tFD.temp)
dm <- dist(scale.traits)
is.euclid(dm)
pcoa <- pcoa(dm)
new.pcoa <- dudi.pco(dm, scannf = FALSE, full = TRUE)
head(new.pcoa$li)
head(pcoa.axes)
head(pcoa$vectors) # all correct! dudi.pco producing same signs as pcoa.
                   # easy plotting. 

# pcoa.test <- pcoa
# pcoa.test$vectors[, 1] <- new.pcoa$li[, 1]
# pcoa.test$vectors[, 2] <- new.pcoa$li[, 2]
# pcoa.test$vectors[, 3] <- new.pcoa$li[, 3]
# pcoa.test$vectors[, 4] <- new.pcoa$li[, 4]
# names(pcoa.test$vectors) <- c("Axis.1", "Axis.2", "Axis.3", "Axis.4")

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
par(cex = 0.7)
biplot(pcoa, scale.traits, plot.axes = c(1, 2))
biplot(pcoa, scale.traits, plot.axes = c(2, 3))
biplot(pcoa, scale.traits, plot.axes = c(3, 4))
par(cex = 1)
biplot(pcoa.test, scale.traits, plot.axes = c(1, 3))

eigens <- pcoa$values
mean(eigens$Eigenvalues)
  
plot.axes = c(1, 2, 3, 4)
n <- nrow(scale.traits)
points.stand <- scale(pcoa$vectors[, plot.axes])
S <- cov(scale.traits, points.stand)
U <- S %*% diag((pcoa$values$Eigenvalues[plot.axes]/(n - 
                                                    1))^(-0.5))
colnames(U) <- colnames(pcoa$vectors[, plot.axes])
U


#####################################

par(mfrow = c(4, 4))
for (i in unique(comm$year)){
  for (j in as.vector(unique(comm$seasn))){
    community.temp <- comm[comm$year == i & comm$seasn == j, ]
    rownames(community.temp) <- paste(community.temp$Site, community.temp$rep, sep = "")
    community.temp <- subset(community.temp, select = - c(s.id, year, seasn, Site, rep))
    community.temp[community.temp > 0] <- 1
    community.temp <- community.temp[rowSums(community.temp) > 0, colSums(community.temp) > 0]
    MDS1 <- metaMDS(community.temp, wascores = FALSE)
    if (j == "J"){
      month <- "January"
    }
    else{
      month <- "September"
    }
    plot(MDS1, type = "t", main = paste(month, i, sep = " "))
  }
}













###############





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
windows(title="Silhouettes - Ward - k = 2 to n-1")
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
windows(title="Optimal number of clusters - Mantel")
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
windows(title="Silhouette plot - Average - K=4")
plot(silo, main="Silhouette plot - Chord - Average", 
     cex.names=0.8, col=cutg+1, nmax.lab=100)



# Final dendrogram with the selected groups
# *****************************************
library(gclus)
spe.chwo <- reorder.hclust(td.clust.aver, td)
summary(spe.chwo)


# Plot reordered dendrogram with group labels
windows(title="Final dendrogram",8,6)
plot(spe.chwo, hang=-1, xlab="5 groups", sub="", 
     ylab="Height", 
     main="Chord - Average (reordered)", 
     labels=alts.s.id$altitude
     #labels=cutree(spe.chwo, k=k)
     )
rect.hclust(spe.chwo, k=k)

# Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
source("hcoplot.R")	       # hcoplot.R must be in the working directory
hcoplot(spe.ch.ward, spe.ch, k=4)


# Heat map
# ********

# Heat map of the distance matrix ordered with the dendrogram
dend <- as.dendrogram(spe.chwo)
windows(title="Heatmap - sites")
heatmap(as.matrix(td), Rowv=dend, symm=TRUE, margin=c(3,3))

# Ordered community table
# Species are ordered by their weighted averages on site scores
or <- vegemite(spe, spe.chwo, scale = "log")

# Heat map of the doubly ordered community table, with dendrogram
library(RColorBrewer)
windows(title="Heatmap - species")
heatmap(t(spe[rev(or$species)]), Rowv=NA, Colv=dend,
        col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4), 
        ylab="Species (weighted averages of sites)", xlab="Sites")




library(indicspecies)

data(wetland) ## Loads species data
## Creates three clusters using kmeans
wetkm = kmeans(wetland, centers=3)
## Determine sensitivity of individual species
B=strassoc(wetland, cluster=wetkm$cluster,func="B")
## Select species with more than 20% of sensitivity for the first group
sel=which(B[,1]>0.2)
## Run indicator analysis with species combinations for the first group
sc= indicators(X=wetland[,sel], cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2)
#Prints the results
print(sc)
## Plots positive predictive power and sensitivity against the order of combinations
plot(sc, type="A")
plot(sc, type="B")

## Run indicator analysis with species combinations for the first group,
## but forcing ’Orysp’ to be in all combinations
sc2= indicators(X=wetland[,sel], cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2, enableFixed=TRUE)



sc= indicators(X=wetland[,sel], cluster=wetkm$cluster, group = 2, verbose=TRUE, At=0.5, Bt=0.1)
sc


ind.group1 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 1, 
                     max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7)
ind.group2 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 2, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7)
ind.group3 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 3, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7)
ind.group4 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 4, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7)
ind.group5 <- indicators(X = spe, cluster = cutree(spe.chwo, k=k), group = 5, 
                         max.order = 3, verbose = TRUE, At = 0.7, Bt = 0.7)

ind.group1
ind.group2
ind.group3
ind.group4
ind.group5




trait[order(trait$Rel.Eye.Width),
      c("Genus", "Species", "Head.Width", "Eye.Width", "Rel.Eye.Width", "Interocular", "Eye.Position", "Rel.Eye.Position") ]