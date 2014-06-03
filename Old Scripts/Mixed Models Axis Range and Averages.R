library(lattice)

xyplot(rangea1 ~ altitude | as.factor(year) + seasn, data = comm.ascores)
xyplot(rangea2 ~ altitude | as.factor(year) + seasn, data = comm.ascores)
xyplot(rangea3 ~ altitude | as.factor(year) + seasn, data = comm.ascores)
xyplot(rangea4 ~ altitude | as.factor(year) + seasn, data = comm.ascores)
xyplot(a1.av ~ altitude | as.factor(year) + seasn, data = comm.ascores)
xyplot(a2.av ~ altitude | as.factor(year) + seasn, data = comm.ascores)
xyplot(a3.av ~ altitude | as.factor(year) + seasn, data = comm.ascores)
xyplot(a4.av ~ altitude | as.factor(year) + seasn, data = comm.ascores)

head(comm.ascores)

fundiv2 <- dbFD(x = tFD.temp, a = community)
fundiv2 <- data.frame(fundiv2)
head(fundiv2)
fundiv2$s.id <- rownames(community)
fundiv2 <- merge(fundiv2, site.info, by = "s.id", all.x = FALSE)
head(fundiv2)

xyplot(FRic ~ altitude | year + seasn, data = fundiv2)
xyplot(CWM.Webers.Length ~ altitude | year + seasn, data = fundiv2)
xyplot(CWM.Rel.Eye.Position ~ altitude | year + seasn, data = fundiv2)
xyplot(CWM.Mand.Index ~ altitude | year + seasn, data = fundiv2)
xyplot(CWM.Rel.Leg.Length ~ altitude | year + seasn, data = fundiv2)



fun.dissim[["2012Sbsne"]]
range.dist <- dist(comm.ascores[comm.ascores$s.id %in% rownames(as.matrix(fun.dissim[["2012Sbsne"]]))
, "a1.av"])


plot(range.dist, fun.dissim[["2012Sbsne"]])

#######
# Mixed Modelling
#######

##################################################
comm.ascores$altitude <- as.numeric(comm.ascores$altitude)
comm.ascores$alti2 <- comm.ascores$altitude^2
comm.ascores$fyear <- as.factor(comm.ascores$year)
comm.ascores$fsite <- as.factor(comm.ascores$Site)
comm.ascores$altiC <- scale(comm.ascores$altitude)
comm.ascores$alti2C <- scale(comm.ascores$alti2)
comm.ascores$seasnB <- comm.ascores$seasn
levels(comm.ascores$seasnB) <- c(0, 1)
comm.ascores$seasnB <- as.numeric(levels(comm.ascores$seasnB))[comm.ascores$seasnB]  

library(lme4)
a1ran.full <- lmer(rangea1 ~ altiC*seasnB + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
summary(a1ran.full)

# test two ways

a1ran.1.1 <- lmer(rangea1 ~ altiC + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
a1ran.1.2 <- lmer(rangea1 ~ altiC*seasnB + alti2C + (1 | fyear / seasn / fsite), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
anova(a1ran.full, a1ran.1.1) # altiC*seasnB = 0.58 remove
anova(a1ran.full, a1ran.1.2) # alti2C*seasnB = 0.3

# 
a1ran.2 <- lmer(rangea1 ~ altiC + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)

a1ran.2.1 <- lmer(rangea1 ~ altiC + alti2C + seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)
anova(a1ran.2, a1ran.2.1) # alti2C*seasnB = 0.0022 KEEP!

# test altitude
a1ran.3 <- lmer(rangea1 ~ altiC + alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)
a1ran.3.1 <- lmer(rangea1 ~ alti2C*seasnB + (1 | fyear / seasn / fsite), data = comm.ascores, 
                control = list(maxIter = 5200), REML = FALSE)
anova(a1ran.3, a1ran.3.1) # altiC = 0.0000 KEEP!




summary(a1ran.3)


a1ran.final <- lmer(rangea1 ~ altiC + alti2C*seasnB + (1 | rep), data = comm.ascores, 
                    control = list(maxIter = 5200), REML = FALSE)
summary(a1ran.final)


a1ran.final.plot <- lmer(rangea1 ~ altiC + alti2C*seasnB + (1 | seasn), data = comm.ascores, 
                    control = list(maxIter = 5200), REML = FALSE)
summary(a1ran.final.plot)








library(MuMIn)
r.squaredGLMM(x = a1ran.final)


plot(fitted(a1ran.final), resid(a1ran.final))


comm.ascores$fit.a1ran <- fitted(a1ran.final.plot)

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

###############################################################################
# a1 average

a1av.full <- lmer(a1.av ~ altiC*seasnB + alti2C*seasnB + (1 | rep), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
summary(a1av.full)

# test two ways

a1av.1.1 <- lmer(a1.av ~ altiC + alti2C*seasnB + (1 | rep), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
a1av.1.2 <- lmer(a1.av ~ altiC*seasnB + alti2C + (1 | rep), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)
anova(a1av.full, a1av.1.1) # altiC*seasnB = 0.02
anova(a1av.full, a1av.1.2) # alti2C*seasnB = 0.3


a1av.2 <- lmer(a1.av ~ altiC*seasnB + alti2C + (1 | rep), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a1av.2.1 <- lmer(a1.av ~ altiC + seasnB + alti2C + (1 | rep), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a1av.2, a1av.2.1) # altiC*seasnB = 0.64 REMOVE


a1av.3 <- lmer(a1.av ~ altiC + seasnB + alti2C + (1 | rep), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a1av.3.1 <- lmer(a1.av ~ seasnB + alti2C + (1 | rep), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a1av.3.2 <- lmer(a1.av ~ altiC +  alti2C + (1 | rep), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a1av.3.3 <- lmer(a1.av ~ altiC + seasnB + (1 | rep), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a1av.3, a1av.3.1) # altiC = 0.05
anova(a1av.3, a1av.3.2) # seasnB = 0.000
anova(a1av.3, a1av.3.3) # alti2C = 0.058 REMOVE


a1av.4 <- lmer(a1.av ~ altiC + seasnB + (1 | rep), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a1av.4.1 <- lmer(a1.av ~ altiC + (1 | rep), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
a1av.4.2 <- lmer(a1.av ~ seasnB + (1 | rep), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a1av.4, a1av.4.1) # seasnB = 0.000
anova(a1av.4, a1av.4.2) # altiC = 0.58 REMOVE

a1av.5 <- lmer(a1.av ~ seasnB + (1 | rep), data = comm.ascores, 
                 control = list(maxIter = 5200), REML = FALSE)
a1av.5.1 <- lmer(a1.av ~ (1 | rep), data = comm.ascores, 
               control = list(maxIter = 5200), REML = FALSE)
anova(a1av.5, a1av.5.1) # seasnB = 0.0000

summary(a1av.5)



a1av.null <- lmer(a1.av ~ (1 | rep), data = comm.ascores, 
                  control = list(maxIter = 5200), REML = FALSE)

anova(a1av.full, a1av.null)

# FINAL
a1av.final <- lmer(a1.av ~ altiC*seasnB + alti2C*seasnB + (1 | rep), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
summary(a1av.final)


a1av.final2.1 <- lmer(a1.av ~ altiC*seasnB + alti2C + seasnB + (1 | rep), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
a1av.final2.2 <- lmer(a1.av ~ altiC + seasnB + alti2C*seasnB + (1 | rep), data = comm.ascores, 
                   control = list(maxIter = 5200), REML = FALSE)
anova(a1av.final, a1av.final2.1)
anova(a1av.final, a1av.final2.2)


a1av.final.plot <- lmer(a1.av ~ seasnB + (1 | rep), data = comm.ascores, 
                                  control = list(maxIter = 5200), REML = FALSE)


library(MuMIn)
r.squaredGLMM(x = a1av.final)

plot(fitted(a1av.final), resid(a1av.final))


comm.ascores$fit.a1av <- fitted(a1av.final.plot)


par(mfrow = c(3, 3), mar = c(5, 5, 4, 2))
for (i in unique(comm.ascores$year)){
  raw <- comm.ascores[comm.ascores$year == i, ]
  plot(fit.a1av ~ seasn, data = raw, pch = 21,
       ylab = "Average A1 Score", main = i)
  
  #raw2 <- comm.ascores[comm.ascores$seasn == "S" & comm.ascores$year == i, ]
  #points(a1.av ~ seasn, data = raw2, pch = 21, col = "white", bg = "lightblue")
  #lines(fit.a1av ~ seasn, data = raw, col = "red")
  #lines(fit.a1av ~ seasn, data = raw2, col = "blue")  
} 