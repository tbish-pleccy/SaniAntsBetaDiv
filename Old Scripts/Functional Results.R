####################################
# Functional Spatial Beta Diversity #
####################################
library(lattice)
# 1. 
# Does beta diversity correlate with elevation?
# 
# Yes in all time slices for Bsor. Positive. 
# Bsim is 4/7 in Jan and 5/7 in Sep. Positive. 
# Bsne is 5/7 positive and 2/7 negative in Jan. No significant results in Sep. 
#

par(mfrow = c(2, 4))
beta.plot(2006, "J", "bsor", "functional")
beta.plot(2007, "J", "bsor", "functional")
beta.plot(2008, "J", "bsor", "functional")
beta.plot(2009, "J", "bsor", "functional")
beta.plot(2010, "J", "bsor", "functional")
beta.plot(2011, "J", "bsor", "functional")
beta.plot(2012, "J", "bsor", "functional")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bsor", "functional")
beta.plot(2007, "S", "bsor", "functional")
beta.plot(2008, "S", "bsor", "functional")
beta.plot(2009, "S", "bsor", "functional")
beta.plot(2010, "S", "bsor", "functional")
beta.plot(2011, "S", "bsor", "functional")
beta.plot(2012, "S", "bsor", "functional")

par(mfrow = c(2, 4))
beta.plot(2006, "J", "bsim", "functional")
beta.plot(2007, "J", "bsim", "functional")
beta.plot(2008, "J", "bsim", "functional")
beta.plot(2009, "J", "bsim", "functional")
beta.plot(2010, "J", "bsim", "functional")
beta.plot(2011, "J", "bsim", "functional")
beta.plot(2012, "J", "bsim", "functional")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bsim", "functional")
beta.plot(2007, "S", "bsim", "functional")
beta.plot(2008, "S", "bsim", "functional")
beta.plot(2009, "S", "bsim", "functional")
beta.plot(2010, "S", "bsim", "functional")
beta.plot(2011, "S", "bsim", "functional")
beta.plot(2012, "S", "bsim", "functional")

par(mfrow = c(2, 4))
beta.plot(2006, "J", "bsne", "functional")
beta.plot(2007, "J", "bsne", "functional")
beta.plot(2008, "J", "bsne", "functional")
beta.plot(2009, "J", "bsne", "functional")
beta.plot(2010, "J", "bsne", "functional")
beta.plot(2011, "J", "bsne", "functional")
beta.plot(2012, "J", "bsne", "functional")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bsne", "functional")
beta.plot(2007, "S", "bsne", "functional")
beta.plot(2008, "S", "bsne", "functional")
beta.plot(2009, "S", "bsne", "functional")
beta.plot(2010, "S", "bsne", "functional")
beta.plot(2011, "S", "bsne", "functional")
beta.plot(2012, "S", "bsne", "functional")


# 2. 
# Are these patterns consistent in their form across time? 
#
# Bsor slope
# Decreases with year


xyplot(slope ~ year | seasn + metric, data = complete[complete$compo == "fun", ])
xyplot(intercept ~ year | seasn + metric, data = complete[complete$compo == "fun", ])

bsor.slope.F1 <- glm(slope ~ year * seasn, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])
bsor.slope.F2 <- glm(slope ~ year + seasn, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])
bsor.slope.F3 <- glm(slope ~ year, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])
bsor.slope.F4 <- glm(slope ~ seasn, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])
bsor.slope.F5 <- glm(slope ~ 1, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])

library(AICcmodavg)
bsor.slope.Fmodels <- list()
for (i in 1:5){
  get(paste("bsor.slope.F", i, sep = ""))
  bsor.slope.Fmodels[[i]] <- get(paste("bsor.slope.F", i, sep = ""))
}

bsor.slope.Fnames <- paste("bsor.slope.F", 1:length(bsor.slope.Fmodels), sep = "")
aictable.bsor.slopeF <- aictab(cand.set = bsor.slope.Fmodels, modnames = bsor.slope.Fnames, sort = TRUE)
aictable.bsor.slopeF 
best.model.bsor.slopeF <- get(as.vector((aictable.bsor.slopeF[1, "Modnames"])))
summary(best.model.bsor.slopeF)

#
# Bsor intercept
# Changes with year and season

bsor.intercept.F1 <- glm(intercept ~ year * seasn, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])
bsor.intercept.F2 <- glm(intercept ~ year + seasn, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])
bsor.intercept.F3 <- glm(intercept ~ year, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])
bsor.intercept.F4 <- glm(intercept ~ seasn, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])
bsor.intercept.F5 <- glm(intercept ~ 1, data = complete[complete$compo == "fun" & complete$metric == "bsor", ])

library(AICcmodavg)
bsor.intercept.Fmodels <- list()
for (i in 1:5){
  get(paste("bsor.intercept.F", i, sep = ""))
  bsor.intercept.Fmodels[[i]] <- get(paste("bsor.intercept.F", i, sep = ""))
}

bsor.intercept.Fnames <- paste("bsor.intercept.F", 1:length(bsor.intercept.Fmodels), sep = "")
aictable.bsor.interceptF <- aictab(cand.set = bsor.intercept.Fmodels, modnames = bsor.intercept.Fnames, sort = TRUE)
aictable.bsor.interceptF 
best.model.bsor.interceptF <- get(as.vector((aictable.bsor.interceptF[1, "Modnames"])))
summary(best.model.bsor.interceptF)

#
# Bsim slope
# Nothing

bsim.slope.F1 <- glm(slope ~ year * seasn, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])
bsim.slope.F2 <- glm(slope ~ year + seasn, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])
bsim.slope.F3 <- glm(slope ~ year, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])
bsim.slope.F4 <- glm(slope ~ seasn, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])
bsim.slope.F5 <- glm(slope ~ 1, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])

library(AICcmodavg)
bsim.slope.Fmodels <- list()
for (i in 1:5){
  get(paste("bsim.slope.F", i, sep = ""))
  bsim.slope.Fmodels[[i]] <- get(paste("bsim.slope.F", i, sep = ""))
}

bsim.slope.Fnames <- paste("bsim.slope.F", 1:length(bsim.slope.Fmodels), sep = "")
aictable.bsim.slopeF <- aictab(cand.set = bsim.slope.Fmodels, modnames = bsim.slope.Fnames, sort = TRUE)
aictable.bsim.slopeF 
best.model.bsim.slopeF <- get(as.vector((aictable.bsim.slopeF[1, "Modnames"])))
summary(best.model.bsim.slopeF)

#
# Bsim intercept
# season

bsim.intercept.F1 <- glm(intercept ~ year * seasn, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])
bsim.intercept.F2 <- glm(intercept ~ year + seasn, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])
bsim.intercept.F3 <- glm(intercept ~ year, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])
bsim.intercept.F4 <- glm(intercept ~ seasn, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])
bsim.intercept.F5 <- glm(intercept ~ 1, data = complete[complete$compo == "fun" & complete$metric == "bsim", ])

library(AICcmodavg)
bsim.intercept.Fmodels <- list()
for (i in 1:5){
  get(paste("bsim.intercept.F", i, sep = ""))
  bsim.intercept.Fmodels[[i]] <- get(paste("bsim.intercept.F", i, sep = ""))
}

bsim.intercept.Fnames <- paste("bsim.intercept.F", 1:length(bsim.intercept.Fmodels), sep = "")
aictable.bsim.interceptF <- aictab(cand.set = bsim.intercept.Fmodels, modnames = bsim.intercept.Fnames, sort = TRUE)
aictable.bsim.interceptF 
best.model.bsim.interceptF <- get(as.vector((aictable.bsim.interceptF[1, "Modnames"])))
summary(best.model.bsim.interceptF)

#
# Bsne slope
# year, season and interaction

bsne.slope.F1 <- glm(slope ~ year * seasn, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])
bsne.slope.F2 <- glm(slope ~ year + seasn, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])
bsne.slope.F3 <- glm(slope ~ year, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])
bsne.slope.F4 <- glm(slope ~ seasn, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])
bsne.slope.F5 <- glm(slope ~ 1, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])

library(AICcmodavg)
bsne.slope.Fmodels <- list()
for (i in 1:5){
  get(paste("bsne.slope.F", i, sep = ""))
  bsne.slope.Fmodels[[i]] <- get(paste("bsne.slope.F", i, sep = ""))
}

bsne.slope.Fnames <- paste("bsne.slope.F", 1:length(bsne.slope.Fmodels), sep = "")
aictable.bsne.slopeF <- aictab(cand.set = bsne.slope.Fmodels, modnames = bsne.slope.Fnames, sort = TRUE)
aictable.bsne.slopeF 
best.model.bsne.slopeF <- get(as.vector((aictable.bsne.slopeF[1, "Modnames"])))
summary(best.model.bsne.slopeF)

#
# Bsne intercept
# Nothing. 

bsne.intercept.F1 <- glm(intercept ~ year * seasn, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])
bsne.intercept.F2 <- glm(intercept ~ year + seasn, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])
bsne.intercept.F3 <- glm(intercept ~ year, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])
bsne.intercept.F4 <- glm(intercept ~ seasn, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])
bsne.intercept.F5 <- glm(intercept ~ 1, data = complete[complete$compo == "fun" & complete$metric == "bsne", ])

library(AICcmodavg)
bsne.intercept.Fmodels <- list()
for (i in 1:5){
  get(paste("bsne.intercept.F", i, sep = ""))
  bsne.intercept.Fmodels[[i]] <- get(paste("bsne.intercept.F", i, sep = ""))
}

bsne.intercept.Fnames <- paste("bsne.intercept.F", 1:length(bsne.intercept.Fmodels), sep = "")
aictable.bsne.interceptF <- aictab(cand.set = bsne.intercept.Fmodels, modnames = bsne.intercept.Fnames, sort = TRUE)
aictable.bsne.interceptF 
best.model.bsne.interceptF <- get(as.vector((aictable.bsne.interceptF[1, "Modnames"])))
summary(best.model.bsne.interceptF)

# 3. 
# Does turnover or nestedness dominate?
# 
# Mantel strength nearly equal in January. Turnover dominates in September.  
# Average pairwise nearly equal in January. Turnover is double in September. 

fun.mantels <- NULL
for (i in unique(site.info$year)){
  for (j in unique(site.info$seasn)){
    for (k in c("bsor", "bsim", "bsne")){
      c <- fun.dissim[[paste(i, j, k, sep = "")]] 
      e <- ele.dissim[[paste(i, j, sep = "")]]
      fun.mantels.t <- data.frame(t(mantel(c ~ e, nperm = 1000, nboot = 1000)))
      fun.mantels.t$year <- i
      fun.mantels.t$seasn <- j
      fun.mantels.t$metric <- k
      fun.mantels.t$av.pairwise <- mean(c)
      fun.mantels <- rbind(fun.mantels, fun.mantels.t)
    }  
  }
}
fun.mantels$metric <- as.factor(fun.mantels$metric)
fun.mantels$seasn <- as.factor(fun.mantels$seasn)

# Average Mantel strength
tapply(fun.mantels$mantelr, list(fun.mantels$seasn, fun.mantels$metric), mean)
xyplot(mantelr ~ metric | seasn, data = fun.mantels)

# Average pairwise matric
tapply(fun.mantels$av.pairwise, list(fun.mantels$seasn, fun.mantels$metric), mean)
xyplot(av.pairwise ~ metric | seasn, data = fun.mantels)
