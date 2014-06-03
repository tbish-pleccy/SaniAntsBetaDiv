####################################
# Taxonomic Spatial Beta Diversity #
####################################
library(lattice)
# 1. 
# Does beta diversity correlate with elevation?
# 
# Yes in all time slices for Bsor and Bsim. 
# Bnes shows mixed reuslts. Some postive in Wet, some negative in Dry. 
# 

par(mfrow = c(2, 4))
beta.plot(2006, "J", "bsor", "taxonomic")
beta.plot(2007, "J", "bsor", "taxonomic")
beta.plot(2008, "J", "bsor", "taxonomic")
beta.plot(2009, "J", "bsor", "taxonomic")
beta.plot(2010, "J", "bsor", "taxonomic")
beta.plot(2011, "J", "bsor", "taxonomic")
beta.plot(2012, "J", "bsor", "taxonomic")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bsor", "taxonomic")
beta.plot(2007, "S", "bsor", "taxonomic")
beta.plot(2008, "S", "bsor", "taxonomic")
beta.plot(2009, "S", "bsor", "taxonomic")
beta.plot(2010, "S", "bsor", "taxonomic")
beta.plot(2011, "S", "bsor", "taxonomic")
beta.plot(2012, "S", "bsor", "taxonomic")

par(mfrow = c(2, 4))
beta.plot(2006, "J", "bsim", "taxonomic")
beta.plot(2007, "J", "bsim", "taxonomic")
beta.plot(2008, "J", "bsim", "taxonomic")
beta.plot(2009, "J", "bsim", "taxonomic")
beta.plot(2010, "J", "bsim", "taxonomic")
beta.plot(2011, "J", "bsim", "taxonomic")
beta.plot(2012, "J", "bsim", "taxonomic")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bsim", "taxonomic")
beta.plot(2007, "S", "bsim", "taxonomic")
beta.plot(2008, "S", "bsim", "taxonomic")
beta.plot(2009, "S", "bsim", "taxonomic")
beta.plot(2010, "S", "bsim", "taxonomic")
beta.plot(2011, "S", "bsim", "taxonomic")
beta.plot(2012, "S", "bsim", "taxonomic")

par(mfrow = c(2, 4))
beta.plot(2006, "J", "bnes", "taxonomic")
beta.plot(2007, "J", "bnes", "taxonomic")
beta.plot(2008, "J", "bnes", "taxonomic")
beta.plot(2009, "J", "bnes", "taxonomic")
beta.plot(2010, "J", "bnes", "taxonomic")
beta.plot(2011, "J", "bnes", "taxonomic")
beta.plot(2012, "J", "bnes", "taxonomic")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bnes", "taxonomic")
beta.plot(2007, "S", "bnes", "taxonomic")
beta.plot(2008, "S", "bnes", "taxonomic")
beta.plot(2009, "S", "bnes", "taxonomic")
beta.plot(2010, "S", "bnes", "taxonomic")
beta.plot(2011, "S", "bnes", "taxonomic")
beta.plot(2012, "S", "bnes", "taxonomic")


# 2. 
# Are these patterns consistent in their form across time? 
#
# Bsor slope
#

complete <- NULL
for (i in unique(site.info$period)){
  for (j in c("tax", "fun")){
    for (k in c("bsor", "bsim", "bnes")){
      if (j == "tax"){
        Y.dissim <- comm.dissim
      }
      else{
        Y.dissim <- fun.dissim
      }
      
      if (j == "fun" & k == "bnes"){
        metric <- "bsne"
      }
      else{
        metric <- k
      }
      
      SI.temp <- site.info[site.info$period == i, ]
      period <- i
      year <- unique(SI.temp$year)
      seasn <- unique(SI.temp$seasn)
      compo <- j
      tax.diss <- Y.dissim[[paste(year, seasn, metric, sep = "")]] 
      ele.diss <- ele.dissim[[paste(year, seasn, sep = "")]]
      slope <- summary(lm(tax.diss ~ ele.diss))[[4]][2]
      intercept <- summary(lm(tax.diss ~ ele.diss))[[4]][1]
      frame.temp <- data.frame(compo, year, seasn, period, slope, intercept, metric)
      complete <- rbind(complete, frame.temp)      
    }  
  }
}


xyplot(slope ~ year | seasn + metric, data = complete[complete$compo == "tax", ])
xyplot(intercept ~ year | seasn + metric, data = complete[complete$compo == "tax", ])

bsor.slope.1 <- glm(slope ~ year * seasn, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])
bsor.slope.2 <- glm(slope ~ year + seasn, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])
bsor.slope.3 <- glm(slope ~ year, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])
bsor.slope.4 <- glm(slope ~ seasn, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])
bsor.slope.5 <- glm(slope ~ 1, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])

library(AICcmodavg)
bsor.slope.models <- list()
for (i in 1:5){
  get(paste("bsor.slope.", i, sep = ""))
  bsor.slope.models[[i]] <- get(paste("bsor.slope.", i, sep = ""))
}

bsor.slope.names <- paste("bsor.slope.", 1:length(bsor.slope.models), sep = "")
aictable.bsor.slope <- aictab(cand.set = bsor.slope.models, modnames = bsor.slope.names, sort = TRUE)
aictable.bsor.slope 
best.model.bsor.slope <- get(as.vector((aictable.bsor.slope[1, "Modnames"])))
summary(best.model.bsor.slope)

#
# Bsor intercept
#

bsor.intercept.1 <- glm(intercept ~ year * seasn, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])
bsor.intercept.2 <- glm(intercept ~ year + seasn, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])
bsor.intercept.3 <- glm(intercept ~ year, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])
bsor.intercept.4 <- glm(intercept ~ seasn, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])
bsor.intercept.5 <- glm(intercept ~ 1, data = complete[complete$compo == "tax" & complete$metric == "bsor", ])

library(AICcmodavg)
bsor.intercept.models <- list()
for (i in 1:5){
  get(paste("bsor.intercept.", i, sep = ""))
  bsor.intercept.models[[i]] <- get(paste("bsor.intercept.", i, sep = ""))
}

bsor.intercept.names <- paste("bsor.intercept.", 1:length(bsor.intercept.models), sep = "")
aictable.bsor.intercept <- aictab(cand.set = bsor.intercept.models, modnames = bsor.intercept.names, sort = TRUE)
aictable.bsor.intercept 
best.model.bsor.intercept <- get(as.vector((aictable.bsor.intercept[1, "Modnames"])))
summary(best.model.bsor.intercept)

#
# Bsim slope
# Bsim slopes do not vary

xyplot(slope ~ year | seasn + metric, data = complete[complete$compo == "tax", ])
xyplot(intercept ~ year | seasn + metric, data = complete[complete$compo == "tax", ])

bsim.slope.1 <- glm(slope ~ year * seasn, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])
bsim.slope.2 <- glm(slope ~ year + seasn, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])
aov(bsim.slope.1, bsim.slope.2)
bsim.slope.3 <- glm(slope ~ year, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])
bsim.slope.4 <- glm(slope ~ seasn, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])
bsim.slope.5 <- glm(slope ~ 1, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])

summary(bsim.slope.1)
summary(bsim.slope.2)
summary(bsim.slope.3)
summary(bsim.slope.4)

library(AICcmodavg)
bsim.slope.models <- list()
for (i in 1:5){
  get(paste("bsim.slope.", i, sep = ""))
  bsim.slope.models[[i]] <- get(paste("bsim.slope.", i, sep = ""))
}

bsim.slope.names <- paste("bsim.slope.", 1:length(bsim.slope.models), sep = "")
aictable.bsim.slope <- aictab(cand.set = bsim.slope.models, modnames = bsim.slope.names, sort = TRUE)
aictable.bsim.slope 
best.model.bsim.slope <- get(as.vector((aictable.bsim.slope[1, "Modnames"])))
summary(best.model.bsim.slope)

#
# Bsim intercept
# Bsim intercepts are higher in September and also increase through the years.

bsim.intercept.1 <- glm(intercept ~ year * seasn, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])
bsim.intercept.2 <- glm(intercept ~ year + seasn, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])
bsim.intercept.3 <- glm(intercept ~ year, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])
bsim.intercept.4 <- glm(intercept ~ seasn, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])
bsim.intercept.5 <- glm(intercept ~ 1, data = complete[complete$compo == "tax" & complete$metric == "bsim", ])

library(AICcmodavg)
bsim.intercept.models <- list()
for (i in 1:5){
  get(paste("bsim.intercept.", i, sep = ""))
  bsim.intercept.models[[i]] <- get(paste("bsim.intercept.", i, sep = ""))
}

bsim.intercept.names <- paste("bsim.intercept.", 1:length(bsim.intercept.models), sep = "")
aictable.bsim.intercept <- aictab(cand.set = bsim.intercept.models, modnames = bsim.intercept.names, sort = TRUE)
aictable.bsim.intercept 
best.model.bsim.intercept <- get(as.vector((aictable.bsim.intercept[1, "Modnames"])))
summary(best.model.bsim.intercept)

#
# Bnes slope
# Bnes slopes vary. Lower slopes in September and through years. Interaction also matters. 

xyplot(slope ~ year | seasn + metric, data = complete[complete$compo == "tax", ])
xyplot(intercept ~ year | seasn + metric, data = complete[complete$compo == "tax", ])

bnes.slope.1 <- glm(slope ~ year * seasn, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])
bnes.slope.2 <- glm(slope ~ year + seasn, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])
bnes.slope.3 <- glm(slope ~ year, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])
bnes.slope.4 <- glm(slope ~ seasn, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])
bnes.slope.5 <- glm(slope ~ 1, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])

summary(bnes.slope.1)
summary(bnes.slope.2)
summary(bnes.slope.3)
summary(bnes.slope.4)

library(AICcmodavg)
bnes.slope.models <- list()
for (i in 1:5){
  get(paste("bnes.slope.", i, sep = ""))
  bnes.slope.models[[i]] <- get(paste("bnes.slope.", i, sep = ""))
}

bnes.slope.names <- paste("bnes.slope.", 1:length(bnes.slope.models), sep = "")
aictable.bnes.slope <- aictab(cand.set = bnes.slope.models, modnames = bnes.slope.names, sort = TRUE)
aictable.bnes.slope 
best.model.bnes.slope <- get(as.vector((aictable.bnes.slope[1, "Modnames"])))
summary(best.model.bnes.slope)

#
# Bnes intercept
# Bnes intercepts do not vary

bnes.intercept.1 <- glm(intercept ~ year * seasn, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])
bnes.intercept.2 <- glm(intercept ~ year + seasn, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])
bnes.intercept.3 <- glm(intercept ~ year, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])
bnes.intercept.4 <- glm(intercept ~ seasn, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])
bnes.intercept.5 <- glm(intercept ~ 1, data = complete[complete$compo == "tax" & complete$metric == "bnes", ])

library(AICcmodavg)
bnes.intercept.models <- list()
for (i in 1:5){
  get(paste("bnes.intercept.", i, sep = ""))
  bnes.intercept.models[[i]] <- get(paste("bnes.intercept.", i, sep = ""))
}

bnes.intercept.names <- paste("bnes.intercept.", 1:length(bnes.intercept.models), sep = "")
aictable.bnes.intercept <- aictab(cand.set = bnes.intercept.models, modnames = bnes.intercept.names, sort = TRUE)
aictable.bnes.intercept 
best.model.bnes.intercept <- get(as.vector((aictable.bnes.intercept[1, "Modnames"])))
summary(best.model.bnes.intercept)

# 3. 
# Does turnover or nestedness dominate?
# 
# Taxonomic turnover. 
# Both the average pairwise metrics and average Mantel strength are greater for
# turnover. 
#

tax.mantels <- NULL
for (i in unique(site.info$year)){
  for (j in unique(site.info$seasn)){
    for (k in c("bsor", "bsim", "bnes")){
     c <- comm.dissim[[paste(i, j, k, sep = "")]] 
     e <- ele.dissim[[paste(i, j, sep = "")]]
     tax.mantels.t <- data.frame(t(mantel(c ~ e, nperm = 1000, nboot = 1000)))
     tax.mantels.t$year <- i
     tax.mantels.t$seasn <- j
     tax.mantels.t$metric <- k
     tax.mantels.t$av.pairwise <- mean(c)
     tax.mantels <- rbind(tax.mantels, tax.mantels.t)
    }  
  }
}
tax.mantels$metric <- as.factor(tax.mantels$metric)
tax.mantels$seasn <- as.factor(tax.mantels$seasn)

# Average Mantel strength
tapply(tax.mantels$mantelr, list(tax.mantels$seasn, tax.mantels$metric), mean)
xyplot(mantelr ~ metric | seasn, data = tax.mantels)

# Average pairwise matric
tapply(tax.mantels$av.pairwise, list(tax.mantels$seasn, tax.mantels$metric), mean)
xyplot(av.pairwise ~ metric | seasn, data = tax.mantels)
