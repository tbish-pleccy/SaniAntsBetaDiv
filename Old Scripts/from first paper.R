########################################
# 0. Set working directory.            #
########################################

#wd <- "C:/Users/Tom Bishop/Documents/Research Projects/PhD/1 First Paper/Paper 1 Restart November 2013/Data"
#setwd(wd)

########################################
# 1. Import and format abundance data. #
########################################
nrow(dd)
dd <- read.csv("Data/Data_for_R_JulyUpdate.csv")

# code to remove duplicates from incorrect data supplied on 27/03/2013
#dd <- (dd[!duplicated(dd), ])  # remove pure duplicates. goes form 6166 to 4694.
#dd$spid <- as.factor(dd$spid)
#campdupes <- as.vector(sort(unique(dd$species)))[c(7:11, 13:17)]
#unique(dd[dd$species %in% campdupes, 7:8])

# levels(dd$spid) <- sub("^15$", "26", levels(dd$spid))
# levels(dd$spid) <- sub("^19$", "30", levels(dd$spid))
# levels(dd$spid) <- sub("^20$", "31", levels(dd$spid))
# levels(dd$spid) <- sub("^16$", "27", levels(dd$spid))
# levels(dd$spid) <- sub("^18$", "29", levels(dd$spid))
# levels(dd$species) <- sub("^Camponatus sp. 1$", "Camponotus sp. 1", levels(dd$species))
# levels(dd$species) <- sub("^Camponatus un01$", "Camponotus un01", levels(dd$species))
# levels(dd$species) <- sub("^Camponatus un10$", "Camponotus un10", levels(dd$species))
# levels(dd$species) <- sub("^Camponatus sp. 3$", "Camponotus sp. 3", levels(dd$species))
# levels(dd$species) <- sub("^Camponatus sp. 8$", "Camponotus sp. 8", levels(dd$species))
# 
# dd <- (dd[!duplicated(dd), ])
# 
# torem <- as.vector(sort(unique(dd$species))[c(70, 79, 102)])
# dd <- dd[!dd$species %in% torem,]  # remove duplicates due to different spelling
# dd$species <- droplevels(dd$species)
# dd$spid <- droplevels(dd$spid)
# temp.remov <- as.vector(sort(unique(dd$species))[c(2, 45, 41, 62, 63, 89)])
# # temporarily remove: aenictus, monomorium 9&10, plagiolepis un02&un03, tetramorium 41
# dd <- dd[!dd$species %in% tempremov,]

# dd[dd$species == "Aenictus sp. 1",]
# dd[dd$species == "Monomorium sp. 10",]
# n<-as.vector(levels(dd$species))
# n2<-c(1:length(n))
# ndf<-data.frame(cbind(n,n2))

# Start of Mark's code.
# Create dataset of zeros (where nothing was collected) and then combine with 
# existing data.

us <- (dd[,7:8])  # unique species identifier and species name
nsp <- nrow(us)  # number of species
s <- sort(rep(1:8, 4))  # site numbers
rp <- rep(letters[1:4], 8)  # letters for reps
rep1 <- paste(s, rp, sep = "")  # replicate names e.g. 1a, 1b, 1c
uoc <- unique(dd[, c(4, 6)])  # unique year and occasions (e.g. J06, S06)
nuoc <- nrow(uoc)  # number of occasions
dat1 <- data.frame()  # create a dataframe
for (j in 1:nuoc){
  yr <- uoc$year[j]  # a particular year
  occ <- uoc$occasion[j]  # a particular occasion
  d1 <- subset(dd, occasion == occ)  # actual data for this occasion
  dff <- data.frame()
  for (i in 1:32){  # for all replicates from 1a to 8d
    site <- s[i]  # site number 
    scode <- paste("s", s[i], sep = "")  # site code e.g. s1
    dy <- data.frame(us, indiv = 0)  # all species with individuals as zero
    # (spid, species, indiv = 0)
    d2 <- subset(d1, rep == rep1[i])  # select data for this 
    # particular replicate
    d3 <- rbind(d2[, 7:9], dy)  # combine actual data with zero individuals data
    d4 <- aggregate(indiv ~ spid + species, data = d3, sum)  # aggregate the 
    # datasets
    df <- data.frame(Site = site, scode, rep = rep1[i], year = yr, 
                     seasn = substr(occ, 1, 1), occasion=occ,d4)
    dff <- rbind(dff, df)
  }
  dat1 <- rbind(dat1, dff)
}
head(dat1)  # this is the full dataset that should be used

# End of Mark's code
# Convert into community data matrix

library(reshape)
comm <- cast(dat1, Site + rep + year + seasn ~ species, value = "indiv")
comm <- as.data.frame(comm)
comm[is.na(comm)] <- 0

comm <- comm[, - which(names(comm) == "Camponatus vestitus")]
comm <- comm[, - which(names(comm) == "Plagiolepis un02")]
comm <- comm[, - which(names(comm) == "Plagiolepis un03")]
head(comm)  # Final community dataframe to use. 
names(comm)

comm[, c(1,2,3,4,11)]
comm[, c(1,2,3,4,6)]

########################################
# 2. Import and format trait data.     #
########################################

########################################
# 3. Error check comm and trait df's   #
########################################
names(comm)[names(comm) == "Tetramorium setiferum"] <- "Tetramorium setuliferum"
names(comm)[names(comm) == "Tetramorium fridgidum"] <- "Tetramorium frigidum"


########################################
# 4. Import area data                  #
########################################

area1 <- read.csv("Data/pixelcounts_v2.csv")
area1 <- area1[2:9, c(2, 3, 6)]
names(area1) <- c("Site", "altitude", "pixels")
area1

area <- read.csv("Data/pixelcounts_v3.csv")
area <- area[2:9, c(2, 3, 6)]
names(area) <- c("Site", "altitude", "pixels")
area

########################################
# 5. Import temperature data           #
########################################

# ############################################################################
# # My own mothertruckin' thing
# 
# require(RODBC)
# con <- odbcDriverConnect("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=C:/Users/Tom Bishop/Documents/Research projects/PhD/1 First Paper/Paper 1 Restart November 2013/Data/SanitemperatureCS.mdb")
# a1<-sqlFetch(con, 'Temperatures')
# close(con)
# head(a1)
# 
# library(lattice)
# xyplot(Temperature ~ year, data = a1[a1$Site == 4 & a1$month == "September", ])
# 
# nrow(a1) # should be 610808
# str(a1)
# hist(a1$Temperature)
# sd(a1$Temperature)
# mean(a1$Temperature)
# 
# str(a1)
# a1$year <- format(a1$Datetime, format = "%Y")
# a1$month <- format(a1$Datetime, format = "%B")
# a1$day <- format(a1$Datetime, format = "%d")
# a1$Site <- substring(a1$SampleCode, 1, 1)
# head(a1)
# 
# a1temp <- a1[- which(a1$Period == 9 & a1$SampleCode == "1c"), ]
# a1 <- a1[- which(a1$Period == 10 & a1$SampleCode == "1c"), ]
# a1 <- a1[- which(a1$Period == 10 & a1$SampleCode == "3d"), ]
# a1 <- a1[- which(a1$Period == 13 & a1$SampleCode == "3a"), ]
# a1 <- a1[- which(a1$Period == 9 & a1$SampleCode == "1c"), ]
# a1 <- a1[- which(a1$Period == 3 & a1$SampleCode == "4d"), ]
# a1 <- a1[- which(a1$Period == 12 & a1$SampleCode == "6d"), ]
# a1 <- a1[- which(a1$Period == 14 & a1$SampleCode == "6a"), ]
# a1 <- a1[- which(a1$Period == 7 & a1$SampleCode == "7d"), ]
# a1 <- a1[- which(a1$Period == 10 & a1$SampleCode == "7d"), ]
# 
# t <- subset(a1, select = Temperature)
# gr <- subset(a1, select = c(Site, year, month))
# mean.temps <- aggregate(t, gr, mean)
# head(mean.temps)
# names(mean.temps)[names(mean.temps) == "Temperature"] <- "meantemp"
# var.temps <- aggregate(t, gr, var)
# names(var.temps)[names(var.temps) == "Temperature"] <- "vartemp"
# max.temps <- aggregate(t, gr, max)
# names(max.temps)[names(max.temps) == "Temperature"] <- "maxtemp"
# min.temps <- aggregate(t, gr, min)
# names(min.temps)[names(min.temps) == "Temperature"] <- "mintemp"
# sd.temps <- aggregate(t, gr, sd)
# names(sd.temps)[names(sd.temps) == "Temperature"] <- "sdtemp"
# 
# Site.temps <- merge(mean.temps, var.temps, all.x = TRUE, all.y = TRUE)
# Site.temps <- merge(Site.temps, max.temps, all.x = TRUE, all.y = TRUE)
# Site.temps <- merge(Site.temps, min.temps, all.x = TRUE, all.y = TRUE)
# Site.temps <- merge(Site.temps, sd.temps, all.x = TRUE, all.y = TRUE)
# Site.temps <- subset(Site.temps, month %in% c("January", "September"))
# Site.temps
# Site.temps$seasn <- ifelse(Site.temps$month == "January", "J", "S")
# Site.temps$rep <- Site.temps$SampleCode
# 
# head(Site.temps[, c("Site", "year", "seasn", "meantemp", "vartemp", "maxtemp", "mintemp", "sdtemp")])
# 
# t1 <- subset(a1, select = Temperature)
# gr2 <- subset(a1, select = c(Site, month))
# t2 <- subset(a1, select = Temperature)
# gr2 <- subset(a1, select = c(Site, month))
# mean.temps2 <- aggregate(t2, gr2, mean)
# head(mean.temps2)
# names(mean.temps)[names(mean.temps) == "Temperature"] <- "meantemp"
# Site.temps2 <- subset(mean.temps2, month %in% c("January", "September"))
# head(Site.temps2)
# Site.temps2
########################################
# 6.Site info frame                    #
########################################

site.info <- comm[, 1:4]
site.info$s.id <- rownames(comm)
site.info <- merge(site.info, area, by = "Site")
head(site.info)
site.info <- merge(site.info, Site.temps[!Site.temps$year %in% c("2005", "2013"), c("Site", "year", "seasn", 
                                                                                    "meantemp", "vartemp", 
                                                                                    "maxtemp", "mintemp", 
                                                                                    "sdtemp")], all.x = TRUE)

table(site.info$Site)


head(site.info)
site.info
nrow(site.info)

par(mfrow = c(2, 2))
plot(meantemp ~ altitude, data = site.info)
plot(sdtemp ~ altitude, data = site.info)
plot(maxtemp ~ altitude, data = site.info)
plot(mintemp ~ altitude, data = site.info)


########################################
# 7. Calculate biodiversity metrics    #
########################################

library(BiodiversityR)
library(vegan)


ric <- cbind(comm[, 1:4], diversityresult(x = comm[, 5:length(comm)], 
                                          method = "s", index = "richness"))
abund <- cbind(comm[, 1:4], diversityresult(x = comm[, 5:length(comm)], 
                                            method = "s", index = "abundance"))
evenness <- cbind(comm[, 1:4], diversityresult(x = comm[, 5:length(comm)], 
                                               method = "s", index = "Jevenness"))
evenness[is.na(evenness)] <- 0

maxab <- NULL  
for (i in 1:nrow(comm)){
  dat <- comm[i,5:length(comm)]
  maxabt <- max(dat)
  maxab <- c(maxab, maxabt)
}
dominance <- maxab/rowSums(comm[,5:length(comm)])
dominance[is.na(dominance)] <- 0
dominance <- cbind(comm[, 1:4], dominance)

PIE <- NULL
for (i in 1:nrow(comm)){
  dat <- comm[i,5:length(comm)]
  dat <- dat[dat >0]
  #dat <- ifelse(sum(dat) == 0, 0, dat)
  PIEt <- (1 - sum((dat/sum(dat))^2))
  PIE <- c(PIE, PIEt)
}
PIE <- cbind(comm[, 1:4], PIE)
nrow(ric)
nrow(dominance)
nrow(evenness)
nrow(PIE)
nrow(abund)

ric <- merge(ric, abund)
ric <- merge(ric, evenness)
ric <- merge(ric, dominance)
ric <- merge(ric, PIE)
ric <- merge(ric, site.info)

plot(richness ~ abundance, data = ric)
plot(richness ~ dominance, data = ric)
plot(richness ~ Jevenness, data = ric)
plot(richness ~ PIE, data = ric)
plot(dominance ~ abundance, data = ric)
plot(dominance ~ Jevenness, data = ric)
plot(dominance ~ PIE, data = ric)
plot(Jevenness ~ PIE, data = ric)
plot(Jevenness ~ abundance, data = ric)
plot(abundance ~ PIE, data = ric)



head(ric)
ric$fyear <- as.factor(ric$year)
ric$fSite <- as.factor(ric$Site)
ric$faltitude <- as.factor(ric$altitude)
ric$seasn2 <- ric$seasn
ric$seasn2 <- c("January", "September")
ric$Nseasn <- as.numeric(ric$seasn)
ric$Nyear <- as.numeric(ric$year)
ric$period <- as.factor(paste(ric$year, ric$seasn, sep = ""))
levels(ric$period)
ric$Nperiod <- as.numeric(ric$period)


ric$altitudepoly <- ric$altitude^2
ric$altitudeC <- scale(ric$altitude)
ric$altpolyC <- scale(ric$altitude^2)
ric$NseasnC <- scale(ric$Nseasn)
ric$pixelsC <- scale(ric$pixels)
ric$yearC <- scale(ric$year)
ric$meantempC <- scale(ric$meantemp)
ric$mintempC <- scale(ric$mintemp)
ric$maxtempC <- scale(ric$maxtemp)
ric$sdtempC <- scale(ric$sdtemp)
ric$vartempC <- scale(ric$vartemp)
ric$NseasnB <- c(0, 1)
library(car)
ric$pielogit <- logit(ric$PIE)
ric$newrep <- substring(ric$rep, 2, 2)
ric.mod <- ric[complete.cases(ric), ]
#ric.mod <- ric.mod[ric.mod$fyear %in% c("2006", "2007", "2008", "2009", "2010"), ]
head(ric.mod)

library(lattice)
xyplot(richness ~ altitude | fyear + seasn, data = ric.mod)
xyplot(abundance ~ altitude | fyear + seasn, data = ric.mod)
xyplot(dominance ~ altitude | fyear + seasn, data = ric.mod)
xyplot(Jevenness ~ altitude | fyear + seasn, data = ric.mod)
xyplot(PIE ~ altitude | fyear + seasn, data = ric.mod) 
# should probably exclude PIE = 1 because this is for communities with no species...

plot(meantemp ~ maxtemp, data = ric)
plot(meantemp ~ mintemp, data = ric)




