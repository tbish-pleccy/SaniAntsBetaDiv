library(FD)
library(betapart)
library(vegan)
library(boot)
library(ecodist)
library(sp)
library(ape)
library(MuMIn)
library(cluster)
library(sm) 

########################################
# 1. Import and format abundance data. #
########################################

dd <- read.csv("Data/Data_for_R_JulyUpdate.csv")

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

#comm <- comm[order(comm$year, comm$seasn, comm$Site, comm$rep), ]

comm$s.id <- seq(1:nrow(comm))
rownames(comm) <- comm$s.id
head(comm)  # Final community dataframe to use. 
names(comm) <- gsub("sp. ", "", x = names(comm))
names(comm) <- gsub(" ", ".", x = names(comm))

#comm$rep <- substring(comm$rep, 2, 2)

########################################
# 2. Import and format trait data.     #
########################################

# Functions for this section

vlookup <- function(val, df, row){
  df[df[1] == val, row][1] }
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Import and create relative traits

raw.traits <- read.csv("Data/Trait Data.csv")
raw.traits$Leg.Length <- raw.traits$Hind.Femur + raw.traits$Hind.Tibia
raw.traits$Eye.Position <- ifelse(raw.traits$Interocular > 0, 
                                  raw.traits$Head.Width - raw.traits$Interocular, 
                                  raw.traits$Head.Width - raw.traits$Head.Width)
  
raw.traits$Mand.Index <- raw.traits$Mandible.Length / raw.traits$Head.Width


raw.traits <- raw.traits[,-16]  # remove pronotum length
raw.traits <- raw.traits[, -26]  # remove notes
rel.traits.list <- as.vector((names(raw.traits[, c(7:15, 17:18, 26:28)])))
for (i in rel.traits.list){
  rel <- raw.traits[,i]/raw.traits[,16]
  name <- paste("Rel.", i, sep = "")
  raw.traits <- cbind(raw.traits, rel)
  colnames(raw.traits)[ncol(raw.traits)] <- name
}

mean.traits.list <- c(7:18, 26:42)  # traits where the mean shoudl be taken
mode.traits.list <- c(19:25)  # traits where the mode should be taken

#subset(raw.traits, select = c("Pilosity", "Alitrunk.Spines", "Petiolar.Spines", 
 #                             "Sculpturing", "Head.Colour", "Mesosoma.Colour", 
  #                            "Gaster.Colour"))


a2 <- as.vector(list(names(raw.traits)))  # df of trait names
b2 <- c(1:length(raw.traits))
c2 <- data.frame(b2,a2)

# Create dummy frame to attach loops to.

morph1 <- data.frame(aggregate(raw.traits$Head.Length,list(raw.traits$Genus,
                                                           raw.traits$Species,
                                                           raw.traits$Sp.Code),
                               mean))
names(morph1) <- c("Genus","Species","Sp.Code","1")
reduced <- morph1[, c(1:3)]

# Loop to generate mean trait values

for (i in mean.traits.list){     
  trait.tab <- raw.traits[, c(1, 2, 3, i)]
{
  names(trait.tab) <- c("Genus", "Species", "Sp.Code", "temp.name")
  trait.tab$temp.name <- as.numeric(trait.tab$temp.name)
  morph1.temp <- data.frame(aggregate(trait.tab$temp.name,
                                      list(trait.tab$Genus, trait.tab$Species,
                                           trait.tab$Sp.Code), 
                                      mean,na.rm = TRUE))
  sn <- toString(vlookup(i, c2,2))                                                                                                                           
  names(morph1.temp) <- c("Genus", "Species", "Sp.Code", sn)
  reduced <- merge(reduced, morph1.temp, all.x = TRUE)
}
}
mean.traits <- reduced
head(mean.traits)

# Loop to generate mode trait values

reduced2 <- morph1[, c(1:3)]
for (i in mode.traits.list){     
  trait.tab <- raw.traits[, c(1, 2, 3, i)]
{
  names(trait.tab) <- c("Genus", "Species", "Sp.Code", "temp.name")
  trait.tab$temp.name <- as.numeric(trait.tab$temp.name)
  morph2.temp <- data.frame(aggregate(trait.tab$temp.name,
                                      list(trait.tab$Genus, trait.tab$Species,
                                           trait.tab$Sp.Code), FUN = Mode))
  sn <- toString(vlookup(i, c2, 2))                                                                                                                           
  names(morph2.temp) <- c("Genus", "Species", "Sp.Code", sn)
  reduced2 <- merge(reduced2, morph2.temp, all.x = TRUE)
}
}
mode.traits <- reduced2
head(mode.traits)

trait <- merge(mean.traits, mode.traits)

# Loop to generate dominant colour variable

sp.code.list<-as.vector(unique(raw.traits[, 3]))
d.col <- NULL
r <- NULL
for (i in sp.code.list){
  sp <- as.matrix(raw.traits[raw.traits$Sp.Code==i, 23:25])
  m <- Mode(sp)
  r <- data.frame(cbind(i, m))
  names(r) <- c("Sp.Code", "Dominant.Colour")
  d.col <- rbind(d.col,r)
}

trait <- merge(d.col, trait)  # Final trait data frame to be used.
head(trait)
trait$name <- paste(trait$Genus, trait$Species, sep = ".")

########################################
# 3. Error check comm and trait df's   #
########################################
names(comm)[names(comm) == "Tetramorium.setiferum"] <- "Tetramorium.setuliferum"
names(comm)[names(comm) == "Tetramorium.fridgidum"] <- "Tetramorium.frigidum"

trait.species <- unique(trait[, c(1, 3, 4)])
nrow(trait.species) == length(names(comm[, 6:ncol(comm)-1]))  # correct number of species...

trait.names <- sort(as.vector(paste(trait.species$Genus, trait.species$Species, 
                                    sep = ".")))
comm.names <- sort(as.vector(names(comm[, 6:ncol(comm)-1])))
name.frame <- data.frame(cbind(trait.names, comm.names))
name.frame$match <- as.vector(ifelse(levels(name.frame$comm.names) == 
                                       levels(name.frame$trait.names), "MATCH", 
                                     "NO"))
name.frame

########################################
# 4. Import area data                  #
########################################

area <- read.csv("Data/pixelcounts_v3.csv")
area <- area[2:9, c(2, 3, 6)]
names(area) <- c("Site", "altitude", "pixels")
area



########################################
# 6.Site info frame                    #
########################################

site.info <- comm[, c("Site", "rep", "year", "seasn", "s.id")]
site.info <- merge(site.info, area, by = "Site")
site.info$rep <- paste(site.info$Site, site.info$rep, sep = "")

spatialcoords <- read.csv("Data/SiteCoordinatesSaniCorrected.csv")
spatialcoords$Site <- substring(spatialcoords$Sitename, 1, 1)
spatialcoords$rep <- as.factor(substring(spatialcoords$Sitename, 2, 2))
levels(spatialcoords$rep) <- c("a", "b", "c", "d")
spatialcoords$rep <- paste(spatialcoords$Site, spatialcoords$rep, sep = "")
spatialcoords <- subset(spatialcoords, select = - Sitename)

site.info <- merge(site.info, spatialcoords)

period.frame <- unique(site.info[, c("year", "seasn")])
period.frame$period <- rownames(period.frame)

site.info <- merge(site.info, period.frame)
rownames(site.info) <- site.info$s.id
site.info <- site.info[order(as.numeric(rownames(site.info))), ]

########################################
# Finalising                           #
########################################

# Community Data
community <- comm
community <- subset(community, select = - c(s.id, Site, rep, year, seasn))
rownames(community) <- comm$s.id
community[community > 0] <- 1
community <- community[, order(colnames(community))] # alphabetically order!
community <- community[as.vector(which(rowSums(community) > 0)), ] # remove empty communities
community <- community[rowSums(community) > 4, ] # only communities with more than 4 species. 
site.info.reps <- site.info[site.info$s.id %in% rownames(community), ]

community.list <- list()
for (i in unique(site.info.reps$year)){
  for (j in unique(site.info.reps$seasn)){
    select <- site.info.reps[site.info.reps$year == i & 
                          site.info.reps$seasn == j, "s.id"]
    temp <- community[rownames(community) %in% select, ]
    temp <- temp[rowSums(temp) > 0, colSums(temp) > 0]
    community.list[[paste(i, j, sep = "")]] <- temp
  }
}



library(BiodiversityR)
library(vegan)


ric <- cbind(comm[, "Site"], diversityresult(x = comm[, 5:length(comm)], 
                                          method = "s", index = "richness"))
ric$s.id <- seq(1:nrow(ric))
names(ric) <- c("Site", "richness", "s.id")
ric <- merge(ric[ric$s.id %in% site.info.reps$s.id, ], site.info.reps)
library(lattice)
xyplot(richness ~ altitude | year + seasn, data = ric)



# Trait Data
tFD <- trait[, c("name", "Webers.Length", "Rel.Eye.Position", "Mand.Index", "Rel.Leg.Length")]
tFD.temp <- tFD
rownames(tFD.temp) <- tFD.temp$name
tFD.temp <- tFD.temp[order(tFD.temp$name), ]
tFD.temp <- subset(tFD.temp, select = - name)

# Calculate PCoA axes
library(FD)
fundiv <- dbFD(x = tFD.temp, a = community, print.pco = TRUE)
pcoa.axes <- as.matrix(fundiv$x.axes)
pcoa.axesDF <- data.frame(pcoa.axes)

testing <- data.frame(dbFD(x = tFD.temp, a = community))
testing <- cbind(site.info.reps, testing)

library(lattice)
xyplot(nbsp ~ altitude | year + seasn, data = testing)


testing<-(data.frame(rownames(community), rowSums(community)))
names(testing) <- c("s.id", "ric")

site.info.test <- site.info
testing <- merge(site.info.test, testing)

xyplot(ric ~ altitude | year + seasn,data = testing)
