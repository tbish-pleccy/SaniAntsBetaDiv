# need Vegan package for this
#wd<-"C:\\projects\\AntsSani"
#wd<-"E:\\AntsSani" # alternative

######################
# data preparation
setwd(wd) # set the working directory
#source(file.path(wd,"functions.r")) # run the functions
dd<-read.csv("AntsSaniJune2013.csv") # read the dataset
head(dd) # view first few lines
# create dataset of zeros (where nothing was collected) and then combine with existing data 
us<-unique(dd[,7:8]) # unique species identifier and species name
nsp<-nrow(us) # number of species
s<-sort(rep(1:8,4)) # site numbers
rp<-rep(letters[1:4],8) # letters for reps
rep1<-paste(s,rp,sep="") # replicate names e.g. 1a, 1b, 1c
uoc<-unique(dd[,c(4,6)]) # unique year and occasions (e.g. J06, S06)
nuoc<-nrow(uoc) # number of occasions
dat1<-data.frame() # create a dataframe
for (j in 1:nuoc){
  yr<-uoc$year[j] # a particular year
  occ<-uoc$occasion[j] # a particular occasion
  d1<-subset(dd,occasion==occ) # actual data for this occasion
  dff<-data.frame()
  for (i in 1:32){ # for all replicates from 1a to 8d
    site<-s[i] # site number 
    scode<-paste("s",s[i],sep="") # site code e.g. s1
    dy<-data.frame(us,indiv=0) # all species with individuals as zero (spid, species, indiv=0)
    d2<-subset(d1,rep==rep1[i]) # select data for this particular replicate
    d3<-rbind(d2[,7:9],dy) # combine actual data with zero individuals data
    d4<-aggregate(indiv~spid+species,data=d3,sum) # aggregate the datasets
    df<-data.frame(Site=site,scode,rep=rep1[i],year=yr,seasn=substr(occ,1,1),occasion=occ,d4)
    dff<-rbind(dff,df)
  }
  dat1<-rbind(dat1,dff)
}
head(dat1)  # this is the full dataset that should be used

names(dat1)
unique(dat1$year)
### MDS for sites only
aa<-aggregate(indiv~rep+species,data=dat1,sum)
#s<-paste(aa$rep,aa$occasion,sep="_")
#aa<-data.frame(Site=s,Species=aa$species,Individiuals=aa$indiv)
names(aa)<-c("Site","Species","Individuals")
us<-sort(unique(aa$Species))  # list of unique species (check these for spelling mistakes and inconsistencies)
nspp<-length(us)    # number of species
reps<-sort(unique(aa$Site))   # replicate names  (these are labeled Traps)
nreps<-length(reps)    # number of replicates (sites)
data<-speciesbysites(aa) # convert into a species by sites matrix
nms<-names(data)
#f1<-match("8D",nms) # find this
#f2<-match("5v",nms) # find this
dat<-data[,-c(1)] # remove the first column (the names)
#dat<-data

nm<-data$Species
coln<-names(dat)
data2<-t(dat)
names(data2)<-nm
rnms<-sort(c(st,st,st,st),decreasing=T)
row.names(data2)=rnms
sol <- metaMDS(data2)
windows(9,9)
plot(sol, type="t")

###########################################################


#an2<-anosim(data2, grouping=as.factor(rnms), permutations = 999, distance = "bray")
#summary(an2)
#plot(an2)

# MDS sites and seasons
aa<-aggregate(indiv~rep+seasn+species,data=dd,sum)
s<-paste(aa$rep,aa$seasn,sep="_")
aa<-data.frame(Site=s,Species=aa$species,Individiuals=aa$indiv)
#names(aa)<-c("Site","Species","Individuals")
us<-sort(unique(aa$Species))  # list of unique species (check these for spelling mistakes and inconsistencies)
nspp<-length(us)    # number of species
reps<-sort(unique(aa$Site))   # replicate names  (these are labeled Traps)
nreps<-length(reps)    # number of replicates (sites)
data<-speciesbysites(aa) # convert into a species by sites matrix
nms<-names(data)
f1<-match("8D_J",nms) # find this
f2<-match("5v_J",nms) # find this
dat<-data[,-c(1,f1,f2)]

nm<-data$Species
coln<-names(dat)
data2<-t(dat)
names(data2)<-nm
rnms<-sort(c(st,st,st,st),decreasing=T)
ss<- rep(c("J","S"),32) # seasons
rnms<-paste(rnms,ss,sep="")
row.names(data2)<-rnms
sol <- metaMDS(data2)
windows(9,9)
plot(sol, type="t")

# anosim
an1<-anosim(data2, grouping=as.factor(rnms), permutations = 999, distance = "bray")
summary(an1)
plot(an1)
