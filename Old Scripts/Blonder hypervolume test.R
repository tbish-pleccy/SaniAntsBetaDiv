library(hypervolume)

data(iris)
hv1 <- hypervolume(subset(iris, Species=="setosa")[,1:4],reps=1000,bandwidth=0.2)
summary(hv1)
hv2 <- hypervolume(subset(iris, Species=="virginica")[,1:4],reps=1000,bandwidth=0.2)


hypervolume_inclusion_test(hv1, points = subset(iris, Species=="virginica")[,1:4])

str(list(hv1, hv2))
# choose fixed axes
plot(hv1, pair=FALSE, npmax=500, varlims=list(c(3,6),c(2,5),c(0,3),c(-1,1)))
plot(hv2, pair=FALSE, npmax=500, varlims=list(c(3,6),c(2,5),c(0,3),c(-1,1)))

hlist <- new("HypervolumeList")
hlist@HVList <- vector(mode = "list",length = 2)
hlist@HVList[[1]] <- hv1
hlist@HVList[[2]] <- hv2


# construct a hypervolume of points in the unit square [0,1] x [0,1]
data = data.frame(x=runif(100,min=0,max=1), y=runif(100,min=0,max=1))
hv = hypervolume(data, reps=1000, bandwidth=0.1)
# test if (0.5,0.5) and (-1,1) are in - should return TRUE FALSE
hypervolume_inclusion_test(hv, points=data.frame(x=c(0.5,-1),y=c(0.5,-1)))

plot(hlist, pair=FALSE, npmax=500, varlims=list(c(3,6),c(2,5),c(0,3),c(-1,1)))



hv_set <- hypervolume_set(hv1, hv2, check_memory = FALSE)
sapply(hv_set@HVList, function(x) {x@Volume})

hv_set2 <- hypervolume_set(hv2, hv1, check_memory = FALSE)
sapply(hv_set2@HVList, function(x) {x@Volume})


demo('finch', package='hypervolume')





sani <- hypervolume(pcoa.axesDF, reps=1000, bandwidth=0.2)
plot(sani, pair=FALSE, npmax=500)

sani1 <- comm[comm$Site == 1, ]
sani1 <- subset(sani1, select = - c(year, Site, rep, seasn, s.id))
sani1 <- names(sani1[, colSums(sani1) > 0])

sani5 <- comm[comm$Site == 5, ]
sani5 <- subset(sani5, select = - c(year, Site, rep, seasn, s.id))
sani5 <- names(sani5[, colSums(sani5) > 0])

pcoa.1 <- pcoa.axesDF[rownames(pcoa.axesDF) %in% sani1, ]
pcoa.5 <- pcoa.axesDF[rownames(pcoa.axesDF) %in% sani5, ]

sanitop <- hypervolume(pcoa.1, reps=1000, bandwidth=0.2)
sanimid <- hypervolume(pcoa.5, reps=1000, bandwidth=0.2)

sani_set <- hypervolume_set(sanitop, sanimid, check_memory = FALSE)
sapply(sani_set@HVList, function(x) {x@Volume})

sanilist <- new("HypervolumeList")
sanilist@HVList <- vector(mode = "list",length = 2)
sanilist@HVList[[1]] <- sanitop
sanilist@HVList[[2]] <- sanimid
plot(sanilist, pair=FALSE, npmax=500, colors = c("blue", "red"))
plot(sanilist, pair=TRUE, npmax=500, colors = c("blue", "red"), varlims = list(c(-3,3),c(-3,3),c(-3,3),c(-3,3)))



inter <- sani_set@HVList$Intersection
union <- sani_set@HVList$Union
uni1.top <- sani_set@HVList$Unique_1
uni2.mid <- sani_set@HVList$Unique_2





plot(inter, varlims = list(c(-3,3),c(-3,3),c(-3,3),c(-3,3)))
plot(union, varlims = list(c(-3,3),c(-3,3),c(-3,3),c(-3,3)))
plot(uni1.top, varlims = list(c(-3,3),c(-3,3),c(-3,3),c(-3,3)))
plot(uni2.mid, varlims = list(c(-3,3),c(-3,3),c(-3,3),c(-3,3)))