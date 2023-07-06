#make sure directory points to:
#oneSimDARcores.txt, oneSimDARprof.txt (simulated data from SHP -- no direct bearing to original as based entirely on simulations)
#multinomDAR.R (code that implements the DAR error to the FLXMRmultinom class in felxmix library)

oneSimDARcores <- as.matrix(read.table(file = "oneSimDARcores.txt", sep = " "))
oneSimDARprof <- as.matrix(read.table(file = "oneSimDARprof.txt", sep = " "))
#data restructuring (wide to long)
df1 <- data.frame(oneSimDARcores)
colnames(df1) <- 1:45 #coding to make things easier later
df1$id <- 1:1650
newdat1 <- reshape::melt(df1, id=c("id"),variable_name="time")
newdat1$time <- as.integer(newdat1$time) #better as a number
colnames(newdat1)[3] <- "coRes"
newDat1 <- newdat1[with(newdat1, order(id, time)), ]  #sort data
#again now for prof
df2 <- data.frame(oneSimDARprof)
colnames(df2) <- 1:45 #coding to make things easier later
df2$id <- 1:1650
newdat2 <- reshape::melt(df2, id=c("id"),variable_name="time")
newdat2$time <- as.integer(newdat2$time) #better as a number
colnames(newdat2)[3] <- "prof"
newDat2 <- newdat2[with(newdat2, order(id, time)), ]  #sort data
newDat <- merge(newDat1,newDat2,by=c("id","time"))
newDat <- newDat[with(newDat, order(id, time)), ]  #re-sort data (needed after merge)

library(flexmix)
library(TraMineR) 
library(TraMineRextras)
#
# reduce nrep, k, if you want it to be faster.
#
set.seed(1001)
bagOfTokens <- stepFlexmix(.~0|id,k=2:4,nrep=5,model=list(FLXMRmultinom(coRes~1), FLXMRmultinom(prof~1)), data=newDat,control=list(iter.max=300,tolerance=1e-05,verbose=1))
#
set.seed(1001)
linearDrift <- stepFlexmix(.~0|id,k=2:4,nrep=5,model=list(FLXMRmultinom(coRes~time), FLXMRmultinom(prof~time)), data=newDat,control=list(iter.max=300,tolerance=1e-05,verbose=1))

# The DAR code uses optim function after an initial guess based on nnet call.
#
# WARNING: optim can be VERY SLOW in this context
#
# Speedup (with greater potential for local maximum) by decreasing maxit in each DAR channel call
# or through usual approaches (iter.max, tolerance), but with usual caveats about local maxima
# greater concern with mixture models is reduction in number of groups due to small size.  Tolerance for small groupsize is a parameter of flexmix, and larger nrep will provide more successful runs with the desired number of groups
#

set.seed(1001)
subset.index <- sort(sample(1:nrow(newDat),15000,replace=F))  # this is still slow, but demonstrates the core capabilities
source("multinomDAR.R")
ids <- newDat$id
times <- as.integer(newDat$time) #just to be sure it is an integer
linearDriftDAR <- stepFlexmix(.~0|id,k=2:4,nrep=5,model=list(FLXMRmultinomDAR(coRes~time,id=ids[subset.index],times=times[subset.index],maxit=300), FLXMRmultinomDAR(prof~time,id=ids[subset.index],times=times[subset.index],maxit=300)), data=newDat[subset.index,],control=list(iter.max=300,tolerance=1e-05,verbose=1))

#extract/assign clusters
mdl <- getModel(linearDriftDAR,"AIC")  #USE OF BIC WILL LIKELY YIELD FEWER CLUSTER WITH THESE DATA.
cluster.tbl <- cbind(ids[subset.index],mdl@cluster)
uniq.id <- unique(cluster.tbl[,1])
clust.asn <- cluster.tbl[match(uniq.id,cluster.tbl[,1]),2]

#plot 2 channels
cores.labels <- c("Two biol. parents", "One biol. parent", "One biol. parent & ptnr", "Alone", "With partner", "W/ Ptnr & biol. child", "W/ ptnr & non-biol. child", "W/ biol. child w/o ptnr", "Friends", "Other")
prof.labels <- c("Full time", "Part time", "Insurances", "Other", "At home", "Retirement", "Education")

seqdplot(seqdef(oneSimDARcores[uniq.id,],labels=cores.labels),group=clust.asn)
seqdplot(seqdef(oneSimDARprof[uniq.id,],labels=prof.labels),group=clust.asn)



