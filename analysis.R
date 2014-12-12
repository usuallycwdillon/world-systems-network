## Beginning 'the project' after having finally found data to build a historic
#  world-systems network to match the COW period 1820--. Of course, there seems
#  to be no (reliable, large) trade data that far back, so this really starts
#  at 1870 and goes through 2009; 140 years. 

library("doParallel") # Optional. Only useful a couple of times
library("foreach")    # Required for plyr to work in parallel
library("plyr")       # Handy way to split, process and recombine data
library("reshape")    # For melting horizontal data
library("mice")       # For imputing missing per capita GDP values
library("igraph")     # For storing and analysing the annual trade graphs

## The data munging activity was getting unhandily long, so I moved it over.
source("munging.R")
## Just to be tidy about it, the code to produce graphs is also moved over.
#  But it gets called part-way through the munging. 


## Network Analysis ############################################################
## This process can also be done in parallel.
# Set up the multi-core processor to act like a small cluster
cnodes <- detectCores()
cl <- makeCluster(cnodes)
registerDoParallel(cl)


## Coreness measures k-cores ##
# world systems analysts abandoned the k-core measure when Bagatti & Everett
# published their algorith to detect core-periphery structures in 1999, but
# it's an available measure and it's fast to compute it. Let's look...
#
# ocoreness <- function(g) {
#   ocores <- graph.coreness(g, mode="out")
#   return(ocores)
# }
# 
# icoreness <- function(g) {
#   icores <- graph.coreness(g, mode="in")
#   return(icores)
# }
# 
# ocorelist <- llply(wgraphs, ocoreness)
# names(ocorelist) <- graph.names
# 
# socorelist <- llply(sgraphs, ocoreness)
# names(socorelist) <- graph.names
# 
# icorelist <- llply(wgraphs, icoreness)
# names(icorelist) <- graph.names
# 
# sicorelist <- llply(sgraphs, icoreness)
# names(sicorelist) <- graph.names


## InfoMap Community ##
# The only algorithm I could find to score modularity and decompose communities
# of directed (and weighted!) graphs is infomap community. The source of the 
# algorithm is by Rosvall and Bergstrom, 2007;
# http://www.pnas.org/content/105/4/1118.abstract
#
ptm <- proc.time() # Just how long *does* it take to do this? 
imclist <- llply(wgraphs, infomap.community, .parallel=TRUE)
names(imclist) <- graph.names
simclist <- llply(sgraphs, infomap.community, .parallel=TRUE)
names(simclist) <- graph.names
proc.time() - ptm # Suprisingly fast, this part


## clusters ##
#
# wwclusters <- llply(wgraphs, clusters, mode="weak", .parallel=TRUE)
# wsclusters <- llply(wgraphs, clusters, mode="strong", .parallel=TRUE)
# swclusters <- llply(sgraphs, clusters, mode="weak", .parallel=TRUE)
# ssclusters <- llply(sgraphs, clusters, mode="strong", .parallel=TRUE)
# names(wwclusters) <- graph.names
# names(wsclusters) <- graph.names
# names(swclusters) <- graph.names
# names(ssclusters) <- graph.names


## cohesive blocks ##
# 
# We can only do cohesive blocking on undirected graphs. Direction may not
# matter as much as it seems like it should. (The core is less obvious but is
# still present, according to __some authors__.)
# wugraphs <- llply(wgraphs, as.undirected, mode="collapse", .parallel=TRUE)
# sugraphs <- llply(sgraphs, as.undirected, mode="collapse", .parallel=TRUE)
# names(wugraphs) <- graph.names
# names(sugraphs) <- graph.names
# 
# wblocks <- llply(wugraphs, cohesive.blocks, labels=TRUE, .parallel=TRUE)
# sblocks <- llply(sugraphs, cohesive.blocks, labels=TRUE, .parallel=TRUE)
# names(wblocks) <- graph.names
# names(sblocks) <- graph.names

## Leading Eigenvector Communities ##
#  Similar to the average nearest-neighbor degree (ANND) proposed by 
# Borgatti and Everett in 1999. This method seems to produce reasonable 
# results.

les <- llply(sugraphs, leading.eigenvector.community, .parallel=TRUE)
names(les) <- graph.names


## Walktrap Communities ##
wcs <- llply(sugraphs, walktrap.community, .parallel=TRUE)
names(wcs) <- graph.names

## Stop the cluster, again
stopCluster(cl)

# Put the data back together again so that we can keep going.
source("reassembling.R")

###############################################################################


## fun with labels
# 
wsa$ccode1 <- wsa$ccode
wsa$StartYear1 <-  wsa$year
wsa$core1 <- wsa$core

c.extra.net <- merge(calf.extra, wsa, by=c("ccode1", "StartYear1"))
c.extra.net$core2 <- FALSE
c.extra.net$ccode2 <- NA
c.extra.net$w.type <- 0

# Exra-state wars are either core-periphery (type 2) or periphery-periphery 
# (type 3), depending on whether SideA (by core1) is core or periphery; the
# country code for SideB is always -8.
for (i in seq_along(c.extra.net$w.type)){
  c.extra.net$w.type[i] <- if(c.extra.net$core1[i]==TRUE) {2} else {3}
}

wsa$CcodeA <- wsa$ccode1
c.intra.net <- merge(calf.intra, wsa, by=c("CcodeA", "StartYear1"))
c.intra.net <- c.intra.net[,c(1,2,3,4,5,6,7,10,11,12,13,22,23,24)]
wsa$CcodeB <- wsa$ccode
c.intra.net2 <- merge(c.intra.net, wsa, by=c("CcodeB", "StartYear1"))
c.intra.net <- c.intra.net2[,c(4,5,2,3,6,1,7,14,27)]
names(c.intra.net) <- c("warnum", "warname", "year", "ccode.A", "side.A",
                        "ccode.B", "side.B", "core.A", "core.B")


#calf.inter$row.names <- row.names(calf.inter)
c.inter.net <- merge(calf.inter, wsa, by=c("ccode", "StartYear1"))
c.inter.net <- c.inter.net[,c(-22,-21,-20, -19)]


warType <- function(d) { #
  warnum <- unique(d$WarNum)
  waryear<- min(d$StartYear1)
  d$lik <- 0
  for(i in seq_along(d$lik)){
    d$lik[i] <- if(d$core[i]==TRUE){1} else {2}
  }
  core.sum <- mean(d$lik)
  w.type <- ifelse(core.sum==1, 1, 
                            ifelse(core.sum==2, 3, 2))
  wt <- data.frame("warnum"=warnum, "waryear"=waryear, "wartype"=w.type)
  return(wt)
}

warType2 <- function(d){
  warnum <- unique(d$warnum)
  waryear<- min(d$year)
  d$lik.A <- 0
  d$lik.B <- 0
  for(i in seq_along(d$lik.A)){
    d$lik.A[i] <- if(d$core.A[i]==TRUE){1} else {2}
  }
  for(i in seq_along(d$lik.B)){
    d$lik.B[i] <- if(d$core.B[i]==TRUE){1} else {2}
  }
  d$lik <- (d$lik.A + d$lik.B) / 2
  core.sum <- mean(d$lik)
  w.type <- ifelse(core.sum==1, 1, 
                   ifelse(core.sum==2, 3, 2))
  wt <- data.frame("warnum"=warnum, "waryear"=waryear, "wartype"=w.type)
  return(wt)
}


n <- c("warnum", "waryear", "wartype")

inter.typifiedWars <- ddply(c.inter.net, .(WarNum), warType)
inter.typifiedWars <- inter.typifiedWars[,c(1,3,4)]
names(inter.typifiedWars) <- n

intra.typifiedWars <- ddply(c.intra.net, .(warnum), warType2)

extra.typifiedWars <- c.extra.net[,c(3,2,25)]
names(extra.typifiedWars) <- n


wars <- rbind(extra.typifiedWars, inter.typifiedWars, intra.typifiedWars)

warsByType <- split(wars, wars$wartype)

# No core-core wars! only between core and periphery or within the periphery.
names(warsByType)  

c.p.wars <- warsByType[[1]]
p.p.wars <- warsByType[[2]]



## Beyond here thar be dragons ################################################# ...not really. This is where the hazard analysis goes. But days and days ago
# (when I started this project) it represented the edge of the project, so I
# wrote that to keep the edge just under the coding horizon. 
###############################################################################

# Sort dataframes by year
c.p.wars <- c.p.wars[order(c.p.wars$waryear),]
p.p.wars <- p.p.wars[order(p.p.wars$waryear),]

c.p.wars$Ot <- 0
p.p.wars$Ot <- 0

for (i in seq_along(c.p.wars$Ot)) {
  j <- i + 1
  if(j == length(c.p.wars$Ot)) {break}
  c.p.wars$Ot[j] <- c.p.wars$waryear[j]-c.p.wars$waryear[i]
}

for (i in seq_along(p.p.wars$Ot)) {
  j <- i + 1
  if(j == length(p.p.wars$Ot)) {break}
  p.p.wars$Ot[j] <- p.p.wars$waryear[j]-p.p.wars$waryear[i]
}

# frequency table of 'time' between war onsets
c.p.obs <- as.data.frame(table(c.p.wars$Ot))
c.p.fre <- summarize(c.p.obs, var=Freq)

p.p.obs <- as.data.frame(table(p.p.wars$Ot))
p.p.fre <- summarize(p.p.obs, var=Freq)

# cumulative sums of frequencies
c.p.cf <- cumsum(c.p.fre)
p.p.cf <- cumsum(p.p.fre)

# cumulative relative frequencies
c.p.rf <- lapply(c.p.cf, function(x) x/max(c.p.cf$var))
p.p.rf <- lapply(p.p.cf, function(x) x/max(p.p.cf$var))

# pull them together
cp.dist <- cbind(c.p.obs, c.p.fre, c.p.cf, c.p.rf)
pp.dist <- cbind(p.p.obs, p.p.fre, p.p.cf, p.p.rf)

n <- c("Ot", "frequency", "freq.int", "cum.freq", "rel.freq" )
names(cp.dist) <- n
names(pp.dist) <- n

# get plottable values from factors
cp.dist$Ot.int <- sapply(cp.dist$Ot, function(x) as.numeric(levels(x)[x]))
pp.dist$Ot.int <- sapply(pp.dist$Ot, function(x) as.numeric(levels(x)[x]))


plot(cp.dist$Ot.int, cp.dist$rel.freq, type="o",
     xlab="Periods between onset", ylab="Relative frequency",
     main="Cumulative frequency of core-periphery wars onsets")

plot(pp.dist$Ot.int, pp.dist$rel.freq, type="o",
     xlab="Periods between onset", ylab="Relative frequency",
     main="Cumulative frequency of war onsets within the periphery")

cp.dist$iv <- cp.dist$Ot.int
l <- length(cp.dist$iv)
k <- l - 1

cp.dist$iv[2:l] <- cp.dist$Ot.int[2:l] - cp.dist$Ot.int[1:k]
cp.dist$mid <- 0
cp.dist$mid[2:l] <- (cp.dist$iv[2:l] / 2) + cp.dist$Ot.int[1:k]

cp.dist$rfi <- cp.dist$rel.freq
cp.dist$rfi[2:l] <- cp.dist$rel.freq[2:l] - cp.dist$rel.freq[1:k]
cp.dist$rfd <- cp.dist$rfi / 2
cp.dist$rfm <- cp.dist$rfd
cp.dist$rfm[2:l] <- cp.dist$rfd[2:l] + cp.dist$rel.freq[1:k]

cp.dist$kme <- cp.dist$rfi / cp.dist$iv
cp.dist


pp.dist$iv <- pp.dist$Ot.int
l <- length(pp.dist$iv)
k <- l - 1

pp.dist$iv[2:l] <- pp.dist$Ot.int[2:l] - pp.dist$Ot.int[1:k]
pp.dist$mid <- 0
pp.dist$mid[2:l] <- (pp.dist$iv[2:l] / 2) + pp.dist$Ot.int[1:k]

pp.dist$rfi <- pp.dist$rel.freq
pp.dist$rfi[2:l] <- pp.dist$rel.freq[2:l] - pp.dist$rel.freq[1:k]
pp.dist$rfd <- pp.dist$rfi / 2
pp.dist$rfm <- pp.dist$rfd
pp.dist$rfm[2:l] <- pp.dist$rfd[2:l] + pp.dist$rel.freq[1:k]

pp.dist$kme <- pp.dist$rfi / pp.dist$iv
pp.dist


plot(cp.dist$iv, cp.dist$kme, type="o",
     xlab="Periods between onset", 
     ylab="Probability density: war onset with core & periphery",
     main="Calculated Kaplan-Meier estimate:  \nprobability density of war onset \nbetween core and periphery countries")


plot(pp.dist$iv, pp.dist$kme, type="o",
     xlab="Periods between onset", 
     ylab="Probability density: war onset within periphery",
     main="Calculated Kaplan-Meier estimate:  \nprobability density of war onset \nwithin the periphery countries")



