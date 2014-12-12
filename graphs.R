## This chunk of code makes graphs from the 'munging.R' -processed data.

library("doParallel") # Optional. Only useful a couple of times
library("foreach")    # Required for plyr to work in parallel
library("plyr")       # Handy way to split, process and recombine data
library("igraph")     # For storing and analysing the annual trade graphs

# This takes a long time, so we do it in parallel to save some time
# Set up the multi-core processor to act like a small cluster
cnodes <- detectCores()
cl <- makeCluster(cnodes)
registerDoParallel(cl)

# This function will take the dataframes from the list and make iGraph objects
# (named by the year), merging the data from the two lists of data frames. 
makeGraphs <- function(x, y) {
  y <- y[order(y$gdppc),]
  g <- graph.data.frame(x,
                        directed = TRUE,
                        vertices = y)
  return(g)
}

graphs <- mapply(makeGraphs, daglist, cya, SIMPLIFY=FALSE)
names(graphs) <- graph.names

# Clean-up some big piles of data not needed anymore. All we need now are the
# two lists: the trade edgelist (daglist) and the country-year attributes (cya).
rm(dag1, dag2, dyadTrade, mpd, gdp, imp.country)


## Relative Openness ###########################################################
# We have graphs but we're still not quite done. Next, we need to calculate
# graph edge weights as relative openness, using Chase-Dunn's formula:
# Openness = Imports / per capita GDP. He does it for the system as a whole, 
# while I'm doing it at the country-level; effectively the _relative_ openness
# to the exporting country.
# 

openness <- function(g) {
  require(igraph)
  # take a graph as input
  # get the edgelist (of names, not indices)
  el <- get.edgelist(g)
  # get a vector of the importer for each edge
  importer <- as.list(el[,1])
  # get the per capita GDP for the importer on each edge
  vgdp <- sapply(importer, function(x) V(g)$gdppc[V(g)$name==x])
  # calculate relative trade openness for each import relationship  
  e.weight <- (E(g)$imports / vgdp)
  v.weight <- V(g)$exports  
  # normalize weight scale [0, 1] and return it to the graph as the edge weight
  e.max <- max(e.weight, na.rm=TRUE)
  v.max <- max(v.weight, na.rm=TRUE)
  rel.e.weight <- e.weight/e.max
  rod.weights  <- E(g)$weight
  rod.weights[is.na(rod.weights)] <- 0
  rel.e.weight[is.na(rel.e.weight)] <- rod.weights[is.na(rel.e.weight)]
  E(g)$weight <- rel.e.weight
  E(g)$weight[is.na(E(g)$weight)] <- 0              # still not there? zero
  V(g)$weight <- v.weight/v.max
  return(g)
}

ptm <- proc.time() # Just how long *does* it take to do this? 
# edge weights as openness
# two <- graphs[1:2]
# wtwo <- llply(two, openness)

wgraphs <- llply(graphs, openness, .parallel=TRUE)
names(wgraphs) <- graph.names

proc.time() - ptm  # A good chunk of time, actually; about 2 min on my laptop

nameGraphs <- function(g, n) {
  name <- paste0("Relative Trade Openness in ", substring(n, 2, 5))
  g <- set.graph.attribute(g, "name", name)
  return(g)
}

wgraphs <- mapply(nameGraphs, wgraphs, graph.names, SIMPLIFY=FALSE)

## I also need a set of subgraphs where edge weights are > 0
getSubGraphs <- function(g) {
  sg <- delete.edges(g, which(E(g)$weight < 0.0001))
  return(sg)
}
sgraphs <- llply(wgraphs, getSubGraphs, .parallel=TRUE)
# names(sgraphs) <- paste0("sg", substring(graph.names,2,5))
names(sgraphs) <- graph.names

stopCluster(cl)
## Now, we have a list of weighted-edge graphs ready for analysis. 
