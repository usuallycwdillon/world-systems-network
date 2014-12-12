## Plot-making. Some of this works better than other parts. ...as in, I have NOT
#  mastered the art of Hive Plots. 

library("doParallel") # Optional. Only useful a couple of times
library("foreach")    # Required for plyr to work in parallel
library("plyr")       # Handy way to split, process and recombine data
library("igraph")     # For storing and analysing the annual trade graphs


## Visualisation ###############################################################

# Visualizations are a critical step in understanding and validating data. 

makePlots <- function(g) {
  require(igraph)
  n <- substring(get.graph.attribute(g, "name"), 28, 32)
  out.file.name <- paste0("plots/subGraphs/tradeNetAnimation_p", n, ".png")
  png(out.file.name, width=480, height=480)  
  #   lrt <- layout.reingold.tilford(g)
  #   lr  <- layout.random(g)
  lc  <- layout.circle(g)
  #   lfr <- layout.fruchterman.reingold(g)
  #   lkk <- layout.kamada.kawai(g)
  #   lgl <- layout.lgl(g)
  #   lg  <- layout.graphopt(g)
  
  plot.igraph(g, layout=lc, main=g$name, vertex.label=V(g)$abb,
              vertex.size=4, vertex.label.cex=.5, vertex.label.dist=.2,
              vertex.frame.color=NA, vertex.color="darkgrey",
              edge.arrow.size=0.05, edge.curved=TRUE, edge.color="darkred",
              edge.width=(E(g)$weight * 10))
  # Let the disk catch-up
  Sys.sleep(1)
  dev.off()
  #return(out.file.name)  
}

cnodes <- detectCores()
cl <- makeCluster(cnodes)
registerDoParallel(cl)

circlePlots <- l_ply(sgraphs, makePlots, .parallel=TRUE)

stopCluster(cl)

## Communtity Plots ##
communityPlots <- function(g) {
  require(igraph)
  n <- substring(get.graph.attribute(g, "name"), 28, 32)
  com <- paste0("g",n)  
  c <- wcs$com
  out.file.name <- paste0("plots/communities/community_p", n, ".png")
  #png(out.file.name, width=640, height=480)  
  #   lrt <- layout.reingold.tilford(g)
  #   lr  <- layout.random(g)
  lc  <- layout.circle(g)
  #   lfr <- layout.fruchterman.reingold(g)
  #   lkk <- layout.kamada.kawai(g)
  #   lgl <- layout.lgl(g)
  #   lg  <- layout.graphopt(g)
  
  plot(c, g, layout=lc, main=g$name, vertex.label=V(g)$abb,
              vertex.size=4, vertex.label.cex=.5, vertex.label.dist=.2,
              vertex.frame.color=NA, vertex.color="darkgrey",
              edge.arrow.size=0.05, edge.curved=TRUE, edge.color="darkred",
              edge.width=(E(g)$weight * 10))
  # Let the disk catch-up
  #Sys.sleep(1)
  #dev.off()
  #return(out.file.name)  
}

cnodes <- detectCores()
cl <- makeCluster(cnodes)
registerDoParallel(cl)

circlePlots <- l_ply(sgraphs, makePlots, .parallel=TRUE)

stopCluster(cl)


g <- sugraphs$g1970
c <- wcs$g1970
lc  <- layout.circle(sugraphs$g2008)
plot(wcs$g2005, g, layout=lc, main=g$name, vertex.label=V(g)$abb,
     vertex.size=4, vertex.label.cex=.5, vertex.label.dist=.2,
     vertex.frame.color=NA, vertex.color="darkgrey",
     edge.arrow.size=0.05, edge.curved=TRUE, edge.color="darkred",
     edge.width=(E(g)$weight * 5 ))



## Hive Plots for Community Detection ##########################################
library(HiveR)
library(RColorBrewer)

ig2hp2d <- function(g) {
  imc <- infomap.community(g)
  mm <- imc$membership
  lew <- as.integer(sort(log10(as.integer(E(g)$weight)+1)))+1
  mew <- max(lew) 
  bc <- brewer.pal(mew, "Pastel2")
  
  ndf <- get.data.frame(g, what="vertices")
  rsi <- max(ndf$gdppc, na.rm=TRUE)
  rpo <- max(ndf$imports, na.rm=TRUE)
  id  <- as.integer(ndf$name)
  lab <- ndf$country
  axis <- as.integer(sapply(mm, function(x) if (x<2) 1 else 2))
  ndf$imports[is.na(ndf$imports)] <- min(ndf$imports, na.rm=TRUE)/2
  radius <- ndf$imports/rpo
  size <- ndf$gdppc/rsi 
  ncolor <- as.character(brewer.pal(max(mm), "Pastel1")[mm])
  nodes <- data.frame(id=id, lab=lab, axis=axis, radius=radius, 
                      size=size, color=ncolor, stringsAsFactors=FALSE)
  
  edf <- get.data.frame(g, what="edges")
  maxi <- max(edf$imports)
  id1 <- as.integer(edf$from)
  id2 <- as.integer(edf$to)
  weight <- edf$weight
  color  <- as.character(bc[lew])
  edges <- data.frame(id1=id1, id2=id2, weight=weight, color=color,
                      stringsAsFactors=FALSE)
  
  htype <- "2D"
  desc <- as.character(ndf$year[1])
  axis.cols <- "black" #brewer.pal(3, "Pastel2")
  h <- list(nodes=nodes, edges=edges, type=htype, desc=desc, axis.cols=axis.cols)
  class(h) <- "HivePlotData"
  
  return(h)
}

hptest <- ig2hp2d(sgraphs$g1919)
chkHPD(hptest)
plotHive(hptest, ch=1, bkgnd="white")
