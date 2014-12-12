## Reassemble the data into the country-year dataframe with the community data
# 
# Recall that cya was a dataframe split by year to become a list of dataframes.
# We need to create a dataframe for each year and pull them together by year.
library("plyr")       # Handy way to split, process and recombine data
library("igraph")     # For pulling info back out of the annual trade graphs

# Here's a list of things we've looked at:
# circlePlots   #  
# icorelist     # 
# ocorelist     # 
# sicorlist     # 
# socorlist     # 
# imclist       # 
# simclist      # 
# sblocks       # 
# wgraphs       # 
# sgraphs       # 
# wugraphs      #  
# sugraphs      # 
# wsclusters    # 
# wwclusters    # 
# les           # 
# wcs           # 

getGraphYear <- function(g) {
  #gn <- graph.attribute("name")
  year <- substring(g,2,5)
  return(year)
}


makeLECs <- function(lec) { 
  lec.df <- data.frame(membership(lec))  
  names(lec.df) <- "lec"
  lec.df$ccode <- row.names(lec.df)
  row.names(lec.df) <- NULL
  return(lec.df)
}

makeWCs <- function(wc) {
  wc.df <- data.frame(membership(wc))
  names(wc.df) <- "wc"
  wc.df$ccode <- row.names(wc.df)
  row.names(wc.df) <- NULL
  return(wc.df)  
}

addYear <- function(d,y) {
  d$year <- y
  df <- data.frame(d)
  return(df)
}

wcs.dfs <- llply(wcs, makeWCs)
les.dfs <- llply(les, makeLECs)

adf <- cbind("g.names"=graph.names, 
             "year"=getGraphYear(graph.names),
             "wcs"=wcs.dfs, 
             "lecs"=les.dfs)

lesdf <- unsplit(les.dfs)

year <- getGraphYear(graph.names)

wcs.y <- mapply(addYear, wcs.dfs, year, SIMPLIFY=FALSE)
les.y <- mapply(addYear, les.dfs, year, SIMPLIFY=FALSE)

wcsdf <- rbind.fill(wcs.y)
lesdf <- rbind.fill(les.y)

adf <- join(wcsdf, lesdf, type="left", match="first")

cy2 <- rbind.fill(cya)

wsa <- join(cy2, adf, type="left", match="first", by=c("ccode", "year"))

## The world systems analysis dataframe (wsa) is nearly complete. Let's use 
# these community measures, along with openness and per capita gdp to make
# a simple core vs periphery (not core) determination. 
# 
# I propose that core countries are in walktrap community 1, are in the 
# leading eigenvector community 1, have above-average per capita GDP and above-
# average trade openness. ...so that new column needs to be calculated first.

wsa$openness <- wsa$imports / wsa$gdppc

wsa$core <- FALSE

makeCore <- function(d){
  med.gdppc <- median(d$gdppc)
  med.open  <- median(d$openness)
  for (i in seq_along(d$core)) {
    d$core[i] <- if(((d$wc[i] + d$lec[i]) < 3) & (d$gdppc[i] > med.gdppc) & (d$openness[i] > med.open)) {TRUE} else {FALSE}
  }
  return(d)
}

wsa <- ddply(wsa, .(year), .fun=makeCore)

# Finally, we need to add rows for -8 "anybody not a country" country code for
# all years, 1870-2009

notaindex <- c(1:140)
notayears <- notaindex + 1869 
notwsa <- data.frame(notaindex, notayears)
names(notwsa) <- c("index", "year")
notwsa$ccode <- -8
notwsa$core <- FALSE
notwsa$country <- "not a country"
notwsa$gdppc <- 0
notwsa$imports <- 0
notwsa$exports <- 0
notwsa$wc <- 0
notwsa$lec <- 0
notwsa$openness <- 0
notwsa$abb   <- "XXX"
notwsa <- notwsa[,-1]
#notwsa <- notwsa[,c(2,1,4,5,6,7,8,9,10,11,3,12,13,14,15,16)]

wsa <- rbind(wsa, notwsa)

## fin ########################################################################
