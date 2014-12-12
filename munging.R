## To be run from the 'analysis' script or on it's own, whatever. 
## Data imports ################################################################
#
#  * Build a directed network of countries in the international system based on
#    trade relationships each year (looks like it'll be from 1870, not 1820)
#  * Import conflict dyads
#
## Trade and Economic data for the world systems ###############################
# banks <- read.csv("~/Code/CSS739_Cioffi_CM/Project/Data/World-Systems/1976_Banks_crossNationalTimeSeries/dataBanksCurrent/CNTSDATA.csv")
#
# There are 41 countries in the COW trade data that Maddison refers to in 
# clumps of countries, _eg_ '14 Carribean' countries. I massaged these into
# the source file (albeit hurriedly) so others will need to contact me for data
# or repeat the process themselves (with different results).
mpd <- read.csv("~/Code/CSS739_Cioffi_CM/Project/Data/World-Systems/Maddison_various/mpd-cwd-country-gdppc.csv")
#
# Barbeiri uses a -9 to indicate NA values. That will cause a problem when 
# calculating centrality, so this needs to be handled on import.
dyadTrade <- read.csv("~/Code/CSS739_Cioffi_CM/Project/Data/COW/COW_Trade_3.0/dyadic_trade_3.0.csv", na.strings="-9")
#
natTrade <- read.csv("~/Code/CSS739_Cioffi_CM/Project/Data/COW/COW_Trade_3.0/national_trade_3.0_cwd.csv", na.strings="-9")
#
# Need to replace missing trade data with Oneal & Russett's data
library(Hmisc) # It's in STATA and their text version is whack; this works
rod <- stata.get("~/Data/World-Systems/2001_Oneal-Russett_triangulatingPeace/TRIANGLE.DTA", convert.factors=FALSE)
#
## War Data ####################################################################
cow.extra <- read.csv("~/Code/CSS739_Cioffi_CM/Project/Data/COW/Extra-StateWarData_v4.0.csv", stringsAsFactors=FALSE)
#
cow.intra <- read.csv("~/Code/CSS739_Cioffi_CM/Project/Data/COW/Intra-StateWarData_v4.1.csv", stringsAsFactors=FALSE)
#
cow.inter <- read.csv("~/Code/CSS739_Cioffi_CM/Project/Data/COW/Inter-StateWarData_v4.0.csv", stringsAsFactors=FALSE)
#
cow.non <- read.csv("~/Code/CSS739_Cioffi_CM/Project/Data/COW/Non-StateWarData_v4.0.csv", stringsAsFactors=FALSE)
#
## No more data ################################################################
################################################################################

library("doParallel") # Optional. Only useful a couple of times
library("foreach")    # Required for plyr to work in parallel
library("plyr")       # Handy way to split, process and recombine data
library("dplyr")      # Functionality I wish were in R already
library("reshape")    # For melting horizontal data
library("mice")       # For imputing missing per capita GDP values


## Begin the data manipulation #################################################
# Melt the columns of Maddison's per capita GDP for years 1870-2010 so that we
# have a long data set.
gdp <- melt(mpd, id=1:3, measure=54:194, variable_name="year")

# Take the 'X' off of the year data
gdp$year <- sapply(gdp$year, function(x) substring(x,2,5))

# Rename the columns f
names(gdp) <- c("abb", "ccode", "mpd_country", "year", "gdppc")

# Function to impute missing values and return the third of five permutations
# I picked the third one because I had to pick one and it seemed lazy to pick
# the first one. 
# This is kind of a big deal because it's really not true that
# a linear relationship exists for per capita gdp; like, when there's a war or 
# something. The trade data determines that, however. The country-year data
# just delivers what the dyadic trade data requires; one way or another.

## This process can be done in parallel.
# Set up the multi-core processor to act like a small cluster
# cnodes <- detectCores()
# cl <- makeCluster(cnodes)
# registerDoParallel(cl)

thirdimputation <- function(x) {
  # Don't impute anything if there are no missing values or no data at all.
  if(any(is.na(x$gdppc)) && !all(complete.cases(x))) {
    ic <- complete(mice(x, m=5), 3)
  } else {
    ic <- x
  }
  return(ic)
}

imp.country <- ddply(gdp, .(ccode), thirdimputation)


# Reshape and reorganize into a list of data frames for each year
cyears <- merge(imp.country, natTrade, by=c("ccode", "year"))
cyears <- cyears[,c(1,3,6,2,5,7,8)]

cya <- split(cyears, cyears$year)
cya <- cya[1:140] #cut out 2010; no trade data
#
# Groupo compacto: country-year attributes (cya) is a list of data frames. 


# The next step is to build networks for each year from the dyadic trade data.
# The dyadic trade data has 14 columns, but we only need 4. 
# Chase-Dunn aruges that trade openness should be calculated as the ratio of
# imports to GDP, and weighted by population. The COW dyadic trade data is 
# ultimately agnostic on this point; preference to importer data but exporter
# and IMF data fill gaps. I'm just calculating as IMPORTS. 
#
# ! Note the two pulls from columns 6 and 7 to build directed edges.

dag1 <- dyadTrade[,c(1,2,6,3)] # Reduce dataframe to IMPORT amount & year
dag2 <- dyadTrade[,c(2,1,7,3)] # Reduce dataframe to IMPORT amount & year

n <- c("to", "from", "imports", "year")        # Update names to allow binding
names(dag1) <- n
names(dag2) <- n

dag <- rbind(dag1, dag2) 
dag$to <- as.integer(dag$to)  #merge Hong Kong back into China
dag$from <- as.integer(dag$from)


## Deal with missing data ######################################################
#

rod1 <- rod[,c(1,2,4,3)]
rod2 <- rod[,c(2,1,5,3)]
n <- c("to", "from", "weight", "year")
names(rod1) <- n
names(rod2) <- n
rods <- rbind(rod1, rod2)

## come up with some reasonable missing values for dead data years
fillVoids <- function(d) {
  if((sum(d$weight, na.rm=TRUE)==0) && is.na(sum(d$weight))){
    return(d)
  }
  if(all(complete.cases(d))){
    return(d)
  } else {
    d <- complete(mice(d, m=5),2) 
    return(d)
  }
}
irods <- ddply(rods, .(to, from), fillVoids)

irods <- irods[irods$year > 1919,]
irods <- irods[irods$year < 1950,]

openNorm <- function(x) {
  # x is a dataframe by year  
  omax <- max(x$weight, na.rm=TRUE)
  x$weight <- x$weight / omax
  return(x)
}

nrods <- ddply(irods, .(year), openNorm)

dag <- merge(dag, nrods, by=c("to", "from", "year"), all.x=TRUE)


# Then I split the dataframe into a list of smaller dataframes for each year.
daglist <- split(dag, dag$year)


# Update the names to match and not be numeric  
graph.names <- paste0("g", names(daglist))
names(daglist) <- graph.names
names(cya) <- graph.names

## Data Validation #############################################################
# Not because I'm that particular, but because there's a discontinuity between
# the annual split dataframes for country-year attributes and the edgelists. 
#
# Uncomment the code chunk below to make certain that all the nodes in the 
# edge lists are described in the annual country data. Extra nodes imply an 
# error in the data but will not stop graphs from being built correctly. The
# other way around breaks things.
#
################################################################################
#
checkVertices <- function (x, y){
  elvs <- unique(c(x$to, x$from))
  clvs <- unique(y$ccode)
  print(paste("Edgelist has: ", setdiff(elvs, clvs)))
  print(paste("Country list has: ", setdiff(clvs, elvs)))
  print(paste("Edgelist is bigger by ",length(elvs) - length(clvs)))
}

vertexValidation <- mapply(checkVertices, daglist, cya)
#
################################################################################

# Now we're ready to actually build the graphs...

################################################################################

source("graphs.R")

## And we'll get the war data ready to start working on when we get back to the 
#  analysis part. 

# small cows are calves, so I'll subset it thataway.
# Slice
calf.intra <- cow.intra[,c(1,2,4,5,6,7,11,14,23,26,27)]
calf.inter <- cow.inter[,c(1,2,4,5,9,12,24)]
calf.extra <- cow.extra[,c(1,2,4,5,6,7,10,13,26)]
# Dice
calf.intra <- calf.intra[calf.intra$StartYear1 >= 1870,]
calf.inter <- calf.inter[calf.inter$StartYear1 >= 1870,]
calf.extra <- calf.extra[calf.extra$StartYear1 >= 1870,]

# calf.intra$core.a <- FALSE
# calf.intra$core.b <- FALSE
# calf.extra$core.a <- FALSE
# calf.extra$core.b <- FALSE
# calf.inter$core   <- FALSE




