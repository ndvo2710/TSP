install.packages("gtools")

#The Assignment P2
# Student: Nguyen Vo
# Student ID: 998153069

ptm <- proc.time() #Set starting time
library(gtools)

#Set working directory
setwd("~/Mat168/allTSP")
directory = "~/Mat168/allTSP"

#read_type function to read type of data (EUC_2D or GEO)
read_type <- function(x){
  type =""
  lines <- readLines(x)
  if (length(grep(pattern = "EDGE_WEIGHT_TYPE : EUC_2D|EDGE_WEIGHT_TYPE: GEO",lines))  == 0) type ="NA" else{
    if( length(grep(pattern = "EDGE_WEIGHT_TYPE : EUC_2D",lines))  == 0) type = "GEO" else type = "EUC_2D"
  } 
  return (type)}

#read_tsp_return_coordinate function to read data and return coordinate
read_tsp_return_coordinate <- function(x){
  lines <- readLines(x)
  data_start_point <- (grep("NODE_COORD_SECTION",lines)+1)
  if (length(grep("EOF", lines)) ==0){
    data_end_point <- length(lines)
  }
  else{
    data_end_point <- (grep("EOF",lines) -1)
  }
  lines<- lines[data_start_point : data_end_point]
  lines <- sub("^[[:space:]]*", "", lines)
  mycoordinate <- strsplit(paste(lines, collapse = " "), "[[:space:]]+")[[1]]
  mycoordinate <- as.numeric(mycoordinate)
  mycoordinate <- matrix( mycoordinate, ncol = 3, byrow= TRUE)
  mycoordinate
  return (mycoordinate)}

# function to read .opt.tour file
read_opt_tour <- function(x){
  lines <- readLines(x)
  lines <- sub("^[[:space:]]*", "", lines)
  lines <- strsplit(paste(lines, collapse = " "), "[[:space:]]+")[[1]]
  if (lines[length(lines)]=="EOF") lines = lines[-length(lines)]
  if (lines[length(lines)]==-1) lines[length(lines)] =1
  data_start_point <- (grep("TOUR_SECTION",lines)+1)
  data_end_point <- length(lines)
  tour <- lines[data_start_point : data_end_point]
  tour <- as.numeric(tour)
  return (tour)
}

# function to calculate distance in EUC_2D
distance_euc2d <-function(i,j){
  xd = x[i] - x[i]
  yd = y[i] - y[j]
  dij = round(sqrt(xd*xd + yd*yd))
  return (dij)}

#function to calculate distance in GEO
distancegeo <- function(i,j){
  q1 <-cos(longitude[i]-longitude[j])
  q2 <-cos(latitude[i]-latitude[j])
  q3 <-cos(latitude[i]+latitude[j])
  
  q12 <- (1+q1)*q2
  q13 <- (1-q1)*q3
  q <- 0.5*(q12 - q13)            
  dij <- as.integer(RRR*acos(q) +1)
  return (dij)}


getHamilton_geo_length <- function(n){
  n = N
  tmp_hamilton_length = 0
  pseudo_index = c(1:n,1)
  for(i in (1:n)){
    tmp_hamilton_length = tmp_hamilton_length + distancegeo(pseudo_index[i],pseudo_index[i+1])
  }
  return (tmp_hamilton_length)}

getHamilton_euc2d_length <- function(n){
  tmp_hamilton_length = 0
  pseudo_index = c(1:n,1)
  for(i in (1:n)){
    tmp_hamilton_length = tmp_hamilton_length + distance_euc2d(pseudo_index[i],pseudo_index[i+1])
  }
  return (tmp_hamilton_length)}



#nearest neighbor function : between one test point and a train set
closest_distance_geo <- function(train, test){
  distanceSet =NA
  for(i in train){
    distanceSet[i] = distancegeo(test,i)
  }
  min_value = min(distanceSet,na.rm = TRUE)
  index = which(distanceSet == min_value)
  min_vector <- c(min_value,index)
  return (min_vector)}

closest_distance_euc2d <- function(train, test){
  distanceSet =NA
  for(i in train){
    distanceSet[i] = distance_euc2d(test,i)
  }
  min_value = min(distanceSet,na.rm = TRUE)
  index = which(distanceSet == min_value)
  min_vector <- c(min_value,index)
  return (min_vector)}





getNearest_Neighbor_tour_geo <-function(n){
  index =""
  test = 1
  neighbor = test
  neighbors_distance = 0
  train = 2:n
  while(length(train) != 0){
    index = closest_distance_geo(train,test)
    test = index[2]
    neighbor <- c(neighbor,test)
    neighbors_distance <- c(neighbors_distance,index[1])
    train <- train[-which(train == test)]
  }
  neighbor <- c(neighbor, 1)
  neighbors_distance <- c(neighbors_distance, distancegeo(index[2],1))
  tour <- data.frame(neighbor,neighbors_distance)
  return (tour)}

getNearest_Neighbor_tour_euc2d <-function(n){
  index =""
  test = 1
  neighbor = test
  neighbors_distance = 0
  train = 2:n
  while(length(train) != 0){
    index = closest_distance_euc2d(train,test)
    test = index[2]
    neighbor <- c(neighbor,test)
    neighbors_distance <- c(neighbors_distance,index[1])
    train <- train[-which(train == test)]
  }
  neighbor <- c(neighbor, 1)
  neighbors_distance <- c(neighbors_distance, distance_euc2d(index[2],1))
  tour <- data.frame(neighbor,neighbors_distance)
  return (tour)}



# Main Program
##################This is Part (1) ###################
#reads a symmetric TSP in the input format TSPLIB
fileList <- list.files(directory)
fileTsp <- list.files(directory, pattern = "\\.tsp$")
fileOptTour <- list.files(directory, pattern = "\\.opt.tour$")
fileOptTour <- gsub(".opt.tour", "", fileOptTour)
fileOptTour <- paste(fileOptTour,".tsp", sep = "")

indexSelectedFiles = rep(0, length(fileTsp))
for (i in (1:length(fileTsp))){
  if (read_type(fileTsp[i]) == "GEO" |read_type(fileTsp[i]) == "EUC_2D") indexSelectedFiles[i] = i
  else indexSelectedFiles[i] = NA}
selectedFiles <- fileTsp[!is.na(indexSelectedFiles)]

resultTable_Hamilton = rep(0,length(selectedFiles))
resultTable_Nearest =  rep(0,length(selectedFiles))

for (i in (1:length(selectedFiles))){
  tmp = selectedFiles[i]
  coordinate <- read_tsp_return_coordinate(tmp)
  #Latitude and Longitude
  x <- coordinate[,2]
  degx <-as.integer(x)
  minx <- (x-as.integer(x))
  latitude <- pi*(degx + 5*minx/3)/180
  
  y <- coordinate[,3]
  degy <-as.integer(y)
  miny <- (y-as.integer(y))
  longitude <- pi*(degy + 5*miny/3)/180
  
  N <- length(x)
  RRR <- 6378.388
  
  if (read_type(tmp) == "GEO"){
    resultTable_Hamilton[i] = getHamilton_geo_length(N)
  }else{
    resultTable_Hamilton[i] = getHamilton_euc2d_length(N)
  }
  
  if (N >5000){
    resultTable_Nearest[i] = -1
  }else{
    if (read_type(tmp) == "GEO"){
      nearest_tour = getNearest_Neighbor_tour_geo(N)
      resultTable_Nearest[i] = sum(nearest_tour[,2])
    }else{
      nearest_tour = getNearest_Neighbor_tour_euc2d(N)
      resultTable_Nearest[i] = sum(nearest_tour[,2])
    }
  }
  rm(x,degx,minx,latitude,y,degy,miny,longitude,N)
}

result <- data.frame(selectedFiles,resultTable_Hamilton, resultTable_Nearest)



indexSelectedFiles = rep(0, length(fileOptTour))
for (i in (1:length(fileOptTour))){
  if (read_type(fileOptTour[i]) == "GEO" |read_type(fileOptTour[i]) == "EUC_2D") indexSelectedFiles[i] = i
  else indexSelectedFiles[i] = NA}
selectedFilesOptTourTsp <- fileOptTour[!is.na(indexSelectedFiles)]
selectedFilesOptTour <- gsub(".tsp","",selectedFilesOptTourTsp)
selectedFilesOptTour <- paste(selectedFilesOptTour,".opt.tour",sep ="")


resultOpt <- rep(0,length(selectedFilesOptTourTsp) )
for (i in (1:length(selectedFilesOptTourTsp))){
  tmpTsp = selectedFilesOptTourTsp[i]
  tmpTour = selectedFilesOptTour[i]
  coordinate <- read_tsp_return_coordinate(tmpTsp)
  tour <- read_opt_tour(tmpTour)
  
  #Latitude and Longitude
  x <- coordinate[,2]
  degx <-as.integer(x)
  minx <- (x-as.integer(x))
  latitude <- pi*(degx + 5*minx/3)/180
  
  y <- coordinate[,3]
  degy <-as.integer(y)
  miny <- (y-as.integer(y))
  longitude <- pi*(degy + 5*miny/3)/180
  
  N <- length(x)
  RRR <- 6378.388
  
  if (read_type(tmpTsp) == "GEO"){
    optimalLength = 0
    for (j in (1:(length(tour)-1))){
      optimalLength = optimalLength + distancegeo(tour[j], tour[j+1])
    } 
  }else{
    optimalLength = 0
    for (j in (1:(length(tour)-1))){
      optimalLength = optimalLength + distance_euc2d(tour[j], tour[j+1])
    } 
  }
  resultOpt[i] = optimalLength
  rm(x,degx,minx,latitude,y,degy,miny,longitude,N,optimalLength)
}

resultOpt <- data.frame(selectedFilesOptTourTsp,resultOpt)

compareTable <- data.frame(result[match(resultOpt[,1],result[,1]),],resultOpt[,2])
names(compareTable) <- c("FILE","Hamilton Length","NN length","Optimal Length")
names(result) <- c("FILE","Hamilton Length","NN length")

#Since there are more .tsp files than .opt.tour files, I create two table
# The first table shows the result of Hamilton length and Nearest Neighbor Length for all .tsp files
# Notice that there are some (-1)s in the value of column "NN length".
# This doesn't mean that my program generate the wrong results, but because those files contain 
# a large number of cities that could take forever for R to process e.g around more than 90 millions of loops
# for about 14000 cities. Therefore I try to choose the cut-off value as 5000 cities, if the number of cities
# is larger than 5000, NN length return value -1.
# Here is my first table:
result

# The second table insert another column for the Optimal Length which get the tour results from .opt.tour file 
# Notice the second table is shorter than the first table. Since there are only 19 .opt.tour files are "GEO" and "EUC_2D"
# Here is my 2nd table:
compareTable


##################################################
##################################################
##################################################
################# BURMA14.TSP ####################

burma = "burma14.tsp" # working file address
coordinate <- read_tsp_return_coordinate(burma) #read coordinates
type <- read_type(burma) # read type

#Latitude and Longitude
x <- coordinate[,2]
degx <-as.integer(x)
minx <- (x-as.integer(x))
latitude <- pi*(degx + 5*minx/3)/180

y <- coordinate[,3]
degy <-as.integer(y)
miny <- (y-as.integer(y))
longitude <- pi*(degy + 5*miny/3)/180

RRR <- 6378.388
#################This is Part (2)###################  
# creates a text file with the computed distances 
index_print = combinations(length(x),2,v = 1:length(x),repeats.allowed=FALSE)
index_print <- rbind(index_print, c(14,1))
if (type == "GEO"){
  txt <- sapply(1:nrow(index_print), function(i) distancegeo(index_print[i,1],index_print[i,2]))
} else {
  txt <- sapply(1:nrow(index_print), function(i) distance_euc_2d(index_print[i,1],index_print[i,2]))
}
txt <- data.frame(index_print, txt)
names(txt) <- c("i","j","dij")
txt
write.table(txt,"S:/burma14cities.txt", sep = "\t")



################This is Part(3)######################
#computes and prints the length of the tour (Hamiltonian cycle) 1-2-3-...-n-1
N <- length(x)
if (type == "GEO"){
  hamilton_length = getHamilton_geo_length(N)
}else{
  hamilton_length = getHamilton_euc2d_length(N)
}
message("The length of the tour(Hamiltonian cycle) from node 1 to ", length(x)," is: ",hamilton_length)


################This is Part(4)######################
# runs the nearest neighbor algorithm, starting at 
#city (vertex) 1; break ties in favor of the smallest-index city
#Print the resulting tour and its length
if (type == "GEO"){
  burma_tour = getNearest_Neighbor_tour_geo(N)
}else{
  burma_tour = getNearest_Neighbor_tour_euc2d(N)
}
names(burma_tour) <- c("Node","Distance of Node[i-1,i]")
burma_tour # Print the tour and distance.
message("The total length of the nearest_neighbor tour is: ",sum(burma_tour[,2]))

##############This is Part 5######################
#prints the total computation time of your program.
proc.time() -ptm
