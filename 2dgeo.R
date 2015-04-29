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
    if (length(grep(pattern = "EDGE_WEIGHT_TYPE : EUC_2D|EDGE_WEIGHT_TYPE: GEO",lines))  == 0) print("This is not type EUC_2D or GEO") else{
      if( length(grep(pattern = "EDGE_WEIGHT_TYPE : EUC_2D",lines))  == 0) type = "GEO" else type = "EUC_2D"
    } 
    return (type)}

#read_tsp_return_coordinate function to read data and return coordinate
  read_tsp_return_coordinate <- function(x){
    lines <- readLines(x)
    data_start_point <- (grep("NODE_COORD_SECTION",lines)+1)
    data_end_point <- (grep("EOF",lines) -1)
    lines<- lines[data_start_point : data_end_point]
    lines <- sub("^[[:space:]]*", "", lines)
    mycoordinate <- strsplit(paste(lines, collapse = " "), "[[:space:]]+")[[1]]
    mycoordinate <- as.numeric(mycoordinate)
    mycoordinate <- matrix( mycoordinate, ncol = 3, byrow= TRUE)
    mycoordinate
    return (mycoordinate)}
  
# function to calculate distance in EUC_2D
  distance_euc_2d <-function(i,j){
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
 

#nearest neighbor function : between one test point and a train set
  closest_distance <- function(train, test){
    distanceSet =NA
    if (type == "GEO"){
      for(i in train){
        distanceSet[i] = distancegeo(test,i)
      }
    }
    else {
      for(i in train){
        distanceSet[i] = distance_euc_2d(test,i)
      }
    }
    min_value = min(distanceSet,na.rm = TRUE)
    index = which(distanceSet == min_value)
    min_vector <- c(min_value,index)
    return (min_vector)}
  
  
  
# Main Program
##################This is Part (1) ###################
#reads a symmetric TSP in the input format TSPLIB

  tmp = "burma14.tsp" # working file address
  coordinate <- read_tsp_return_coordinate(tmp) #read coordinates
  type <- read_type(tmp) # read type
  
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

  if (type == "GEO"){
    table <- sapply(1:nrow(index_print), function(i) distancegeo(index_print[i,1],index_print[i,2]))
  } else {
      table <- sapply(1:nrow(index_print), function(i) distance_euc_2d(index_print[i,1],index_print[i,2]))
  }
  table <- data.frame(index_print, table)
  names(table) <- c("i","j","dij")
  table
  write.table(table,"S:/burma14cities.txt", sep = "\t")
    


################This is Part(3)######################
#computes and prints the length of the tour (Hamiltonian cycle) 1-2-3-...-n-1
  hamilton_length = 0
  for(i in 1:(length(latitude)-1)){
    hamilton_length = hamilton_length + distancegeo(i,i+1)
  }
  message("The length of the tour(Hamiltonian cycle) from node 1 to ", length(x)," is: ",hamilton_length)


################This is Part(4)######################
# runs the nearest neighbor algorithm, starting at 
#city (vertex) 1; break ties in favor of the smallest-index city
#Print the resulting tour and its length
    index =""
    test = 1
    neighbor = test
    neighbors_distance = 0
    train = 2:length(latitude)
    while(length(train) != 0){
      index = closest_distance(train,test)
      test = index[2]
      neighbor <- c(neighbor,test)
      neighbors_distance <- c(neighbors_distance,index[1])
      train <- train[-which(train == test)]
    }
    tour <- data.frame(neighbor,neighbors_distance)
  
  tour # Print the tour and distance.
  message("The total length of the nearest_neighbor tour is: ",sum(tour[,2]))

##############This is Part 5######################
#prints the total computation time of your program.
proc.time() -ptm
 