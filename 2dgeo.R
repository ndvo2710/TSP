  
  setwd("~/Mat168/allTSP")
  directory = "~/Mat168/allTSP"
  fileList = list.files(path = directory, pattern = ".tsp")
  data <- sapply(fileList,readLines)
  index <- 1:111
  for(i in 1:length(data)){
    if ( length(grep(pattern = "EDGE_WEIGHT_TYPE : EUC_2D|EDGE_WEIGHT_TYPE: GEO",data[[i]]))  == 0) index[i]= 0 else index[i]=i
  }
  selectedIndex <- index[index>0] 
  selectedData <- data[selectedIndex]
  
  
  
  
  sapply(1:length(selectedData),function(i) grep(pattern = "NODE_COORD_SECTION", selectedData[[i]]))        
                  
  
  
  
  lines <- readLines("ulysses22.tsp")
  metadata <- c(grep(pattern = "[[:alpha:]]" ,lines), which(lines==""))
  lines <- lines[-metadata]
  lines <- sub("^[[:space:]]*", "", lines)
  mycoordinate <- strsplit(paste(lines, collapse = " "), "[[:space:]]+")[[1]]
  mycoordinate <- as.numeric(mycoordinate)
  mycoordinate <- matrix( mycoordinate, ncol = 3, byrow= TRUE)
  return (mycoordinate)}

calculateGeoDistance <- function(mycoordinate){
  x <- mycoordinate[,2]
  deg <-round(x)
  min <- (x-round(x))
  latitude <- pi*(deg + 5*min/3)/180
  
  y <- mycoordinate[,3]
  deg <-round(y)
  min <- (y-round(y))
  longitude <- pi*(deg + 5*min/3)/180
  
  RRR <- 6378.388
  
  q1 <-sapply(1:(length(longitude)-1), function(i) cos(longitude[i]-longitude[i+1]))
  q2 <-sapply(1:(length(latitude)-1), function(i) cos(latitude[i]-latitude[i+1]))
  q3 <-sapply(1:(length(latitude)-1), function(i) cos(latitude[i]+latitude[i+1]))
  
  q <- sapply(1:length(q1), function(i) (0.5*((1+q1[i])*q2[i] - (1-q1[i])*q3[i])))
  
  dij <- round(sapply(1:length(q),function(i) (RRR*acos(q[i]) +1)))
  return (sum(dij))}
