read_tsp_return_coordinate <- function(x){
  lines <- readLines(x])
  data_start_point <- (grep("NODE_COORD_SECTION",lines)+1)
  data_end_point <- (grep("EOF",lines) -1)
  lines<- lines[data_start_point : data_end_point]
  lines <- sub("^[[:space:]]*", "", lines)
  mycoordinate <- strsplit(paste(lines, collapse = " "), "[[:space:]]+")[[1]]
  mycoordinate <- as.numeric(mycoordinate)
  mycoordinate <- matrix( mycoordinate, ncol = 3, byrow= TRUE)
  mycoordinate
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

tmp = "~/Downloads/ulysses22.tsp"
coordinate <- read_tsp_return_coordinate(tmp)
tour_length <- calculateGeoDistance(coordinate)
message("Tour length of the Hamilton graph in your file is: ",tour_length)