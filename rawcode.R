
lines <- readLines("burma14.tsp")

metadata <- c(grep(pattern ="[[:alpha:]]" ,lines), which(lines==""))
index = intersect(metadata,e_lines)
lines <- lines[-metadata]
lines <- sub("^[[:space:]]*", "", lines)
coordinate <- strsplit(paste(lines, collapse = " "), "[[:space:]]+")[[1]]
coordinate <- as.numeric(coordinate)
coordinate <- matrix(coordinate, ncol =3, byrow= TRUE)

x <- coordinate[,2]
degx <-as.integer(x)
minx <- (x-as.integer(x))
latitude <- pi*(degx + 5*minx/3)/180

y <- coordinate[,3]
degy <-as.integer(y)
miny <- (y-as.integer(y))
longitude <- pi*(degy + 5*miny/3)/180

RRR <- 6378.388



neighbor_tour <- function(){
  
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


}

test = 1
neighbor = test
neighbors_distance = 0
train = 2:length(latitude)
while(length(train) > 1){
  index = closest_distance(train,test)
  test = index[2]
  neighbor <- c(neighbor,test)
  neighbors_distance <- c(neighbors_distance,neighbor[1])
  train <- train[-which(train = test)]
}




distance <- function(i,j){
  q1 <-cos(longitude[i]-longitude[j])
  q2 <-cos(latitude[i]-latitude[j])
  q3 <-cos(latitude[i]+latitude[j])
  
  q12 <- (1+q1)*q2
  q13 <- (1-q1)*q3
  q <- 0.5*(q12 - q13)            
  dij <- as.integer(RRR*acos(q) +1)
  return (dij)}

closest_distance <- function(train, test){
  distanceSet =NA
  for(i in train){
    distanceSet[i] = distance(test,i)
  }
  min_value = min(distanceSet,na.rm = TRUE)
  index = which(distanceSet == min_value)
  min_vector <- c(min_value,index)
  return (min_vector)}


