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
  return sum(dij)
}