
lines <- readLines("~/Downloads/ulysses22.tsp")

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

q1 <-sapply(1:(length(longitude)-1), function(i) cos(longitude[i]-longitude[i+1]))
q2 <-sapply(1:(length(latitude)-1), function(i) cos(latitude[i]-latitude[i+1]))
q3 <-sapply(1:(length(latitude)-1), function(i) cos(latitude[i]+latitude[i+1]))

q12 <- sapply(1:length(q1), function(i) (1+q1[i])*q2[i])
q13 <- sapply(1:length(q1), function(i) (1-q1[i])*q3[i])
q <- sapply(1:length(q1), function(i) 0.5*(q12[i] - q13[i]))            
dij <- as.integer(sapply(1:length(q),function(i) (RRR*acos(q[i]) +1)))

sum(dij)


