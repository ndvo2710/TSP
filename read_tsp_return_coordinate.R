read_tsp_return_coordinate <- function(x){
  lines <- readLines(x)
  
  metadata <- c(grep(pattern = "[[:alpha:]]" ,lines), which(lines==""))
  lines <- lines[-metadata]
  lines <- sub("^[[:space:]]*", "", lines)
  coordinate <- strsplit(paste(lines, collapse = " "), "[[:space:]]+")[[1]]
  coordinate <- as.numeric(coordinate)
  coordinate <- matrix(coordinate, ncol =3, byrow= TRUE)
  return coordinate
}