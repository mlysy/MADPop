#Internal function
# sorts a matrix by first column, breaking ties with second column, breaking those ties with 3rd, etc.
colSort <- function(X) {
  if(!is.matrix(X)) return(sort(X))
  X[do.call(order, lapply(data.frame(X), c)),]
}
