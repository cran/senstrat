evall <-
function(sc, z, g, method){
  #Computes ev for all m
  bigN<-length(z)
  stopifnot(g>=1)
  stopifnot(length(sc)==bigN)
  o<-matrix(NA,bigN-1,3)
  for (i in 1:(bigN-1)){
    res<-ev(sc,z,i,g, method)
    o[i,1]<-i
    o[i,2]<-res$expect
    o[i,3]<-res$vari
  }
  colnames(o)<-c("m","expect","var")
  as.data.frame(o)
}
