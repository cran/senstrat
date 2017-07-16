zeta <-
function(bigN,n,m,g){
  #The zeta function in expression (8) of Rosenbaum and Krieger (1995), JASA, 85, 493-498.
  low<-max(0,(m+n)-bigN)
  if ((m<0)|(n<0)) o<-0
  else {
    high<-min(n,m)
    a<-low:high
    term<-choose(m,a)*choose(bigN-m,n-a)*(g^a)
    o<-sum(term)
  }
  o
}
