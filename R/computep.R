computep <-
function(bigN,n,m,g){
  #Expressions (10) and (11) of Rosenbaum and Krieger (1995), JASA, 85, 493-498
  denom<-zeta(bigN,n,m,g)
  p1<-g*zeta(bigN-1,n-1,m-1,g)/denom
  p0<-zeta(bigN-1,n-1,m,g)/denom
  p11<-g*g*zeta(bigN-2,n-2,m-2,g)/denom
  p10<-g*zeta(bigN-2,n-2,m-1,g)/denom
  p01<-p10
  p00<-zeta(bigN-2,n-2,m,g)/denom
  list(p1=p1,p0=p0,p11=p11,p10=p10,p01=p01,p00=p00)
}
