sen2sample <-
function(sc,z,gamma=1,alternative="greater", method="BU"){
  stopifnot((method=="RK")|(method=="BU")|(method=="LS")|(method=="AD"))
  stopifnot(length(gamma)==1)
  stopifnot(gamma>=1)
  stopifnot((alternative=="greater")|(alternative=="less"))
  stopifnot(all((z==1)|(z==0)))
  Nbig<-length(sc)
  stopifnot(length(z)==Nbig)
  if (alternative=="less") sc<- (-sc)
  ts<-sum(z*sc)
  ev<-evall(sc,z,gamma,method)
  devs<-(ts-ev$expect)/sqrt(ev$var)
  m<-which.min(devs)
  dev<-min(devs)
  pval<-1-pnorm(dev)
  ex<-ev$expect[m]
  if (alternative=="less"){
    dev<-(-dev)
    ex<-(-ex)
    ts<-(-ts)
  }
  o2<-c(dev,ts,ex,ev$var[m],gamma)
  if (gamma==1) m<-NA
  o<-c(dev,m,sum(z),sum(1-z))
  names(o)<-c("deviate","m","n.treated","n.control")
  names(o2)<-c("Deviate","Statistic","Expectation","Variance","Gamma")
  list(pval=pval,description=o,detail=o2)
}
