mscores<-function(y,z,st=NULL,inner=0,trim=3,lambda=.5,tau=0){
  stopifnot(length(tau)==1)
  stopifnot(inner>=0)
  stopifnot(trim>inner)
  stopifnot((lambda>0)&(lambda<1))
  stopifnot(length(y)==length(z))
  if (is.null(st)) st<-rep(1,length(y))
  stopifnot(length(st)==length(y))
  ust<-unique(st)
  nst<-length(ust)
  if (tau!=0) y <- y - z*tau
  adif<-NULL
  for (i in 1:nst){
    who<-st==ust[i]
    yi<-y[who]
    ni<-length(yi)
    if (ni>=2) {
      adifi<-as.vector(abs(outer(yi,yi,"-")[outer(1:ni,1:ni,"<")]))
      adif<-c(adif,adifi)
    }
  }
  sig<-as.vector(stats::quantile(adif,lambda))
  sc<-rep(NA,length(y))
  for (i in 1:nst){
    who<-st==ust[i]
    ni<-sum(who)
    if (ni==1){
      sc[who]<-0
    }
    else if (ni>1) {
      yi<-y[who]
      difi<-outer(yi,yi,"-")
      if ((trim<Inf)) sci<-sign(difi)*trim*pmin(trim-inner,pmax(0,abs(difi/sig)-inner))/(trim-inner)
      else if ((trim==Inf)) sci<- difi
      sc[who]<-apply(sci,1,sum)/(ni-1)
    }
  }
  sc
}
