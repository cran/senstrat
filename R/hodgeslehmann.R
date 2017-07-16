hodgeslehmann<-function(y,z,st,align="median",tau=0){
  stopifnot((align=="huber")|(align=="median")|(align=="mean")|(align=="hl"))
  stopifnot(length(tau)==1)
  stopifnot(length(y)==length(st))
  stopifnot(length(z)==length(y))
  stopifnot(all(!is.na(y)))
  stopifnot(all(!is.na(st)))
  stopifnot(all(!is.na(z)))
  if (is.factor(st)) st<-as.integer(st)
  ust<-sort(unique(st))
  nst<-length(ust)
  if (tau!=0) y<-y-tau*z
  sc<-rep(NA,length(y))

  for (i in 1:nst){
    who<-(st==ust[i])
    yi<-y[who]
    if (length(yi)==1) {
      vi<-0
      sc[who]<-vi
    }
    else{
      if (align=="huber") ctr<-MASS::huber(yi)$mu
      else if (align=="median") ctr<-stats::median(yi)
      else if (align=="hl") suppressWarnings(ctr<-stats::wilcox.test(yi,conf.int=TRUE)$estimate)
      else ctr<-mean(yi)
      vi<-yi-ctr
      sc[who]<-vi
    }
  }
  sc<-rank(sc)
  sc
}
