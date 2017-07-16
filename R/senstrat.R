senstrat <-
function(sc,z,st,gamma=1,alternative="greater",level=0.05,method="BU",detail=FALSE){
  stopifnot((method=="RK")|(method=="BU")|(method=="LS")|(method=="AD"))
  stopifnot(length(level)==1)
  stopifnot((level>0)&(level<1))
  kappa<-stats::qnorm(1-level)
  stopifnot(length(gamma)==1)
  stopifnot(gamma>=1)
  stopifnot((alternative=="greater")|(alternative=="less"))
  stopifnot(all((z==1)|(z==0)))
  stopifnot((min(z)==0)&(max(z)==1))
  stopifnot(length(z)==length(sc))
  stopifnot(length(z)==length(st))
  stopifnot(all(!is.na(sc)))
  stopifnot(all(!is.na(z)))
  stopifnot(all(!is.na(st)))
  if (1==length(unique(as.integer(st)))) warning("You have only one stratum with data.  Use sen2sample instead.")
  stopifnot((detail==TRUE)|(detail==FALSE))
  ustf<-NULL #stratum names, if any
  if (is.factor(st)) {
    ustf<-sort(unique(st))
    st<-as.integer(st)
    names(ustf)<-sort(as.integer(ustf))
  }
  tabst<-table(z,st)
  ust<-as.numeric(colnames(tabst))
  use<-(tabst[1,]>0)&(tabst[2,]>0)
  nstd<-length(use) #number of strata including degenerate strata
  ust<-ust[use]
  nst<-length(ust)
  tabst<-tabst[,use,drop=FALSE]
  skipsome<-sum(use==FALSE)>0 #Are some strata degenerate?
  notused<-length(z)-sum(as.vector(tabst)) #number of observations in degenerate strata
  if (skipsome){
    mess1<-paste("Of ",nstd," strata, ",nstd-nst," strata with ",notused,
                 " individuals contained only treated subjects or only controls.")
    mess1<-c(mess1," These strata do not affect permutation inferences and were removed.")
  }
  else mess1<-paste("All ",nst," strata were used.")
  maxn<-max(tabst[1,]+tabst[2,])
  mu<-matrix(NA,nst,maxn)
  nu<-matrix(NA,nst,maxn)
  if (alternative=="less") sc<- (-sc)
  ts<-rep(NA,nst) #stratum specific test statistics
  ex<-rep(NA,nst) #stratum specific max expectations
  vr<-rep(NA,nst) #stratum specific variances for separable approximation
  exL<-rep(NA,nst) #stratum specific expectations at the linear adjustment
  vrL<-rep(NA,nst) #stratum specific variances at the linear adjustment
  nt<-rep(NA,nst) #stratum specific number treated
  nc<-rep(NA,nst) #stratum specific number control
  m<-rep(NA,nst) #stratum specific u=1's at separable
  mL<-rep(NA,nst) #stratum specific u=1's at linear
  lt<-rep(NA,nst) #linear adjustment
  multiple<-rep(FALSE,nst)

  #Computations for the separable approximation
  for (i in 1:nst){
    who<-st==ust[i]
    sci<-sc[who]
    zi<-z[who]
    BigNi<-length(zi)
    nt[i]<-sum(zi)
    nc[i]<-BigNi-nt[i]
    ts[i]<-sum(zi*sci)
    res<-evall(sci,zi,gamma,method)
    nterms<-dim(res)[1]
    mu[i,1:nterms]<-res$expect
    nu[i,1:nterms]<-res$var
    maxi<-max(res$expect)
    ex[i]<-maxi
    maxs<-res$expect==maxi
    if (sum(maxs)>1) multiple[i]<-TRUE
    vari<-max(res$var[maxs])
    vr[i]<-vari
    m[i]<-(1:(BigNi-1))[(res$expect==maxi)&(res$var==vari)][1]
  }

  teststat<-sum(ts)
  expected<-sum(ex)
  variance<-sum(vr)
  dev<-(teststat-expected)/sqrt(variance)

  #Computations for the Taylor adjustment to the separable approximation
  for (i in 1:nst){
    use<-!is.na(as.vector(mu[i,]))
    mui<-as.vector(mu[i,use])
    nui<-as.vector(nu[i,use])
    adj<-(mui-ex[i])+(kappa*nui/(2*sqrt(variance)))-(kappa*vr[i]/(2*sqrt(variance)))
    mL[i]<-which.max(adj)[1]
    exL[i]<-mui[mL[i]]
    vrL[i]<-nui[mL[i]]
    lt[i]<-adj[mL[i]]
  }


  expectedL<-sum(exL)
  varianceL<-sum(vrL)
  devL<-(teststat-expectedL)/sqrt(varianceL)
  pval<-1-stats::pnorm(dev)
  pvalL<-1-stats::pnorm(devL)
  cvS<-(expected-teststat)+kappa*sqrt(variance)
  cvA<-cvS+sum(lt)
  approxS<-c(cvS,cvA)
  names(approxS)<-c("Separable","Linear Taylor bound")
  if (alternative=="less") {
    teststat<-(-teststat)
    expected<-(-expected)
    dev<-(-dev)
    expectedL<-(-expectedL)
    devL<-(-devL)
  }

  result<-c(pval,dev,teststat,expected,variance,gamma)
  names(result)<-c("P-value","Deviate","Statistic","Expected","Variance","Gamma")

  resultL<-c(pvalL,devL,teststat,expectedL,varianceL,gamma)
  names(resultL)<-c("P-value","Deviate","Statistic","Expected","Variance","Gamma")

  description<-c(nst,sum(nt),sum(nc))
  names(description)<-c("Strata","Treated","Control")

  if (approxS[2]<=0) messmain<-paste("The null hypothesis is rejected at level ",level,
                                     " in the presence of a bias of at most Gamma = ",gamma,".")
  else messmain<-paste("The null hypothesis is not rejected at level ",level,
                       " in the presence of a bias of at most Gamma = ",gamma,".")

  #Do separable approximation and Taylor adjustment agree about rejection?
  if ((cvS<=0)&(cvA>0)) remark<-paste("In testing at level ",level,
    " in the presence of a bias of at most Gamma = ",gamma,
    ", the linear bound and the separable approximation disagree; see their P-values. ",
    " In all cases, it is safe to report the stated conclusion based on the linear bound. ",
    " However, that conclusion may be slightly conservative when the two methods disagree. ")
  else remark<-paste("In testing at level ",level,
                     " in the presence of a bias of at most Gamma = ",gamma,
                     ", the linear bound and the separable approximation agree. ")

  #Report Taylor adjustment alone unless detail is requested.
  if (!detail) list(Conclusion=messmain,Result=resultL,Description=description,StrataUse=mess1)
  else list(Conclusion=messmain,LinearBoundResult=resultL,Separable=result,
            Description=description,Remark=remark,StrataUse=mess1,lambda=approxS)
  }
