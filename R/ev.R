ev <-
  function(sc,z,m,g,method){
    #Computes the expectation and variance for one stratum with m 1's in u
    #sc are scores, z is treated/control
    #large determines whether exact moments are computed or a large
    #sample approximation is used.  Large sample moments are discussed
    #on page 143 of Rosenbaum (2002) Observational Studies (2nd ed.)
    bigN<-length(z)
    n<-sum(z)
    stopifnot((bigN>n)&(bigN>m)&(1<=n)&(1<=m))
    q<-sort(sc)
    u<-c(rep(0,bigN-m),rep(1,m)) #worst case u
    if (method=="AD"){
      if ((bigN>200)&(n>50)&((bigN-n)>50)&(m>50)&((bigN-m)>50)) method<-"LS"
      else method<-"BU"
    }

    if ((method=="LS")|(method=="BU")){
      if (method=="LS"){ #Large sample, Section 4.6.4 of Rosenbaum (2002) Observational Studies
        quad<-function(ets){
          (ets*ets*(g-1))+(g*n*m)-ets*((g-1)*(n+m)+bigN)
        }
        et<-stats::uniroot(quad,lower=max(0,m+n-bigN),upper=min(n,m))$root
        #expression (4.30), Rosenbaum (2002)
        vt<-1/((1/et)+(1/(m-et))+(1/(n-et))+(1/((bigN+et)-(m+n))))
      }
      else{#Exact moments, Proposition 20, page 155 of Rosenbaum (2002) Observational Studies
        et<-BiasedUrn::meanFNCHypergeo(m,bigN-m,n,g)
        vt<-BiasedUrn::varFNCHypergeo(m,bigN-m,n,g)
      }
      qb1<-sum(q[u==1])/m
      qb0<-sum(q[u==0])/(bigN-m)
      if (m==1) w1<-0
      else w1<-stats::var(q[u==1])
      if ((bigN-m)==1) w0<-0
      else w0<-stats::var(q[u==0])
      expect<-et*qb1+(n-et)*qb0 #expression (4.31), page 143 of Rosenbaum (2002)
      term<-(w1-w0)*et
      term<-term-(et*et+vt)*((w1/m)+(w0/(bigN-m)))
      term<-term+n*(bigN+2*et-(m+n))*w0/(bigN-m)
      vari<-term+vt*((qb1-qb0)^2) #expression (4.32), page 143 of Rosenbaum (2002)
    }
    else{ #Use exact moments calculated as in Rosenbaum and Krieger (1990)
      o<-computep(bigN,n,m,g)
      expect<-sum(((o$p1*u)+(o$p0*(1-u)))*q)
      cv<-(outer(u,u,"*")*(o$p11-o$p1*o$p1))+(outer(u,1-u,"*")*(o$p10-o$p1*o$p0))+
        (outer(1-u,u,"*")*(o$p10-o$p0*o$p1))+(outer(1-u,1-u,"*")*(o$p00-o$p0*o$p0))
      diag(cv)<-((o$p1-o$p1^2)*u)+((o$p0-o$p0^2)*(1-u)) #pi=pii
      vari<-0
      for (i in 1:bigN){
        for (j in 1:bigN){
          vari<-vari+q[i]*cv[i,j]*q[j]
        }
      }
    }
    list(expect=expect,vari=vari)
  }


