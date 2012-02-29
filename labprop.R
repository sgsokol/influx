# 2011-03-16
# Label propagation in linear chain.
# All turnover rates are supposed to be pairwise different

# Copyright 2011, INRA France
# Author: Serguei SOKOL
# License: GPL2 + citation of a possible corresponding paper.
source("/home/sokol/sysbio/dev/ftbl2sys/opt_cumo_tools.R")
source("/home/sokol/R/Rtools/tools_ssg.R")
lab_chain=function(ti, lambda, crodi=outer(lambda, lambda, "-")) {
   # return vector of the same length as time vector ti
   # with n labeled metabolites where n is the lentgh of
   # turnover rates lambda.
   # crodi is the optional matrix of crossdifferences
   n=length(lambda)
   if (n < 1) {
      return(rep(1., length(ti)))
   }
   mexp=exp(-ti%o%lambda)
   lab=apply(t(1:n), 2, function(i){(1-apply(apply(t(1:i), 2, function(j)
      {prod(lambda[1:i][-j]/crodi[1:i,1:i, drop=F][,j][-j])*mexp[,j]}), 1, sum))})
   return(lab)
   #return(1-apply(apply(t(1:n), 2, function(j)
   #   {prod(lambda[-j]/crodi[,j][-j])*exp(-lambda[j]*ti)}), 1, sum))
}
lab_rev2=function(ti, f, m, r) {
   # return matrix of two columns
   # first column is label1 and second label2
   # as time function at moments ti
   # f: input flux (scalar)
   # m: vector of two metabolite pools
   # r: f2/f1 ratio of revers and forward fluxes [0; 1) (f1-f2=f)
   m1=m[1]; m2=m[2]
   f1=f/(1.-r); f2=r*f1
   
   lam_11=f1/m1; lam_12=f1/m2; lam_21=f2/m1; lam1=f/m1
   delta=lam_11-lam_12; somme=lam_11+lam_12
   racine=sqrt(delta*delta+4*lam_12*lam_21)
   lambda1=(-somme+racine)*0.5; lambda2=(-somme-racine)*0.5
   alpha=(delta+racine)*0.5; bet=(delta-racine)*0.5
   return(lam1/(lambda1-lambda2) * ((1.-exp(lambda1*ti))%o%c(bet,-lam_12)/lambda1 - (1.-exp(lambda2*ti))%o%c(alpha, -lam_12)/lambda2))
}

lab_rev=function(ti, f, m, r) {
   # linear chain of reversible reactions
   # ti: nt time vector
   # f: net flux (scalar)
   # m: n vector of metabolite pools
   # r: n-1 vector of rev/fwd ratio (in [0,1[). fwd-rev=f
   # r[n] is assumed be 0 (final non reversible output reaction)
   #
   # returns a matrix of label fractions as function of time (nt*n)
   
   # fluxes constructions
   n=length(m)
   if (length(r)==n-1) {
      r=c(r, 0.)
   } else if (length(r)==n) {
      r[n]=0.
   } else {
      stop("Wrong r vector dimension")
   }
   # forward output fluxes
   f1=f/(1.-r)
   # revers
   f2=r*f1
   # matrix construction
   a=matrix(0., n, n)
   # main diag
   diag(a)=-f1
   diag(a)[-1]=diag(a)[-1]-f2[-n]
   # subdiag (precedent metabs)
   a[row(a)-col(a)==1]=f1[-n]
   # overdiag (next metabs)
   a[row(a)-col(a)==-1]=f2[-n]
   # scaling by pool size
   a=a/m
   
   # add vector
   b=rep(0., n)
   b[1]=f/m[1]
   
   # eigen decomposition
   ea=eigen(a)
   uinv=solve(ea$vectors)
   
   # stationary solution
   xst=solve(a, -b)
   
   # matrix exponent
   lab=t(apply(as.matrix(ti), 1, function(tit) {ea$vectors%*%((uinv*(1.-exp(ea$values*tit)))%*%xst)}))
   return(lab)
}
lab_netw=function(ti, a, b) {
   # label propagation in a network described by matrix a and
   # constant vector b: dl/dt=a*l+b
   ea=eigen(a)
   uinv=solve(ea$vectors)
   
   # stationary solution
   xst=solve(a, -b)
   
   # matrix exponent
   lab=t(apply(as.matrix(ti), 1,
      function(tit) ea$vectors%*%((uinv*(1.-exp(ea$values*tit)))%*%xst)
   ))
   return(lab)
}
lab_cumo=function(spAbr, li_m2x, muamp_xinp, f, m, nm_incu) {
   # labeling propagation in a cumomer network.
   # The network is determined by a matrix a (a function of fluxes f
   # and metabolite pool vector m) and source term s corresponding to
   # the input and lighter weight cumomers).
   # The list spAbr allows a passage from
   # vectors f and x to matrix a and source term s
   # for each cumomer weight (one weight per list item)
   # The items of the list li_m2x provide a mapping
   # of metabolite pool vector m on cumomer vectors x.
   # The source term is supposed to be causal (0 before time 0)
   # and composed exclusively
   # of exponential terms (possibly complexe in case of
   # sinusoidal input or in a presence of loops in the network)
   # Output:
   # a list with $mu and $amp fields which can be passed to plot_muamp()
   # Data structures:
   # spAbr structure is the same as in influx_s project.
   # muamp${amp,mu} - vectors of amplitudes and exponential mu
   # 2012-02-23 sokol
   
   weights=1:length(spAbr)
   muamp_xinp$amp=as.matrix(muamp_xinp$amp)
   muamp_x=muamp_xinp # we start with input exp's and extend it with each weight
   nb_xw=0
   nb_xl=nrow(muamp_xinp$amp)+1 # +1 is for "one" in incu names
   for (w in weights) {
      cat("w=", w, "dim(xmap)=", dim(muamp_x$amp), "\n", sep=" ")
      if (w == 3) {
         break
      }
      invm=(1./m)[li_m2x[[w]]] # result has length of xw
      sp=spAbr[[w]]
      nb_xl=nb_xl+nb_xw # add precedet weights to get ligther cumo size
      nb_xw=dim(sp$tA)[1]
      
      # construct the a matrix of this weight
      a=fwrv2Abr(f, sp, NULL, nm_incu[nb_xl+(1:nb_xw)], getb=F)$A*invm

      # construct the sources for this weight
      ## product terms, i.e. condensation of lighter cumomers or
      ## product with 1 when mono input cumomer
      ## each exponent is added to each other
      crosum=outer(muamp_x$mu, muamp_x$mu, "+")
      # ilow is indexes of underdiagonal terms
      ilow=row(crosum)>col(crosum)
      muamp_s=list()
      muamp_s$mu=c(diag(crosum), crosum[ilow])
      nb_smu=length(muamp_s$mu)
      ## nmu_in_lmu defines where go old (narrow) mu in new (large) mu
      nmu_in_lmu=apply(t(muamp_x$mu), 2, function(mu)which(mu==muamp_s$mu))
      ## source amplitude matrix
      zmu=rep(0., nb_smu)
      sp_amp=t(as.matrix(apply(t(1:length(sp$ind_x1)), 2, function(i) {
         f[sp$ind_fb[i]]*prodlc(muamp_x$amp, sp$ind_x1[i], sp$ind_x2[i], ilow, nmu_in_lmu)
      })))
      ## prune superflucious mu
      muamp_s$amp=sp_amp # temporary storage
      muamp_s=shrink_muamp(muamp_s)
      ## sum according to b_pre@p
      p=sp$b_pre@p
      sp_amp=muamp_s$amp
      sp_amp=t(apply(t(1:dim(sp$b_pre)[2]), 2, function(i) {
         if (p[i+1]==p[i]) {return(zmu)}
         else if (p[i+1]==p[i]+1) {return(sp_amp[p[i+1],])}
         return(apply(sp_amp[(p[i]+1):p[i+1],], 2, sum))
      }))
      muamp_s$amp=matrix(0., nrow=dim(sp$b)[1], ncol=ncol(sp_amp))
      muamp_s$amp[sp$b@i+1,]=sp_amp
      muamp_s$amp=muamp_s$amp*invm
      # construct muamp_xw for this weight
      muamp_xw=list()
      ## e.v of a
      eva=eigen(a)
      u=eva$vectors
      #uinv=solve(u)
      bu=solve(u, muamp_s$amp) # uinv%*%muamp_s$amp
      lam=eva$values
      muamp_xw$mu=c(muamp_s$mu, lam)
      #muamp_xw$mu=unique(muamp_s$mu, lam)
      #if (!all(table(muamp_xw$mu)==1)) { # all mu values are supposed to be unique
      #   warning("Not all mu values are unique in xw")
      #   browser()
      #   return(muamp_x)
      #}
      muamp_xw$amp=bu/outer(-lam, muamp_s$mu, "+")
      # horizontal sum of bu to construct lam part of new amp matrix
      buh=-apply(muamp_xw$amp, 1, sum)
      # dispatch buh in matrix
      bulam=diag(buh, nrow=length(buh))
      muamp_xw$amp=cbind(muamp_xw$amp, bulam)
      # shrink xw
      muamp_xw=shrink_muamp(muamp_xw)
      # go back to natural coords
      muamp_xw$amp=u%*%muamp_xw$amp
      # extend previous weight x@amp to new lam and reorder it acordingly to the new mu vector
      omu2nmu=apply(t(muamp_x$mu), 2, function(mu)which.min(abs(mu-muamp_xw$mu)))
      tmp=matrix(0., nrow=nrow(muamp_x$amp), ncol=ncol(muamp_xw$amp))
      tmp[,omu2nmu]=muamp_x$amp
      muamp_x$amp=rbind(tmp, muamp_xw$amp)
      muamp_x$mu=muamp_xw$mu
      rownames(muamp_x$amp)=nm_incu[-1][1:nrow(muamp_x$amp)]
   }
   return(muamp_x)
}
prodlc=function(amp, i, j, ilow, nmu_in_lmu) {
   # function doing a product of linear combinations of rows i and j
   # of the matrix and groupping the terms according to ilow indications
   if (j > 1) {
      # i cannot be equal to 1 by construction of b
      p=outer(amp[i-1,], amp[j-1,], "*")
      res=c(diag(p), p[ilow]+t(p)[ilow])
   } else {
      # no product, just input cumomer
      n=ncol(amp)
      res=numeric(n*(n+1)/2)
      res[nmu_in_lmu]=amp[i-1,]
   }
   return(res)
}
shrink_muamp=function(muamp, relim=-1000, nearz=1.e-20) {
   # the following actions are performed on muamp
   # - mu having near 0 amplitudes are removed
   # - Re(mu) < relim are constrained to have Re(mu)=relim
   # - repeated mu are grouped together
   return(muamp) # for debugging only
   ikeep=apply(muamp$amp, 2, function(v) !all(abs(v) < nearz))
   #ikeep=ikeep & Re(muamp$mu)>=relim
   muamp$amp=muamp$amp[,ikeep,drop=F]
   muamp$mu=muamp$mu[ikeep]
   ilim=Re(muamp$mu)<relim
   #muamp$mu[ilim]=relim+Im(muamp$mu[ilim])
   # new mu
   nmu=as.complex(names(table(muamp$mu)))
   # sum up old mus
   namp=apply(t(nmu), 2, function(mu) {
      i=which.min(abs(mu-muamp$mu));
      if (length(i) == 1) {
         return(muamp$amp[,i])
      } else {
         return(apply(muamp$amp[,i],1,sum))
      }
   })
   namp=matrix(namp, ncol=length(nmu))
   rownames(namp)=rownames(muamp$amp)
   return(list(mu=nmu, amp=namp))
}
lab_cubs=function(ti, spAbr, li_m2x, muamp_xinp, f, m, nm_incu, measmat=NULL) {
   # labeling propagation in a cumomer network by cubic splines.
   # The network is determined by a matrix a (a function of fluxes f
   # and metabolite pool vector m) and source term s corresponding to
   # the input and lighter weight cumomers).
   # The list spAbr allows a passage from
   # vectors f and x to matrix a and source term s
   # for each cumomer weight (one weight per list item)
   # The items of the list li_m2x provide a mapping
   # of metabolite pool vector m on cumomer vectors x.
   # The source term is supposed to be causal (0 before time 0)
   # and composed exclusively
   # of exponential terms (possibly complexe in case of
   # sinusoidal input??? see in detail later)
   # Output:
   # labeling matrix where each measurement (measmat%*%x) run in time
   # takes one column. If measmat is NULL, just x is returned
   # 
   # Data structures:
   # muamp${amp,mu} - vectors of amplitudes and exponential mu
   # x[icumo, ti] - the running index is cumomer number not the time.
   # xp is first derivative of x in time
   # 2012-02-29 sokol
   
   weights=1:length(spAbr)
   muamp_xinp$amp=as.matrix(muamp_xinp$amp)
   x=t(muamp2ti(ti, muamp_xinp))
   muamp_xinp$amp=t(t(muamp_xinp$amp)*muamp_xinp$mu)
   xp=t(muamp2ti(ti, muamp_xinp))
   # cumomer number in a cumomer weigth
   nb_xw=0
   # cumomer number in all lighter weights
   nb_xl=nrow(muamp_xinp$amp)+1 # +1 is for "one" in incu names
   rownames(x)=nm_incu[2:nb_xl]
   
   nb_ti=length(ti)
   dt=diff(ti)
   for (w in weights) {
      invm=(1./m)[li_m2x[[w]]] # result has length of xw
      sp=spAbr[[w]]
      nb_xl=nb_xl+nb_xw # add precedet weights to get ligther cumo size
      nb_xw=dim(sp$tA)[1]
      nm_cuw=nm_incu[nb_xl+(1:nb_xw)]
      
      # construct the matrix a of this weight
      a=fwrv2Abr(f, sp, NULL, nm_incu[nb_xl+(1:nb_xw)], getb=F)$A*invm

      # construct the sources s and its derivatives for this weight
      s_pre=f[sp$ind_fb]*x[sp$ind_x1,]*x[sp$ind_x2,]
      p=sp$b_pre@p
      # sum up what should be summed to give the non zero part of b
      s_pre=apply(t(1:dim(sp$b_pre)[2]), 2, function(i) {
         if (p[i+1]==p[i]+1) {return(s_pre[p[i+1],])}
         return(apply(s_pre[(p[i]+1):p[i+1],], 2, sum))
      })
      s=matrix(0., nrow=nb_xw, ncol=nb_ti)
      s[sp$b@i+1,]=s_pre
      s_pre=f[sp$ind_fb]*(xp[sp$ind_x1,]*x[sp$ind_x2,]+x[sp$ind_x1,]*xp[sp$ind_x2,])
      # sum up what should be summed to give the non zero part of b
      s_pre=apply(t(1:dim(sp$b_pre)[2]), 2, function(i) {
         if (p[i+1]==p[i]+1) {return(s_pre[p[i+1],])}
         return(apply(s_pre[(p[i]+1):p[i+1],], 2, sum))
      })
      sd=matrix(0., nrow=nb_xw, ncol=nb_ti)
      sd[sp$b@i+1,]=s_pre
      # decompose the matrix a
      eva=eigen(a)
      u=eva$vectors
      uinv=solve(u)
      lam=eva$values
      # go to eigen cumomers
      s=uinv%*%s
      sd=uinv%*%sd
      # auxiliary terms for integration
      tmp1=6*(s[,-1]-s[,-nb_ti])
      tmp2=sd[,-1]+sd[,-nb_ti]
      explam=exp(outer(lam, dt, "*"))
      xw=matrix(0., nrow=nb_xw, ncol=nb_ti)
      linv=1./lam
      ltinv=outer(linv, dt, "/")
      # second derivative on the left
      sddl=tmp1-(tmp2+sd[,-nb_ti])%mrv%(2*dt)
      # second derivative on the right
      dsdd=-2*tmp1+tmp2%mrv%(6*dt)
      # third derivative
      sddd=-tmp1+tmp2%mrv%(6*dt)
      tmp1=(s+ltinv*(sd+ltinv*(sddl+ltinv*sddd)))
      tmp2=tmp1+dsdd*ltinv*ltinv
      # main part: integration of ode
      xw[,-1]=-linv*(tmp2-tmp1*explam)
      # go back to natural cumomers
      xw=u%*%xw
      rownames(xw)=nm_cuw
      # calculate derivatives
      xwp=a%*%xw+s
      # extend x and xp matrices by xw, xwp
      x=rbind(x, xw)
      xp=rbind(xp, xwp)
   }
   if (is.null(measmat)) {
      return(t(x))
   } else {
      return(t(measmat%*%x))
   }
}
eitrid=function(a,b,c,maxit=1) {
   # reccursive algorithm for eigenvalues of tridiagonal system
   # the system is supposed to have all distinct real eigenvalues
   # matrix has diagonals: a_n, b_n, c_n (b is central)
   # a[1] and c[n] are not consulted
   #
   # NOT WORKING.
   #
   A=abc2trid(a,b,c)
   n=length(b)
   if (n<0) {
      return(c())
   } else if (n==1) {
      return(b)
   } else if (n==2) {
      qb=-sum(b[1:2])/2
      qc=prod(b[1:2])-a[2]*c[1]
      di=sqrt(qb*qb-qc)
      return(c(-qb-di,-qb+di))
   }
   lam_m2=b[1]
   qb=-sum(b[1:2])/2
   qc=prod(b[1:2])-a[2]*c[1]
   di=sqrt(qb*qb-qc)
   lam_m1=c(-qb-di,-qb+di)
   for (nord in 3:n) {
      nord1=nord-1
      lam=lam_m1
      lam_m1=c(lam_m1,b[nord])
      dlam=numeric(nord-1)
      ei=eigen(A[1:nord,1:nord])
      print(ei$values)
      # newton iterations
      for (it in 1:maxit) {
         # crossdifferences of short (-last term) lambda
         crodi1=outer(lam_m1, lam, "-")
         crodi2=outer(lam_m2, lam, "-")
         if (it > 1) {
            pm1=apply(t(1:nord1),2,function(i)prod(crodi1[,i]))
            dpm1=-pm1*apply(t(1:nord1),2,function(i)sum(1./crodi1[,i]))
         } else {
            pm1=0.
            dpm1=-apply(t(1:nord1),2,function(i)prod(crodi1[-i,i]))
         }
         ca=(a[nord]*c[nord-1])*apply(t(1:nord1),2,function(i)prod(crodi2[,i]))
         # derive term
         dlam=dpm1+ca*apply(t(1:nord1),2,function(i)sum(1./crodi2[,i]))
         # take care of deriv near zero
         iz=which(abs(dlam)<=1.e-10)
         dlam[iz]=sqrt(abs((pm1-cb[iz])*prod(lam_m1[c(-iz,-nord)])/2.))
         # ratio -P/P'
         dlam=-(pm1-ca)/dlam
         lam=lam+dlam
         if (sqrt(crossprod(dlam))<=1.e-14) {
            cat("converged at it=",it,"\n")
            break
         }
         print(sort(lam))
      }
      # find last eivalue
      lam[nord]=(prod(lam_m1)-b[nord]*c[nord-1]*prod(lam_m2))/prod(lam[-nord])
      print(sort(lam))
      lam_m2=lam_m1[-nord]
      lam_m1=lam
   }
   return(lam)
}
abc2trid=function(a,b,c) {
   # matrix has diagonals: a_n, b_n, c_n (b is central)
   # a[1] and c[n] are not consulted
   n=length(b)
   A=diag(b)
   A[col(A)-row(A)==1]=c[-n]
   A[col(A)-row(A)==-1]=a[-1]
   return(A)
}
plot_tilab=function(ti, lab) {
   ncurv=ncol(lab)
   matplot(ti, lab, t="l", ylab="Labeling", xlab="Time")
   legend("bottomright", legend=colnames(lab), lty=1:ncurv, col=1:ncurv, lwd=2)
}
muamp2ti=function(ti, muamp) {
   # produce time curse from muamp structure
   ex=exp(outer(ti, muamp$mu, "*"))
   res=ex%mmt%muamp$amp
   colnames(res)=rownames(muamp$amp)
   return(res)
}

n=5 # metab number in the chain
lambda=runif(n)

# matrix of propagation, a column corresponds to one metab
nti=100;
ti=seq(0, 30, len=nti)
lab=cbind(1., lab_chain(ti, lambda))
# check that lab is the solution of equadiff
print(apply(t(1:n), 2, function(i)range(diff(lab[,i+1])/diff(ti)-lambda[i]*0.5*
      ((lab[,i]-lab[,i+1])[-1]+(lab[,i]-lab[,i+1])[-nti]))))

matplot(ti, lab, t="l", main="Not reversible linear chain", ylab="Labeling", xlab="Time (s)", lwd=2)
legend("bottomright", legend=c("input", paste("lam=", round(lambda, 2), sep="")), lty=1:(n+1), col=1:(n+1), lwd=2)

#-------------
# reversible reaction, two metabolites
nti=100;
ti=seq(0, 30., len=nti)
f=1.
m=2*c(1., 3.)
r=0.75

l_rev=lab_rev2(ti, f, m, r)
l_ref=lab_chain(ti, f/m)
l_ref=cbind(l_ref, lab_chain(ti, f/(m[1]+m[2])))

matplot(ti, cbind(l_ref, l_rev, apply(l_rev, 1, function(v)sum(v*m)/sum(m))), main=paste("Reversible reaction, r=", r, sep=""), t="l", ylab="Labeling", xlab="Time (s)", lwd=2, lty=1:6, col=1:6)
legend("bottomright", legend=c("m1 not rev", "m2 not rev", "lump(m1,m2)", "m1 rev", "m2 rev", "mean(m1,m2)"), lty=1:6, col=1:6, lwd=2)

#-------------
# reversible reactions, many metabolites
nti=100;
ti=seq(0, 10., len=nti)
f=1.
n=5
m=2*runif(n)
r=runif(n-1)

l_revm=lab_rev(ti, f, m, r)
if (n==2) {
   matplot(ti, cbind(l_rev, l_revm), main=paste("reversible reactions, r=", paste(r, collapse=",", sep=", "), sep=""), t="l")
   legend("bottomright", legend=c("m1 rev", "m2 rev", "m1", "m2"), lty=1:4, col=1:4)
} else if (all(r==0.)) {
   matplot(ti, cbind(lab_chain(ti, f/m), l_revm), main="not reversible reactions", t="l")
   legend("bottomright", legend=c(paste("lam=", round(f/m, 2), sep=""), round(m, 2)), lty=1:(2*n), col=1:(2*n))
} else {
   matplot(ti, l_revm, main=paste("Reversible reactions, r=", paste(round(r,2), collapse=",", sep=", "), sep=""), t="l", lty=1:n, col=1:n, lwd=2)
   legend("bottomright", legend=paste("m=", round(m, 2), sep=""), lty=1:n, col=1:n, lwd=2)
}

#------------
# loop in a non reversible network
nti=100;
ti=seq(0, 20., len=nti)
n=10 # total metabolite number
k=round(n-2) # k metabolites in the upper branch of the loop
f=1 # input flux
r=0.2 # rev/fwd
f1=f/(1.-r)
f2=r*f1
m=0.5+runif(n); #m=1+(1:n)/n
a=f1/m
a[-(1:k)]=f2/m[-(1:k)]
a[1]=f2/m[1]
b=-f1/m
b[-(1:k)]=-f2/m[-(1:k)]
c=rep(0., n)
A=abc2trid(a,b,c)
A[1,n]=a[1]
v=rep(0., n)
v[1]=f/m[1]

l_loop=lab_netw(ti, A, v)
matplot(ti, l_loop, main=paste("Non reversible loop, r=", r, sep=""), t="l", lty=1:n, col=1:n, lwd=2, ylab="Labeling", xlab="Time (s)")
legend("bottomright", legend=paste("m", 1:n, "=", round(m*c(rep(1, k), rep(-1, n-k)), 2), sep=""), lty=1:n, col=1:n, lwd=2)

#------------
# semi-analytical pulse for 1 box (ex_i_1box.ftbl)
nti=100;
ti=seq(0, 10., len=nti)
# read spAbr
# save(spAbr, fwrv, nm_xi, nm_incu, file="1box.Rdata")
load("1box.Rdata") # creates spAbr, fwrv, nm_xi, nm_incu for one box
nb_xi=length(nm_xi)

# vector of metabolite pools
m=c(1)
# mapping of metabolites on cumomer vector
li_m2x=list(rep(1., length(nm_incu)-nb_xi-1))
# input muamp
# all are just at 1 (step pulse)
muamp_xinp=list(mu=0., amp=matrix(1., nrow=nb_xi, ncol=1))

muamp_x=lab_cumo(spAbr, li_m2x, muamp_xinp, fwrv, m, nm_incu)
lab=muamp2ti(ti, muamp_x)
plot_tilab(ti, lab)
# check diff with exacte solution. Must be 0
range(lab[,-1]-lab_chain(ti, 1./m))

#------------
# semi-analytical sinusoidal wave 1+exp(i*t)
muamp_xinp=list(mu=c(0., 1i), amp=matrix(0.5, nrow=nb_xi, ncol=2))
muamp_x=lab_cumo(spAbr, li_m2x, muamp_xinp, fwrv, m, nm_incu)
lab=muamp2ti(ti, muamp_x)
plot_tilab(ti, lab)

#------------
# semi-analytical pulse 2 box (ex_i_2box.ftbl)
nti=100;
ti=seq(0, 20., len=nti)
# read spAbr
# save(spAbr, fwrv, nm_xi, nm_incu, nb_rw, file="2box.Rdata")
load("2box.Rdata") # creates spAbr, fwrv, nm_xi, nm_incu, nb_rw for one box
nb_xi=length(nm_xi)

# vector of metabolite pools
m=c(1., 0.5)
# mapping of metabolites on cumomer vector
li_m2x=lapply(1:nb_rw, function(i) {
   v=numeric(dim(spAbr[[i]]$tA)[1]);
   v[grep("B", nm_incu)-nb_xi-1]=1
   v[grep("C", nm_incu)-nb_xi-1]=2
   return(v)
})
# input muamp
# all are just at 1 (step pulse)
muamp_xinp=list(mu=0., amp=matrix(1., nrow=nb_xi, ncol=1))

muamp_x=lab_cumo(spAbr, li_m2x, muamp_xinp, fwrv, m, nm_incu)
lab=muamp2ti(ti, muamp_x)
# sort in alphabetical order
lab=lab[, order(colnames(lab))]
plot_tilab(ti, lab)
# check diff with exacte solution. Must be 0
range(lab[,-1]-lab_chain(ti, 1./m))
#plot_tilab(ti, cbind(lab, lab_chain(ti, 1./m)))

# sinusoidal wave 1+exp(i*t)
muamp_xinp=list(mu=c(0., 1i), amp=matrix(c(1.-0.5, 0.5), nrow=nb_xi, ncol=2))
muamp_x=lab_cumo(spAbr, li_m2x, muamp_xinp, fwrv, m, nm_incu)
labs=muamp2ti(ti, muamp_x)
# sort in alphabetical order
labs=labs[, order(colnames(labs))]
plot_tilab(ti, labs)

#------------
# semi-analytical pulse 2 reversible box (ex_i_2box_rev.ftbl)
nti=100;
ti=seq(0, 20., len=nti)
# read spAbr
# save(spAbr, fwrv, nm_xi, nm_incu, nb_rw, file="2boxr.Rdata")
load("2boxr.Rdata") # creates spAbr, fwrv, nm_xi, nm_incu, nb_rw for one box
nb_xi=length(nm_xi)

# vector of metabolite pools
m=c(1., 0.5)
# mapping of metabolites on cumomer vector
li_m2x=lapply(1:nb_rw, function(i) {
   v=numeric(dim(spAbr[[i]]$tA)[1]);
   v[grep("B", nm_incu)-nb_xi-1]=1
   v[grep("C", nm_incu)-nb_xi-1]=2
   return(v)
})
# input muamp
# all are just at 1 (step pulse)
muamp_xinp=list(mu=0., amp=matrix(1., nrow=nb_xi, ncol=1))

muamp_x=lab_cumo(spAbr, li_m2x, muamp_xinp, fwrv, m, nm_incu)
labr=muamp2ti(ti, muamp_x)
# sort in alphabetical order
labr=labr[, order(colnames(labr))]
plot_tilab(ti, labr)
# check diff with exacte solution. Must be 0
r=fwrv["rev.v2"]/fwrv["fwd.v2"]
range(labr[,-1]-lab_rev2(ti, 1, m, r))
#plot_tilab(ti, labr[,-1]-lab_rev2(ti, 1, m, r))
#plot_tilab(ti, cbind(lab, lab_chain(ti, 1./m)))

# sinusoidal wave 1+exp(i*t)
muamp_xinp=list(mu=c(0., 1i), amp=matrix(c(1.-0.5, 0.5), nrow=nb_xi, ncol=2))
muamp_x=lab_cumo(spAbr, li_m2x, muamp_xinp, fwrv, m, nm_incu)
labrs=muamp2ti(ti, muamp_x)
# sort in alphabetical order
labrs=labrs[, order(colnames(labrs))]
# see the difference with non reversible sinusoidal wave
plot_tilab(ti, labrs-labs)
plot_tilab(ti, cbind(labrs,labs))

#------------
setwd("algomics")
# semi-analyticalreal world problem (chlamy_i.ftbl)
nti=100;
ti=seq(0, 6., len=nti)
# read spAbr
# save(spAbr, fwrv, nm_xi, nm_incu, nb_rw, nb_rcumos, file="chlamy.Rdata")
load("chlamy.Rdata") # creates spAbr, fwrv, nm_xi, nm_incu, nb_rw, nb_rcumos for one box
nb_xi=length(nm_xi)

# get metab pools form ftbl
mp=system("fgrep -n METABOLITE_POOLS chlamy_i.ftbl", intern=T)
nskip=as.integer(strsplit(mp, ":")[[1]][1])
mdata=read.table("chlamy_i.ftbl", header=F, skip=nskip+1, sep="\t", quote="", strip.white=T, comment.char="/")
m=mdata[,"V3"]
nm_metab=as.character(mdata[,"V2"])
nb_metab=length(m)
names(m)=nm_metab
# mapping of metabolites on cumomer vector
li_m2x=lapply(1:nb_rw, function(i) {
   nbase=1+nb_xi+c(0, cumsum(nb_rcumos))[i]
   v=numeric(nb_rcumos[i]);
   apply(t(1:nb_metab), 2, function(j) {
      metab=nm_metab[j]
      metab=gsub("[", "\\[", metab, fixed=T)
      metab=gsub("]", "\\]", metab, fixed=T)
      icumo=grep(paste("^", metab,":", collapse="", sep=""), nm_incu[nbase+(1:nb_rcumos[i])], fixed=F)
      v[icumo] <<- j
      #cat(metab, "grep=", grep(paste(metab,":", collapse="", sep=""), nm_incu[nbase+(1:nb_rcumos[i])], fixed=T), "i=", metab2i[metab], "\n")
   })
   return(v)
})
# input muamp
# all are just at 1 (step pulse)
muamp_xinp=list(mu=0., amp=matrix(1., nrow=nb_xi, ncol=1))

#muamp_x=lab_cumo(spAbr, li_m2x, muamp_xinp, fwrv, m, nm_incu)
lab=Re(muamp2ti(ti, muamp_x))
# sort in alphabetical order
#lab=lab[, order(colnames(lab))]
plot_tilab(ti, lab)
