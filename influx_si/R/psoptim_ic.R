# This work is based on a CRAN package pso (https://CRAN.R-project.org/package=pso ) by Claus Bendtsen <papyrus.bendtsen at gmail.com>.
# As the original license requires, this modified file is released under terms of LGPL-3.
# modification's author: Serguei Sokol (sokol <at> insa-toulouse <dot> fr)
# modification's Copyright: INRAE/CNRS/INSA, 2020.

aap=arrApply::arrApply

psoptim_ic <- function (par, fn, ..., magnitude=1, mean=0,
                     control = list(), ui=NULL, ci=NULL) {
  fn1 <- function(par) fn(par, ...)/p.fnscale
  mrunif <- function(n,m,magnitude,mean) {
    return(matrix(runif(n*m,-1,1),nrow=n,ncol=m)*magnitude+mean)
  }
  norm <- function(x) sqrt(sum(x*x))
  rsphere.unif <- function(n,r) {
    temp <- runif(n, 0., 1.)
    return((runif(1,min=0,max=r)/norm(temp))*temp)
  }
  svect <- function(a,b,n,k) {
    temp <- rep(a,n)
    temp[k] <- b
    return(temp)
  }
  normcol=function(mat)
    aap(mat, 1, "norm")
  mrsphere.unif <- function(n,r) {
    m <- length(r)
    temp <- matrix(runif(n*m, 0, 1),n,m)
    return(aap(temp, 2, "multv", v=runif(m,min=0,max=r)/normcol(temp)))
  }
  almost1=1.-1.e5*.Machine$double.eps
  npar <- length(par)
  nm_par=names(par)
  con <- list(trace = 0, fnscale = 1, maxit = 1000L, maxf = Inf,
              abstol = -Inf, reltol = 0, REPORT = 10,
              s = NA, k = 3, p = NA, w = 1/(2*log(2)),
              c.p = .5+log(2), c.g = .5+log(2), d = NA,
              v.max = NA, rand.order = TRUE, max.restart=Inf,
              maxit.stagnate = Inf,
              trace.stats = FALSE, type = "SPSO2007",
              tolineq=1.e-10)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  ## Argument error checks
  if (any(!is.finite(magnitude)))
    stop("magnitude must be finite value")
  if (any(!is.finite(mean)))
    stop("magnitude must be finite value")
  p.type <- pmatch(con[["type"]],c("SPSO2007","SPSO2011"))-1
  if (is.na(p.type)) stop("type should be one of \"SPSO2007\", \"SPSO2011\"")
  if (!is.null(ui)) {
    if (is.null(ci))
      stop("both 'ui' and ci' must be NULL or not NULL simultaneously")
    if (nrow(ui) != NROW(ci))
      stop("row number in 'ui' and 'ci' must be the same")
    if (ncol(ui) != length(par))
      stop("column number in 'ui' must be equal to length(par)")
  }
  if (!is.null(ci) && is.null(ui))
    stop("both 'ui' and ci' must be NULL or not NULL simultaneously")
  
  p.trace <- con[["trace"]]>0L # provide output on progress?
  p.fnscale <- con[["fnscale"]] # scale funcion by 1/fnscale
  p.maxit <- con[["maxit"]] # maximal number of iterations
  p.maxf <- con[["maxf"]] # maximal number of function evaluations
  p.abstol <- con[["abstol"]] # absolute tolerance for convergence
  p.reltol <- con[["reltol"]] # relative minimal tolerance for restarting
  p.report <- as.integer(con[["REPORT"]]) # output every REPORT iterations
  p.s <- ifelse(is.na(con[["s"]]),ifelse(p.type==0,floor(10+2*sqrt(npar)),40),
                con[["s"]]) # swarm size
  p.p <- ifelse(is.na(con[["p"]]),1-(1-1/p.s)^con[["k"]],con[["p"]]) # average % of informants
  p.w0 <- con[["w"]] # exploitation constant
  if (length(p.w0)>1) {
    p.w1 <- p.w0[2]
    p.w0 <- p.w0[1]
  } else {
    p.w1 <- p.w0
  }
  p.c.p <- con[["c.p"]] # local exploration constant
  p.c.g <- con[["c.g"]] # global exploration constant
  p.d <- ifelse(is.na(con[["d"]]),norm(magnitude),con[["d"]]) # domain diameter
  p.vmax <- con[["v.max"]]*p.d # maximal velocity
  p.randorder <- as.logical(con[["rand.order"]]) # process particles in random order?
  p.maxrestart <- con[["max.restart"]] # maximal number of restarts
  p.maxstagnate <- con[["maxit.stagnate"]] # maximal number of iterations without improvement
  p.trace.stats <- as.logical(con[["trace.stats"]]) # collect detailed stats?
  p.has_ic <- !is.null(ui)
  p.tolineq <- con[["tolineq"]]
  
  if (p.trace) {
    cat("S=",p.s,", K=",con[["k"]],", p=",signif(p.p,4),", w0=",
            signif(p.w0,4),", w1=",
            signif(p.w1,4),", c.p=",signif(p.c.p,4),
            ", c.g=",signif(p.c.g,4), "\n")
    cat("v.max=",signif(con[["v.max"]],4),
            ", d=",signif(p.d,4), "\n")
    if (p.trace.stats) {
      stats.trace.it <- c()
      stats.trace.error <- c()
      stats.trace.f <- NULL
      stats.trace.x <- NULL
    }
  }
  ## Initialization
  if (p.reltol!=0) p.reltol <- p.reltol*p.d
  X <- mrunif(npar,p.s,magnitude,mean)
  if (!anyNA(par)) X[,1] <- par
  # put inside
  if (p.has_ic) {
    X <- vapply(seq(ncol(X)), function(i) put_inside(X[,i], ui, ci), X[,1L])
    mean=rowMeans(X)
  }
#browser()
  if (p.type==0) {
    V <- mrunif(npar,p.s,magnitude,0)
  } else { ## p.type==1
    V <- matrix(runif(npar*p.s,min=-magnitude,max=magnitude),npar,p.s)
    p.c.p2 <- p.c.p/2 # precompute constants
    p.c.p3 <- p.c.p/3
    p.c.g3 <- p.c.g/3
    p.c.pg3 <- p.c.p3+p.c.g3
  }
  if (!is.na(p.vmax)) { # scale to maximal velocity
    temp <- normcol(V)
    temp <- pmin.int(temp,p.vmax)/temp
    V <- aap(V, 2, "multv", v=temp)
  }
  f.x <- apply(X,2,fn1) # first evaluations
  stats.feval <- p.s
  P <- X
  f.p <- f.x
  P.improved <- rep(FALSE,p.s)
  i.best <- which.min(f.p)
  error <- f.p[i.best]
  init.links <- TRUE
  if (p.trace && p.report==1) {
    cat("It 1: fitness=",signif(error,4), "\n")
    if (p.trace.stats) {
      stats.trace.it <- c(stats.trace.it,1)
      stats.trace.error <- c(stats.trace.error,error)
      stats.trace.f <- c(stats.trace.f,list(f.x))
      stats.trace.x <- c(stats.trace.x,list(X))
    }
  }
  ## Iterations
  msgcode=0
  stats.iter <- 1
  stats.restart <- 0
  stats.stagnate <- 0
  while (stats.iter<p.maxit && stats.feval<p.maxf && error>p.abstol &&
         stats.restart<p.maxrestart && stats.stagnate<p.maxstagnate) {
    stats.iter <- stats.iter+1
    if (p.p!=1 && init.links) {
      links <- matrix(runif(p.s*p.s,0,1)<=p.p,p.s,p.s)
      diag(links) <- TRUE
    }
    ## The swarm moves
    if (p.p==1)
      j <- rep(i.best,p.s)
    else # best informant
      j <- sapply(1:p.s,function(i)
                  which(links[,i])[which.min(f.p[links[,i]])]) 
    temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
    V <- temp*V # exploration tendency
    if (p.type==0) {
      V <- V+mrunif(npar,p.s,p.c.p*0.5,p.c.p*0.5)*(P-X) # exploitation
      temp <- j!=(1:p.s)
      V[,temp] <- V[,temp]+mrunif(npar,sum(temp),p.c.p*0.5,p.c.p*0.5)*(P[,j[temp]]-X[,temp])
    } else { # SPSO 2011
      temp <- j==(1:p.s)
      temp <- aap(P, 2, "multv", v=svect(p.c.p3,p.c.p2,p.s,temp))+
        aap(P[,j], 2, "multv", v=svect(p.c.g3,0,p.s,temp))-
        aap(X, 2, "multv", v=svect(p.c.pg3,p.c.p2,p.s,temp)) # G-X
      V <- V+temp+mrsphere.unif(npar,normcol(temp))
    }
    if (!is.na(p.vmax)) {
      temp <- normcol(V)
      temp <- pmin.int(temp,p.vmax)/temp
      V <- aap(V, 2, "multv", v=temp)
    }
#print(stats.iter)
#print(P-X)
#print("V generated")
#print(V)
    if (p.has_ic) {
      # enforce constraints
      X1 <- X+V
      ine1 <- ui%*%X1 - ci
      ine <- ui%*%X - ci
      if (any(ine < -p.tolineq)) {
        iw=which(ine<0, arr.ind=TRUE)
        nm_i=if (is.null(rownames(ui))) iw[,1L] else rownames(ui)[iw[,1L]]
        mes=format(rbind(c("violation", "ineq", "iparticle"), cbind(ine[ine<0], nm_i, iw[,2L])))
        mes=paste0(apply(mes, 1, paste0, collapse=" "), collapse="\n\t")
#browser()
        msg=paste0("not all particles satisfy inequality constraints (iter=", stats.iter, "):\n\t", mes)
        msgcode=5
        break
      }
      for (j in seq(ncol(X))) {
        if (!any(ibad <- (ine1[,j] < -p.tolineq))) {
          X[,j] <- X[,j]+V[,j]
          next
        }
        alpha=min((ine[ibad,j]-p.tolineq)/(ine[ibad,j]-ine1[ibad,j])) # aim at a value slightly over 0
        if (!is.finite(alpha))
          alpha <- 0.
        alpha <- (if (alpha > 0.) almost1 else 1./almost1)*alpha # alpha can be < 0 if X was already slightly outside of the domain
        V[,j] <- alpha*V[,j]
#browser()
        if (any(ui%*%(X[,j]+V[,j])-ci < -p.tolineq))
          V[,j] <- (if (alpha  > 0.) almost1 else 1./almost1)*V[,j]
        X[,j] <- X[,j]+V[,j]
        # remove the normal component of V_j looking outside of the domain (normal to the constraint hyperplane)
        ibad <- which.min(ui%*%X[,j] - ci)
        nrm=ui[ibad,]%*%V[,j]
        if (nrm < 0.)
          V[,j] <- V[,j]-ui[ibad,]*c(nrm/sum(ui[ibad,]**2))
      }
    } else {
      X <- X+V
    }
#print("V used")
#print(V)
    ## Evaluate function
    f.x <- apply(X,2,fn1)
    stats.feval <- stats.feval+p.s
    temp <- f.x<f.p
    if (any(temp)) { # improvement
      P[,temp] <- X[,temp]
      f.p[temp] <- f.x[temp]
      i.best <- which.min(f.p)
    }
    if (stats.feval>=p.maxf) break
    if (p.reltol!=0) {
      d <- X-P[,i.best]
      d <- max(normcol(d))
      if (d<p.reltol) {
        # restart
        X <- mrunif(npar,p.s,magnitude,mean)
        # put inside
        if (p.has_ic) {
          X <- vapply(seq(ncol(X)), function(i) put_inside(X[,i], ui, ci), X[,1L])
          mean=rowMeans(X)
        }

        V <- mrunif(npar,p.s,magnitude,0)
        if (!is.na(p.vmax)) {
          temp <- normcol(V)
          temp <- pmin.int(temp,p.vmax)/temp
          V <- aap(V, 2, "multv", v=temp)
        }
        stats.restart <- stats.restart+1
        if (p.trace) cat("It ", stats.iter, ": restarting", "\n", sep="")
      }
    }
    init.links <- f.p[i.best]==error # if no overall improvement
    stats.stagnate <- ifelse(init.links,stats.stagnate+1,0)
    error <- f.p[i.best]
    if (p.trace && stats.iter%%p.report==0) {
      if (p.reltol!=0) 
        cat("It ",stats.iter,": fitness=",signif(error,4),
                ", swarm diam.=",signif(d,4), "\n", sep="")
      else
        cat("It ",stats.iter,": fitness=",signif(error,4), "\n", sep="")
      if (p.trace.stats) {
        stats.trace.it <- c(stats.trace.it,stats.iter)
        stats.trace.error <- c(stats.trace.error,error)
        stats.trace.f <- c(stats.trace.f,list(f.x))
        stats.trace.x <- c(stats.trace.x,list(X))
      }
    }
  }
  if (error<=p.abstol) {
    msg <- "Converged"
    msgcode <- 0
  } else if (stats.feval>=p.maxf) {
    msg <- "Maximal number of function evaluations reached"
    msgcode <- 1
  } else if (stats.iter>=p.maxit) {
    msg <- "Maximal number of iterations reached"
    msgcode <- 2
  } else if (stats.restart>=p.maxrestart) {
    msg <- "Maximal number of restarts reached"
    msgcode <- 3
  } else if (msgcode==5) {
    ; # inequalities failed
  } else {
    msg <- "Maximal number of iterations without improvement reached"
    msgcode <- 4
  }
  if (p.trace) cat(msg, "\n")
  o <- list(par=structure(P[,i.best], names=nm_par),value=f.p[i.best],
            counts=c("function"=stats.feval,"iteration"=stats.iter,
              "restarts"=stats.restart),
            convergence=msgcode,message=msg)
  if (p.trace && p.trace.stats) o <- c(o,list(stats=list(it=stats.trace.it,
                                                error=stats.trace.error,
                                                f=stats.trace.f,
                                                x=stats.trace.x)))
  return(o)
}
