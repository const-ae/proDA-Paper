# This file was extracted from the supplementary materials of 
# Koopmans, F., Cornelisse, L. N., Heskes, T. & Dijkstra, T. M. H. 
#   Empirical bayesian random censoring threshold model improves detection of differentially abundant proteins. 
#   J. Proteome Res. 13, 3871–3880 (2014).

# It contains an additional bug fix where a call to ifelse in line 130
# sometimes failed to produce the expected numeric matrix.

RanCenHier2Fun_mod <- function(Y1, Y2, frac.prior.mu.x=0.1, frac.prior.sigma.x=1.0,
                           frac.prior.mu.c=0.1, frac.prior.sigma.c=0.01,
                           n.burnin=200, n.post.samp=1000, n.chain=3, sleep.time=0) {
  # Fit Bayesian random censoring threshold model with Gibbs sampling
  # to a pair of protein abundance data sets (biosample 1 and 2)
  # Input is two data frames Y1, Y2 with N.prot rows and N.rep[1,2] columns
  # Output n.chain, n.burnin, n.iter, 
  #  mu.x.prior, sigma.x.prior, mu.c.prior, sigma.c.prior,
  #  frac.prior.mu.x, frac.prior.sigma.x, frac.prior.mu.c, frac.prior.sigma.c,
  #  post.mu.x1.mean, post.mu.x2.mean, post.sigma.x.mean, post.mu.c.mean, post.sigma.c.mean,
  #  post.mu.x1.sd, post.mu.x2.sd, post.sigma.x.sd, post.mu.c.sd, post.sigma.c.sd,
  #   post.p.cens1, post.p.cens2, post.x1, post.x2, post.c1, post.c2
  #
  # Derived from version 1.9 of GibbsRanCenHier2.R
  # Version 0.3, 10 January 2014 by Tjeerd Dijkstra
  # Tested with R 3.0.1 under OS X 10.7.5 on an i7-2635QM@2GHz CPU
  
  # random samples from a truncated normal variable via cdf inversion method
  # with safeguard for high truncation point by mapping
  # the upper tail pnorm ~ 1 to the lower tail pnorm ~ 0
  # this code will work up to 37 SD's from the mean
  # note that this code is vectorized, par.mu, par.sigma, lo and hi
  # should either be of length 1 (scalar) or length n
  rnorm.trunc2 <- function(n, par.mu, par.sigma, lo=-Inf, hi=Inf)
  { p.lo <- pnorm(lo, par.mu, par.sigma)
  p.hi <- pnorm(hi, par.mu, par.sigma)
  u <- runif(n, p.lo, p.hi)
  out <- qnorm(u, par.mu, par.sigma)
  if(any(p.lo > 1 - sqrt(.Machine$double.eps)))
  { p.lo.flag <- (p.lo > 1 - sqrt(.Machine$double.eps))
  p.lo <- pnorm(-lo, -par.mu, par.sigma)
  p.hi <- pnorm(-hi, -par.mu, par.sigma)
  u <- runif(n, p.hi, p.lo)
  if(length(par.mu) == 1) par.mu <- rep(par.mu, n)
  if(length(par.sigma) == 1) par.sigma <- rep(par.sigma, n)
  out[p.lo.flag] <- -qnorm(u[p.lo.flag], -par.mu[p.lo.flag], par.sigma[p.lo.flag])
  }
  return(out)
  }
  
  # two-dimensional normal random variable with truncation such that x1 > x2
  # n: number of samples
  # m1, m2: means
  # s1, s2: standard deviations
  normtrunc2drnd <- function(n, m1, s1, m2, s2)
  { mmin <- m1-m2
  mplus <- m1+m2
  vplus <- s1^2+s2^2
  ymin <- rnorm.trunc2(n, mmin, sqrt(vplus), 0, Inf)
  
  rho <- (s1^2-s2^2)/vplus
  vcond <- (1-rho^2)*vplus
  mcond <- mplus + rho*(ymin-mmin)
  yplus <- rnorm(n, mcond, sqrt(vcond))
  
  return(list(x1=(yplus+ymin)/2, x2=(yplus-ymin)/2))
  }
  
  
  ### preprocessing
  # treat proteins with all observations missing in both biosample 1 or 2
  all.mis.flag <- (rowSums(!is.na(Y1)) == 0) & (rowSums(!is.na(Y2)) == 0)
  N.all.mis <- sum(all.mis.flag)
  if (N.all.mis > 0)
  { cat(sprintf("%d proteins out of %d have all observations missing\n", N.all.mis, nrow(Y1)))
    Y1.orig <- Y1; Y2.orig <- Y2
    Y1 <- rbind(Y1[!all.mis.flag, ], rep(NA, ncol(Y1)))
    Y2 <- rbind(Y2[!all.mis.flag, ], rep(NA, ncol(Y2)))
  }
  
  # convenient derived quantities and constants
  N.rep1 <- ncol(Y1); N.rep2 <- ncol(Y2); N.rep12 <- N.rep1 + N.rep2
  N.prot <- nrow(Y1); N <- N.prot*N.rep12
  # censor flag is called "w" in the Tech Note
  Cen.flag1 <- is.na(Y1); Obs.flag1 <- !Cen.flag1
  N.cens1 <- rowSums(Cen.flag1); P.cens1 <- N.cens1/N.rep1
  N.obs1 <- rowSums(Obs.flag1); P.obs1 <- N.obs1/N.rep1
  Cen.flag2 <- is.na(Y2); Obs.flag2 <- !Cen.flag2
  N.cens2 <- rowSums(Cen.flag2); P.cens2 <- N.cens2/N.rep2
  N.obs2 <- rowSums(Obs.flag2); P.obs2 <- N.obs2/N.rep2
  
  
  ### specification of prior
  # prior means from rough ad-hoc estimates
  tmp.sd <- apply(cbind(Y1, Y2), 1, sd, na.rm=TRUE)
  min.prot.abun1 <- apply(Y1[(N.cens1 > 0), ], 1, min, na.rm=TRUE)
  min.prot.abun2 <- apply(Y2[(N.cens2 > 0), ], 1, min, na.rm=TRUE)
  min.prot.abun <- c(min.prot.abun1, min.prot.abun2)
  
  mu.x.prior <- mean(rowMeans(cbind(Y1, Y2), na.rm=TRUE), na.rm=TRUE)
  sigma.x.prior <- mean(tmp.sd, na.rm=TRUE)
  mu.c.prior <- mean(min.prot.abun)
  sigma.c.prior <- sd(min.prot.abun)
  
  # prior weights
  n.prior.pseudo.obs.mu.x <- frac.prior.mu.x*N.rep12/2
  n.prior.pseudo.obs.sigma.x <- frac.prior.sigma.x*N.rep12
  n.prior.pseudo.obs.mu.c <- frac.prior.mu.c*N.rep12
  n.prior.pseudo.obs.sigma.c <- frac.prior.sigma.c*N
  
  # md0, kd0, ad0, bd0: parameters of normal-Gamma prior for data model
  md0 <- mu.x.prior # prior data mean
  kd0 <- n.prior.pseudo.obs.mu.x # number of pseudo observations of prior mean
  bd0 <- (n.prior.pseudo.obs.sigma.x/2)*sigma.x.prior^2 # prior data precision
  ad0 <- n.prior.pseudo.obs.sigma.x/2 # twice number of pseudo observations of prior precision
  # mc0, kc0, ac0, bc0: parameters of normal-Gamma prior for censor model
  mc0 <- mu.c.prior # prior censor mean
  kc0 <- n.prior.pseudo.obs.mu.c # number of pseudo observations of prior mean
  bc0 <- (n.prior.pseudo.obs.sigma.c/2)*sigma.c.prior^2 # prior censor precision
  ac0 <- n.prior.pseudo.obs.sigma.c/2 # twice number of pseudo observations of prior precision
  
  
  ### run Gibbs sampler
  # MCMC parameters
  t1 <- proc.time()
  n.iter <- n.post.samp + n.burnin
  
  # inital values, same parametrization as prior
  n.init.pseudo.obs.x <- 0.5*N.rep12
  md.init <- mu.x.prior;                              kd.init <- n.init.pseudo.obs.x
  bd.init <- (n.init.pseudo.obs.x/2)*sigma.x.prior^2; ad.init <- n.init.pseudo.obs.x/2
  n.init.pseudo.obs.c <- 0.1*N
  mc.init <- mu.c.prior;                              kc.init <- n.init.pseudo.obs.x
  bc.init <- (n.init.pseudo.obs.c/2)*sigma.c.prior^2; ac.init <- n.init.pseudo.obs.c/2
  
  # Xc[1,2] and Xd[1,2] are of the same dimensions as Y1 and Y2. Xc[1,2] contains
  # samples from the censoring threshold distribution and Xd[1,2] contains the data
  # (if record is observed) and otherwise samples from the data distribution
  # in the Tech Note Xc is called c and Xd is called x
  Xd1 <- as.matrix(Y1)
  Xd1[Cen.flag1] <- 0
  Xd2 <- as.matrix(Y2)
  Xd2[Cen.flag2] <- 0
  Xc1 <- matrix(0, N.prot, N.rep1); Xc2 <- matrix(0, N.prot, N.rep2)
  # arrays to hold the Gibbs samples from each Markov chain
  chains.mu.x1 <- chains.mu.x2 <- array(0, c(N.prot, n.iter-n.burnin, n.chain))
  chains.mu.c <- chains.sigma.x <- array(0, c(N.prot, n.iter-n.burnin, n.chain))
  chains.sigma.c <- matrix(0, n.iter-n.burnin, n.chain)
  
  # storing posterior parameter means and complete data for each chain
  post.mu.x1 <- matrix(0, N.prot, n.iter); post.mu.x2 <- matrix(0, N.prot, n.iter)
  post.sigma.x <- matrix(0, N.prot, n.iter); post.mu.c <- matrix(0, N.prot, n.iter)
  post.sigma.c <- rep(0, n.iter)
  
  # for debugging only store data from last chain
  debug.mean.x1 <- matrix(0, N.prot, n.iter)
  debug.mean.x2 <- matrix(0, N.prot, n.iter)
  debug.var.x12 <- matrix(0, N.prot, n.iter)
  # note that storage for Xd1 and Xd2 is inefficient as the observed Y1 and Y2
  # are stored for every Gibbs sample
  Xd1.all <- array(0, c(N.prot, N.rep1, n.iter))
  Xd2.all <- array(0, c(N.prot, N.rep2, n.iter))
  Xc1.all <- array(0, c(N.prot, N.rep1, n.iter))
  Xc2.all <- array(0, c(N.prot, N.rep2, n.iter))
  
  for (i.chain in 1:n.chain)
  { # draw first sample from inital value distribution  
    post.sigma.x[, 1] <- 1/sqrt(rgamma(N.prot, ad.init, rate=bd.init))
    post.mu.x1[, 1] <- rnorm(N.prot, md.init, post.sigma.x[, 1]/sqrt(kd.init))
    post.mu.x2[, 1] <- rnorm(N.prot, md.init, post.sigma.x[, 1]/sqrt(kd.init))
    post.sigma.c[1] <- 1/sqrt(rgamma(1, ac.init, rate=bc.init))
    post.mu.c[, 1] <- rnorm(N.prot, mc.init, post.sigma.c[1]/sqrt(kc.init))
    for (j.samp in 2:n.iter)
    { # vectorized code to avoid looping over proteins for biosample 1
      # uncensored data: draw normal random variable for censoring value truncated above by data value
      Post.mu.c1 <- matrix(post.mu.c[, j.samp-1], N.prot, N.rep1, byrow=FALSE)
      Xc1[Obs.flag1] <- rnorm.trunc2(sum(N.obs1), Post.mu.c1[Obs.flag1],
                                     post.sigma.c[j.samp-1], -Inf, Xd1[Obs.flag1])
      # censored data: draw two normal random variables with
      # constraint that censoring value is larger than data value
      Post.mu.x1 <- matrix(post.mu.x1[, j.samp-1], N.prot, N.rep1, byrow=FALSE)
      Post.sigma.x1 <- matrix(post.sigma.x[, j.samp-1], N.prot, N.rep1, byrow=FALSE)
      res <- normtrunc2drnd(sum(N.cens1), Post.mu.c1[Cen.flag1], post.sigma.c[j.samp-1],
                            Post.mu.x1[Cen.flag1], Post.sigma.x1[Cen.flag1])
      Xc1[Cen.flag1] <- res$x1; Xd1[Cen.flag1] <- res$x2
      
      # same for biosample 2
      Post.mu.c2 <- matrix(post.mu.c[, j.samp-1], N.prot, N.rep2, byrow=FALSE)
      Xc2[Obs.flag2] <- rnorm.trunc2(sum(N.obs2), Post.mu.c2[Obs.flag2],
                                     post.sigma.c[j.samp-1], -Inf, Xd2[Obs.flag2])
      Post.mu.x2 <- matrix(post.mu.x2[, j.samp-1], N.prot, N.rep2, byrow=FALSE)
      Post.sigma.x2 <- matrix(post.sigma.x[, j.samp-1], N.prot, N.rep2, byrow=FALSE)
      res <- normtrunc2drnd(sum(N.cens2), Post.mu.c2[Cen.flag2], post.sigma.c[j.samp-1],
                            Post.mu.x2[Cen.flag2], Post.sigma.x2[Cen.flag2])
      Xc2[Cen.flag2] <- res$x1; Xd2[Cen.flag2] <- res$x2
      
      
      # update normal-gamma distributions given data samples
      # See Murphy 2007 eqs 21-24
      mean.x1 <- rowMeans(Xd1); mean.x2 <- rowMeans(Xd2); mean.x12 <- (mean.x1 + mean.x2)/2
      var.x <- rowMeans((cbind(Xd1, Xd2) - mean.x12)^2) # note vectorized code
      
      kappa.x1.next <- kd0 + N.rep1
      mu.x1.next <- (kd0*md0 + N.rep1*mean.x1)/kappa.x1.next
      
      kappa.x2.next <- kd0 + N.rep2
      mu.x2.next <- (kd0*md0 + N.rep2*mean.x2)/kappa.x2.next
      
      alpha.x.next <- ad0 + N.rep12/2
      beta.x.next <- bd0 + (N.rep12/2)*(var.x + kd0*(mean.x12-md0)^2/(kd0 + N.rep12))
      
      # draw new means and standard deviations for data mean and sd
      post.sigma.x[, j.samp] <- 1/sqrt(rgamma(N.prot, alpha.x.next, rate=beta.x.next))
      post.mu.x1[, j.samp] <- rnorm(N.prot, mu.x1.next, post.sigma.x[, j.samp]/sqrt(kappa.x1.next))
      post.mu.x2[, j.samp] <- rnorm(N.prot, mu.x2.next, post.sigma.x[, j.samp]/sqrt(kappa.x2.next))
      
      # update normal-gamma distributions given censoring threshold samples
      mean.c <- rowMeans((cbind(Xc1, Xc2))); var.c <- mean((cbind(Xc1, Xc2) - mean.c)^2)
      kappa.c.next <- kc0 + N.rep12
      mu.c.next <- (kc0*mc0 + N.rep12*mean.c)/kappa.c.next
      alpha.c.next <- ac0 + N/2
      beta.c.next <- bc0 + (N/2)*(var.c + kc0*(mean(mean.c)-mc0)^2/(kc0 + N))
      # draw new means and standard deviations for censoring threshold mean and sd
      post.sigma.c[j.samp] <- 1/sqrt(rgamma(1, alpha.c.next, rate=beta.c.next))
      post.mu.c[, j.samp] <- rnorm(N.prot, mu.c.next, post.sigma.c[j.samp]/sqrt(kappa.c.next))
      
      # save samples
      debug.mean.x1[, j.samp] <- mean.x1; debug.mean.x2[, j.samp] <- mean.x2
      debug.var.x12[, j.samp] <- var.x
      Xd1.all[, , j.samp] <- Xd1; Xc1.all[, , j.samp] <- Xc1
      Xd2.all[, , j.samp] <- Xd2; Xc2.all[, , j.samp] <- Xc2
      
      # sanity check
      stopifnot(!is.na(post.sigma.c[j.samp]), is.finite(post.sigma.c[j.samp]))
    }
    # discard burnin samples and store each Markov chain
    chains.mu.x1[, , i.chain] <- post.mu.x1[, (1+n.burnin):n.iter];
    chains.mu.x2[, , i.chain] <- post.mu.x2[, (1+n.burnin):n.iter];
    chains.sigma.x[, , i.chain] <- post.sigma.x[, (1+n.burnin):n.iter]
    chains.mu.c[, , i.chain] <- post.mu.c[, (1+n.burnin):n.iter];
    chains.sigma.c[, i.chain] <- post.sigma.c[(1+n.burnin):n.iter]
    
    Sys.sleep(sleep.time) # to keep laptop cool
  }
  t1 <- (proc.time() - t1)[3] - sleep.time*n.chain
  
  
  # flatten all chains 
  post.mu.x1 <- chains.mu.x1;     dim(post.mu.x1) <- c(N.prot, (n.iter-n.burnin)*n.chain)
  post.mu.x2 <- chains.mu.x2;     dim(post.mu.x2) <- c(N.prot, (n.iter-n.burnin)*n.chain)
  post.sigma.x <- chains.sigma.x; dim(post.sigma.x) <- c(N.prot, (n.iter-n.burnin)*n.chain)
  post.mu.c <- chains.mu.c;       dim(post.mu.c) <- c(N.prot, (n.iter-n.burnin)*n.chain)
  post.sigma.c <- c(chains.sigma.c)
  # sample N.cens[1,2] from posterior x distribution and N.rep[1,2] from posterior c distribution
  # for the x this is work as we need to sample a variable number depending on N.cens[1,2]
  post.x1 <- vector("list", N.prot); post.x2 <- vector("list", N.prot)
  post.c1 <- matrix(0, N.prot, N.rep1); post.c2 <- matrix(0, N.prot, N.rep2)
  for (k.prot in 1:N.prot) 
  { if (N.cens1[k.prot] > 0) {
    Xd1.tmp <- Xd1.all[k.prot, Cen.flag1[k.prot, ], (1+n.burnin):n.iter]
    post.x1[[k.prot]] <- sample(c(Xd1.tmp), N.cens1[k.prot])
  }
    post.c1[k.prot, ] <- sample(c(Xc1.all[k.prot, , (1+n.burnin):n.iter]), N.rep1)
    
    if (N.cens2[k.prot] > 0) {
      Xd2.tmp <- Xd2.all[k.prot, Cen.flag2[k.prot, ], (1+n.burnin):n.iter]
      post.x2[[k.prot]] <- sample(c(Xd2.tmp), N.cens2[k.prot])
    }
    post.c2[k.prot, ] <- sample(c(Xc2.all[k.prot, , (1+n.burnin):n.iter]), N.rep2)
  }
  
  # summary statistics of each parameter
  post.mu.x1.mean <- rowMeans(post.mu.x1);      post.mu.x1.sd <- apply(post.mu.x1, 1, sd)
  post.mu.x2.mean <- rowMeans(post.mu.x2);      post.mu.x2.sd <- apply(post.mu.x2, 1, sd)
  post.sigma.x.mean <- rowMeans(post.sigma.x);  post.sigma.x.sd <- apply(post.sigma.x, 1, sd)
  post.mu.c.mean <- rowMeans(post.mu.c);        post.mu.c.sd <- apply(post.mu.c, 1, sd)
  post.sigma.c.mean <- mean(post.sigma.c);      post.sigma.c.sd <- sd(post.sigma.c)
  post.p.cens1 <- pnorm((post.mu.c.mean - post.mu.x1.mean)/sqrt(post.sigma.c.mean^2 + post.sigma.x.mean^2))
  post.p.cens2 <- pnorm((post.mu.c.mean - post.mu.x2.mean)/sqrt(post.sigma.c.mean^2 + post.sigma.x.mean^2))
  
  # print results to console
  cat(sprintf("RanCenHier2Fun n.prot %d  n.rep1 %d  n.rep2 %d  ", N.prot, N.rep1, N.rep2))
  cat(sprintf("n.chain %d  n.burnin %d  n.iter %d took %4.1f s  %3d samp/s\n",
              n.chain, n.burnin, n.iter-n.burnin, t1, round(n.chain*n.iter/t1)))
  cat("variable  frac.pri   mean prior  Gibbs samp\n")
  cat(sprintf("mu.x1      %5.4f     %5.4f   %5.4f+/-%5.4f\n",
              frac.prior.mu.x, mu.x.prior, mean(post.mu.x1.mean), mean(post.mu.x1.sd)))
  cat(sprintf("mu.x2      %5.4f     %5.4f   %5.4f+/-%5.4f\n",
              frac.prior.mu.x, mu.x.prior, mean(post.mu.x2.mean), mean(post.mu.x2.sd)))
  cat(sprintf("sigma.x    %5.4f     %5.4f   %5.4f+/-%5.4f\n",
              frac.prior.sigma.x, sigma.x.prior, mean(post.sigma.x.mean), mean(post.sigma.x.sd)))
  cat(sprintf("mu.c       %5.4f     %5.4f   %5.4f+/-%5.4f\n",
              frac.prior.mu.c, mu.c.prior, mean(post.mu.c.mean), mean(post.mu.c.sd)))
  cat(sprintf("sigma.c    %5.4f     %5.4f   %5.4f+/-%5.4f\n",
              frac.prior.sigma.c, sigma.c.prior, post.sigma.c.mean, post.sigma.c.sd))
  
  ## to be done in version > 0.1
  # # treat proteins with all observations missing
  # if (N.all.mis > 0)
  # { N.prot.orig <- nrow(Y.orig)
  #   post.mu.x.mean.tmp <- post.mu.x.mean
  #   post.mu.x.mean <- vector("numeric", N.prot.orig)
  #   post.mu.x.mean[!all.mis.flag] <- post.mu.x.mean.tmp[-N.prot]
  #   post.mu.x.mean[all.mis.flag] <- post.mu.x.mean.tmp[N.prot]
  # 
  #   post.sigma.x.mean.tmp <- post.sigma.x.mean
  #   post.sigma.x.mean <- vector("numeric", N.prot.orig)
  #   post.sigma.x.mean[!all.mis.flag] <- post.sigma.x.mean.tmp[-N.prot]
  #   post.sigma.x.mean[all.mis.flag] <- post.sigma.x.mean.tmp[N.prot]
  #   post.p.cens <- pnorm((post.mu.c.mean - post.mu.x.mean)/sqrt(post.sigma.c.mean^2 + post.sigma.x.mean^2))
  # 
  #   post.mu.x.sd.tmp <- post.mu.x.sd
  #   post.mu.x.sd <- vector("numeric", N.prot.orig)
  #   post.mu.x.sd[!all.mis.flag] <- post.mu.x.sd.tmp[-N.prot]
  #   post.mu.x.sd[all.mis.flag] <- post.mu.x.sd.tmp[N.prot]
  # 
  #   post.sigma.x.sd.tmp <- post.sigma.x.sd
  #   post.sigma.x.sd <- vector("numeric", N.prot.orig)
  #   post.sigma.x.sd[!all.mis.flag] <- post.sigma.x.sd.tmp[-N.prot]
  #   post.sigma.x.sd[all.mis.flag] <- post.sigma.x.sd.tmp[N.prot]
  # 
  #   post.x.tmp <- post.x
  #   post.x <- vector("list", N.prot.orig)
  #   post.x[!all.mis.flag] <- post.x.tmp[-N.prot]
  #   post.x[all.mis.flag] <- post.x.tmp[N.prot]
  # 
  #   post.c.tmp <- post.c
  #   post.c <- matrix(0, N.prot.orig, N.rep)
  #   post.c[!all.mis.flag, ] <- post.c.tmp[-N.prot, ]
  #   post.c[all.mis.flag, ] <- post.c.tmp[N.prot, ]
  # }
  
  return(list(n.chain=n.chain, n.burnin=n.burnin, n.iter=n.iter, run.time=t1,
              mu.x.prior=mu.x.prior, sigma.x.prior=sigma.x.prior,
              mu.c.prior=mu.c.prior, sigma.c.prior=sigma.c.prior,
              frac.prior.mu.x=frac.prior.mu.x, frac.prior.sigma.x=frac.prior.sigma.x,
              frac.prior.mu.c=frac.prior.mu.c, frac.prior.sigma.c=frac.prior.sigma.c,
              post.mu.x1.mean=post.mu.x1.mean, post.mu.x1.sd=post.mu.x1.sd,
              post.mu.x2.mean=post.mu.x2.mean, post.mu.x2.sd=post.mu.x2.sd,
              post.sigma.x.mean=post.sigma.x.mean, post.sigma.x.sd=post.sigma.x.sd,
              post.mu.c.mean=post.mu.c.mean, post.mu.c.sd=post.mu.c.sd,
              post.sigma.c.mean=post.sigma.c.mean, post.sigma.c.sd=post.sigma.c.sd,
              post.p.cens1=post.p.cens1, post.p.cens2=post.p.cens2,
              post.x1=post.x1, post.x2=post.x2, post.c1=post.c1, post.c2=post.c2))
}