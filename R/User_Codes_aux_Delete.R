###############################################################################
###############################################################################
######################## ALL CODES FOR Weibull  MODEL #########################
###############################################################################
###############################################################################

# DIC
DIC.WEI=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; n=length(Time); L<-rep(0,times=N)
  # LOG-LIKELIHOOD FOR EACH DRAW
  for(iter in 1:N)
  {
    L[iter]=log.lik.WEI(Time,Cens,X,beta=as.vector(chain[iter,1:k]),gam=chain[iter,k+1],EXP)
  }
  aux=apply(chain[,1:(k+1)],2,"median"); pd=-2*mean(L)+2*log.lik.WEI(Time,Cens,X,beta=aux[1:k],gam=aux[k+1],EXP)
  if(EXP==FALSE){pd.aux=k+1}
  else(pd.aux=k)
  DIC=-2*mean(L)+pd

  print(pd); print(pd.aux); return(DIC)
}

# LOG-MARGINAL LIKELIHOOD ESTIMATOR
LML.WEI<-function(thin,Time,Cens,X,chain,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN
  #######################################################################
  chain=as.matrix(chain)
  # SAMPLE SIZE, NUMBER OF DRAWS AND NUMBER OF REGRESSORS
  n=length(t); N=dim(chain)[1]; k=dim(X)[2]
  # SETTING PROPOSAL VARIANCES AS THE MEDIAN VALUE (BASED ON ADAPTIVE VERSION)
  if(EXP==FALSE) {omega2.gam=exp(median(chain[,2*k+2]))}
  else {omega2.gam=1}
  if(k>1)	{omega2.beta=exp(apply(chain[,(k+2):(2*k+1)],2,"median"))}
  else		{omega2.beta=exp(median(chain[,3]))}

  # NON-ADAPTIVE RUN OF THE CHAIN
  chain.nonadapt=MCMC.WEI.NonAdapt(N=0.5*N*thin,thin=thin,beta0=t(chain[N,1:k]),gam0=chain[N,k+1],Time=Time,Cens=Cens,X=X,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,EXP)
  chain.nonadapt=chain.nonadapt[-1,]
  gam.star=median(chain.nonadapt[,k+1])
  if(k>1)	{beta.star=apply(chain.nonadapt[,1:k],2,"median")}
  else		{beta.star=median(chain.nonadapt[,1])}
  N.aux=dim(chain.nonadapt)[1]

  # LIKELIHOOD ORDINATE
  LL.ord=log.lik.WEI(Time,Cens,X,beta=beta.star,gam=gam.star,EXP)
  print("Likelihood ordinate ready!")

  # PRIOR ORDINATE
  if(EXP==FALSE) {P.ord=dgamma(gam.star,shape=hyp1.gam,rate=hyp2.gam)}
  else {P.ord=1}
  print("Prior ordinate ready!")

  if(EXP==FALSE)
  {
    # POSTERIOR ORDINATE - gam
    chain.gam=MCMCR.WEI.gam(N=0.5*N*thin,thin=thin,beta0=t(chain.nonadapt[N.aux,1:k]),gam0=gam.star,Time=Time,Cens=Cens,X=X,ar=0.44,omega2.beta)
    chain.gam=chain.gam[-1,]
    po1.gam=rep(0,times=N.aux); po2.gam=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      po1.gam[i]=alpha.gam(gam.0=chain.nonadapt[i,k+1],gam.1=gam.star,Time,Cens,X,beta=as.vector(chain.nonadapt[i,1:k]),theta=NA,lambda=rep(1,times=n),typ.theta=NA,hyp.theta=NA,hyp1.gam,hyp2.gam,mixing="None",lower.bound=0.06)*dnorm(x=gam.star,mean=chain.nonadapt[i,k+1],sd=sqrt(omega2.gam))
      gam.aux=rnorm(n=1, mean=gam.star,sd=sqrt(omega2.gam))
      if(gam.aux<=0.06) {po2.gam[i]=0}
      else              {po2.gam[i]=alpha.gam(gam.0=gam.star,gam.1=gam.aux,Time,Cens,X,beta=as.vector(chain.gam[i,1:k]),theta=NA,lambda=rep(1,times=n),typ.theta=NA,hyp.theta=NA,hyp1.gam,hyp2.gam,mixing="None",lower.bound=0.06)}
    }
    PO.gam=mean(po1.gam)/mean(po2.gam)
    print("Posterior ordinate gamma ready!")
  }

  # POSTERIOR ORDINATE - beta
  if(EXP==FALSE) {chain.prev=chain.gam}
  else {chain.prev=chain.nonadapt}
  PO.beta=rep(0,times=k)
  for(j.beta in 0:(k-1))
  {
    print(j.beta)
    beta0=t(chain.prev[N.aux,1:k]); beta0[j.beta+1]=beta.star[j.beta+1]
    chain.next=MCMCR.WEI.betaJ.gam(N=0.5*N*thin,thin=thin,beta0=beta0,gam0=gam.star,Time,Cens,X,J=j.beta+1,ar=0.44,omega2.beta)
    chain.next=chain.next[-1,]
    po1.beta=rep(0,times=N.aux); po2.beta=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      beta.0=as.vector(t(chain.prev[i,1:k]))
      beta.1=beta.0; beta.1[j.beta+1]=beta.star[j.beta+1]
      po1.beta[i]=alpha.beta.j(beta.0=beta.0,beta.1=beta.1,gam=chain.prev[i,k+1],Time=Time,Cens=Cens,X=X,lambda=rep(1,times=n),j=j.beta+1)*dnorm(x=beta.star[j.beta+1],mean=as.numeric(chain.prev[i,j.beta+1]),sd=sqrt(omega2.beta)[j.beta+1])
      betaj.aux=rnorm(n=1, mean=beta.star[j.beta+1],sd=sqrt(omega2.beta)[j.beta+1])
      beta.2=beta.star; beta.2[j.beta+1]=betaj.aux
      po2.beta[i]=alpha.beta.j(beta.0=beta.star,beta.1=beta.2,gam=chain.next[i,k+1],Time=Time,Cens=Cens,X=X,lambda=rep(1,times=n),j=j.beta+1)
    }
    PO.beta[j.beta+1]=mean(po1.beta)/mean(po2.beta)
    chain.prev=chain.next
  }
  print("Posterior ordinate beta ready!")

  # TAKING LOGARITHM
  if(EXP==FALSE) {LPO.gam=log(PO.gam)}
  else {LPO.gam=0}
  LP.ord=log(P.ord); LPO.beta=log(PO.beta)

  # LOG-MARGINAL LIKELIHOOD
  LML=LL.ord+LP.ord-LPO.gam-sum(LPO.beta)

  list("LL.ord"=LL.ord, "LP.ord"=LP.ord, "LPO.gam"=LPO.gam, "LPO.beta"=sum(LPO.beta), "LML"=LML)
}

# CASE DELETION ANALYSIS
CaseDeletion.WEI=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; n=dim(X)[1]; k=dim(X)[2]; logCPO=rep(0,times=n); KL.aux=rep(0,times=n)
  for(s in 1:n)
  {
    aux1=rep(0,times=N); aux2=rep(0,times=N)
    # OFSET USED TO AVOID OVER/UNDERFLOWS
    aux0=log.lik.WEI(Time=Time[s],Cens=Cens[s],X=X[s,],beta=as.vector(chain[1,1:k]),gam=chain[1,k+1],EXP)
    for(ITER in 1:N)
    {
      aux2[ITER]=log.lik.WEI(Time=Time[s],Cens=Cens[s],X=X[s,],beta=as.vector(chain[ITER,1:k]),gam=chain[ITER,k+1],EXP)
      aux1[ITER]=exp(-aux2[ITER]+aux0)
    }
    logCPO[s]=-log(mean(aux1))+aux0
    KL.aux[s]=mean(aux2)
  }
  KL=KL.aux-logCPO
  CALIBRATION=0.5*(1+sqrt(1-exp(-2*KL)))
  return(cbind(logCPO,KL,CALIBRATION))
}

###############################################################################
###############################################################################
######################## ALL CODES FOR RMWEXP   MODEL #########################
###############################################################################
###############################################################################

# DIC
DIC.RMWEXP=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; n=length(Time); L<-rep(0,times=N)
  # LOG-LIKELIHOOD FOR EACH DRAW
  for(iter in 1:N)
  {
    L[iter]=log.lik.RMWEXP(Time,Cens,X,beta=as.vector(chain[iter,1:k]),gam=chain[iter,k+1],EXP)
  }
  aux=apply(chain[,1:(k+1)],2,"median"); pd=-2*mean(L)+2*log.lik.RMWEXP(Time,Cens,X,beta=aux[1:k],gam=aux[k+1],EXP)
  if(EXP==FALSE){pd.aux=k+1}
  else(pd.aux=k)
  DIC=-2*mean(L)+pd
  print(pd); print(pd.aux); return(DIC)
}

# LOG-MARGINAL LIKELIHOOD ESTIMATOR
LML.RMWEXP<-function(thin,Q,Time,Cens,X,chain,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  # SAMPLE SIZE, NUMBER OF DRAWS AND NUMBER OF REGRESSORS
  n=length(Time); N=dim(chain)[1]; k=dim(X)[2]
  # SETTING PROPOSAL VARIANCES AS THE MEDIAN VALUE (BASED ON ADAPTIVE VERSION)
  if(EXP==FALSE) {omega2.gam=exp(median(chain[,n+2*k+2]))}
  else {omega2.gam=1}
  if(k>1)	{omega2.beta=exp(apply(chain[,(n+k+2):(n+2*k+1)],2,"median"))}
  else		{omega2.beta=exp(median(chain[,n+3]))}

  # NON-ADAPTIVE RUN OF THE CHAIN
  chain.nonadapt=MCMC.RMWEXP.NonAdapt(N=0.5*N*thin,thin=thin,Q,beta0=t(chain[N,1:k]),gam0=chain[N,k+1],Time,Cens,X,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,EXP)
  chain.nonadapt=chain.nonadapt[-1,]
  gam.star=median(chain.nonadapt[,k+1])
  if(k>1)	{beta.star=apply(chain.nonadapt[,1:k],2,"median")}
  else		{beta.star=median(chain.nonadapt[,1])}
  N.aux=dim(chain.nonadapt)[1]

  # LIKELIHOOD ORDINATE
  LL.ord=log.lik.RMWEXP(Time,Cens,X,beta=beta.star,gam=gam.star,EXP)
  print("Likelihood ordinate ready!")

  # PRIOR ORDINATE
  if(EXP==FALSE) {P.ord=dgamma(gam.star,shape=hyp1.gam,rate=hyp2.gam)}
  else {P.ord=1}
  print("Prior ordinate ready!")

  if(EXP==FALSE)
  {
    # POSTERIOR ORDINATE - gam
    chain.gam=MCMCR.RMWEXP.gam(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.nonadapt[N.aux,1:k]),gam0=gam.star,lambda0=t(chain.nonadapt[N.aux,(k+2):(k+n+1)]),Time,Cens,X,ar=0.44,omega2.beta)
    chain.gam=chain.gam[-1,]
    po1.gam=rep(0,times=N.aux); po2.gam=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      po1.gam[i]=alpha.gam(gam.0=chain.nonadapt[i,k+1],gam.1=gam.star,Time,Cens,X,beta=as.vector(chain.nonadapt[i,1:k]),theta=NA,lambda=t(chain.nonadapt[i,(k+2):(k+n+1)]),typ.theta=NA,hyp.theta=NA,hyp1.gam,hyp2.gam,mixing="None",lower.bound=0.06)*dnorm(x=gam.star,mean=chain.nonadapt[i,k+1],sd=sqrt(omega2.gam))
      gam.aux=rnorm(n=1, mean=gam.star,sd=sqrt(omega2.gam))
      if(gam.aux<=0.06) {po2.gam[i]=0}
      else              {po2.gam[i]=alpha.gam(gam.0=gam.star,gam.1=gam.aux,Time,Cens,X,beta=as.vector(chain.gam[i,1:k]),theta=NA,lambda=t(chain.gam[i,(k+2):(k+n+1)]),typ.theta=NA,hyp.theta=NA,hyp1.gam,hyp2.gam,mixing="None",lower.bound=0.06)}
    }
    PO.gam=mean(po1.gam)/mean(po2.gam)
    print("Posterior ordinate gamma ready!")
  }

  # POSTERIOR ORDINATE - beta
  if(EXP==FALSE) {chain.prev=chain.gam}
  else {chain.prev=chain.nonadapt}
  PO.beta=rep(0,times=k)

  for(j.beta in 0:(k-1))
  {
    print(j.beta)
    beta0=t(chain.prev[N.aux,1:k]); beta0[j.beta+1]=beta.star[j.beta+1]
    chain.next=MCMCR.RMWEXP.betaJ.gam(N=0.5*N*thin,thin=thin,Q,beta0=beta0,gam0=gam.star,lambda0=t(chain.prev[N.aux,(k+2):(k+n+1)]),Time,Cens,X,J=j.beta+1,ar=0.44,omega2.beta)
    chain.next=chain.next[-1,]

    po1.beta=rep(0,times=N.aux); po2.beta=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      beta.0=as.vector(t(chain.prev[i,1:k]))
      beta.1=beta.0; beta.1[j.beta+1]=beta.star[j.beta+1]
      po1.beta[i]=alpha.beta.j(beta.0=beta.0,beta.1=beta.1,gam=chain.prev[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.prev[i,(k+2):(k+n+1)]),j=j.beta+1)*dnorm(x=beta.star[j.beta+1],mean=as.numeric(chain.prev[i,j.beta+1]),sd=sqrt(omega2.beta)[j.beta+1])
      betaj.aux=rnorm(n=1, mean=beta.star[j.beta+1],sd=sqrt(omega2.beta)[j.beta+1])
      beta.2=beta.star; beta.2[j.beta+1]=betaj.aux
      po2.beta[i]=alpha.beta.j(beta.0=beta.star,beta.1=beta.2,gam=chain.next[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.next[i,(k+2):(k+n+1)]),j=j.beta+1)
    }
    PO.beta[j.beta+1]=mean(po1.beta)/mean(po2.beta)

    chain.prev=chain.next
  }
  print("Posterior ordinate beta ready!")

  # TAKING LOGARITHM
  if(EXP==FALSE) {LPO.gam=log(PO.gam)}
  else {LPO.gam=0}
  LP.ord=log(P.ord); LPO.beta=log(PO.beta)

  # MARGINAL LOG-LIKELIHOOD
  LML=LL.ord+LP.ord-LPO.gam-sum(LPO.beta)

  list("LL.ord"=LL.ord, "LP.ord"=LP.ord, "LPO.gam"=LPO.gam, "LPO.beta"=sum(LPO.beta), "LML"=LML)

}

# CASE DELETION ANALYSIS
CaseDeletion.RMWEXP=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; logCPO=rep(0,times=n); KL.aux=rep(0,times=n)
  for(s in 1:n)
  {
    aux1=rep(0,times=N); aux2=rep(0,times=N)
    for(ITER in 1:N)
    {
      aux2[ITER]=log.lik.RMWEXP(Time=Time[s],Cens=Cens[s],X=X[s,],beta=as.vector(chain[ITER,1:k]),gam=chain[ITER,k+1],EXP)
      aux1[ITER]=exp(-aux2[ITER])
    }
    logCPO[s]=-log(mean(aux1))
    KL.aux[s]=mean(aux2)
  }
  KL=KL.aux-logCPO
  CALIBRATION=0.5*(1+sqrt(1-exp(-2*KL)))
  return(cbind(logCPO,KL,CALIBRATION))
}


# OUTLIER DETECTION FOR A SPECIFIC OBSERVATION I
BF.lambda.obs.RMWEXP<-function(ref,obs,Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N<-dim(chain)[1]; k<-dim(X)[2]; aux1<-rep(0,times=N)
  for(j in 1:N)
  {
    RATE.ref=as.numeric(exp(-chain[j,k+1]*X%*%as.vector(chain[j,1:k])))*ref
    aux1[j]=exp(log.lik.WEI.ref1(Time=Time[obs],Cens=Cens[obs],RATE=RATE.ref[obs],gam=chain[j,k+1],EXP)-log.lik.RMWEXP(Time[obs],Cens[obs],X[obs,],beta=as.vector(chain[j,1:k]),gam=chain[j,k+1],EXP))
  }
  aux=mean(aux1)
  return(aux)
}


###############################################################################
###############################################################################
######################## ALL CODES FOR RMWGAM   MODEL #########################
###############################################################################
###############################################################################

# DIC
DIC.RMWGAM=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; n=length(Time); L<-rep(0,times=N)
  # LOG-LIKELIHOOD FOR EACH DRAW
  for(iter in 1:N)
  {
    L[iter]=log.lik.RMWGAM(Time,Cens,X,beta=as.vector(chain[iter,1:k]),gam=chain[iter,k+1],theta=chain[iter,k+2],EXP)
  }

  aux=apply(chain[,1:(k+2)],2,"median"); pd=-2*mean(L)+2*log.lik.RMWGAM(Time,Cens,X,beta=aux[1:k],gam=aux[k+1],theta=aux[k+2],EXP)
  if(EXP==FALSE){pd.aux=k+2}
  else(pd.aux=k+1)
  DIC=-2*mean(L)+pd

  print(pd); print(pd.aux); return(DIC)
}

LML.RMWGAM<-function(thin,Q,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  # SAMPLE SIZE, NUMBER OF DRAWS AND NUMBER OF REGRESSORS
  n=length(Time); N=dim(chain)[1]; k=dim(X)[2]
  # SETTING PROPOSAL VARIANCES AS THE MEDIAN VALUE (BASED ON ADAPTIVE VERSION)
  if(EXP==FALSE) {omega2.gam=exp(median(chain[,n+2*k+3]))}
  else {omega2.gam=1}
  if(k>1)	{omega2.beta=exp(apply(chain[,(k+n+3):(2*k+n+2)],2,"median"))}
  else		{omega2.beta=exp(median(chain[,n+4]))}
  omega2.theta=exp(median(chain[,2*k+n+4]))

  # NON-ADAPTIVE RUN OF THE CHAIN
  chain.nonadapt=MCMC.RMWGAM.NonAdapt(N=0.5*N*thin,thin=thin,Q,beta0=t(chain[N,1:k]),gam0=chain[N,k+1],theta0=chain[N,k+2],Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,omega2.theta,EXP)
  chain.nonadapt=chain.nonadapt[-1,]
  gam.star=median(chain.nonadapt[,k+1]); theta.star=median(chain.nonadapt[,k+2])
  if(k>1)	{beta.star=apply(chain.nonadapt[,1:k],2,"median")}
  else		{beta.star=median(chain.nonadapt[,1])}
  N.aux=dim(chain.nonadapt)[1]

  # LIKELIHOOD ORDINATE
  LL.ord=log.lik.RMWGAM(Time,Cens,X,beta=beta.star,gam=gam.star,theta=theta.star,EXP)
  print("Likelihood ordinate ready!")

  # PRIOR ORDINATE
  if(EXP==FALSE) {LP.ord=dgamma(gam.star,shape=hyp1.gam,rate=hyp2.gam,log=TRUE)+log.prior.theta(theta=theta.star,gam=gam.star,a=hyp.theta,type=typ.theta,mixing="Gamma")}
  else {LP.ord=log.prior.theta(theta=theta.star,gam=gam.star,a=hyp.theta,type=typ.theta,mixing="Gamma")}
  print("Prior ordinate ready!")

  # POSTERIOR ORDINATE - theta
  chain.theta=MCMCR.RMWGAM.theta(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.nonadapt[N.aux,1:k]),gam0=chain.nonadapt[N.aux,k+1],theta0=theta.star,lambda0=t(chain.nonadapt[N.aux,(k+3):(k+n+2)]),Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,EXP)
  chain.theta=chain.theta[-1,]
  po1.theta=rep(0,times=N.aux); po2.theta=rep(0,times=N.aux)
  for(i in 1:N.aux)
  {
    po1.theta[i]=alphaRMWGAM.theta(theta.0=chain.nonadapt[i,k+2],theta.1=theta.star,gam=chain.nonadapt[i,k+1],lambda=t(chain.nonadapt[i,(k+3):(k+n+2)]),typ.theta,hyp.theta)*dnorm(x=theta.star,mean=chain.nonadapt[i,k+2],sd=sqrt(omega2.theta))
    theta.aux=rnorm(n=1, mean=theta.star,sd=sqrt(omega2.theta))
    if(theta.aux<=2/gam.star) {po2.theta[i]=0}
    else             {po2.theta[i]=alphaRMWGAM.theta(theta.0=theta.star,theta.1=theta.aux,gam=chain.theta[i,k+1],lambda=t(chain.theta[i,(k+3):(k+n+2)]),typ.theta,hyp.theta)}
  }
  PO.theta=mean(po1.theta)/mean(po2.theta)
  print("Posterior ordinate theta ready!")

  if(EXP==FALSE)
  {
    # POSTERIOR ORDINATE - gam
    chain.gam=MCMCR.RMWGAM.gam.theta(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.theta[N.aux,1:k]),gam0=gam.star,theta0=theta.star,lambda0=t(chain.theta[N.aux,(k+3):(k+n+2)]),Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
    chain.gam=chain.gam[-1,]
    po1.gam=rep(0,times=N.aux); po2.gam=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      po1.gam[i]=alpha.gam(gam.0=chain.theta[i,k+1],gam.1=gam.star,Time,Cens,X,beta=as.vector(chain.theta[i,1:k]),theta=theta.star,lambda=t(chain.theta[i,(k+3):(k+n+2)]),typ.theta,hyp.theta,hyp1.gam,hyp2.gam,mixing="Gamma",lower.bound=0.06)*dnorm(x=gam.star,mean=chain.theta[i,k+1],sd=sqrt(omega2.gam))
      gam.aux=rnorm(n=1, mean=gam.star,sd=sqrt(omega2.gam))
      if(gam.aux<=0.06 | gam.aux<=2/theta.star) {po2.gam[i]=0} # FIXED? FIX THIS, IT MUST BE ZERO FOR LESS THAN 2/THETA AS WELL.
      else              {po2.gam[i]=alpha.gam(gam.0=gam.star,gam.1=gam.aux,Time,Cens,X,beta=as.vector(chain.gam[i,1:k]),theta=theta.star,lambda=t(chain.gam[i,(k+3):(k+n+2)]),typ.theta,hyp.theta,hyp1.gam,hyp2.gam,mixing="Gamma",lower.bound=0.06)}
    }
    PO.gam=mean(po1.gam)/mean(po2.gam)
    print("Posterior ordinate gamma ready!")
  }

  # POSTERIOR ORDINATE - beta
  if(EXP==FALSE) {chain.prev=chain.gam}
  else {chain.prev=chain.theta}
  PO.beta=rep(0,times=k)

  for(j.beta in 0:(k-1))
  {
    print(j.beta)
    beta0=t(chain.prev[N.aux,1:k]); beta0[j.beta+1]=beta.star[j.beta+1]
    chain.next=MCMCR.RMWGAM.betaJ.gam.theta(N=0.5*N*thin,thin=thin,Q,beta0=beta0,gam0=gam.star,theta0=theta.star,lambda0=t(chain.prev[N.aux,(k+3):(k+n+2)]),Time,Cens,X,J=j.beta+1,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
    chain.next=chain.next[-1,]

    po1.beta=rep(0,times=N.aux); po2.beta=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      beta.0=as.vector(t(chain.prev[i,1:k]))
      beta.1=beta.0; beta.1[j.beta+1]=beta.star[j.beta+1]
      po1.beta[i]=alpha.beta.j(beta.0=beta.0,beta.1=beta.1,gam=chain.prev[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.prev[i,(k+3):(k+n+2)]),j=j.beta+1)*dnorm(x=beta.star[j.beta+1],mean=as.numeric(chain.prev[i,j.beta+1]),sd=sqrt(omega2.beta)[j.beta+1])
      betaj.aux=rnorm(n=1, mean=beta.star[j.beta+1],sd=sqrt(omega2.beta)[j.beta+1])
      beta.2=beta.star; beta.2[j.beta+1]=betaj.aux
      po2.beta[i]=alpha.beta.j(beta.0=beta.star,beta.1=beta.2,gam=chain.next[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.next[i,(k+3):(k+n+2)]),j=j.beta+1)
    }
    PO.beta[j.beta+1]=mean(po1.beta)/mean(po2.beta)

    chain.prev=chain.next
  }
  print("Posterior ordinate beta ready!")

  # TAKING LOGARITHM
  if(EXP==FALSE) {LPO.gam=log(PO.gam)}
  else {LPO.gam=0}
  LPO.theta=log(PO.theta); LPO.beta=log(PO.beta)

  # MARGINAL LOG-LIKELIHOOD
  LML=LL.ord+LP.ord-LPO.theta-LPO.gam-sum(LPO.beta)

  list("LL.ord"=LL.ord, "LP.ord"=LP.ord, "LPO.theta"=LPO.theta, "LPO.gam"=LPO.gam, "LPO.beta"=sum(LPO.beta), "LML"=LML)

}

CaseDeletion.RMWGAM=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; logCPO=rep(0,times=n); KL.aux=rep(0,times=n)

  for(s in 1:n)
  {
    aux1=rep(0,times=N); aux2=rep(0,times=N)
    for(ITER in 1:N)
    {
      aux2[ITER]=log.lik.RMWGAM(Time=Time[s],Cens=Cens[s],X=X[s,],beta=as.vector(chain[ITER,1:k]),gam=chain[ITER,k+1],theta=chain[ITER,k+2],EXP)
      aux1[ITER]=exp(-aux2[ITER])
    }
    logCPO[s]=-log(mean(aux1))
    KL.aux[s]=mean(aux2)
  }
  KL=KL.aux-logCPO
  CALIBRATION=0.5*(1+sqrt(1-exp(-2*KL)))
  return(cbind(logCPO,KL,CALIBRATION))
}

################################################################################
# OUTLIER DETECTION
BF.lambda.obs.RMWGAM<-function(ref,obs,N,thin,Q,burn,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  aux1=Post.lambda.obs.RMWGAM(ref,obs,Time,Cens,X,chain)
  aux2=CFP.obs.RMWGAM(ref,obs,N,thin,Q,burn,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar,EXP)
  return(aux1*aux2)
}


###############################################################################
###############################################################################
######################## ALL CODES FOR RMWIGAM  MODEL #########################
###############################################################################
###############################################################################

MCMC.RMWIGAM<-function(N,thin,Q,beta0,gam0=1,theta0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  # NUMBER OF REGRESSORS, SAMPLE SIZE AND NUMBER OF STORED DRAWS
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  # INITIALIZATION OF ELEMENTS WHERE DRAWS ARE STORED
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1)
  # WHERE EMPIRICAL ACCEPTANCE PROBABILITIES ARE STORED
  accept.beta=rep(0,times=k); pbeta.aux=rep(0,times=k); accept.theta=0; ptheta.aux=0; accept.gam=0; pgam.aux=0
  # WHERE ADAPTIVE PROPOSAL VARIANCES ARE STORED (LOG-SCALE)
  ls.beta=matrix(rep(0,times=(N.aux+1)*k),ncol=k); ls.gam=rep(0,times=N.aux+1); ls.theta=rep(0,times=N.aux+1)

  # INITIALIZATION OF AUXILIARY PARAMETERS
  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.beta.aux=ls.beta[1,]; ls.gam.aux=ls.gam[1]; ls.theta.aux=ls.theta[1]

  # BATCH INDICATOR FOR ADAPTIVE STEPS (WHEN TO ADAPT)
  i_batch=0;

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    # DRAWS FROM BETA
    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=exp(ls.beta.aux)[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1; pbeta.aux[j]=pbeta.aux[j]+1}
    }

    # DRAWS FROM GAMMA (IF NOT FIXED)
    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMW.gam(N=1,omega2=exp(ls.gam.aux),gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="InvGamma",lower.bound=0.06)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1; pgam.aux=pgam.aux+1}
    }

    # DRAWS FROM THETA
    MH.theta=GRWMH.RMWIGAM.theta(N=1,omega2=exp(ls.theta.aux),theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
    theta.aux<-MH.theta$theta
    if(MH.theta$ind==1) {accept.theta=accept.theta+1; ptheta.aux=ptheta.aux+1}

    # DRAWS FROM MIXING PARAMETERS
    if((iter-1)%%Q==0)
    {
      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      for(ind in 1:n)
      {
        lambda.aux[ind]=rgig(n=1,lambda=-theta.aux+Cens[ind],chi=2,psi=psi.aux[ind])
      }
    }

    # ADAPTIVE STEP
    if(i_batch==50)
    {
      pbeta.aux=pbeta.aux/50; Pbeta.aux=as.numeric(pbeta.aux<rep(ar,times=k))
      ls.beta.aux=ls.beta.aux+((-1)^Pbeta.aux)*min(0.01,1/sqrt(iter))
      pgam.aux=pgam.aux/50; Pgam.aux=as.numeric(pgam.aux<ar)
      ls.gam.aux=ls.gam.aux+((-1)^Pgam.aux)*min(0.01,1/sqrt(iter))
      ptheta.aux=ptheta.aux/50; Ptheta.aux=as.numeric(ptheta.aux<ar)
      ls.theta.aux=ls.theta.aux+((-1)^Ptheta.aux)*min(0.01,1/sqrt(iter))
      i_batch=0; pbeta.aux=rep(0,times=k); pgam.aux=0; ptheta.aux=0;
    }

    # DRAWS STORAGE
    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; theta[iter/thin+1]<-theta.aux; lambda[iter/thin+1,]<-lambda.aux
      ls.beta[iter/thin+1,]=ls.beta.aux; ls.gam[iter/thin+1]=ls.gam.aux; ls.theta[iter/thin+1]=ls.theta.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))
  print(paste("AR gamma :",round(accept.gam/N,2)))
  print(paste("AR theta :",round(accept.theta/N,2)))

  chain=cbind(beta,gam,theta,lambda,ls.beta,ls.gam,ls.theta)
  return(chain)
}

DIC.RMWIGAM=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; n=length(Time); L<-rep(0,times=N)
  # LOG-LIKELIHOOD FOR EACH DRAW
  for(iter in 1:N)
  {
    L[iter]=log.lik.RMWIGAM(Time,Cens,X,beta=as.vector(chain[iter,1:k]),gam=chain[iter,k+1],theta=chain[iter,k+2],EXP)
  }

  aux=apply(chain[,1:(k+2)],2,"median"); pd=-2*mean(L)+2*log.lik.RMWIGAM(Time,Cens,X,beta=aux[1:k],gam=aux[k+1],theta=aux[k+2],EXP)
  if(EXP==FALSE){pd.aux=k+2}
  else(pd.aux=k+1)
  DIC=-2*mean(L)+pd

  print(pd); print(pd.aux); return(DIC)
}

LML.RMWIGAM<-function(thin,Q,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  # SAMPLE SIZE, NUMBER OF DRAWS AND NUMBER OF REGRESSORS
  n=length(Time); N=dim(chain)[1]; k=dim(X)[2]
  # SETTING PROPOSAL VARIANCES AS THE MEDIAN VALUE (BASED ON ADAPTIVE VERSION)
  if(EXP==FALSE) {omega2.gam=exp(median(chain[,n+2*k+3]))}
  else {omega2.gam=1}
  if(k>1)	{omega2.beta=exp(apply(chain[,(k+n+3):(2*k+n+2)],2,"median"))}
  else		{omega2.beta=exp(median(chain[,n+4]))}
  omega2.theta=exp(median(chain[,2*k+n+4]))

  # NON-ADAPTIVE RUN OF THE CHAIN
  chain.nonadapt=MCMC.RMWIGAM.NonAdapt(N=0.5*N*thin,thin=thin,Q,beta0=t(chain[N,1:k]),gam0=chain[N,k+1],theta0=chain[N,k+2],Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,omega2.theta,EXP)
  chain.nonadapt=chain.nonadapt[-1,]
  gam.star=median(chain.nonadapt[,k+1]); theta.star=median(chain.nonadapt[,k+2])
  if(k>1)	{beta.star=apply(chain.nonadapt[,1:k],2,"median")}
  else		{beta.star=median(chain.nonadapt[,1])}
  N.aux=dim(chain.nonadapt)[1]

  # LIKELIHOOD ORDINATE
  LL.ord=log.lik.RMWIGAM(Time,Cens,X,beta=beta.star,gam=gam.star,theta=theta.star,EXP)
  print("Likelihood ordinate ready!")

  # PRIOR ORDINATE
  if(EXP==FALSE) {LP.ord=dgamma(gam.star,shape=hyp1.gam,rate=hyp2.gam,log=TRUE)+log.prior.theta(theta=theta.star,gam=gam.star,a=hyp.theta,type=typ.theta,mixing="InvGamma")}
  else {LP.ord=log.prior.theta(theta=theta.star,gam=gam.star,a=hyp.theta,type=typ.theta,mixing="InvGamma")}
  print("Prior ordinate ready!")

  # POSTERIOR ORDINATE - theta
  chain.theta=MCMCR.theta.RMWIGAM(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.nonadapt[N.aux,1:k]),gam0=chain.nonadapt[N.aux,k+1],theta0=theta.star,lambda0=t(chain.nonadapt[N.aux,(k+3):(k+n+2)]),Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,EXP)
  chain.theta=chain.theta[-1,]
  po1.theta=rep(0,times=N.aux); po2.theta=rep(0,times=N.aux)
  for(i in 1:N.aux)
  {
    po1.theta[i]=alphaRMWIGAM.theta(theta.0=chain.nonadapt[i,k+2],theta.1=theta.star,gam=chain.nonadapt[i,k+1],lambda=t(chain.nonadapt[i,(k+3):(k+n+2)]),typ.theta,hyp.theta)*dnorm(x=theta.star,mean=chain.nonadapt[i,k+2],sd=sqrt(omega2.theta))
    theta.aux=rnorm(n=1, mean=theta.star,sd=sqrt(omega2.theta))
    if(theta.aux<=1) {po2.theta[i]=0}
    else             {po2.theta[i]=alphaRMWIGAM.theta(theta.0=theta.star,theta.1=theta.aux,gam=chain.theta[i,k+1],lambda=t(chain.theta[i,(k+3):(k+n+2)]),typ.theta,hyp.theta)}
  }
  PO.theta=mean(po1.theta)/mean(po2.theta)
  print("Posterior ordinate theta ready!")

  if(EXP==FALSE)
  {
    # POSTERIOR ORDINATE - gam
    chain.gam=MCMCR.gam.theta.RMWIGAM(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.theta[N.aux,1:k]),gam0=gam.star,theta0=theta.star,lambda0=t(chain.theta[N.aux,(k+3):(k+n+2)]),Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
    chain.gam=chain.gam[-1,]
    po1.gam=rep(0,times=N.aux); po2.gam=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      po1.gam[i]=alpha.gam(gam.0=chain.theta[i,k+1],gam.1=gam.star,Time,Cens,X,beta=as.vector(chain.theta[i,1:k]),theta=theta.star,lambda=t(chain.theta[i,(k+3):(k+n+2)]),typ.theta,hyp.theta,hyp1.gam,hyp2.gam,mixing="InvGamma",lower.bound=0.06)*dnorm(x=gam.star,mean=chain.theta[i,k+1],sd=sqrt(omega2.gam))
      gam.aux=rnorm(n=1, mean=gam.star,sd=sqrt(omega2.gam))
      if(gam.aux<=0.06) {po2.gam[i]=0}
      else              {po2.gam[i]=alpha.gam(gam.0=gam.star,gam.1=gam.aux,Time,Cens,X,beta=as.vector(chain.gam[i,1:k]),theta=theta.star,lambda=t(chain.gam[i,(k+3):(k+n+2)]),typ.theta,hyp.theta,hyp1.gam,hyp2.gam,mixing="InvGamma",lower.bound=0.06)}
    }
    PO.gam=mean(po1.gam)/mean(po2.gam)
    print("Posterior ordinate gamma ready!")
  }

  # POSTERIOR ORDINATE - beta
  if(EXP==FALSE) {chain.prev=chain.gam}
  else {chain.prev=chain.theta}
  PO.beta=rep(0,times=k)
  for(j.beta in 0:(k-1))
  {
    print(j.beta)
    beta0=t(chain.prev[N.aux,1:k]); beta0[j.beta+1]=beta.star[j.beta+1]
    chain.next=MCMCR.betaJ.gam.theta.RMWIGAM(N=0.5*N*thin,thin=thin,Q,beta0=beta0,gam0=gam.star,theta0=theta.star,lambda0=t(chain.prev[N.aux,(k+3):(k+n+2)]),Time,Cens,X,J=j.beta+1,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
    chain.next=chain.next[-1,]

    po1.beta=rep(0,times=N.aux); po2.beta=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      beta.0=as.vector(t(chain.prev[i,1:k]))
      beta.1=beta.0; beta.1[j.beta+1]=beta.star[j.beta+1]
      po1.beta[i]=alpha.beta.j(beta.0=beta.0,beta.1=beta.1,gam=chain.prev[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.prev[i,(k+3):(k+n+2)]),j=j.beta+1)*dnorm(x=beta.star[j.beta+1],mean=as.numeric(chain.prev[i,j.beta+1]),sd=sqrt(omega2.beta)[j.beta+1])
      betaj.aux=rnorm(n=1, mean=beta.star[j.beta+1],sd=sqrt(omega2.beta)[j.beta+1])
      beta.2=beta.star; beta.2[j.beta+1]=betaj.aux
      po2.beta[i]=alpha.beta.j(beta.0=beta.star,beta.1=beta.2,gam=chain.next[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.next[i,(k+3):(k+n+2)]),j=j.beta+1)
    }
    PO.beta[j.beta+1]=mean(po1.beta)/mean(po2.beta)

    chain.prev=chain.next
  }
  print("Posterior ordinate beta ready!")

  # TAKING LOGARITHM
  if(EXP==FALSE) {LPO.gam=log(PO.gam)}
  else {LPO.gam=0}
  LPO.theta=log(PO.theta); LPO.beta=log(PO.beta)

  # MARGINAL LOG-LIKELIHOOD
  MLL=LL.ord+LP.ord-LPO.theta-LPO.gam-sum(LPO.beta)

  list("LL.ord"=LL.ord, "LP.ord"=LP.ord, "LPO.theta"=LPO.theta, "LPO.gam"=LPO.gam, "LPO.beta"=sum(LPO.beta), "MLL"=MLL)

}

CaseDeletion.RMWIGAM=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; logCPO=rep(0,times=n); KL.aux=rep(0,times=n)

  for(s in 1:n)
  {
    aux1=rep(0,times=N); aux2=rep(0,times=N)
    for(ITER in 1:N)
    {
      aux2[ITER]=log.lik.RMWIGAM(Time=Time[s],Cens=Cens[s],X=X[s,],beta=as.vector(chain[ITER,1:k]),gam=chain[ITER,k+1],theta=chain[ITER,k+2],EXP)
      aux1[ITER]=exp(-aux2[ITER])
    }
    logCPO[s]=-log(mean(aux1))
    KL.aux[s]=mean(aux2)
  }
  KL=KL.aux-logCPO
  CALIBRATION=0.5*(1+sqrt(1-exp(-2*KL)))
  return(cbind(logCPO,KL,CALIBRATION))
}

################################################################################
# CORRECTION FACTOR FOR OUTLIER DETECTION RMWIGAM MODELS
OD.Correction.RMWIGAM<-function(chain,Time,X)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  k=dim(X)[2]
  if(k>1)  {beta=apply(chain[,1:k],2,"median")}
  else		{beta=median(chain[,1])}
  gam=median(chain[,k+1]); theta=median(chain[,k+2])
  aux0=2*sqrt((exp(-X%*%beta)*Time)^gam)
  aux=(besselK(x=aux0,nu=-theta+1,expon.scaled = FALSE)^2)/(besselK(x=aux0,nu=-theta+2,expon.scaled = FALSE)*besselK(x=aux0,nu=-theta,expon.scaled = FALSE))
  return(aux)
}


################################################################################
# OUTLIER DETECTION
BF.lambda.obs.RMWIGAM<-function(ref,obs,N,thin,Q,burn,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  aux1=Post.lambda.obs.RMWIGAM(ref,obs,Time,Cens,X,chain)
  aux2=CFP.obs.RMWIGAM(ref,obs,N,thin,Q,burn,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar,EXP)
  return(aux1*aux2)
}


###############################################################################
###############################################################################
######################## ALL CODES FOR RMWIGAUSSMODEL #########################
###############################################################################
###############################################################################

MCMC.RMWIGAUSS<-function(N,thin,Q,beta0,gam0=1,theta0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  # NUMBER OF REGRESSORS, SAMPLE SIZE AND NUMBER OF STORED DRAWS
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  # INITIALIZATION OF ELEMENTS WHERE DRAWS ARE STORED
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1)
  # WHERE EMPIRICAL ACCEPTANCE PROBABILITIES ARE STORED
  accept.beta=rep(0,times=k); pbeta.aux=rep(0,times=k); accept.theta=0; ptheta.aux=0; accept.gam=0; pgam.aux=0
  # WHERE ADAPTIVE PROPOSAL VARIANCES ARE STORED (LOG-SCALE)
  ls.beta=matrix(rep(0,times=(N.aux+1)*k),ncol=k); ls.gam=rep(0,times=N.aux+1); ls.theta=rep(0,times=N.aux+1)

  # INITIALIZATION OF AUXILIARY PARAMETERS
  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.beta.aux=ls.beta[1,]; ls.gam.aux=ls.gam[1]; ls.theta.aux=ls.theta[1]

  # BATCH INDICATOR FOR ADAPTIVE STEPS (WHEN TO ADAPT)
  i_batch=0;

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    # DRAWS FROM BETA
    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=exp(ls.beta.aux)[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1; pbeta.aux[j]=pbeta.aux[j]+1}
    }

    # DRAWS FROM GAMMA (IF NOT FIXED)
    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMW.gam(N=1,omega2=exp(ls.gam.aux),gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="InvGauss",lower.bound=0.06)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1; pgam.aux=pgam.aux+1}
    }

    # DRAWS FROM THETA
    MH.theta=GRWMH.RMWIGAUSS.theta(N=1,omega2=exp(ls.theta.aux),theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
    theta.aux<-MH.theta$theta
    if(MH.theta$ind==1) {accept.theta=accept.theta+1; ptheta.aux=ptheta.aux+1}

    # DRAWS FROM MIXING PARAMETERS
    if((iter-1)%%Q==0)
    {
      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux+theta.aux^(-2)
      for(ind in 1:n)
      {
        lambda.aux[ind]=rgig(n=1,lambda=Cens[ind]-1/2,chi=1,psi=psi.aux[ind])
      }
    }

    # ADAPTIVE STEP
    if(i_batch==50)
    {
      pbeta.aux=pbeta.aux/50; Pbeta.aux=as.numeric(pbeta.aux<rep(ar,times=k))
      ls.beta.aux=ls.beta.aux+((-1)^Pbeta.aux)*min(0.01,1/sqrt(iter))
      pgam.aux=pgam.aux/50; Pgam.aux=as.numeric(pgam.aux<ar)
      ls.gam.aux=ls.gam.aux+((-1)^Pgam.aux)*min(0.01,1/sqrt(iter))
      ptheta.aux=ptheta.aux/50; Ptheta.aux=as.numeric(ptheta.aux<ar)
      ls.theta.aux=ls.theta.aux+((-1)^Ptheta.aux)*min(0.01,1/sqrt(iter))
      i_batch=0; pbeta.aux=rep(0,times=k); pgam.aux=0; ptheta.aux=0;
    }

    # DRAWS STORAGE
    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; theta[iter/thin+1]<-theta.aux; lambda[iter/thin+1,]<-lambda.aux
      ls.beta[iter/thin+1,]=ls.beta.aux; ls.gam[iter/thin+1]=ls.gam.aux; ls.theta[iter/thin+1]=ls.theta.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}
  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))
  print(paste("AR gamma :",round(accept.gam/N,2)))
  print(paste("AR theta :",round(accept.theta/N,2)))

  chain=cbind(beta,gam,theta,lambda,ls.beta,ls.gam,ls.theta)
  return(chain)
}

DIC.RMWIGAUSS=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; n=length(Time); L<-rep(0,times=N)
  # LOG-LIKELIHOOD FOR EACH DRAW
  for(iter in 1:N)
  {
    L[iter]=log.lik.RMWIGAUSS(Time,Cens,X,beta=as.vector(chain[iter,1:k]),gam=chain[iter,k+1],theta=chain[iter,k+2],EXP)
  }

  aux=apply(chain[,1:(k+2)],2,"median"); pd=-2*mean(L)+2*log.lik.RMWIGAUSS(Time,Cens,X,beta=aux[1:k],gam=aux[k+1],theta=aux[k+2],EXP)
  if(EXP==FALSE){pd.aux=k+2}
  else(pd.aux=k+1)
  DIC=-2*mean(L)+pd

  print(pd); print(pd.aux); return(DIC)
}

LML.RMWIGAUSS<-function(thin,Q,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  # SAMPLE SIZE, NUMBER OF DRAWS AND NUMBER OF REGRESSORS
  n=length(Time); N=dim(chain)[1]; k=dim(X)[2]
  # SETTING PROPOSAL VARIANCES AS THE MEDIAN VALUE (BASED ON ADAPTIVE VERSION)
  if(EXP==FALSE) {omega2.gam=exp(median(chain[,n+2*k+3]))}
  else {omega2.gam=1}
  if(k>1)	{omega2.beta=exp(apply(chain[,(k+n+3):(2*k+n+2)],2,"median"))}
  else		{omega2.beta=exp(median(chain[,n+4]))}
  omega2.theta=exp(median(chain[,2*k+n+4]))

  # NON-ADAPTIVE RUN OF THE CHAIN
  chain.nonadapt=MCMC.RMWIGAUSS.NonAdapt(N=0.5*N*thin,thin=thin,Q,beta0=t(chain[N,1:k]),gam0=chain[N,k+1],theta0=chain[N,k+2],Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,omega2.theta,EXP)
  chain.nonadapt=chain.nonadapt[-1,]
  gam.star=median(chain.nonadapt[,k+1]); theta.star=median(chain.nonadapt[,k+2])
  if(k>1)	{beta.star=apply(chain.nonadapt[,1:k],2,"median")}
  else		{beta.star=median(chain.nonadapt[,1])}
  N.aux=dim(chain.nonadapt)[1]

  # LIKELIHOOD ORDINATE
  LL.ord=log.lik.RMWIGAUSS(Time,Cens,X,beta=beta.star,gam=gam.star,theta=theta.star,EXP)
  print("Likelihood ordinate ready!")

  # PRIOR ORDINATE
  if(EXP==FALSE) {LP.ord=dgamma(gam.star,shape=hyp1.gam,rate=hyp2.gam,log=TRUE)+log.prior.theta(theta=theta.star,gam=gam.star,a=hyp.theta,type=typ.theta,mixing="InvGauss")}
  else {LP.ord=log.prior.theta(theta=theta.star,gam=gam.star,a=hyp.theta,type=typ.theta,mixing="InvGauss")}
  print("Prior ordinate ready!")

  # POSTERIOR ORDINATE - theta
  chain.theta=MCMCR.theta.RMWIGAUSS(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.nonadapt[N.aux,1:k]),gam0=chain.nonadapt[N.aux,k+1],theta0=theta.star,lambda0=t(chain.nonadapt[N.aux,(k+3):(k+n+2)]),Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,EXP)
  chain.theta=chain.theta[-1,]
  po1.theta=rep(0,times=N.aux); po2.theta=rep(0,times=N.aux)
  for(i in 1:N.aux)
  {
    po1.theta[i]=alphaRMWIGAUSS.theta(theta.0=chain.nonadapt[i,k+2],theta.1=theta.star,gam=chain.nonadapt[i,k+1],lambda=t(chain.nonadapt[i,(k+3):(k+n+2)]),typ.theta,hyp.theta)*dnorm(x=theta.star,mean=chain.nonadapt[i,k+2],sd=sqrt(omega2.theta))
    theta.aux=rnorm(n=1, mean=theta.star,sd=sqrt(omega2.theta))
    if(theta.aux<=0) {po2.theta[i]=0}
    else             {po2.theta[i]=alphaRMWIGAUSS.theta(theta.0=theta.star,theta.1=theta.aux,gam=chain.theta[i,k+1],lambda=t(chain.theta[i,(k+3):(k+n+2)]),typ.theta,hyp.theta)}
  }
  PO.theta=mean(po1.theta)/mean(po2.theta)
  print("Posterior ordinate theta ready!")

  if(EXP==FALSE)
  {
    # POSTERIOR ORDINATE - gam
    chain.gam=MCMCR.gam.theta.RMWIGAUSS(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.theta[N.aux,1:k]),gam0=gam.star,theta0=theta.star,lambda0=t(chain.theta[N.aux,(k+3):(k+n+2)]),Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
    chain.gam=chain.gam[-1,]
    po1.gam=rep(0,times=N.aux); po2.gam=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      po1.gam[i]=alpha.gam(gam.0=chain.theta[i,k+1],gam.1=gam.star,Time,Cens,X,beta=as.vector(chain.theta[i,1:k]),theta=theta.star,lambda=t(chain.theta[i,(k+3):(k+n+2)]),typ.theta,hyp.theta,hyp1.gam,hyp2.gam,mixing="InvGauss",lower.bound=0.06)*dnorm(x=gam.star,mean=chain.theta[i,k+1],sd=sqrt(omega2.gam))
      gam.aux=rnorm(n=1, mean=gam.star,sd=sqrt(omega2.gam))
      if(gam.aux<=0.06) {po2.gam[i]=0}
      else              {po2.gam[i]=alpha.gam(gam.0=gam.star,gam.1=gam.aux,Time,Cens,X,beta=as.vector(chain.gam[i,1:k]),theta=theta.star,lambda=t(chain.gam[i,(k+3):(k+n+2)]),typ.theta,hyp.theta,hyp1.gam,hyp2.gam,mixing="InvGauss",lower.bound=0.06)}
    }
    PO.gam=mean(po1.gam)/mean(po2.gam)
    print("Posterior ordinate gamma ready!")
  }

  # POSTERIOR ORDINATE - beta
  if(EXP==FALSE) {chain.prev=chain.gam}
  else {chain.prev=chain.theta}
  PO.beta=rep(0,times=k)
  for(j.beta in 0:(k-1))
  {
    print(j.beta)
    beta0=t(chain.prev[N.aux,1:k]); beta0[j.beta+1]=beta.star[j.beta+1]
    chain.next=MCMCR.betaJ.gam.theta.RMWIGAUSS(N=0.5*N*thin,thin=thin,Q,beta0=beta0,gam0=gam.star,theta0=theta.star,lambda0=t(chain.prev[N.aux,(k+3):(k+n+2)]),Time,Cens,X,J=j.beta+1,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
    chain.next=chain.next[-1,]

    po1.beta=rep(0,times=N.aux); po2.beta=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      beta.0=as.vector(t(chain.prev[i,1:k]))
      beta.1=beta.0; beta.1[j.beta+1]=beta.star[j.beta+1]
      po1.beta[i]=alpha.beta.j(beta.0=beta.0,beta.1=beta.1,gam=chain.prev[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.prev[i,(k+3):(k+n+2)]),j=j.beta+1)*dnorm(x=beta.star[j.beta+1],mean=as.numeric(chain.prev[i,j.beta+1]),sd=sqrt(omega2.beta)[j.beta+1])
      betaj.aux=rnorm(n=1, mean=beta.star[j.beta+1],sd=sqrt(omega2.beta)[j.beta+1])
      beta.2=beta.star; beta.2[j.beta+1]=betaj.aux
      po2.beta[i]=alpha.beta.j(beta.0=beta.star,beta.1=beta.2,gam=chain.next[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.next[i,(k+3):(k+n+2)]),j=j.beta+1)
    }
    PO.beta[j.beta+1]=mean(po1.beta)/mean(po2.beta)

    chain.prev=chain.next
  }
  print("Posterior ordinate beta ready!")

  # TAKING LOGARITHM
  if(EXP==FALSE) {LPO.gam=log(PO.gam)}
  else {LPO.gam=0}
  LPO.theta=log(PO.theta); LPO.beta=log(PO.beta)

  # MARGINAL LOG-LIKELIHOOD
  LML=LL.ord+LP.ord-LPO.theta-LPO.gam-sum(LPO.beta)

  list("LL.ord"=LL.ord, "LP.ord"=LP.ord, "LPO.theta"=LPO.theta, "LPO.gam"=LPO.gam, "LPO.beta"=sum(LPO.beta), "LML"=LML)

}

CaseDeletion.RMWIGAUSS=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; logCPO=rep(0,times=n); KL.aux=rep(0,times=n)

  for(s in 1:n)
  {
    aux1=rep(0,times=N); aux2=rep(0,times=N)
    for(ITER in 1:N)
    {
      aux2[ITER]=log.lik.RMWIGAUSS(Time=Time[s],Cens=Cens[s],X=X[s,],beta=as.vector(chain[ITER,1:k]),gam=chain[ITER,k+1],theta=chain[ITER,k+2],EXP)
      aux1[ITER]=exp(-aux2[ITER])
    }
    logCPO[s]=-log(mean(aux1))
    KL.aux[s]=mean(aux2)
  }
  KL=KL.aux-logCPO
  CALIBRATION=0.5*(1+sqrt(1-exp(-2*KL)))
  return(cbind(logCPO,KL,CALIBRATION))
}

################################################################################
# CORRECTION FACTOR FOR OUTLIER DETECTION RMWIGAUSS MODELS
OD.Correction.RMWIGAUSS<-function(chain,Time,X)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  k=dim(X)[2]
  if(k>1)  {beta=apply(chain[,1:k],2,"median")}
  else  	{beta=median(chain[,1])}
  gam=median(chain[,k+1]); theta=median(chain[,k+2])
  aux0=sqrt(2*(exp(-X%*%beta)*Time)^gam+theta^(-2))
  aux=(besselK(x=aux0,nu=1/2,expon.scaled = FALSE)^2)/(besselK(x=aux0,nu=-1/2,expon.scaled = FALSE)*besselK(x=aux0,nu=3/2,expon.scaled = FALSE))
  return(aux)
}

################################################################################
# OUTLIER DETECTION
BF.lambda.obs.RMWIGAUSS<-function(ref,obs,N,thin,Q,burn,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  aux1=Post.lambda.obs.RMWIGAUSS(ref,obs,Time,Cens,X,chain)
  aux2=CFP.obs.RMWIGAUSS(ref,obs,N,thin,Q,burn,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar,EXP)
  return(aux1*aux2)
}

###############################################################################
###############################################################################
######################## ALL CODES FOR RMWLN    MODEL #########################
###############################################################################
###############################################################################

MCMC.RMWLN<-function(N,thin,Q,beta0,gam0=1,theta0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE,FIX.THETA=FALSE)
{
  # NUMBER OF REGRESSORS, SAMPLE SIZE AND NUMBER OF STORED DRAWS
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  # INITIALIZATION OF ELEMENTS WHERE DRAWS ARE STORED
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1)
  # WHERE EMPIRICAL ACCEPTANCE PROBABILITIES ARE STORED
  accept.beta=rep(0,times=k); pbeta.aux=rep(0,times=k); accept.theta=0; ptheta.aux=0; accept.gam=0; pgam.aux=0
  # WHERE ADAPTIVE PROPOSAL VARIANCES ARE STORED (LOG-SCALE)
  ls.beta=matrix(rep(0,times=(N.aux+1)*k),ncol=k); ls.gam=rep(0,times=N.aux+1); ls.theta=rep(0,times=N.aux+1)
  accept.lambda=rep(0,times=n); plambda.aux=rep(0,times=n)

  # INITIALIZATION OF AUXILIARY PARAMETERS
  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.beta.aux=ls.beta[1,]; ls.gam.aux=ls.gam[1]; ls.theta.aux=ls.theta[1]
  ls.lambda.aux=rep(0,times=n)

  # BATCH INDICATOR FOR ADAPTIVE STEPS (WHEN TO ADAPT)
  i_batch=0;

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    # DRAWS FROM BETA
    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=exp(ls.beta.aux)[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1; pbeta.aux[j]=pbeta.aux[j]+1}
    }

    # DRAWS FROM GAMMA (IF NOT FIXED)
    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMWLN.gam(N=1,omega2=exp(ls.gam.aux),gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="LogNormal",lower.bound=0.06,FIX.THETA=FIX.THETA)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1; pgam.aux=pgam.aux+1}
    }

    # DRAWS FROM THETA
    if(FIX.THETA==FALSE)
    {
      MH.theta=GRWMH.RMWLN.theta(N=1,omega2=exp(ls.theta.aux),theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
      theta.aux<-MH.theta$theta
      if(MH.theta$ind==1) {accept.theta=accept.theta+1; ptheta.aux=ptheta.aux+1}
    }

    # DRAWS FROM MIXING PARAMETERS
    if((iter-1)%%Q==0)
    {
      RATE.aux=(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      MH.lambda=GRWMH.RMWLN.lambda(N=1,omega2=exp(ls.lambda.aux),lambda0=lambda.aux,RATE=RATE.aux,theta=theta.aux,Cens=Cens)
      lambda.aux=MH.lambda$lambda
      accept.lambda=accept.lambda+MH.lambda$ind; plambda.aux=plambda.aux+MH.lambda$ind
    }

    # ADAPTIVE STEP
    if(i_batch==50)
    {
      pbeta.aux=pbeta.aux/50; Pbeta.aux=as.numeric(pbeta.aux<rep(ar,times=k))
      ls.beta.aux=ls.beta.aux+((-1)^Pbeta.aux)*min(0.01,1/sqrt(iter))
      pgam.aux=pgam.aux/50; Pgam.aux=as.numeric(pgam.aux<ar)
      ls.gam.aux=ls.gam.aux+((-1)^Pgam.aux)*min(0.01,1/sqrt(iter))
      ptheta.aux=ptheta.aux/50; Ptheta.aux=as.numeric(ptheta.aux<ar)
      ls.theta.aux=ls.theta.aux+((-1)^Ptheta.aux)*min(0.01,1/sqrt(iter))
      plambda.aux=Q*plambda.aux/50; Plambda.aux=as.numeric(plambda.aux<rep(ar,times=n))
      ls.lambda.aux=ls.lambda.aux+((-1)^Plambda.aux)*min(0.01,1/sqrt(iter))

      i_batch=0; pbeta.aux=rep(0,times=k); pgam.aux=0; ptheta.aux=0; plambda.aux=rep(0,times=n)
    }

    # DRAWS STORAGE
    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; theta[iter/thin+1]<-theta.aux; lambda[iter/thin+1,]<-lambda.aux
      ls.beta[iter/thin+1,]=ls.beta.aux; ls.gam[iter/thin+1]=ls.gam.aux; ls.theta[iter/thin+1]=ls.theta.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))
  print(paste("AR gamma :",round(accept.gam/N,2)))
  print(paste("AR theta :",round(accept.theta/N,2)))
  print(paste("Min AR lambda :",round(min(Q*accept.lambda/N),2)))
  print(paste("Max AR lambda :",round(max(Q*accept.lambda/N),2)))
  print(paste("Mean AR lambda :",round(mean(Q*accept.lambda/N),2)))
  print(paste("Median AR lambda :",round(median(Q*accept.lambda/N),2)))

  chain=cbind(beta,gam,theta,lambda,ls.beta,ls.gam,ls.theta)
  return(chain)
}

DIC.RMWLN=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; n=length(Time); L<-rep(0,times=N)
  # LOG-LIKELIHOOD FOR EACH DRAW
  for(iter in 1:N)
  {
    L[iter]=log.lik.RMWLN(Time,Cens,X,beta=as.vector(chain[iter,1:k]),gam=chain[iter,k+1],theta=chain[iter,k+2],EXP)
  }

  aux=apply(chain[,1:(k+2)],2,"median"); pd=-2*mean(L)+2*log.lik.RMWLN(Time,Cens,X,beta=aux[1:k],gam=aux[k+1],theta=aux[k+2],EXP)
  if(EXP==FALSE){pd.aux=k+2}
  else(pd.aux=k+1)
  DIC=-2*mean(L)+pd

  print(pd); print(pd.aux); return(DIC)
}

LML.RMWLN<-function(thin,Q,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE,FIX.THETA=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  # SAMPLE SIZE, NUMBER OF DRAWS AND NUMBER OF REGRESSORS
  n=length(Time); N=dim(chain)[1]; k=dim(X)[2]
  # SETTING PROPOSAL VARIANCES AS THE MEDIAN VALUE (BASED ON ADAPTIVE VERSION)
  if(EXP==FALSE) {omega2.gam=exp(median(chain[,n+2*k+3]))}
  else {omega2.gam=1}
  if(k>1)	{omega2.beta=exp(apply(chain[,(k+n+3):(2*k+n+2)],2,"median"))}
  else		{omega2.beta=exp(median(chain[,n+4]))}
  omega2.theta=exp(median(chain[,2*k+n+4]))

  # NON-ADAPTIVE RUN OF THE CHAIN
  chain.nonadapt=MCMC.RMWLN.NonAdapt(N=0.5*N*thin,thin=thin,Q,beta0=t(chain[N,1:k]),gam0=chain[N,k+1],theta0=chain[N,k+2],Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,omega2.theta,EXP,FIX.THETA)
  chain.nonadapt=chain.nonadapt[-1,]
  gam.star=median(chain.nonadapt[,k+1]); theta.star=median(chain.nonadapt[,k+2])
  if(k>1)	{beta.star=apply(chain.nonadapt[,1:k],2,"median")}
  else		{beta.star=median(chain.nonadapt[,1])}
  N.aux=dim(chain.nonadapt)[1]

  # LIKELIHOOD ORDINATE
  LL.ord=log.lik.RMWLN(Time,Cens,X,beta=beta.star,gam=gam.star,theta=theta.star,EXP)
  print("Likelihood ordinate ready!")

  # PRIOR ORDINATE
  aux=0
  if(FIX.THETA==FALSE){aux=log.prior.theta(theta=theta.star,gam=gam.star,a=hyp.theta,type=typ.theta,mixing="LogNormal")}
  if(EXP==FALSE) {LP.ord=dgamma(gam.star,shape=hyp1.gam,rate=hyp2.gam,log=TRUE)+aux}
  else {LP.ord=aux}
  print("Prior ordinate ready!")

  PO.theta=1
  if(FIX.THETA==FALSE)
  {
    # POSTERIOR ORDINATE - theta
    chain.theta=MCMCR.theta.RMWLN(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.nonadapt[N.aux,1:k]),gam0=chain.nonadapt[N.aux,k+1],theta0=theta.star,lambda0=t(chain.nonadapt[N.aux,(k+3):(k+n+2)]),Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,EXP,FIX.THETA)
    chain.theta=chain.theta[-1,]
    po1.theta=rep(0,times=N.aux); po2.theta=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      po1.theta[i]=alphaRMWLN.theta(theta.0=chain.nonadapt[i,k+2],theta.1=theta.star,gam=chain.nonadapt[i,k+1],lambda=t(chain.nonadapt[i,(k+3):(k+n+2)]),typ.theta,hyp.theta)*dnorm(x=theta.star,mean=chain.nonadapt[i,k+2],sd=sqrt(omega2.theta))
      theta.aux=rnorm(n=1, mean=theta.star,sd=sqrt(omega2.theta))
      if(theta.aux<=0) {po2.theta[i]=0}
      else             {po2.theta[i]=alphaRMWLN.theta(theta.0=theta.star,theta.1=theta.aux,gam=chain.theta[i,k+1],lambda=t(chain.theta[i,(k+3):(k+n+2)]),typ.theta,hyp.theta)}
    }
    PO.theta=mean(po1.theta)/mean(po2.theta)
    print("Posterior ordinate theta ready!")
  }
  if(FIX.THETA==TRUE){chain.theta=chain.nonadapt}

  if(EXP==FALSE)
  {
    # POSTERIOR ORDINATE - gam
    chain.gam=MCMCR.gam.theta.RMWLN(N=0.5*N*thin,thin=thin,Q,beta0=t(chain.theta[N.aux,1:k]),gam0=gam.star,theta0=theta.star,lambda0=t(chain.theta[N.aux,(k+3):(k+n+2)]),Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
    chain.gam=chain.gam[-1,]
    po1.gam=rep(0,times=N.aux); po2.gam=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      po1.gam[i]=alphaRMWLN.gam(gam.0=chain.theta[i,k+1],gam.1=gam.star,Time,Cens,X,beta=as.vector(chain.theta[i,1:k]),theta=theta.star,lambda=t(chain.theta[i,(k+3):(k+n+2)]),typ.theta,hyp.theta,hyp1.gam,hyp2.gam,mixing="LogNormal",lower.bound=0.06,FIX.THETA)*dnorm(x=gam.star,mean=chain.theta[i,k+1],sd=sqrt(omega2.gam))
      gam.aux=rnorm(n=1, mean=gam.star,sd=sqrt(omega2.gam))
      if(gam.aux<=0.06) {po2.gam[i]=0}
      else              {po2.gam[i]=alphaRMWLN.gam(gam.0=gam.star,gam.1=gam.aux,Time,Cens,X,beta=as.vector(chain.gam[i,1:k]),theta=theta.star,lambda=t(chain.gam[i,(k+3):(k+n+2)]),typ.theta,hyp.theta,hyp1.gam,hyp2.gam,mixing="LogNormal",lower.bound=0.06,FIX.THETA)}
    }
    PO.gam=mean(po1.gam)/mean(po2.gam)
    print("Posterior ordinate gamma ready!")
  }

  # POSTERIOR ORDINATE - beta
  if(EXP==FALSE) {chain.prev=chain.gam}
  else {chain.prev=chain.theta}
  PO.beta=rep(0,times=k)
  for(j.beta in 0:(k-1))
  {
    print(j.beta)
    beta0=t(chain.prev[N.aux,1:k]); beta0[j.beta+1]=beta.star[j.beta+1]
    chain.next=MCMCR.betaJ.gam.theta.RMWLN(N=0.5*N*thin,thin=thin,Q,beta0=beta0,gam0=gam.star,theta0=theta.star,lambda0=t(chain.prev[N.aux,(k+3):(k+n+2)]),Time,Cens,X,J=j.beta+1,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
    chain.next=chain.next[-1,]

    po1.beta=rep(0,times=N.aux); po2.beta=rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      beta.0=as.vector(t(chain.prev[i,1:k]))
      beta.1=beta.0; beta.1[j.beta+1]=beta.star[j.beta+1]
      po1.beta[i]=alpha.beta.j(beta.0=beta.0,beta.1=beta.1,gam=chain.prev[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.prev[i,(k+3):(k+n+2)]),j=j.beta+1)*dnorm(x=beta.star[j.beta+1],mean=as.numeric(chain.prev[i,j.beta+1]),sd=sqrt(omega2.beta)[j.beta+1])
      betaj.aux=rnorm(n=1, mean=beta.star[j.beta+1],sd=sqrt(omega2.beta)[j.beta+1])
      beta.2=beta.star; beta.2[j.beta+1]=betaj.aux
      po2.beta[i]=alpha.beta.j(beta.0=beta.star,beta.1=beta.2,gam=chain.next[i,k+1],Time=Time,Cens=Cens,X=X,lambda=t(chain.next[i,(k+3):(k+n+2)]),j=j.beta+1)
    }
    PO.beta[j.beta+1]=mean(po1.beta)/mean(po2.beta)

    chain.prev=chain.next
  }
  print("Posterior ordinate beta ready!")

  # TAKING LOGARITHM
  if(EXP==FALSE) {LPO.gam=log(PO.gam)}
  else {LPO.gam=0}
  LPO.theta=log(PO.theta); LPO.beta=log(PO.beta)

  # MARGINAL LOG-LIKELIHOOD
  LML=LL.ord+LP.ord-LPO.theta-LPO.gam-sum(LPO.beta)

  list("LL.ord"=LL.ord, "LP.ord"=LP.ord, "LPO.theta"=LPO.theta, "LPO.gam"=LPO.gam, "LPO.beta"=sum(LPO.beta), "LML"=LML)

}

CaseDeletion.RMWLN=function(Time,Cens,X,chain,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  N=dim(chain)[1]; k=dim(X)[2]; logCPO=rep(0,times=n); KL.aux=rep(0,times=n)

  for(s in 1:n)
  {
    aux1=rep(0,times=N); aux2=rep(0,times=N)
    for(ITER in 1:N)
    {
      aux2[ITER]=log.lik.RMWLN(Time=Time[s],Cens=Cens[s],X=X[s,],beta=as.vector(chain[ITER,1:k]),gam=chain[ITER,k+1],theta=chain[ITER,k+2],EXP)
      aux1[ITER]=exp(-aux2[ITER])
    }
    logCPO[s]=-log(mean(aux1))
    KL.aux[s]=mean(aux2)
  }
  KL=KL.aux-logCPO
  CALIBRATION=0.5*(1+sqrt(1-exp(-2*KL)))
  return(cbind(logCPO,KL,CALIBRATION))
}

################################################################################
# CORRECTION FACTOR FOR OUTLIER DETECTION RMWLN MODEL
OD.Correction.RMWLN<-function(chain,Time,X)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  n=length(Time); k=dim(X)[2]
  if(k>1)	{beta=apply(chain[,1:k],2,"median")}
  else		{beta=median(chain[,1])}
  gam=median(chain[,k+1]); theta=median(chain[,k+2])
  R0=(exp(-X%*%beta)*Time)^gam
  norm0=rep(0,times=n); norm1=rep(0,times=n)
  exp0=rep(0,times=n); exp1=rep(0,times=n)
  for(i in 1:n)
  {
    norm0[i]=integrate(f.lambda.aux,lower=0,upper=Inf,RATE=R0[i],theta=theta,Cens=0)$value
    norm1[i]=integrate(f.lambda.aux,lower=0,upper=Inf,RATE=R0[i],theta=theta,Cens=1)$value
    exp0[i]=integrate(E.lambda.aux,lower=0,upper=Inf,RATE=R0[i],theta=theta,Cens=0)$value
    exp1[i]=integrate(E.lambda.aux,lower=0,upper=Inf,RATE=R0[i],theta=theta,Cens=1)$value
  }
  aux=(norm1/norm0)*(exp0/exp1)
  return(aux)
}

################################################################################
# OUTLIER DETECTION
BF.lambda.obs.RMWLN<-function(ref,obs,N,thin,Q,burn,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  #######################################################################
  # IMPORTANT: THIS TAKES AS ARGUMENT THE CHAIN AFTER BURN.
  #######################################################################
  chain=as.matrix(chain)
  aux1=Post.lambda.obs.RMWLN(ref,obs,Time,Cens,X,chain)
  aux2=CFP.obs.RMWLN(ref,obs,N,thin,Q,burn,Time,Cens,X,chain,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar,EXP)
  return(aux1*aux2)
}
