################################################################################
################################################################################
# FUNCTIONS SPECIFIC FOR RMWEXP MODEL
################################################################################
################################################################################

# LOG-LIKELIHOOD
log.lik.RMWEXP<-function(Time,Cens,X,beta,gam=1,EXP=FALSE)
{
  RATE=exp(-gam*X%*%beta)
  aux=Cens*(log(RATE)-2*log(RATE*Time^gam+1))+(1-Cens)*(-log(RATE*Time^gam+1))+I(EXP==FALSE)*(Cens*log(gam)+Cens*(gam-1)*log(Time))
  return(sum(aux))
}

# AUXILIARY FUNCTION FOR OUTLIER DETECTION
log.lik.WEI.ref1<-function(Time,Cens,RATE,gam=1,EXP=FALSE)
{
  n=length(Time); aux<-rep(0,n)
  aux=Cens*(log(RATE)-RATE*Time^gam)+(1-Cens)*(-RATE*Time^gam)+I(EXP==FALSE)*(Cens*log(gam)+Cens*(gam-1)*log(Time))
  return(sum(aux))
}

################################################################################
################################################################################
# FUNCTIONS SPECIFIC FOR RMWGAM MODEL
################################################################################
################################################################################

# ACCEPTANCE PROBABILITY THETA
alphaRMWGAM.theta<-function(theta.0,theta.1,gam,lambda,typ.theta,hyp.theta)
{
  l1=n*(theta.1*log(theta.1)-theta.0*log(theta.0))-n*(lgamma(theta.1)-lgamma(theta.0))+(theta.1-theta.0)*sum(log(lambda))-(theta.1-theta.0)*sum(lambda)
  l2=log(prior.theta(theta=theta.1,gam=gam,a=hyp.theta,type=typ.theta,mixing="Gamma")/prior.theta(theta=theta.0,gam=gam,a=hyp.theta,type=typ.theta,mixing="Gamma"))
  aux=min(1,exp(l1+l2))

  return(aux)
}

# LOG-LIKELIHOOD
log.lik.RMWGAM<-function(Time,Cens,X,beta,gam=1,theta,EXP=FALSE)
{
  #n=length(T);
  RATE=exp(-gam*X%*%beta)
  aux=Cens*(log(RATE)-(theta+1)*log((RATE/theta)*Time^gam+1))+(1-Cens)*(-theta*log((RATE/theta)*Time^gam+1))+I(EXP==FALSE)*(Cens*log(gam)+Cens*(gam-1)*log(Time))
  return(sum(aux))
}

# REDUCED CHAIN GIVEN A FIXED VALUE OF LAMBDA[i] (REQUIRED FOR BF.lambda.obs.RMWGAM ONLY)
MCMCR.RMWGAM.lambda.obs<-function(ref,obs,N,thin,Q,beta0,gam0=1,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1); lambda[1,obs]=ref
  accept.beta=rep(0,times=k); pbeta.aux=rep(0,times=k); accept.theta=0; ptheta.aux=0; accept.gam=0; pgam.aux=0
  ls.beta=matrix(rep(0,times=(N.aux+1)*k),ncol=k); ls.gam=rep(0,times=N.aux+1); ls.theta=rep(0,times=N.aux+1)

  i_batch=0;

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.beta.aux=ls.beta[1,]; ls.gam.aux=ls.gam[1]; ls.theta.aux=ls.theta[1]

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=exp(ls.beta.aux)[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1; pbeta.aux[j]=pbeta.aux[j]+1}
    }

    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMW.gam(N=1,omega2=exp(ls.gam.aux),gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="Gamma",lower.bound=0.06)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1; pgam.aux=pgam.aux+1}
    }

    MH.theta=GRWMH.RMWGAM.theta(N=1,omega2=exp(ls.theta.aux),theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
    theta.aux<-MH.theta$theta
    if(MH.theta$ind==1) {accept.theta=accept.theta+1; ptheta.aux=ptheta.aux+1}

    if((iter-1)%%Q==0)
    {
      lambda.aux=rgamma(n=n,shape=theta.aux+Cens,rate=theta.aux+(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux)
      lambda.aux[obs]=ref
    }

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

################################################################################
# MARGINAL POSTERIOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.RMWGAM ONLY)
Post.lambda.obs.RMWGAM<-function(ref,obs,Time,Cens,X,chain)
{
  N<-dim(chain)[1]; k=dim(X)[2]; n<-length(Time); aux1<-rep(0,times=N); aux2<-rep(0,times=N)
  for(iter in 1:N)
  {
    aux1[iter]<-chain[iter,k+2]+(exp(-as.numeric(X[obs,]%*%as.vector(chain[iter,1:k])))*Time[obs])^(chain[iter,k+1])
    aux2[iter]<-dgamma(x=ref,shape=chain[iter,k+2]+Cens[obs],rate=aux1[iter])
  }
  aux=mean(aux2)
  return(aux)
}

################################################################################
# CORRECTION FACTOR/PRIOR FOR BAYES FACTOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.RMWGAM ONLY)
CFP.obs.RMWGAM<-function(ref,obs,N,thin,Q,burn,Time,Cens,X,chain0,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  n=length(Time); k=dim(X)[2]
  beta0=apply(chain0[,1:k],2,"median"); gam0=median(chain0[,k+1]); theta0=median(chain0[,k+2]); lambda0=apply(chain0[,(k+3):(k+2+n)],2,"median")

  chain=MCMCR.RMWGAM.lambda.obs(ref,obs,N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar,EXP)
  chain=chain[-(1:burn),]; N.aux=dim(chain)[1]; aux1=rep(0,times=N.aux)

  for(iter in 1:N.aux)
  {
    aux1[iter]=1/dgamma(x=ref,shape=chain[iter,k+2],rate=chain[iter,k+2])
  }

  aux=mean(aux1)
  return(aux)
}



################################################################################
################################################################################
# FUNCTIONS SPECIFIC FOR RMWIGAM MODEL
################################################################################
################################################################################



# ACCEPTANCE PROBABILITY THETA
alphaRMWIGAM.theta<-function(theta.0,theta.1,gam,lambda,typ.theta,hyp.theta)
{
  l1=n*(lgamma(theta.0)-lgamma(theta.1))+(theta.0-theta.1)*sum(log(lambda))
  l2=log.prior.theta(theta=theta.1,gam=gam,a=hyp.theta,type=typ.theta,mixing="InvGamma")-log.prior.theta(theta=theta.0,gam=gam,a=hyp.theta,type=typ.theta,mixing="InvGamma")
  aux=min(1,exp(l1+l2))

  return(aux)
}

# LOG-LIKELIHOOD
log.f.RME.IGAM<-function(Time,alpha,theta)
{
  aux1=log(2)+((theta+1)/2)*log(alpha)+((theta-1)/2)*log(Time)-lgamma(theta)
  aux2=log(besselK(x=2*sqrt(alpha*Time), nu=-(theta-1), expon.scaled=FALSE))
  aux=aux1+aux2
  return(aux)
}

log.S.RME.IGAM<-function(Time,alpha,theta)
{
  aux1=log(2)+(theta/2)*log(alpha*Time)-lgamma(theta)
  aux2=log(besselK(x=2*sqrt(alpha*Time), nu=-theta, expon.scaled=FALSE))
  aux=aux1+aux2
  return(aux)
}

log.lik.RMWIGAM<-function(Time,Cens,X,beta,gam=1,theta,EXP=FALSE)
{
  RATE=exp(-gam*X%*%beta)
  aux=Cens*log.f.RME.IGAM(Time^gam,RATE,theta)+(1-Cens)*log.S.RME.IGAM(Time^gam,RATE,theta)+I(EXP==FALSE)*(Cens*log(gam)+Cens*(gam-1)*log(Time))
  return(sum(aux))
}



# AUXILIARY FUNCTION FOR OUTLIER DETECTION
# REDUCED CHAIN GIVEN A FIXED VALUE OF LAMBDA[i] (REQUIRED FOR BF.lambda.obs.RMWIGAM ONLY)
MCMCR.RMWIGAM.lambda.obs<-function(ref,obs,N,thin,Q,beta0,gam0=1,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1); lambda[1,obs]=ref
  accept.beta=rep(0,times=k); pbeta.aux=rep(0,times=k); accept.theta=0; ptheta.aux=0; accept.gam=0; pgam.aux=0
  ls.beta=matrix(rep(0,times=(N.aux+1)*k),ncol=k); ls.gam=rep(0,times=N.aux+1); ls.theta=rep(0,times=N.aux+1)

  i_batch=0;

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.beta.aux=ls.beta[1,]; ls.gam.aux=ls.gam[1]; ls.theta.aux=ls.theta[1]

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=exp(ls.beta.aux)[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1; pbeta.aux[j]=pbeta.aux[j]+1}
    }

    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMW.gam(N=1,omega2=exp(ls.gam.aux),gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="InvGamma",lower.bound=0.06)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1; pgam.aux=pgam.aux+1}
    }

    MH.theta=GRWMH.RMWIGAM.theta(N=1,omega2=exp(ls.theta.aux),theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
    theta.aux<-MH.theta$theta
    if(MH.theta$ind==1) {accept.theta=accept.theta+1; ptheta.aux=ptheta.aux+1}

    if((iter-1)%%Q==0)
    {
      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      for(ind in 1:n){lambda.aux[ind]=rgig(n=1,lambda=-theta.aux+Cens[ind],chi=2,psi=psi.aux[ind])}
      lambda.aux[obs]=ref

      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      for(ind in 1:n)
      {
        lambda.aux[ind]=rgig(n=1,lambda=-theta.aux+Cens[ind],chi=2,psi=psi.aux[ind])
      }
    }

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

################################################################################
# MARGINAL POSTERIOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.RMWIGAM ONLY)
Post.lambda.obs.RMWIGAM<-function(ref,obs,Time,Cens,X,chain)
{
  N<-dim(chain)[1]; k=dim(X)[2]; n<-length(Time); aux1<-rep(0,times=N); aux2<-rep(0,times=N)
  for(iter in 1:N)
  {
    aux1[iter]<-2*(exp(-as.numeric(X[obs,]%*%as.vector(chain[iter,1:k])))*Time[obs])^(chain[iter,k+1])
    aux2[iter]<-dgig(x=ref,lambda=-chain[iter,k+2]+Cens[obs],chi=2,psi=aux1[iter])
  }
  aux=mean(aux2)
  return(aux)
}

################################################################################
# CORRECTION FACTOR/PRIOR FOR BAYES FACTOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.RMWIGAM ONLY)
CFP.obs.RMWIGAM<-function(ref,obs,N,thin,Q,burn,Time,Cens,X,chain0,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  n=length(Time); k=dim(X)[2]
  beta0=apply(chain0[,1:k],2,"median"); gam0=median(chain0[,k+1]); theta0=median(chain0[,k+2]); lambda0=apply(chain0[,(k+3):(k+2+n)],2,"median")

  chain=MCMCR.RMWIGAM.lambda.obs(ref,obs,N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar,EXP)
  chain=chain[-(1:burn),]; N.aux=dim(chain)[1]; aux1=rep(0,times=N.aux)

  for(iter in 1:N.aux)
  {
    aux1[iter]=(ref^2)/dgamma(x=1/ref,shape=chain[iter,k+2],rate=1)
  }

  aux=mean(aux1)
  return(aux)
}

################################################################################
################################################################################
# FUNCTIONS SPECIFIC FOR RMWIGAUSS MODEL
################################################################################
################################################################################

# GAUSSIAN RANDOM WALK METROPOLIS STEP FOR THETA
GRWMH.RMWIGAUSS.theta<-function(N=1,omega2,theta0,gam,lambda,type.prior,hyper)
{
  n<-length(lambda); theta<-rep(0,times=N+1); theta[1]<-theta0; ind<-rep(0,N+1)

  for(l in 1:N)
  {
    y<-rnorm(n=1,mean=theta[l],sd=sqrt(omega2))
    if(y<=0) {theta[l+1]<-theta[l]; ind[l]<-2}
    else
    {
      u.aux<-runif(1,min=0,max=1)

      aux1=0.5*(theta[l]^(-2)-y^(-2))*sum(lambda)-n*(theta[l]^(-1)-y^(-1))
      aux2=log.prior.theta(theta=y,gam=gam,a=hyper,type=type.prior,mixing="InvGauss")-log.prior.theta(theta=theta[l],gam=gam,a=hyper,type=type.prior,mixing="InvGauss")
      log.aux<-aux1+aux2
      if(log(u.aux)<log.aux) {theta[l+1]<-y; ind[l+1]<-1}
      else {theta[l+1]<-theta[l]; ind[l+1]<-0}
    }
  }
  theta<-theta[-1]; ind=ind[-1]
  list("theta"=theta,"ind"=ind)
}

# ACCEPTANCE PROBABILITY THETA
alphaRMWIGAUSS.theta<-function(theta.0,theta.1,gam,lambda,typ.theta,hyp.theta)
{
  l1=0.5*(theta.0^(-2)-theta.1^(-2))*sum(lambda)-n*(theta.0^(-1)-theta.1^(-1))
  l2=log.prior.theta(theta=theta.1,gam=gam,a=hyp.theta,type=typ.theta,mixing="InvGauss")-log.prior.theta(theta=theta.0,gam=gam,a=hyp.theta,type=typ.theta,mixing="InvGauss")
  aux=min(1,exp(l1+l2))

  return(aux)
}

# LOG-LIKELIHOOD
log.lik.RMWIGAUSS<-function(Time,Cens,X,beta,gam=1,theta,EXP=FALSE)
{
  RATE=exp(-gam*X%*%beta); AUX=2*RATE*(Time^gam)+theta^(-2)
  aux=Cens*(log(RATE)+1/theta-0.5*log(AUX)-sqrt(AUX))+(1-Cens)*(1/theta-sqrt(AUX))+I(EXP==FALSE)*(Cens*log(gam)+Cens*(gam-1)*log(Time))
  return(sum(aux))
}

# NON-ADAPTIVE MCMC (MARGINAL LIKELIHOOD CALCULATION)
MCMC.RMWIGAUSS.NonAdapt<-function(N,thin,Q,beta0,gam0,theta0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,omega2.beta,omega2.gam=1,omega2.theta,EXP=FALSE)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1)
  accept.beta=rep(0,times=k); accept.theta=0; accept.gam=0

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]

  for(iter in 2:(N+1))
  {
    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=omega2.beta[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1}
    }

    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMW.gam(N=1,omega2=omega2.gam,gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="InvGauss",lower.bound=0.06)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1}
    }

    MH.theta=GRWMH.RMWIGAUSS.theta(N=1,omega2=omega2.theta,theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
    theta.aux<-MH.theta$theta
    if(MH.theta$ind==1) {accept.theta=accept.theta+1}

    if((iter-1)%%Q==0)
    {
      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux+theta.aux^(-2)
      for(ind in 1:n)
      {
        lambda.aux[ind]=rgig(n=1,lambda=Cens[ind]-1/2,chi=1,psi=psi.aux[ind])
      }
    }

    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; theta[iter/thin+1]<-theta.aux; lambda[iter/thin+1,]<-lambda.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))
  print(paste("AR gamma :",round(accept.gam/N,2)))
  print(paste("AR theta :",round(accept.theta/N,2)))

  chain=cbind(beta,gam,theta,lambda)
  return(chain)
}


# REDUCED MCMC RUN WITH FIXED THETA (MARGINAL LIKELIHOOD CALCULATION)
MCMCR.theta.RMWIGAUSS<-function(N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,omega2.beta,omega2.gam=1,EXP=FALSE)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  theta<-rep(theta0,times=N.aux+1)
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-lambda0
  accept.beta=rep(0,times=k); accept.gam=0

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]

  for(iter in 2:(N+1))
  {
    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=omega2.beta[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1}
    }

    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMW.gam(N=1,omega2=omega2.gam,gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="InvGauss",lower.bound=0.06)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1}
    }

    if((iter-1)%%Q==0)
    {
      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux+theta.aux^(-2)
      for(ind in 1:n)
      {
        lambda.aux[ind]=rgig(n=1,lambda=Cens[ind]-1/2,chi=1,psi=psi.aux[ind])
      }
    }

    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; lambda[iter/thin+1,]<-lambda.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))
  print(paste("AR gamma :",round(accept.gam/N,2)))

  chain=cbind(beta,gam,theta,lambda)
  return(chain)
}


# REDUCED MCMC RUN WITH FIXED THETA AND GAMMA (MARGINAL LIKELIHOOD CALCULATION)
MCMCR.gam.theta.RMWIGAUSS<-function(N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  gam<-rep(gam0,times=N.aux+1)
  theta<-rep(theta0,times=N.aux+1)
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-lambda0
  accept.beta=rep(0,times=k)

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]

  for(iter in 2:(N+1))
  {
    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=omega2.beta[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1}
    }

    if((iter-1)%%Q==0)
    {
      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux+theta.aux^(-2)
      for(ind in 1:n)
      {
        lambda.aux[ind]=rgig(n=1,lambda=Cens[ind]-1/2,chi=1,psi=psi.aux[ind])
      }
    }

    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; lambda[iter/thin+1,]<-lambda.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))

  chain=cbind(beta,gam,theta,lambda)
  return(chain)
}

# REDUCED RUN MCMC WITH FIXED THETA, GAMMA AND beta[1],...,beta[J] (MARGINAL LIKELIHOOD CALCULATION)
MCMCR.betaJ.gam.theta.RMWIGAUSS<-function(N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,J,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  gam<-rep(gam0,times=N.aux+1)
  theta<-rep(theta0,times=N.aux+1)
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-lambda0
  accept.beta=rep(0,times=k)

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]

  for(iter in 2:(N+1))
  {
    if(J<k)
    {
      for(ind.b in (J+1):k)
      {
        MH.beta=GRWMH.RMW.beta.j(N=1,omega2=omega2.beta[ind.b],j=ind.b,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
        beta.aux<-MH.beta$beta
        if(MH.beta$ind==1) {accept.beta[ind.b]=accept.beta[ind.b]+1}
      }
    }

    if((iter-1)%%Q==0)
    {
      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux+theta.aux^(-2)
      for(ind in 1:n)
      {
        lambda.aux[ind]=rgig(n=1,lambda=Cens[ind]-1/2,chi=1,psi=psi.aux[ind])
      }
    }

    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; lambda[iter/thin+1,]<-lambda.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))

  chain=cbind(beta,gam,theta,lambda)
  return(chain)
}

# AUXILIARY FUNCTION FOR OUTLIER DETECTION
# REDUCED CHAIN GIVEN A FIXED VALUE OF LAMBDA[i] (REQUIRED FOR BF.lambda.obs.RMWIGAM ONLY)
MCMCR.RMWIGAUSS.lambda.obs<-function(ref,obs,N,thin,Q,beta0,gam0=1,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1); lambda[1,obs]=ref
  accept.beta=rep(0,times=k); pbeta.aux=rep(0,times=k); accept.theta=0; ptheta.aux=0; accept.gam=0; pgam.aux=0
  ls.beta=matrix(rep(0,times=(N.aux+1)*k),ncol=k); ls.gam=rep(0,times=N.aux+1); ls.theta=rep(0,times=N.aux+1)

  i_batch=0;

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.beta.aux=ls.beta[1,]; ls.gam.aux=ls.gam[1]; ls.theta.aux=ls.theta[1]

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=exp(ls.beta.aux)[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1; pbeta.aux[j]=pbeta.aux[j]+1}
    }

    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMW.gam(N=1,omega2=exp(ls.gam.aux),gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="InvGauss",lower.bound=0.06)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1; pgam.aux=pgam.aux+1}
    }

    MH.theta=GRWMH.RMWIGAM.theta(N=1,omega2=exp(ls.theta.aux),theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
    theta.aux<-MH.theta$theta
    if(MH.theta$ind==1) {accept.theta=accept.theta+1; ptheta.aux=ptheta.aux+1}

    if((iter-1)%%Q==0)
    {
      psi.aux=2*(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux+theta.aux^(-2)
      for(ind in 1:n){lambda.aux[ind]=rgig(n=1,lambda=Cens[ind]-1/2,chi=1,psi=psi.aux[ind])}
      lambda.aux[obs]=ref
    }

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

################################################################################
# MARGINAL POSTERIOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.RMWIGAM ONLY)
Post.lambda.obs.RMWIGAUSS<-function(ref,obs,Time,Cens,X,chain)
{
  N<-dim(chain)[1]; k=dim(X)[2]; n<-length(Time); aux1<-rep(0,times=N); aux2<-rep(0,times=N)
  for(iter in 1:N)
  {
    aux1[iter]<-2*(exp(-as.numeric(X[obs,]%*%as.vector(chain[iter,1:k])))*Time[obs])^(chain[iter,k+1])
    aux2[iter]<-dgig(x=ref,lambda=Cens[obs]-1/2,chi=1,psi=aux1[iter])
  }
  aux=mean(aux2)
  return(aux)
}

################################################################################
# CORRECTION FACTOR/PRIOR FOR BAYES FACTOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.RMWIGAM ONLY)
CFP.obs.RMWIGAUSS<-function(ref,obs,N,thin,Q,burn,Time,Cens,X,chain0,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  n=length(Time); k=dim(X)[2]
  beta0=apply(chain0[,1:k],2,"median"); gam0=median(chain0[,k+1]); theta0=median(chain0[,k+2]); lambda0=apply(chain0[,(k+3):(k+2+n)],2,"median")

  chain=MCMCR.RMWIGAUSS.lambda.obs(ref,obs,N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar,EXP)
  chain=chain[-(1:burn),]; N.aux=dim(chain)[1]; aux1=rep(0,times=N.aux)
  for(iter in 1:N.aux)
  {
    aux1[iter]=1/dinvgauss(x=ref,mean=chain[iter,k+2],dispersion=1)
  }

  aux=mean(aux1)
  return(aux)
}

################################################################################
################################################################################
# FUNCTIONS SPECIFIC FOR RMWLN MODEL
################################################################################
################################################################################

# GAUSSIAN RANDOM WALK METROPOLIS STEP FOR THETA
GRWMH.RMWLN.theta<-function(N=1,omega2,theta0,gam,lambda,type.prior,hyper)
{
  n<-length(lambda); theta<-rep(0,times=N+1); theta[1]<-theta0; ind<-rep(0,N+1)
  for(l in 1:N)
  {
    y<-rnorm(n=1,mean=theta[l],sd=sqrt(omega2))
    if(y<=0) {theta[l+1]<-theta[l]; ind[l+1]<-2}
    else
    {
      u.aux<-runif(1,min=0,max=1)
      aux1=(n/2)*log(theta[l]/y)-0.5*(y^(-1)-theta[l]^(-1))*sum((log(lambda))^2)
      aux2=log(prior.theta(theta=y,gam=gam,a=hyper,type=type.prior,mixing="LogNormal")/prior.theta(theta=theta[l],gam=gam,a=hyper,type=type.prior,mixing="LogNormal"))
      log.aux<-aux1+aux2
      if(log(u.aux)<log.aux) {theta[l+1]<-y; ind[l+1]<-1}
      else {theta[l+1]<-theta[l]; ind[l+1]<-0}
    }
  }
  theta<-theta[-1]; ind=ind[-1]
  list("theta"=theta,"ind"=ind)
}

# ACCEPTANCE PROBABILITY THETA
alphaRMWLN.theta<-function(theta.0,theta.1,gam,lambda,typ.theta,hyp.theta)
{
  l1=(n/2)*log(theta.0/theta.1)-0.5*(theta.1^(-1)-theta.0^(-1))*sum((log(lambda))^2)
  l2=log(prior.theta(theta=theta.1,gam=gam,a=hyp.theta,type=typ.theta,mixing="LogNormal")/prior.theta(theta=theta.0,gam=gam,a=hyp.theta,type=typ.theta,mixing="LogNormal"))
  aux=min(1,exp(l1+l2))
  return(aux)
}

# GAUSSIAN RANDOM WALK METROPOLIS STEP FOR LAMBDA
GRWMH.RMWLN.lambda<-function(N=1,omega2,lambda0,RATE,theta,Cens)
{
  n<-length(lambda0); lambda<-matrix(rep(0,times=n*(N+1)),ncol=n); lambda[1,]<-lambda0; ind<-matrix(rep(0,times=n*(N+1)),ncol=n)

  for(l in 1:N)
  {
    y<-rnorm(n=n,mean=lambda[l,],sd=sqrt(omega2))
    ind.aux=I(y>=0); y=abs(y)
    log.aux=(Cens-1)*log(y/lambda[l,])-RATE*(y-lambda[l,])-(1/(2*theta))*((log(y))^2-(log(lambda[l,]))^2)
    u.aux<-runif(n=n,min=0,max=1)
    aux=I(log(u.aux)<log.aux)
    aux=as.numeric(aux)*as.numeric(ind.aux)
    lambda[l+1,]=aux*y+(1-aux)*lambda[l,]; ind[l+1,]=as.numeric(aux)
  }
  lambda<-lambda[-1,]; ind=ind[-1,]
  list("lambda"=lambda,"ind"=ind)
}


# LOG-LIKELIHOOD
f.joint.RME.LN=function(lambda,Time,alpha,theta)
{
  aux1=exp(-alpha*lambda*Time)
  aux2=exp(-(1/(2*theta))*(log(lambda))^2)
  aux=aux1*aux2

  return(aux)
}

f.RME.LN=function(Time,alpha,theta)
{
  aux1=alpha*(2*pi*theta)^(-1/2)
  aux2=integrate(f.joint.RME.LN,lower=0,upper=Inf,Time=Time,alpha=alpha,theta=theta)$value
  aux=aux1*aux2

  return(aux)
}

S.joint.RME.LN=function(lambda,Time,alpha,theta)
{
  aux1=(lambda^(-1))*exp(-alpha*lambda*Time)
  aux2=exp(-(1/(2*theta))*(log(lambda))^2)
  aux=aux1*aux2

  return(aux)
}

S.RME.LN=function(Time,alpha,theta)
{
  aux1=(2*pi*theta)^(-1/2)
  aux2=integrate(S.joint.RME.LN,lower=0.00001,upper=Inf,Time=Time,alpha=alpha,theta=theta)$value
  aux=aux1*aux2

  return(aux)
}

# LOG-LIKELIHOOD FUNCTION
log.lik.RMW.LN.aux<-function(Time,Cens,X,beta,gam=1,theta)
{
  n=length(Time); aux<-rep(0,n); RATE=exp(-gam*X%*%beta)
  for(i in 1:n)
  {
    if(Cens[i]==1) {aux[i]=log(f.RME.LN(Time=Time[i],alpha=RATE[i],theta=theta))}
    if(Cens[i]==0) {aux[i]=log(S.RME.LN(Time=Time[i],alpha=RATE[i],theta=theta))}
  }
  return(aux)
}

log.lik.RMWLN<-function(Time,Cens,X,beta,gam=1,theta,EXP=FALSE)
{
  aux=sum(log.lik.RMW.LN.aux(Time^gam,Cens,X,beta,gam,theta)+I(EXP==FALSE)*(Cens*log(gam)+Cens*(gam-1)*log(Time)))
  return(aux)
}

# NON-ADAPTIVE MCMC (MARGINAL LIKELIHOOD CALCULATION)
MCMC.RMWLN.NonAdapt<-function(N,thin,Q,beta0,gam0,theta0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta,omega2.gam,omega2.theta,EXP=FALSE,FIX.THETA=FALSE)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1)
  accept.beta=rep(0,times=k); accept.theta=0; accept.gam=0
  accept.lambda=rep(0,times=n); plambda.aux=rep(0,times=n)

  i_batch=0;

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.lambda.aux=rep(0,times=n)

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=omega2.beta[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1}
    }

    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMWLN.gam(N=1,omega2=omega2.gam,gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="LogNormal",lower.bound=0.06,FIX.THETA=FIX.THETA)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1}
    }

    if(FIX.THETA==FALSE)
    {
      MH.theta=GRWMH.RMWLN.theta(N=1,omega2=omega2.theta,theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
      theta.aux<-MH.theta$theta
      if(MH.theta$ind==1) {accept.theta=accept.theta+1}
    }

    if((iter-1)%%Q==0)
    {
      RATE.aux=(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      MH.lambda=GRWMH.RMWLN.lambda(N=1,omega2=exp(ls.lambda.aux),lambda0=lambda.aux,RATE=RATE.aux,theta=theta.aux,Cens=Cens)
      lambda.aux=MH.lambda$lambda
      accept.lambda=accept.lambda+MH.lambda$ind; plambda.aux=plambda.aux+MH.lambda$ind
    }

    if(i_batch==50)
    {
      plambda.aux=Q*plambda.aux/50; Plambda.aux=as.numeric(plambda.aux<rep(ar,times=n))
      ls.lambda.aux=ls.lambda.aux+((-1)^Plambda.aux)*min(0.01,1/sqrt(iter))

      i_batch=0; plambda.aux=rep(0,times=n)
    }

    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; theta[iter/thin+1]<-theta.aux; lambda[iter/thin+1,]<-lambda.aux
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

  chain=cbind(beta,gam,theta,lambda)
  return(chain)
}


# REDUCED MCMC RUN WITH FIXED THETA (MARGINAL LIKELIHOOD CALCULATION)
MCMCR.theta.RMWLN<-function(N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,omega2.beta,omega2.gam=1,EXP=FALSE,FIX.THETA=FALSE)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1)
  accept.beta=rep(0,times=k); accept.theta=0; accept.gam=0
  accept.lambda=rep(0,times=n); plambda.aux=rep(0,times=n)

  i_batch=0;

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.lambda.aux=rep(0,times=n)

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=omega2.beta[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1}
    }

    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMWLN.gam(N=1,omega2=omega2.gam,gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="LogNormal",lower.bound=0.06,FIX.THETA=FIX.THETA)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1}
    }

    if((iter-1)%%Q==0)
    {
      RATE.aux=(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      MH.lambda=GRWMH.RMWLN.lambda(N=1,omega2=exp(ls.lambda.aux),lambda0=lambda.aux,RATE=RATE.aux,theta=theta.aux,Cens=Cens)
      lambda.aux=MH.lambda$lambda
      accept.lambda=accept.lambda+MH.lambda$ind; plambda.aux=plambda.aux+MH.lambda$ind
    }

    if(i_batch==50)
    {
      plambda.aux=Q*plambda.aux/50; Plambda.aux=as.numeric(plambda.aux<rep(ar,times=n))
      ls.lambda.aux=ls.lambda.aux+((-1)^Plambda.aux)*min(0.01,1/sqrt(iter))

      i_batch=0; plambda.aux=rep(0,times=n)
    }

    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; theta[iter/thin+1]<-theta.aux; lambda[iter/thin+1,]<-lambda.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))
  print(paste("AR gamma :",round(accept.gam/N,2)))
  print(paste("Min AR lambda :",round(min(Q*accept.lambda/N),2)))
  print(paste("Max AR lambda :",round(max(Q*accept.lambda/N),2)))
  print(paste("Mean AR lambda :",round(mean(Q*accept.lambda/N),2)))
  print(paste("Median AR lambda :",round(median(Q*accept.lambda/N),2)))

  chain=cbind(beta,gam,theta,lambda)
  return(chain)
}

# REDUCED MCMC RUN WITH FIXED THETA AND GAMMA (MARGINAL LIKELIHOOD CALCULATION)
MCMCR.gam.theta.RMWLN<-function(N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1)
  accept.beta=rep(0,times=k); accept.theta=0; accept.gam=0
  accept.lambda=rep(0,times=n); plambda.aux=rep(0,times=n)

  i_batch=0;

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.lambda.aux=rep(0,times=n)

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=omega2.beta[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1}
    }

    if((iter-1)%%Q==0)
    {
      RATE.aux=(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      MH.lambda=GRWMH.RMWLN.lambda(N=1,omega2=exp(ls.lambda.aux),lambda0=lambda.aux,RATE=RATE.aux,theta=theta.aux,Cens=Cens)
      lambda.aux=MH.lambda$lambda
      accept.lambda=accept.lambda+MH.lambda$ind; plambda.aux=plambda.aux+MH.lambda$ind
    }

    if(i_batch==50)
    {
      plambda.aux=Q*plambda.aux/50; Plambda.aux=as.numeric(plambda.aux<rep(ar,times=n))
      ls.lambda.aux=ls.lambda.aux+((-1)^Plambda.aux)*min(0.01,1/sqrt(iter))

      i_batch=0; plambda.aux=rep(0,times=n)
    }

    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; theta[iter/thin+1]<-theta.aux; lambda[iter/thin+1,]<-lambda.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))
  print(paste("AR gamma :",round(accept.gam/N,2)))
  print(paste("Min AR lambda :",round(min(Q*accept.lambda/N),2)))
  print(paste("Max AR lambda :",round(max(Q*accept.lambda/N),2)))
  print(paste("Mean AR lambda :",round(mean(Q*accept.lambda/N),2)))
  print(paste("Median AR lambda :",round(median(Q*accept.lambda/N),2)))

  chain=cbind(beta,gam,theta,lambda)
  return(chain)
}

# REDUCED RUN MCMC WITH FIXED THETA, GAMMA AND beta[1],...,beta[J] (MARGINAL LIKELIHOOD CALCULATION)
MCMCR.betaJ.gam.theta.RMWLN<-function(N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,J,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar=0.44,omega2.beta)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1)
  accept.beta=rep(0,times=k); accept.theta=0; accept.gam=0
  accept.lambda=rep(0,times=n); plambda.aux=rep(0,times=n)

  i_batch=0;

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.lambda.aux=rep(0,times=n)

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    if(J<k)
    {
      for(ind.b in (J+1):k)
      {
        MH.beta=GRWMH.RMW.beta.j(N=1,omega2=omega2.beta[ind.b],j=ind.b,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
        beta.aux<-MH.beta$beta
        if(MH.beta$ind==1) {accept.beta[ind.b]=accept.beta[ind.b]+1}
      }
    }

    if((iter-1)%%Q==0)
    {
      RATE.aux=(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      MH.lambda=GRWMH.RMWLN.lambda(N=1,omega2=exp(ls.lambda.aux),lambda0=lambda.aux,RATE=RATE.aux,theta=theta.aux,Cens=Cens)
      lambda.aux=MH.lambda$lambda
      accept.lambda=accept.lambda+MH.lambda$ind; plambda.aux=plambda.aux+MH.lambda$ind
    }

    if(i_batch==50)
    {
      plambda.aux=Q*plambda.aux/50; Plambda.aux=as.numeric(plambda.aux<rep(ar,times=n))
      ls.lambda.aux=ls.lambda.aux+((-1)^Plambda.aux)*min(0.01,1/sqrt(iter))

      i_batch=0; plambda.aux=rep(0,times=n)
    }

    if(iter%%thin==0)
    {
      beta[iter/thin+1,]<-beta.aux; gam[iter/thin+1]<-gam.aux; theta[iter/thin+1]<-theta.aux; lambda[iter/thin+1,]<-lambda.aux
    }
    if((iter-1)%%100000==0) {print(iter-1)}

  }

  print(paste("AR beta",1:k,":",round(accept.beta/N,2)))
  print(paste("AR gamma :",round(accept.gam/N,2)))
  print(paste("Min AR lambda :",round(min(Q*accept.lambda/N),2)))
  print(paste("Max AR lambda :",round(max(Q*accept.lambda/N),2)))
  print(paste("Mean AR lambda :",round(mean(Q*accept.lambda/N),2)))
  print(paste("Median AR lambda :",round(median(Q*accept.lambda/N),2)))

  chain=cbind(beta,gam,theta,lambda)
  return(chain)
}

# AUXILIARY FUNCTION FOR OUTLIER DETECTION
f.lambda.aux<-function(lambda,RATE,theta,Cens)
{
  aux=(lambda^(Cens-1))*exp(-RATE*lambda)*exp(-((log(lambda))^2)/(2*theta))
  return(aux)
}

E.lambda.aux<-function(lambda,RATE,theta,Cens)
{
  aux=(lambda^Cens)*exp(-RATE*lambda)*exp(-((log(lambda))^2)/(2*theta))
  return(aux)
}

# REDUCED CHAIN GIVEN A FIXED VALUE OF LAMBDA[i] (REQUIRED FOR BF.lambda.obs.RMWLN ONLY)
MCMCR.RMWLN.lambda.obs<-function(ref,obs,N,thin,Q,beta0,gam0=1,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  k<-length(beta0); n<-length(Time); N.aux<-N/thin
  beta<-matrix(rep(0,times=(N.aux+1)*k),ncol=k); beta[1,]<-beta0
  theta<-rep(0,times=N.aux+1); theta[1]<-theta0
  gam<-rep(0,times=N.aux+1); gam[1]<-gam0
  lambda<-matrix(rep(0,times=(N.aux+1)*n),ncol=n); lambda[1,]<-rgamma(n,shape=theta0,rate=1); lambda[1,obs]=ref
  accept.beta=rep(0,times=k); pbeta.aux=rep(0,times=k); accept.theta=0; ptheta.aux=0; accept.gam=0; pgam.aux=0
  accept.lambda=rep(0,times=n); plambda.aux=rep(0,times=n)
  ls.beta=matrix(rep(0,times=(N.aux+1)*k),ncol=k); ls.gam=rep(0,times=N.aux+1); ls.theta=rep(0,times=N.aux+1)

  i_batch=0;

  beta.aux=beta[1,]; gam.aux=gam[1]; theta.aux=theta[1]; lambda.aux=lambda[1,]
  ls.beta.aux=ls.beta[1,]; ls.gam.aux=ls.gam[1]; ls.theta.aux=ls.theta[1]
  ls.lambda.aux=rep(0,times=n)

  for(iter in 2:(N+1))
  {
    i_batch=i_batch+1;

    for(j in 1:k)
    {
      MH.beta=GRWMH.RMW.beta.j(N=1,omega2=exp(ls.beta.aux)[j],j=j,beta0=beta.aux,Time=Time,Cens=Cens,X=X,gam=gam.aux,lambda=lambda.aux)
      beta.aux<-MH.beta$beta
      if(MH.beta$ind==1) {accept.beta[j]=accept.beta[j]+1; pbeta.aux[j]=pbeta.aux[j]+1}
    }

    if(EXP==FALSE)
    {
      MH.gam=GRWMH.RMWLN.gam(N=1,omega2=exp(ls.gam.aux),gam0=gam.aux,Time=Time,Cens=Cens,X=X,beta=beta.aux,theta=theta.aux,lambda=lambda.aux,typ.theta=typ.theta,hyp.theta=hyp.theta,hyp1.gam=hyp1.gam,hyp2.gam=hyp2.gam,mixing="InvGauss",lower.bound=0.06)
      gam.aux<-MH.gam$gam
      if(MH.gam$ind==1) {accept.gam=accept.gam+1; pgam.aux=pgam.aux+1}
    }

    MH.theta=GRWMH.RMWIGAM.theta(N=1,omega2=exp(ls.theta.aux),theta0=theta.aux,gam=gam.aux,lambda=lambda.aux,type.prior=typ.theta,hyper=hyp.theta)
    theta.aux<-MH.theta$theta
    if(MH.theta$ind==1) {accept.theta=accept.theta+1; ptheta.aux=ptheta.aux+1}

    if((iter-1)%%Q==0)
    {
      RATE.aux=(exp(-as.numeric(X%*%(beta.aux)))*Time)^gam.aux
      MH.lambda=GRWMH.RMWLN.lambda(N=1,omega2=exp(ls.lambda.aux),lambda0=lambda.aux,RATE=RATE.aux,theta=theta.aux,Cens=Cens)
      lambda.aux=MH.lambda$lambda
      accept.lambda=accept.lambda+MH.lambda$ind; plambda.aux=plambda.aux+MH.lambda$ind
      lambda.aux[obs]=ref
    }

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

################################################################################
# MARGINAL POSTERIOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.RMWLN ONLY)
Post.lambda.obs.RMWLN<-function(ref,obs,Time,Cens,X,chain)
{
  N<-dim(chain)[1]; k=dim(X)[2]; n<-length(Time); aux1<-rep(0,times=N); aux2<-rep(0,times=N); aux3<-rep(0,times=N)
  for(iter in 1:N)
  {
    aux1[iter]<-(exp(-as.numeric(X[obs,]%*%as.vector(chain[iter,1:k])))*Time[obs])^(chain[iter,k+1])
    aux2[iter]<-f.lambda.aux(lambda=ref,RATE=aux1[iter],theta=chain[iter,k+3],Cens=Cens[obs])
    aux3[iter]<-integrate(f.lambda.aux,lower=0,upper=Inf,RATE=aux1[iter],theta=chain[iter,k+3],Cens=Cens[obs])$value
  }
  aux=mean(aux2/aux3)
  return(aux)
}

################################################################################
# CORRECTION FACTOR/PRIOR FOR BAYES FACTOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.RMWLN ONLY)
CFP.obs.RMWLN<-function(ref,obs,N,thin,Q,burn,Time,Cens,X,chain0,typ.theta,hyp.theta,hyp1.gam=1,hyp2.gam=1,ar=0.44,EXP=FALSE)
{
  n=length(Time); k=dim(X)[2]
  beta0=apply(chain0[,1:k],2,"median"); gam0=median(chain0[,k+1]); theta0=median(chain0[,k+2]); lambda0=apply(chain0[,(k+3):(k+2+n)],2,"median")

  chain=MCMCR.RMWLN.lambda.obs(ref,obs,N,thin,Q,beta0,gam0,theta0,lambda0,Time,Cens,X,typ.theta,hyp.theta,hyp1.gam,hyp2.gam,ar,EXP)
  chain=chain[-(1:burn),]; N.aux=dim(chain)[1]; aux1=rep(0,times=N.aux)

  for(iter in 1:N.aux)
  {
    aux1[iter]=1/dlnorm(x=ref,meanlog=0,sdlog=(chain[iter,k+2])^2)
  }

  aux=mean(aux1)
  return(aux)
}

