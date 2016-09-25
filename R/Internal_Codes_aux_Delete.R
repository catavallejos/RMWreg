################################################################################
################################################################################
# FUNCTIONS SPECIFIC FOR RMWEXP MODEL
################################################################################
################################################################################

################################################################################
################################################################################
# FUNCTIONS SPECIFIC FOR RMWGAM MODEL
################################################################################
################################################################################

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

