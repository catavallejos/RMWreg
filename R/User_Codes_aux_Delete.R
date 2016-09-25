###############################################################################
###############################################################################
######################## ALL CODES FOR Weibull  MODEL #########################
###############################################################################
###############################################################################


###############################################################################
###############################################################################
######################## ALL CODES FOR RMWEXP   MODEL #########################
###############################################################################
###############################################################################



###############################################################################
###############################################################################
######################## ALL CODES FOR RMWGAM   MODEL #########################
###############################################################################
###############################################################################

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
