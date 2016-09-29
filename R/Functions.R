
#' @title MCMC sampler to fit a RMW AFT regression model
#'
#' @description MCMC sampler to fit a RMW AFT regression model
#'
#' @param N Total number of iterations for the MCMC sampler. Use \code{N>=max(4,Thin)}, \code{N} being a multiple of \code{Thin}.
#' @param Thin Thining period for the MCMC sampler. Use \code{Thin>=2}.
#' @param Burn Burn-in period for the MCMC sampler. Use \code{Burn>=1}, \code{Burn<N}, \code{Burn} being a multiple of \code{Thin}.
#' @param Time Vector of length \code{n} containing survival times
#' @param Event Vector of length \code{n} containing event indicators (\code{TRUE} / \code{1} is the event is observed, \code{FALSE} / \code{0} is the event is Eventored)
#' @param Mixing Mixing distribution assigned to the (frailty) random effects. Possible values are
#' \code{"None", "Exponential", "Gamma", "InvGamma", "InvGauss", "LogNormal"}
#' @param BaseModel If \code{BaseModel = "Weibull"}, a RMW regression is used. If \code{BaseModel = "Exponential"}, a RME regression is used.
#' @param PriorCV Type of prior assigned to the coefficient of variation of the survival times.
#' Possible values are \code{"Pareto", "TruncExp"}
#' @param PriorMeanCV Ellicited prior mean of the coefficient of variation (\code{PriorMeanCV > 1}).
#' If \code{Mixing = "InvGamma"}, \code{PriorMeanCV} must be below \code{sqrt{3}}.
#' If \code{Mixing = "InvGauss"}, \code{PriorMeanCV} must be below \code{sqrt{5}}.
#' Default: \code{PriorMeanCV = 1.5}
#' @param Hyp1Gam Shape hyper-parameter for the Gamma(\code{Hyp1Gam}, \code{Hyp2Gam}) assigned to \code{gam}.
#' @param Hyp2Gam Rate hyper-parameter for the Gamma(\code{Hyp1Gam}, \code{Hyp2Gam}) assigned to \code{gam}.
#' @param ... Optional parameters.
#' \describe{
#'
#' \item{\code{AR}}{Optimal acceptance rate for adaptive Metropolis Hastings updates.
#' It must be a positive number between 0 and 1.
#' Default (and recommended): \code{ar = 0.44}}.
#' \item{\code{StopAdapt}}{Iteration at which adaptive proposals are not longer adapted.
#' Use \code{stopAdapt>=1}. Default: \code{StopAdapt = Burn}.}
#' \item{\code{StoreChains}}{If \code{StoreChains = TRUE}, MCMC chains of each parameter are stored in separate .txt files.
#' (\code{RunName} argument used for file names). Default: \code{StoreChains = FALSE}.}
#' \item{\code{StoreAdapt}}{If \code{StoreAdapt = TRUE}, trajectory of adaptive proposal variances (log scale)
#' for each parameter are stored in separate .txt files. (\code{RunName} argument used for file names).
#' Default: \code{StoreAdapt = FALSE}.}
#' \item{\code{StoreDir}}{Directory where MCMC chain will be stored (only required if \code{Store = TRUE}).
#' Default: \code{StoreDir = getwd()}.}
#' \item{\code{RunName}}{Run-name to be used when storing chains and/or adaptive proposal variances in .txt files.}
#' \item{\code{PrintProgress}}{If \code{PrintProgress = TRUE}, intermediate output is displayed in the console. }
#' }
#'
#' @return A \code{list} containing MCMC draws for all parameters.
#'
#' @examples
#'
#' library(KMsurv)
#' data(alloauto)
#' n=dim(alloauto)[1]; k=2
#' Intercept=rep(1,times=n); x1=alloauto$type-1
#' DesignMat=cbind(Intercept,x1); rm(Intercept)
#' Time=alloauto$time; Event=alloauto$delta
#'
#' Chain <- RMWreg_MCMC(N = 100, Thin = 2, Burn = 50,
#'                      Time, Event, DesignMat,
#'                      Mixing = "None", BaseModel = "Weibull",
#'                      PriorCV = "Pareto", PriorMeanCV = 1.5,
#'                      Hyp1Gam = 1, Hyp2Gam = 1)
#'
#' @author Catalina A. Vallejos \email{cvallejos@@turing.ac.uk}
RMWreg_MCMC <- function(N, Thin, Burn,
                        Time, Event, DesignMat,
                        Mixing = "None",
                        BaseModel = "Weibull",
                        PriorCV = "Pareto", PriorMeanCV = 1.5,
                        Hyp1Gam = 1, Hyp2Gam = 1,
                        ... )
{
  # No of samples and covariates
  n = nrow(DesignMat); k = ncol(DesignMat)

  # Validity checks for parameter values
  if (!(length(N) == 1 | length(Thin) == 1 | length(Burn) == 1))
    stop("Invalid parameter values.")
  if (!(N%%Thin == 0 & N >= max(4, Thin)))
    stop("Please use an integer value for N. It must also be a multiple of thin (N>=4)).")
  if (!(Thin%%1 == 0 & Thin >= 2))
    stop("Please use an integer value for Thin (Thin>=2).")
  if (!(Burn%%Thin == 0 & Burn < N & Burn >= 1))
    stop("Please use an integer value for Burn. It must also be lower than N and a multiple of thin (Burn>=1).")

  if(length(Time) != n | length(Event) != n)
    stop("The dimensions of the input dataset are not compatible")

  if(!(Mixing %in% c("None", "Exponential", "Gamma", "InvGamma", "InvGauss", "LogNormal")))
    stop("Invalid value for 'Mixing'")
  if(!(BaseModel %in% c("Exponential", "Weibull")))
    stop("Invalid value for 'BaseModel'")
  if(!(PriorCV %in% c("Pareto", "TruncExp")))
    stop("Invalid value for 'PriorCV'")
  if(PriorMeanCV <= 1)
    stop("Invalid value for 'PriorMeanCV'")
  if(PriorMeanCV >= sqrt(3) & Mixing == "InvGamma" & BaseModel == "Exponential")
    stop("Invalid value for 'PriorMeanCV' (must be below sqrt(3) when Mixing == 'InvGamma')")
  if(PriorMeanCV >= sqrt(5) & Mixing == "InvGauss" & BaseModel == "Exponential")
    stop("Invalid value for 'PriorMeanCV' (must be below sqrt(5) when Mixing == 'InvGauss')")

  args <- list(...)

  # Starting values
  if ("Start" %in% names(args)) { Start = args$Start }
  else Start = list("beta0" = rnorm(k),
                    "gam0" = rexp(1,1)+1,
                    "theta0" = rexp(1,1)+1)
  if(Mixing == "Gamma") {Start$theta0 = max(Start$theta0, 2/Start$gam0 + 1)}
  if(Mixing == "InvGamma") {Start$theta0 = max(Start$theta0, 1 + 1)}

  if ("StartAdapt" %in% names(args)) { StartAdapt = args$StartAdapt }
  else StartAdapt = list("LSbeta0" = rep(0, times = k),
                         "LSgam0" = 0,
                         "LStheta0" = 0)

  if("lambdaPeriod" %in% names(args)) { lambdaPeriod = args$lambdaPeriod }
  else lambdaPeriod = 5

  # Additional parameters related to adaptive proposals
  if ("Adapt" %in% names(args)) { Adapt = args$Adapt }
  else Adapt = TRUE
  StopAdapt = ifelse("StopAdapt" %in% names(args), args$StopAdapt, Burn)
  AR = ifelse("AR" %in% names(args), args$AR, 0.44)

  # Extra parameters to allow fixed values of selected parameters
  FixBetaJ = ifelse("FixBetaJ" %in% names(args), args$FixBetaJ, 0)
  FixGam = ifelse("FixGam" %in% names(args), args$FixGam, FALSE)
  FixTheta = ifelse("FixTheta" %in% names(args), args$FixTheta, FALSE)
  FixLambdaI = ifelse("FixLambdaI" %in% names(args), args$FixLambdaI, 0)
  RefLambda = ifelse("RefLambda" %in% names(args), args$RefLambda, 1)

  # Storage/console parameters
  PrintProgress = ifelse("PrintProgress" %in% names(args), args$PrintProgress, TRUE)
  StoreChains = ifelse("StoreChains" %in% names(args), args$StoreChains, FALSE)
  StoreAdapt = ifelse("StoreAdapt" %in% names(args), args$StoreAdapt, FALSE)
  StoreDir = ifelse("StoreDir" %in% names(args), args$StoreDir, getwd())
  RunName = ifelse("RunName" %in% names(args), args$RunName, "")

  # RME family
  if(BaseModel == "Exponential") { FixGam = TRUE; Start$gam0 = 1 }

  # Models for which theta is not required
  if(Mixing %in% c("None", "Exponential")) {Start$theta0 = 0; FixTheta = TRUE}

  # Hyper-parameter for theta
  if(PriorCV == "TruncExp") { HypTheta = 1 / (PriorMeanCV - 1)} # E(cv) = 1 + 1/a
  if(PriorCV == "Pareto") { HypTheta = PriorMeanCV / (PriorMeanCV - 1)} # E(cv) = b/(b − 1)

  Time = system.time(Chain <- HiddenRMWreg_MCMC(N, Thin, Burn, Time, Event, DesignMat,
                                                Mixing, Hyp1Gam, Hyp2Gam, PriorCV, HypTheta,
                                                Start$beta0, Start$gam0, Start$theta0,
                                                as.numeric(Adapt), AR, as.numeric(StoreAdapt), StopAdapt,
                                                StartAdapt$LSbeta0, StartAdapt$LSgam0, StartAdapt$LStheta0,
                                                FixBetaJ, as.numeric(FixGam), as.numeric(FixTheta),
                                                as.numeric(PrintProgress), lambdaPeriod,
                                                FixLambdaI, RefLambda))

  cat("-------------------------------------------------------------------- \n")
  cat("MCMC running time \n")
  cat("-------------------------------------------------------------------- \n")
  print(Time)
  cat("\n")

  OldDir = getwd()

  if (StoreChains)
  {
    setwd(StoreDir)
    cat("-------------------------------------------------------------------- \n")
    cat("Storing MCMC chains of model parameters as .txt files in \n")
    cat(paste0("'", StoreDir, "' directory ... \n"))
    cat(paste0("Files are indexed using the name '", RunName, "'... \n"))
    cat("-------------------------------------------------------------------- \n")
    write.table(Chain$beta, paste0("chain_beta_", RunName, ".txt"),
                col.names = FALSE, row.names = FALSE)
    if(BaseModel == "Weibull")
    {
      write.table(Chain$gam, paste0("chain_gam_", RunName, ".txt"),
                  col.names = FALSE, row.names = FALSE)
    }
    if(Mixing %in% c("Gamma", "InvGamma", "InvGauss", "LogNormal"))
    {
      write.table(Chain$theta, paste0("chain_theta_", RunName, ".txt"),
                  col.names = FALSE, row.names = FALSE)
    }
    setwd(OldDir)
  }

  if (StoreAdapt & StoreChains)
  {
    setwd(StoreDir)
    cat("-------------------------------------------------------------------- \n")
    cat("Storing trajectories of adaptive proposal variances (log-scale) as .txt files in \n")
    cat(paste0("'", StoreDir, "' directory ... \n"))
    cat(paste0("Files are indexed using the name '", RunName, "'... \n"))
    cat("-------------------------------------------------------------------- \n")

    write.table(Chain$ls.beta, paste0("chain_ls.beta_", RunName, ".txt"),
                col.names = FALSE, row.names = FALSE)
    if(BaseModel == "Weibull")
    {
      write.table(Chain$ls.gam, paste0("chain_ls.gam_", RunName, ".txt"),
                  col.names = FALSE, row.names = FALSE)
    }
    if(Mixing %in% c("Gamma", "InvGamma", "InvGauss", "LogNormal"))
    {
      write.table(Chain$ls.theta, paste0("chain_ls.theta_", RunName, ".txt"),
                  col.names = FALSE, row.names = FALSE)
    }
    setwd(OldDir)
  }

  cat("-------------------------------------------------------------------- \n")
  cat("Output \n")
  cat("-------------------------------------------------------------------- \n")

  return(Chain)

}

# LOG-LIKELIHOOD
Hiddenf.joint.RME.LN=function(lambda,Time,alpha,theta)
{
  aux = exp(-alpha*lambda*Time) * exp(-(1/(2*theta))*(log(lambda))^2)
  return(aux)
}

Hiddenf.RME.LN=function(Time,alpha,theta)
{
  aux = alpha*(2*pi*theta)^(-1/2) * integrate(Hiddenf.joint.RME.LN,lower=0,upper=Inf,Time=Time,alpha=alpha,theta=theta)$value
  return(aux)
}

HiddenS.joint.RME.LN=function(lambda,Time,alpha,theta)
{
  aux = (lambda^(-1))*exp(-alpha*lambda*Time) * exp(-(1/(2*theta))*(log(lambda))^2)
  return(aux)
}

HiddenS.RME.LN=function(Time,alpha,theta)
{
  aux = (2*pi*theta)^(-1/2) * integrate(HiddenS.joint.RME.LN,lower=0.00001,upper=Inf,Time=Time,alpha=alpha,theta=theta)$value
  return(aux)
}

Hiddenlog.lik.RMW.LN.aux<-function(Time, Event, DesignMat, beta, gam = 1, theta)
{
  RATE = exp(-gam * DesignMat%*%beta)
  f.aux <- function(i, Event, Time, RATE, theta)
  {
    if(Event[i]==1) {out = log(Hiddenf.RME.LN(Time=Time[i],alpha=RATE[i],theta=theta))}
    if(Event[i]==0) {out = log(HiddenS.RME.LN(Time=Time[i],alpha=RATE[i],theta=theta))}
    return(out)
  }
  aux = sapply(as.list(1:length(Time)), FUN = f.aux,
               Event = Event, Time = Time, RATE = RATE, theta = theta, simplify = TRUE)
  return(aux)
}

Hiddenlog.lik.RMWLN<-function(Time, Event, DesignMat, beta, gam=1, theta, BaseModel)
{
  aux = sum(Hiddenlog.lik.RMW.LN.aux(Time^gam,Event,DesignMat,beta,gam,theta) + I(BaseModel=="Weibull")*(Event*log(gam)+Event*(gam-1)*log(Time)))
  return(aux)
}

HiddenRMWreg_DIC_LN <- function(Chain, Time, Event, DesignMat, Mixing, BaseModel)
{
  beta = Chain$beta; gam = Chain$gam; theta = Chain$theta
  beta_hat = apply(beta,2,"median");
  gam_hat = median(gam); theta_hat = median(theta)

  # LOG-LIKELIHOOD FOR EACH DRAW
  f.aux <- function(iter, Time, Event, DesignMat,
                    beta, gam, theta, BaseModel)
  {
    out = Hiddenlog.lik.RMWLN(Time, Event, DesignMat,
                              beta=beta[iter,], gam = gam[iter,], theta = theta[iter,], BaseModel)
    return(out)
  }
  L = sapply(as.list(1:nrow(beta)), FUN = f.aux,
             Time = Time, Event = Event, DesignMat = DesignMat,
             beta = beta, gam = gam, theta = theta, BaseModel = BaseModel, simplify = TRUE)

#  for(iter in 1:nrow(beta))
#  {
#    L[iter] = Hiddenlog.lik.RMWLN(Time, Event, DesignMat,
#                                  beta=beta[iter,], gam = gam[iter,], theta = theta[iter,], BaseModel)
#  }
  pd = -2*mean(L) + 2*Hiddenlog.lik.RMWLN(Time, Event, DesignMat,
                                          beta = beta_hat, gam = gam_hat, theta = theta_hat, BaseModel)
  return(-2*mean(L) + pd)
}

RMWreg_DIC <- function(Chain, Time, Event, DesignMat, Mixing, BaseModel)
{
  if(Mixing %in% c("None", "Exponential", "Gamma","InvGamma", "InvGauss"))
  {
    DIC = HiddenRMWreg_DIC(Chain, Time, Event, DesignMat, Mixing, BaseModel)
  }
  if(Mixing == "LogNormal")
  {
    DIC = HiddenRMWreg_DIC_LN(Chain, Time, Event, DesignMat, Mixing, BaseModel)
  }
  return(DIC)
}

HiddenRMWreg_CaseDeletion_LN <- function(Chain, Time, Event, DesignMat, Mixing, BaseModel)
{
  beta = Chain$beta; gam = Chain$gam; theta = Chain$theta

  N = nrow(beta); n = length(Time);
  logCPO=rep(0,times=n); KL.aux=rep(0,times=n)

  f.aux <- function(iter, Time, Event, DesignMat,
                    beta, gam, theta, BaseModel)
  {
    out = Hiddenlog.lik.RMWLN(Time, Event, DesignMat,
                              beta = beta[iter,], gam = gam[iter,], theta = theta[iter,],
                              BaseModel)
    return(out)
  }

  for(i in 1:n)
  {
    aux2 = sapply(as.list(1:N), FUN = f.aux,
                  Time = Time[i], Event = Event[i], DesignMat = DesignMat[i,],
                  beta = beta, gam = gam, theta = theta, BaseModel = BaseModel, simplify = TRUE)
    aux1 = exp(-aux2)

#    aux1 = rep(0, times = N); aux2 = rep(0, times = N)
#    for(ITER in 1:N)
#    {
#      aux2[ITER] = Hiddenlog.lik.RMWLN(Time[i], Event[i], DesignMat[i,],
#                                       beta = beta[ITER,], gam = gam[ITER,], theta = theta[ITER,],
#                                       BaseModel)
#      aux1[ITER] = exp(-aux2[ITER])
#    }
    logCPO[i] = -log(mean(aux1))
    KL.aux[i] = mean(aux2)
  }
  KL=KL.aux-logCPO
  CALIBRATION=0.5*(1+sqrt(1-exp(-2*KL)))
  return(cbind(logCPO,KL,CALIBRATION))
}

RMWreg_CaseDeletion <- function(Chain, Time, Event, DesignMat, Mixing, BaseModel)
{
  if(Mixing %in% c("None", "Exponential", "Gamma","InvGamma", "InvGauss"))
  {
    CD = HiddenRMWreg_CaseDeletion(Chain, Time, Event, DesignMat, Mixing, BaseModel)
  }
  if(Mixing == "LogNormal")
  {
    CD = HiddenRMWreg_CaseDeletion_LN(Chain, Time, Event, DesignMat, Mixing, BaseModel)
  }
  return(CD)
}

# ACCEPTANCY PROBABILITY FOR GRWMH beta[j]
HiddenAcceptProb_betaJ <- function(beta0, beta1, gam, Time, Event, DesignMat, lambda, j)
{
  # AUXILIARY QUANTITIES
  Xj <- as.matrix(DesignMat[,-j]); bj <- beta0[-j]; aux0 <- exp(-gam*as.vector(Xj%*%bj))
  # ACCEPTANCE PROBABILITY
  prob = (beta0[j] - beta1[j]) * gam * sum(DesignMat[,j]*Event)
  prob = prob + sum(aux0 * (Time^gam) * lambda * (exp(-gam*beta0[j]*DesignMat[,j])-exp(-gam*beta1[j]*DesignMat[,j])))

  prob=min(1,exp(prob))

  return(prob)
}

# ACCEPTANCY PROBABILITY FOR GRWMH GAMMA
HiddenAcceptProb_gam <- function(gam0, gam1, Time, Event, DesignMat, beta, theta, lambda,
                                 PriorCV, HypTheta, Hyp1Gam, Hyp2Gam, Mixing,
                                 lower.bound = 0.06)
{
  if(gam1 < lower.bound) { prob = 0 }
  else
  {
    prob = sum(Event) * log(gam1/gam0) + (gam1-gam0) * sum(Event*(log(Time)-as.numeric(DesignMat%*%beta)))
    prob = prob - sum(lambda*((exp(-as.numeric(DesignMat%*%beta))*Time)^gam1 - (exp(-as.numeric(DesignMat%*%beta))*Time)^gam0))
    prob = prob + dgamma(gam1, shape = Hyp1Gam, rate = Hyp2Gam, log = TRUE)
    prob = prob - dgamma(gam0, shape = Hyp1Gam, rate = Hyp2Gam, log = TRUE)

    if(!(Mixing %in% c("None", "Exponential")))
    {
      prob = prob + HiddenLogPriorTheta(theta, gam1, HypTheta, PriorCV, Mixing)
      prob = prob - HiddenLogPriorTheta(theta, gam0, HypTheta, PriorCV, Mixing)
    }

    prob = min(1,exp(prob))
  }

  return(prob)
}


# ACCEPTANCE PROBABILITY THETA
HiddenAcceptProb_theta <- function(theta0, theta1, gam, lambda, PriorCV, HypTheta, Mixing)
{
  n = length(lambda)

  if(Mixing == "Gamma")
  {
    prob = n * (theta1*log(theta1) - theta0*log(theta0)) -n*(lgamma(theta1) - lgamma(theta0))
    prob = prob + (theta1 - theta0) * sum(log(lambda)) - (theta1 - theta0) * sum(lambda)
  }
  if(Mixing == "InvGamma")
  {
    prob = n * (lgamma(theta0) - lgamma(theta1)) + (theta0 - theta1) * sum(log(lambda))
  }
  if(Mixing == "InvGauss")
  {
    prob = 0.5 * (theta0^(-2) - theta1^(-2)) * sum(lambda) - n * (theta0^(-1)-theta1^(-1))
  }
  if(Mixing == "LogNormal")
  {
    prob = (n/2) * log(theta0/theta1) - 0.5 * (theta1^(-1) - theta0^(-1)) * sum((log(lambda))^2)
  }

  prob = prob + HiddenLogPriorTheta(theta1, gam, HypTheta, PriorCV, Mixing)
  prob = prob - HiddenLogPriorTheta(theta0, gam, HypTheta, PriorCV, Mixing)

  prob = min(1,exp(prob))

  return(prob)
}



# LOG-MARGINAL LIKELIHOOD ESTIMATOR
RMWreg_logML <- function(Chain,
                         Time, Event, DesignMat,
                         PriorCV = "Pareto", PriorMeanCV = 1.5,
                         Hyp1Gam = 1, Hyp2Gam = 1,
                         Thin = 10, lambdaPeriod = 5, AR = 0.44, Mixing, BaseModel)
{
  # EXTRACTING MCMC CHAINS
  beta = Chain$beta; gam = Chain$gam; theta = Chain$theta
  ls.beta = Chain$ls.beta; ls.gam = Chain$ls.gam; ls.theta = Chain$ls.theta

  # SAMPLE SIZE, NUMBER OF DRAWS AND NUMBER OF REGRESSORS
  n = length(Time); N = nrow(beta); k = ncol(beta)

  # Hyper-parameter for theta
  if(PriorCV == "TruncExp") { HypTheta = 1 / (PriorMeanCV - 1)} # E(cv) = 1 + 1/a
  if(PriorCV == "Pareto") { HypTheta = PriorMeanCV / (PriorMeanCV - 1)} # E(cv) = b/(b − 1)

  # Starting values
  if(k>1)	{ beta0 = apply(beta, 2, "median") }
  else		{ beta0 = median(beta) }
  gam0 = median(gam)
  theta0 = median(theta)
  if(k>1)	{ ls.beta0 = apply(ls.beta, 2, "median") }
  else		{ ls.beta0 = median(ls.beta) }
  ls.gam0 = median(ls.gam)
  ls.theta0 = median(ls.theta)

  # NON-ADAPTIVE RUN OF THE CHAIN
  chain.nonadapt = RMWreg_MCMC(N*Thin, Thin, Burn = round(0.25*N)*Thin, Time, Event, DesignMat,
                               Mixing = Mixing, BaseModel = BaseModel,
                               PriorCV = PriorCV, PriorMeanCV = PriorMeanCV,
                               Hyp1Gam = Hyp1Gam, Hyp2Gam = Hyp2Gam, AR = AR,
                               lambdaPeriod = lambdaPeriod,
                               PrintProgress = FALSE, Adapt = FALSE,
                               Start = list("beta0" = beta0,
                                            "gam0" = gam0,
                                            "theta0" = theta0),
                               StartAdapt = list("LSbeta0" = ls.beta0,
                                                 "LSgam0" = ls.gam0,
                                                 "LStheta0" = ls.theta0))
  if(k>1)	{ beta.star = apply(chain.nonadapt$beta, 2, "median") }
  else		{ beta.star = median(chain.nonadapt$beta) }
  gam.star = median(chain.nonadapt$gam)
  theta.star = median(chain.nonadapt$theta)
  N.aux = nrow(chain.nonadapt$beta)

  # LIKELIHOOD ORDINATE
  if(Mixing %in% c("None", "Exponential", "Gamma", "InvGamma", "InvGauss"))
  {
    LL.ord = HiddenLogLik(Time, Event, DesignMat,
                          beta.star, gam.star, theta.star, Mixing, BaseModel)
  }
  if(Mixing == "LogNormal")
  {
    LL.ord = Hiddenlog.lik.RMWLN(Time, Event, DesignMat,
                                 beta.star, gam.star, theta.star, BaseModel)
  }
  print("Likelihood ordinate ready!")

  # PRIOR ORDINATE
  LP.ord = 0
  if(BaseModel == "Weibull")
  {
    LP.ord = LP.ord + dgamma(gam.star, shape = Hyp1Gam, rate = Hyp2Gam, log = TRUE)
  }
  if(!(Mixing %in% c("None", "Exponential")))
  {
    LP.ord = LP.ord + HiddenLogPriorTheta(theta.star, gam.star, HypTheta, PriorCV, Mixing)
  }
  print("Prior ordinate ready!")


  if(!(Mixing %in% c("None", "Exponential")))
  {
    # REDUCED CHAIN WITH FIXED THETA
    chain.theta = RMWreg_MCMC(N*Thin, Thin, Burn = round(0.25*N)*Thin, Time, Event, DesignMat,
                              Mixing = Mixing, BaseModel = BaseModel,
                              PriorCV = PriorCV, PriorMeanCV = PriorMeanCV,
                              Hyp1Gam = Hyp1Gam, Hyp2Gam = Hyp2Gam, AR = AR,
                              lambdaPeriod = lambdaPeriod,
                              PrintProgress = FALSE, Adapt = FALSE,
                              Start = list("beta0" = beta0,
                                           "gam0" = gam0,
                                           "theta0" = theta0),
                              StartAdapt = list("LSbeta0" = ls.beta0,
                                                "LSgam0" = ls.gam0,
                                                "LStheta0" = ls.theta0),
                              FixTheta = TRUE)

    # POSTERIOR ORDINATE - theta
    po1.theta = rep(0,times=N.aux); po2.theta = rep(0,times=N.aux)
    for(i in 1:N.aux)
    {
      po1.theta[i] = HiddenAcceptProb_theta(theta0 = chain.nonadapt$theta[i], theta1 = theta.star,
                                            gam = chain.nonadapt$gam[i], lambda = t(chain.nonadapt$lambda[i,]),
                                            PriorCV, HypTheta, Mixing) * dnorm(x = theta.star,
                                                                               mean = chain.nonadapt$theta[i],
                                                                               sd = sqrt(exp(ls.theta0)))
      theta.aux = rnorm(n = 1, mean = theta.star, sd = sqrt(exp(ls.theta0)) )
      if(theta.aux <= 2/gam.star & Mixing == "Gamma") { po2.theta[i] = 0 }
      else
      {
        if(theta.aux <= 1 & Mixing == "InvGamma") { po2.theta[i] = 0 }
        else
        {
          if(theta.aux <= 0) { po2.theta[i] = 0 }
          else
          {
            po2.theta[i] = HiddenAcceptProb_theta(theta0 = theta.star, theta1 = theta.aux,
                                                  gam = chain.theta$gam[i],
                                                  lambda = t(chain.theta$lambda[i,]),
                                                  PriorCV, HypTheta, Mixing)
          }
        }
      }
    }
    LPO.theta = log(mean(po1.theta)/mean(po2.theta))
    print("Posterior ordinate theta ready!")
  }
  else { chain.theta = chain.nonadapt; LPO.theta = 0 }

  if(BaseModel == "Weibull")
  {
    # REDUCED CHAIN WITH FIXED THETA + GAMMA
    chain.gam = RMWreg_MCMC(N*Thin, Thin, Burn = round(0.25*N)*Thin, Time, Event, DesignMat,
                            Mixing = Mixing, BaseModel = BaseModel,
                            PriorCV = PriorCV, PriorMeanCV = PriorMeanCV,
                            Hyp1Gam = Hyp1Gam, Hyp2Gam = Hyp2Gam, AR = AR,
                            lambdaPeriod = lambdaPeriod,
                            PrintProgress = FALSE, Adapt = FALSE,
                            Start = list("beta0" = beta0,
                                         "gam0" = gam0,
                                         "theta0" = theta0),
                            StartAdapt = list("LSbeta0" = ls.beta0,
                                              "LSgam0" = ls.gam0,
                                              "LStheta0" = ls.theta0),
                            FixTheta = TRUE, FixGam = TRUE)

    # POSTERIOR ORDINATE - gam
    po1.gam = rep(0, times = N.aux); po2.gam = rep(0, times = N.aux)
    for(i in 1:N.aux)
    {
      po1.gam[i] = HiddenAcceptProb_gam(gam0 = chain.theta$gam[i], gam1 = gam.star,
                                        Time, Event, DesignMat,
                                        beta = as.vector(chain.theta$beta[i,]), theta = theta.star,
                                        lambda = t(chain.theta$lambda[i,]),
                                        PriorCV, HypTheta, Hyp1Gam, Hyp2Gam,
                                        Mixing, lower.bound = 0.06) * dnorm(x = gam.star,
                                                                            mean = chain.theta$gam[i],
                                                                            sd = sqrt(exp(ls.gam0)))
      gam.aux = rnorm(n = 1, mean = gam.star, sd = sqrt(exp(ls.gam0)))
      if(gam.aux <= 0.06 | gam.aux <= 2/theta.star & Mixing == "Gamma") { po2.gam[i]=0 }
      else
      {
        if(gam.aux <= 0.06) { po2.gam[i] = 0 }
        else
        {
          po2.gam[i] = HiddenAcceptProb_gam(gam0 = gam.star, gam1 = gam.aux,
                                            Time, Event, DesignMat,
                                            beta = as.vector(chain.gam$beta[i,]), theta = theta.star,
                                            lambda = t(chain.gam$lambda[i,]),
                                            PriorCV, HypTheta, Hyp1Gam, Hyp2Gam,
                                            Mixing, lower.bound = 0.06)
          if(is.nan(po2.gam[i])) {print(gam.aux)}
        }
      }
    }
    LPO.gam = log(mean(po1.gam) / mean(po2.gam))
    print("Posterior ordinate gamma ready!")
  }
  else { LPO.gam = 0 }

  # POSTERIOR ORDINATE - beta
  if(BaseModel == "Weibull") { chain.prev = chain.gam }
  else { chain.prev = chain.theta }

  LPO.beta = rep(0, times = k)

  for(j.beta in 0:(k-1))
  {
    print(j.beta)
    beta0 = t(chain.prev$beta[N.aux,]); beta0[j.beta+1] = beta.star[j.beta+1]
    chain.next = RMWreg_MCMC(N*Thin, Thin, Burn = round(0.25*N)*Thin, Time, Event, DesignMat,
                             Mixing = Mixing, BaseModel = BaseModel,
                             PriorCV = PriorCV, PriorMeanCV = PriorMeanCV,
                             Hyp1Gam = Hyp1Gam, Hyp2Gam = Hyp2Gam, AR = AR,
                             lambdaPeriod = lambdaPeriod,
                             PrintProgress = FALSE, Adapt = FALSE,
                             Start = list("beta0" = beta0,
                                          "gam0" = gam0,
                                          "theta0" = theta0),
                             StartAdapt = list("LSbeta0" = ls.beta0,
                                               "LSgam0" = ls.gam0,
                                               "LStheta0" = ls.theta0),
                             FixTheta = TRUE, FixGam = TRUE, FixBetaJ = j.beta + 1)

    po1.beta = rep(0, times = N.aux); po2.beta = rep(0, times = N.aux)
    for(i in 1:N.aux)
    {
      beta0 = as.vector(t(chain.prev$beta[i,]))
      beta1 = beta0; beta1[j.beta+1] = beta.star[j.beta+1]
      po1.beta[i] = HiddenAcceptProb_betaJ(beta0, beta1, gam = chain.prev$gam[i],
                                           Time, Event, DesignMat,
                                           lambda = t(chain.prev$lambda[i,]), j = j.beta+1) *dnorm(x = beta.star[j.beta+1],
                                                                                                   mean = as.numeric(chain.prev$beta[i,j.beta+1]),
                                                                                                   sd = sqrt(exp(ls.beta0[j.beta+1])) )
      betaj.aux = rnorm(n = 1, mean = beta.star[j.beta+1], sd = sqrt(exp(ls.beta0[j.beta+1])) )
      beta2 = beta.star; beta2[j.beta+1] = betaj.aux
      po2.beta[i] = HiddenAcceptProb_betaJ(beta0 = beta.star, beta1 = beta2, gam = chain.next$gam[i],
                                           Time, Event, DesignMat,
                                           lambda = t(chain.next$lambda[i,]), j = j.beta+1)
    }
    LPO.beta[j.beta+1] = log(mean(po1.beta) / mean(po2.beta))
    chain.prev=chain.next
  }
  print("Posterior ordinate beta ready!")

  # LOG-MARGINAL LIKELIHOOD
  print(paste("LL.ord:", LL.ord))
  print(paste("LP.ord:", LP.ord))
  print(paste("LPO.theta:", LPO.theta))
  print(paste("LPO.gam:", LPO.gam))
  print(paste("sum(LPO.beta):", (LPO.beta)))

  LML = LL.ord + LP.ord - LPO.theta - LPO.gam - sum(LPO.beta)

  return(LML)
}


RMWreg_BFoutlier <- function(Chain, RefLambda,
                             Time, Event, DesignMat,
                             Mixing, BaseModel, ...)
{
  args <- list(...)

  if(Mixing == "None") stop("This function cannot be applied to 'Mixing = 'None''")
  if(Mixing %in% c("Gamma", "InvGamma", "InvGauss", "LogNormal"))
  {
    if(!("PriorCV" %in% names(args))) stop("Value of 'PriorCV' is missing")
    else { PriorCV = args$PriorCV }
    if(!("HypTheta" %in% names(args))) stop("Value of 'HypTheta' is missing")
    else { HypTheta = args$HypTheta }
    if(!("Hyp1Gam" %in% names(args))) stop("Value of 'Hyp1Gam' is missing")
    else { Hyp1Gam = args$Hyp1Gam }
    if(!("Hyp2Gam" %in% names(args))) stop("Value of 'Hyp2Gam' is missing")
    else { Hyp2Gam = args$Hyp2Gam }
    if(!("Thin" %in% names(args))) stop("Value of 'Thin' is missing")
    else { Thin = args$Thin }
    if(!("Burn" %in% names(args))) stop("Value of 'Burn' is missing")
    else { Burn = args$Burn }
    if(!("lambdaPeriod" %in% names(args))) stop("Value of 'lambdaPeriod' is missing")
    else { lambdaPeriod = args$lambdaPeriod }
    if(!("AR" %in% names(args))) stop("Value of 'AR' is missing")
    else { AR = args$AR }

    BF = 1
  }
  if(Mixing == "Exponential")
  {
    BF = HiddenRMWreg_BFoutlier(Chain, RefLambda,
                                Time, Event, DesignMat,
                                PriorCV = "None", HypTheta = 0,
                                Hyp1Gam = 0, Hyp2Gam = 0,
                                Mixing, BaseModel,
                                thin = 0, lambdaPeriod = 0, ar = 0)
  }
  return(BF)
}

