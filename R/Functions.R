
#' @title MCMC sampler to fit a RMW AFT regression model
#'
#' @description MCMC sampler to fit a RMW AFT regression model
#'
#' @param N Total number of iterations for the MCMC sampler. Use \code{N>=max(4,Thin)}, \code{N} being a multiple of \code{Thin}.
#' @param Thin Thining period for the MCMC sampler. Use \code{Thin>=2}.
#' @param Burn Burn-in period for the MCMC sampler. Use \code{Burn>=1}, \code{Burn<N}, \code{Burn} being a multiple of \code{Thin}.
#' @param Time Vector of length \code{n} containing survival times
#' @param Event Vector of length \code{n} containing event indicators (\code{TRUE} / \code{1} is the event is observed, \code{FALSE} / \code{0} is the event is censored)
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
#' X=cbind(Intercept,x1); rm(Intercept)
#' Time=alloauto$time; Event=alloauto$delta
#'
#' RMWreg_MCMC(N = 100, Thin = 2, Burn = 50,
#'             Time, Event, X,
#'             Mixing = "None", BaseModel = "Weibull",
#'             PriorCV = "Pareto", PriorMeanCV = 1.5,
#'             Hyp1Gam = 1, Hyp2Gam = 1)
#'
#' @author Catalina A. Vallejos \email{cvallejos@@turing.ac.uk}
RMWreg_MCMC <- function(N, Thin, Burn,
                        Time, Event, X,
                        Mixing = "None",
                        BaseModel = "Weibull",
                        PriorCV = "Pareto", PriorMeanCV = 1.5,
                        Hyp1Gam = 1, Hyp2Gam = 1,
                        ... )
{
  # No of samples and covariates
  n = nrow(X); k = ncol(X)

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
  if(PriorMeanCV >= sqrt(3) & Mixing == "InvGamma")
    stop("Invalid value for 'PriorMeanCV' (must be below sqrt(3) when Mixing == 'InvGamma')")
  if(PriorMeanCV >= sqrt(5) & Mixing == "InvGauss")
    stop("Invalid value for 'PriorMeanCV' (must be below sqrt(5) when Mixing == 'InvGauss')")

  args <- list(...)

  # Starting values
  if ("Start" %in% names(args)) { Start = args$Start }
  else Start = list("beta0" = rnorm(k),
                    "gam0" = rexp(1,1)+1,
                    "theta0" = rexp(1,1)+1)
  if ("StartAdapt" %in% names(args)) { StartAdapt = args$StartAdapt }
  else StartAdapt = list("LSbeta0" = rep(0, times = k),
                         "LSgam0" = 0,
                         "LStheta0" = 0)

  # Additional parameters related to adaptive proposals
  if ("Adapt" %in% names(args)) { Adapt = args$Adapt }
  else Adapt = TRUE
  StopAdapt = ifelse("StopAdapt" %in% names(args), args$StopAdapt, Burn)
  AR = ifelse("AR" %in% names(args), args$AR, 0.44)

  # Extra parameters to allow fixed values of selected parameters
  FixBetaJ = ifelse("FixBetaJ" %in% names(args), args$FixBetaJ, 0)
  FixGam = ifelse("FixGam" %in% names(args), args$FixGam, FALSE)
  FixTheta = ifelse("FixTheta" %in% names(args), args$FixTheta, FALSE)

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
  if(PriorCV == "Pareto") { HypTheta = 1 / (PriorMeanCV - 1)} # E(cv) = 1 + 1/a
  if(PriorCV == "TruncExp") { HypTheta = PriorMeanCV / (PriorMeanCV - 1)} # E(cv) = b/(b − 1)

  Time = system.time(Chain <- HiddenRMWreg_MCMC(N, Thin, Burn,
                            Time, Event, X,
                            Mixing,
                            Hyp1Gam, Hyp2Gam,
                            PriorCV, HypTheta,
                            Start$beta0, Start$gam0, Start$theta0,
                            as.numeric(Adapt), AR,
                            as.numeric(StoreAdapt), StopAdapt,
                            StartAdapt$LSbeta0, StartAdapt$LSgam0, StartAdapt$LStheta0,
                            FixBetaJ, as.numeric(FixGam), as.numeric(FixTheta),
                            as.numeric(PrintProgress)))

  cat("-------------------------------------------------------------------- \n")
  cat("MCMC running time")
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

  if (StoreAdapt)
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




