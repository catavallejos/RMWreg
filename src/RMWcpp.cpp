#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
//File: matprod_arma.cpp
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

#define ZTOL sqrt(DOUBLE_EPS)

/*
* Rgig is an adaptation of the rgig.c function implemented by
* Ester Pantaleo and Robert B. Gramacy, 2010
* (which in turn, was adapted from the C code in the ghyp v_1.5.2 package for R, rgig function)
* Our adaptation makes the function compatible with the Rcpp classes
*/

/* Auxiliary function taken from rgig.c */
double gig_y_gfn(double y, double m, double beta, double lambda)
{
  /*
  * gig_y_gfn:
  *
  * evaluate the function that we need to find the root
  * of in order to construct the optimal rejection
  * sampler to obtain GIG samples
  */
  double y2, g;
  y2 = y * y;
  g = 0.5 * beta * y2 * y;
  g -= y2 * (0.5 * beta * m + lambda + 1.0);
  g += y * ((lambda - 1.0) * m - 0.5 * beta) + 0.5 * beta * m;
  return(g);
}

/* Auxiliary function taken from rgig.c */
double rinvgauss(const double mu, const double lambda)
{
  /*
  * rinvgauss:
  *
  * Michael/Schucany/Haas Method for generating Inverse Gaussian
  * random variable with mean mu and shape parameter lambda,
  * as given in Gentle's book on page 193
  */
  double u, y, x1, mu2, l2;

  y = norm_rand(); /// ******
  y *= y;
  mu2 = mu * mu;
  l2 = 2.0*lambda;
  x1 = mu + mu2*y/l2 - (mu/l2)* sqrt(4.0*mu*lambda*y + mu2*y*y);

  u = unif_rand(); // ******
  if(u <= mu/(mu + x1)) return x1;
  else return mu2/x1;
}

/* Auxiliary function taken from rgig.c */
double zeroin_gig(double ax,double bx,double f(double x, double m, double beta, double lambda),double tol,double m,double beta,double lambda)  /* An estimate to the root  */
//double ax;				/* Left border | of the range	*/
//double bx;  				/* Right border| the root is seeked*/
/* Function under investigation	*/
//double (*f)(double x, double m, double beta, double lambda);
//double tol;				/* Acceptable tolerance	*/
//double m;                               /* specific to gig_y_gfn */
//double beta;                            /* specific to gig_y_gfn */
//double lambda;                          /* specific to gig_y_gfn */
{

  double a,b,c;				/* Abscissae, descr. see above	*/
double fa;				/* f(a)				*/
double fb;				/* f(b)				*/
double fc;				/* f(c)				*/

a = ax;  b = bx;  fa = (f)(a, m, beta, lambda);  fb = (f)(b, m, beta, lambda);
c = a;   fc = fa;

for(;;)		/* Main iteration loop	*/
{
  double prev_step = b-a;		/* Distance from the last but one*/
/* to the last approximation	*/
double tol_act;			/* Actual tolerance		*/
double p;      			/* Interpolation step is calcu- */
double q;      			/* lated in the form p/q; divi- */
/* sion operations is delayed   */
/* until the last moment	*/
double new_step;      		/* Step at this iteration       */

if( fabs(fc) < fabs(fb) )
{                         		/* Swap data for b to be the 	*/
a = b;  b = c;  c = a;          /* best approximation		*/
fa=fb;  fb=fc;  fc=fa;
}
tol_act = 2.0*DOUBLE_EPS*fabs(b) + tol/2.0;
new_step = (c-b)/2.0;

if( fabs(new_step) <= tol_act || fb == (double)0 )
  return b;				/* Acceptable approx. is found	*/

/* Decide if the interpolation can be tried	*/
if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
&& fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
{					/* Interpolatiom may be tried	*/
register double t1,cb,t2;
  cb = c-b;
  if( a==c )			/* If we have only two distinct	*/
  {				/* points linear interpolation 	*/
t1 = fb/fa;			/* can only be applied		*/
p = cb*t1;
q = 1.0 - t1;
  }
  else				/* Quadric inverse interpolation*/
  {
    q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
    p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
    q = (q-1.0) * (t1-1.0) * (t2-1.0);
  }
  if( p>(double)0 )		/* p was calculated with the op-*/
q = -q;			/* posite sign; make p positive	*/
else				/* and assign possible minus to	*/
p = -p;			/* q				*/

if( p < (0.75*cb*q-fabs(tol_act*q)/2.0)	/* If b+p/q falls in [b,c]*/
&& p < fabs(prev_step*q/2.0) )	/* and isn't too large	*/
new_step = p/q;			/* it is accepted	*/
/* If p/q is too large then the	*/
/* bissection procedure can 	*/
/* reduce [b,c] range to more	*/
/* extent			*/
}

if( fabs(new_step) < tol_act ) {	/* Adjust the step to be not less*/
if( new_step > (double)0 )	/* than tolerance		*/
new_step = tol_act;
else
  new_step = -tol_act;
}

a = b;  fa = fb;			/* Save the previous approx.	*/
b += new_step;  fb = (*f)(b, m, beta, lambda);  /* Do step to a new approxim. */
if( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) )
{                 			/* Adjust c for it to have a sign*/
c = a;  fc = fa;                  /* opposite to that of b	*/
}
}

}

/* Random sample generator from a Generalized Inverse Gaussian distribution (modified version of rgig.c using Rcpp classes)*/
NumericVector Rgig(const int n, const double lambda, const double chi, const double psi)
{
  NumericVector samps(n);
  /* special case which is basically a gamma distribution */
  if((chi < ZTOL) & (lambda > 0.0)) {
    int i;
    for(i=0; i<n; i++) samps(i) = R::rgamma(lambda, 2.0/psi);
    return samps;
  }

  /* special cases which is basically an inverse gamma distribution */
  if((psi < ZTOL) & (lambda < 0.0)) {
    int i;
    for(i=0; i<n; i++) samps(i) = 1.0/R::rgamma(0.0-lambda, 2.0/chi);
    return samps;
  }

  /* special case which is basically an inverse gaussian distribution */
  if(lambda == -0.5) {
    double alpha;
    int i;
    alpha = sqrt(chi/psi);
    for(i=0; i<n; i++) samps(i) = rinvgauss(alpha, chi);
    return samps;
  }

  /*
  * begin general purpose rgig code, which was basically
  * translated from the R function rgig in the ghyp package v_1.5.2
  */

  double alpha, beta, beta2, m, m1, lm1, lm12, upper, yM, yP, a, b, c, R1, R2, Y;
  int i, need;

  alpha = sqrt(chi/psi);
  beta2 = psi*chi;
  beta = sqrt(psi*chi);
  lm1 = lambda - 1.0;
  lm12 = lm1*lm1;
  m = (lm1 + sqrt(lm12 + beta2))/beta;
  m1 = m + 1.0/m;

  upper = m;
  while (gig_y_gfn(upper, m, beta, lambda) <= 0) { upper *= 2.0; }

  yM = zeroin_gig(0.0, m, gig_y_gfn, ZTOL, m, beta, lambda);
  yP = zeroin_gig(m, upper, gig_y_gfn, ZTOL, m, beta, lambda);

  a = (yP - m) * pow(yP/m, 0.5 * lm1);
  a *=  exp(-0.25 * beta * (yP + 1.0/yP - m1));
  b = (yM - m) * pow(yM/m, 0.5 * lm1);
  b *= exp(-0.25 * beta * (yM + 1/yM - m1));
  c = -0.25 * beta * m1 + 0.5 * lm1 * log(m);

  for (i=0; i<n; i++) {
    need = 1;
    while (need) {
      R1 = unif_rand();
      R2 = unif_rand();

      Y = m + a * R2/R1 + b * (1.0 - R2)/R1;
      if (Y > 0.0) {
        if (-log(R1) >= - 0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0/Y) + c) {
          need = 0;
        }
      }
    }
    samps[i] = Y*alpha;
  }
  return(samps);
}

/* Random sample generator from a Generalized Inverse Gaussian distribution (modified version of rgig.c using Rcpp classes)*/
double RgigDouble(const double lambda, const double chi, const double psi)
{
  double samps;
  /* special case which is basically a gamma distribution */
  if((chi < ZTOL) & (lambda > 0.0)) {
    samps = R::rgamma(lambda, 2.0/psi);
    return samps;
  }

  /* special cases which is basically an inverse gamma distribution */
  if((psi < ZTOL) & (lambda < 0.0)) {
    samps = 1.0/R::rgamma(0.0-lambda, 2.0/chi);
    return samps;
  }

  /* special case which is basically an inverse gaussian distribution */
  if(lambda == -0.5) {
    double alpha;
    alpha = sqrt(chi/psi);
    samps = rinvgauss(alpha, chi);
    return samps;
  }

  /*
  * begin general purpose rgig code, which was basically
  * translated from the R function rgig in the ghyp package v_1.5.2
  */

  double alpha, beta, beta2, m, m1, lm1, lm12, upper, yM, yP, a, b, c, R1, R2, Y;
  int need;

  alpha = sqrt(chi/psi);
  beta2 = psi*chi;
  beta = sqrt(psi*chi);
  lm1 = lambda - 1.0;
  lm12 = lm1*lm1;
  m = (lm1 + sqrt(lm12 + beta2))/beta;
  m1 = m + 1.0/m;

  upper = m;
  while (gig_y_gfn(upper, m, beta, lambda) <= 0) { upper *= 2.0; }

  yM = zeroin_gig(0.0, m, gig_y_gfn, ZTOL, m, beta, lambda);
  yP = zeroin_gig(m, upper, gig_y_gfn, ZTOL, m, beta, lambda);

  a = (yP - m) * pow(yP/m, 0.5 * lm1);
  a *=  exp(-0.25 * beta * (yP + 1.0/yP - m1));
  b = (yM - m) * pow(yM/m, 0.5 * lm1);
  b *= exp(-0.25 * beta * (yM + 1/yM - m1));
  c = -0.25 * beta * m1 + 0.5 * lm1 * log(m);

  need = 1;
  while (need) {
    R1 = unif_rand();
    R2 = unif_rand();

    Y = m + a * R2/R1 + b * (1.0 - R2)/R1;
    if (Y > 0.0) {
      if (-log(R1) >= - 0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0/Y) + c) {
        need = 0;
      }
    }
  }
  samps = Y*alpha;
  return(samps);
}

/* END OF BORROWED CODE */

/* Auxiliary function that converts Rcpp::NumericVector elements into arma::vec elements */
arma::vec as_arma(
    NumericVector& x) /* Vector to be converted into an arma::vec class */
{
  return arma::vec(x.begin(), x.length(), false);
}

/* Auxiliary function that converts Rcpp::NumericMatrix elements into arma::mac elements */
arma::mat as_arma(NumericMatrix& x) /* Matrix to be converted into an arma::mat class */
{
  return arma::mat(x.begin(), x.nrow(), x.ncol(), false);
}

/* Auxiliary function: logarithm of a Gamma function applied element-wise to an arma::mat element */
arma::mat lgamma_cpp(arma::mat const& x)
{
  /*
  * NumericMatrix x: matrix to which the logarithm of the Gamma function will be applied element-wise
  */
  arma::mat output = x;
  for (unsigned int i=0; i<arma::size(x,0); i++)
  {
    for(unsigned int j=0; j<arma::size(x,1); j++)
    {
      output(i,j) = R::lgammafn(x(i,j));
    }
  }
  return output;
}

/* Auxiliary function: logarithm of a Gamma function applied element-wise to an arma::mat element */
arma::vec lgamma_cpp_vec(arma::vec const& x)
{
  /*
  * NumericMatrix x: matrix to which the logarithm of the Gamma function will be applied element-wise
  */
  arma::vec output = x;
  for (unsigned int i=0; i < output.size(); i++)
  {
    output(i) = R::lgammafn(x(i));
  }
  return output;
}

/* Auxiliary function: to avoid numerical overflows/underflows when computing a sum of exponentiated values */
double log_sum_exp_cpp(arma::vec const& x_arma)
{
  double offset;
  if ( max(abs(x_arma)) > max(x_arma) ) { offset = min(x_arma);}
  else { offset = max(x_arma);}
  return log(sum(exp(x_arma - offset))) + offset;
}



// PRIOR FOR THE COEFFICIENT OF VARIATION
// Debug: function tested to match original R implementation
double Prior(double const& x,
             double const& a,
             double const& trunc,
             double const& lower,
             double const& upper,
             String const& type) /* 1: "TruncExp"; 2: "Pareto" */
{
  double aux;

  if(type == "TruncExp")
  {
    if(trunc == 1)
    {
      if(x <= lower){ aux = 0; }
      else{ aux = a*exp(-a*(x-lower)); }
    }
    if(trunc==2)
    {
      if(x<=lower || x>=upper){ aux = 0;}
      else{ aux = a * exp(-a*(x-lower)) / (1 - exp(-a*(upper-lower))); }
    }
  }
  if(type == "Pareto")
  {
    if(trunc == 1){ aux = a * pow(lower,a) * pow(x,-(a+1)); }
    if(trunc == 2){ aux = (a * pow(lower,a) * pow(x,-(a+1))) / (1 - pow(lower/upper,a)); }
  }
  return aux;
}


// PRIOR FOR THE COEFFICIENT OF VARIATION (log-scale)
// Debug: function tested to match original R implementation
double LogPrior(double const& x,
                double const& a,
                double const& trunc,
                double const& lower,
                double const& upper,
                String const& type) /* 1: "TruncExp"; 2: "Pareto" */
{
  double aux;

  if(type == "TruncExp")
  {
    if(trunc == 1)
    {
      if(x <= lower){ aux = R_NegInf; }
      else{ aux = log(a) - a*(x-lower); }
    }
    if(trunc == 2)
    {
      if(x <= lower || x >= upper){ aux = R_NegInf;}
      else{ aux = log(a) - a * (x-lower) - log(1-exp(-a*(upper-lower))); }
    }
  }
  if(type == "Pareto")
  {
    if(trunc == 1){ aux = log(a) + a * log(lower) - (a+1) * log(x); }
    if(trunc == 2){ aux = log(a) + a * log(lower) - (a+1) * log(x) - log(1 - pow(lower/upper,a));}
  }
  return aux;
}

// CV AS A FUNCTION OF GAMMA AND THETA
// Debug: function tested to match original R implementation
arma::vec cvStarSquare(double const& gam,
                       double const& theta,
                       String const& mixing)
{
  // 1st element: cv.star2; 2nd element: dcv.star2
  arma::vec aux = arma::ones(2);

  if(mixing == "Gamma")
  {
    aux(0) = exp( lgamma(theta) + lgamma(theta-2/gam) - 2 * lgamma(theta-1/gam) ) - 1;
    aux(1) = exp( lgamma(theta) + lgamma(theta-2/gam) - 2 * lgamma(theta-1/gam)) * (Rf_digamma(theta) + Rf_digamma(theta-2/gam) - 2 * Rf_digamma(theta-1/gam));
  }
  if(mixing == "InvGamma")
  {
    aux(0) = exp( lgamma(theta) + lgamma(theta+2/gam) - 2 * lgamma(theta+1/gam)) - 1;
    aux(1) = exp( lgamma(theta) + lgamma(theta+2/gam) - 2 * lgamma(theta+1/gam)) * (Rf_digamma(theta) + Rf_digamma(theta+2/gam) - 2 * Rf_digamma(theta+1/gam));
  }
  if(mixing=="InvGauss")
  {

    //Using the asymptotic expansion as theta->0
    if(theta<0.0015) {aux(0) = 0; aux(1) = 0;}
    else
    {
      double K1 = R::bessel_k(1/theta, -((1/gam) + 0.5), 1);
      double K2 = R::bessel_k(1/theta, -((2/gam) + 0.5), 1);
      aux(0) = sqrt(M_PI/2) * (exp(-1/theta) * pow(theta, 1.0/2.0)) * (K2/K1) * (1/K1) - 1;
      aux(1) = sqrt(M_PI/2) * pow(theta,-3.0/2.0) * exp(-1/theta) / K1;
      aux(1) *= (K2/K1) + R::bessel_k(1/theta, -((2/gam)-0.5), 1) / K1 - 2*(K2/K1)*(1/K1)*R::bessel_k(1/theta,-((1/gam)-0.5),1);
    }
  }
  if(mixing == "LogNormal")
  {
    aux(0) = exp(theta * pow(gam, -2)) - 1;
    aux(1) = exp(theta * pow(gam, -2)) * pow(gam, -2);
  }

  return aux;
}

// PRIOR FOR THETA GIVEN GAMMA
// Debug: function tested to match original R implementation
double PriorTheta(double const& theta,
                  double const& gam,
                  double const& a,
                  String const& type,
                  String const& mixing)
{
  // Auxiliary variables
  double out;
  double lower; double upper; double trunc;
  double aux0; double aux1;
  double cv; double dcv;

  if(mixing == "None"){ out = 1; }
  else
  {
    lower = sqrt( exp( lgamma(1 + (2/gam)) - 2 * lgamma(1 + (1/gam)) ) - 1);

    if(lower == R_PosInf) { out = 0; }
    else
    {
      if(mixing == "Gamma") { trunc = 1; upper = R_PosInf; }
      if(mixing == "InvGamma")
      {
        trunc = 2;
        upper = sqrt(exp(2*lgamma(1 + (2/gam)) - 4*lgamma(1 + (1/gam))) - 1);
        if(upper == R_PosInf) { trunc = 1; }
      }
      if(mixing == "InvGauss")
      {
        trunc = 2;
        upper = sqrt(sqrt(M_PI) * exp(lgamma(1 + (2/gam))- 2 * lgamma(1 + (1/gam)) + lgamma((2/gam) + 0.5) - 2*lgamma((1/gam) + 0.5))-1);
        if(upper == R_PosInf){ trunc = 1; }
      }
      if(mixing == "LogNormal")	{ trunc = 1; upper = R_PosInf;}

      aux0 = (lgamma(1 + (2/gam)) - 2 * lgamma(1 + (1/gam)));
      aux1 = exp(aux0) - 1;
      arma::vec CVS = cvStarSquare(gam, theta, mixing);
      cv = sqrt(exp(aux0 + log(CVS(0))) + aux1);
      dcv = std::abs((exp(aux0)/(2*cv)) * CVS(1));

      if(R_IsNA(cv) == TRUE){ out = NA_REAL; }
      else
      {
        out = Prior(cv, a, trunc, lower, upper, type) * dcv;
      }
    }
  }
  return out;
}

// PRIOR FOR THETA GIVEN GAMMA (LOGARITHMIC SCALE)
// Debug: function tested to match original R implementation
// [[Rcpp::export]]
double HiddenLogPriorTheta(double const& theta,
                     double const& gam,
                     double const& a,
                     String const& type,
                     String const& mixing)
{
  double out;
  double lower; double upper; double trunc;
  double aux0; double aux1;
  double cv; double dcv;

  if(mixing == "None"){ out = 1; }
  else
  {
    lower = sqrt(exp(lgamma(1 + (2/gam)) - 2*lgamma(1 + (1/gam))) - 1);

    if(lower == R_PosInf) { out = 0; }
    else
    {
      if(mixing == "Gamma") { trunc = 1; upper = R_PosInf; }
      if(mixing == "InvGamma")
      {
        trunc = 2;
        upper = sqrt(exp(2*lgamma(1 + (2/gam)) - 4*lgamma(1 + (1/gam))) - 1);
        if(upper == R_PosInf) { trunc = 1; }
      }
      if(mixing == "InvGauss")
      {
        trunc = 2;
        upper = sqrt(sqrt(M_PI) * exp(lgamma(1 + (2/gam))- 2 * lgamma(1 + (1/gam)) + lgamma((2/gam) + 0.5) - 2*lgamma((1/gam) + 0.5))-1);
        if(upper == R_PosInf){ trunc = 1; }
      }
      if(mixing == "LogNormal")	{ trunc = 1; upper = R_PosInf;}

      aux0 = (lgamma(1 + (2/gam)) - 2 * lgamma(1 + (1/gam)));
      aux1 = exp(aux0) - 1;
      arma::vec CVS = cvStarSquare(gam, theta, mixing);
      cv = sqrt(exp(aux0 + log(CVS(0))) + aux1);
      dcv = std::abs((exp(aux0)/(2*cv)) * CVS(1));

      if(R_IsNA(cv) == TRUE){ out = NA_REAL; }
      else
      {
        out = LogPrior(cv, a, trunc, lower, upper, type) + log(dcv);
      }
    }
  }
  return out;
}

arma::mat RemoveColMat(arma::mat const& X,
                       int const& j)
{
  arma::mat aux;
  int k = X.n_cols;

  if(j == 0) {aux = X.cols(1, k-1);}
  else
  {
    if(j == k-1) {aux = X.cols(0, k-2);}
    else {aux = join_rows(X.cols(0, j-1), X.cols(j+1, k-1));}
  }
  return aux;
}

arma::vec RemoveElemVec(arma::vec const& beta,
                        int const& j)
{
  using arma::span;

  arma::vec aux;
  int k = beta.n_elem;

  if(j == 0) {aux = beta(span(1, k-1));}
  else
  {
    if(j == k-1) {aux = beta(span(0, k-2));}
    else {aux = join_cols(beta(span(0, j-1)), beta(span(j+1, k-1)));}
  }
  return aux;
}

// ################################################################################
// ################################################################################
// # FUNCTIONS REQUIRED TO RUN THE MCMC SAMPLER
// ################################################################################
// ################################################################################

// GAUSSIAN RANDOM WALK METROPOLIS HASTINGS FOR beta[j]
// USE lambda=rep(1,times=n) FOR WEIBULL MODEL WITH MO MIXTURE
// Debug: function tested to match original R implementation
arma::vec betaUpdate(double const& omega2,
                     int const& j,
                     arma::vec const& beta0,
                     arma::vec const& Time,
                     arma::vec const& Event,
                     arma::mat const& X,
                     double const& gam,
                     arma::vec const& lambda,
                     int const& k,
                     int const& n)
{
  // INITIALIZING THE VECTOR WHERE DRAWS ARE STORED
  arma::vec out(2);

  // AUXILIARY QUANTITIES
  arma::mat Xj = RemoveColMat(X, j); arma::vec bj = RemoveElemVec(beta0, j);
  arma::vec SumOtherEffects = exp( - gam * (Xj * bj) );

  // PROPOSAL STEP
  double y = R::rnorm(0,1) * sqrt(omega2) + beta0(j);
  double u = R::runif(0,1);
  // ACCEPTANCE STEP
  double log_aux = (beta0(j) - y) * gam * sum(X.col(j) % Event);
  log_aux += sum( SumOtherEffects % pow(Time,gam) % lambda % (exp(-gam * beta0(j) * X.col(j)) - exp(-gam * y * X.col(j))));
  if(log(u) < log_aux) { out(0) = y; out(1) = 1;}
  if(log(u) >= log_aux) { out(0) = beta0(j); out(1) = 0;}

  // OUTPUT
  return out;
}

// GAUSSIAN RANDOM WALK METROPOLIS HASTINGS FOR GAMMA
// USE lambda=rep(1,times=n) FOR WEIBULL MODEL WITH MIXTURE
// Debug: function tested to match original R implementation
// Debug: includes GRWMH.RMWLN.gam function as a special case (rather than a separate function)
arma::vec gamUpdate(double const& omega2,
                    double const& gam0,
                    arma::vec const& Time,
                    arma::vec const& Event,
                    arma::mat const& X,
                    arma::vec const& beta,
                    double const& theta,
                    arma::vec const& lambda,
                    String const& PriorCV,
                    double const& HypTheta,
                    double const& Hyp1Gam,
                    double const& Hyp2Gam,
                    String const& mixing,
                    double const& lower_bound,
                    int FIX_THETA = 0)
{
  // INITIALIZATION OF OUTPUT VARIABLES
  arma::vec out(2);

  // PROPOSAL STEP
  double y = R::rnorm(0,1) * sqrt(omega2) + gam0;
  double u = R::runif(0,1);

  // VERIFYING RESTRICTION IMPOSED FOR NUMERICAL ACCURACY (REJEC IF NOT SATISFIED)
  if(y <= lower_bound) { out(0) = gam0; out(1) = 0;}
  else
  {
    // VERIFYING OTHER RESTRICTIONS IN THE PARAMETER SPACE (REJEC IF NOT SATISFIED)
    if((mixing == "Gamma") & (y <= 2/theta)) {out(0) = gam0; out(1) = 0;}
    // ACCEPTANCE STEP
    else
    {
      double log_aux = sum(Event) * log(y / gam0);
      log_aux += (y-gam0) * sum( Event % (log(Time) - X * beta));
      log_aux -= sum(lambda % (pow(exp(-X*beta) % Time, y) - pow(exp(-X*beta) % Time, gam0)));
      if(FIX_THETA == 0)
      {
        log_aux += HiddenLogPriorTheta(theta, y, HypTheta, PriorCV, mixing);
        log_aux -= HiddenLogPriorTheta(theta, gam0, HypTheta, PriorCV, mixing);
      }
      log_aux += R::dgamma(y, Hyp1Gam, Hyp2Gam, 1);
      log_aux -= R::dgamma(gam0, Hyp1Gam, Hyp2Gam, 1);
      if(mixing == "LogNormal")
      {
        log_aux -= Time.n_elem * log(y / gam0);
        log_aux -= (pow(y, -2.0) - pow(gam0, -2.0)) * (0.5/theta) * sum(pow(log(lambda), 2.0));
      }
      if(log(u) < log_aux) { out(0) = y; out(1) = 1; }
      else { out(0) = gam0; out(1) = 0; }
    }
  }

  // OUTPUT
  return(out);
}

// GAUSSIAN RANDOM WALK METROPOLIS STEP FOR THETA
// USE lambda=rep(1,times=n) FOR WEIBULL MODEL WITH MIXTURE
// Debug: function tested to match original R implementation
// Debug: includes GRWMH.RMWLN.gam function as a special case (rather than a separate function)
arma::vec thetaUpdate(double const& omega2,
                      double const& theta0,
                      double const& gam,
                      arma::vec const& lambda,
                      String const& PriorCV,
                      double const& HypTheta,
                      String const& mixing,
                      double const& n)
{
  // INITIALIZATION OF OUTPUT VARIABLES
  arma::vec out(2);

  // PROPOSAL STEP
  double y = R::rnorm(0,1) * sqrt(omega2) + theta0;
  double u = R::runif(0,1);

  double log_aux;
  if(mixing == "Gamma")
  {
    if(y <= 2/gam) {out(0) = theta0; out(1) = 0;}
    else
    {
      log_aux = n * (y*log(y) - R::lgammafn(y) - theta0*log(theta0) + R::lgammafn(theta0));
      log_aux += (y-theta0) * sum(log(lambda)) - (y-theta0) * sum(lambda);
    }
  }
  if(mixing == "InvGamma")
  {
    if(y <= 1) {out(0) = theta0; out(1) = 0;}
    else
    {
      log_aux = n * (R::lgammafn(theta0) - R::lgammafn(y)) + (theta0-y) * sum(log(lambda));
    }
  }
  if(mixing == "InvGauss")
  {
    log_aux = 0.5 * (pow(theta0,-2) - pow(y,-2)) * sum(lambda);
    log_aux -= n * (pow(theta0,-1) - pow(y,-1));
  }
  if(mixing == "LogNormal")
  {
    log_aux = (n/2) * log(theta0/y);
    log_aux -= 0.5 * (pow(y,-1) - pow(theta0,-1)) * sum(pow(log(lambda),2));
  }
  // PRIOR FACTOR
  log_aux += HiddenLogPriorTheta(y, gam, HypTheta, PriorCV, mixing);
  log_aux -= HiddenLogPriorTheta(theta0, gam, HypTheta, PriorCV, mixing);

  // ACCEPT/REJECT
  if((arma::is_finite(log_aux)) & (log(u) < log_aux)) { out(0) = y; out(1) = 1; }
  else { out(0) = theta0; out(1) = 0; }

  // REJECTING INVALID VALUES
  if((mixing == "Gamma") & (y <= 2/gam)) { out(0) = theta0; out(1) = 0; }
  if((mixing == "InvGamma") & (y <= 1)) { out(0) = theta0; out(1) = 0; }

  // OUTPUT
  return(out);
}

// Sampling of mixing parameters
// Debug: function tested to match original R implementation
arma::vec lambdaUpdate(String const& mixing,
                       arma::vec const& Time,
                       arma::vec const& Event,
                       arma::mat const& X,
                       arma::vec const& beta,
                       double const& gam,
                       double const& theta,
                       int const& n)
{
  arma::vec lambda = arma::ones(n);
  arma::vec ScaleAux;

  if(mixing == "Exponential")
  {
    ScaleAux = 1 / (1 + pow(exp(-X*beta) % Time, gam));
    for(int i=0; i < n; i++) { lambda(i) = R::rgamma(1 + Event(i), ScaleAux(i) );}
  }
  if(mixing == "Gamma")
  {
    ScaleAux = 1 / (theta + pow(exp(-X*beta) % Time, gam));
    for(int i=0; i < n; i++) { lambda(i) = R::rgamma(theta + Event(i), ScaleAux(i) );}
  }
  if(mixing == "InvGamma")
  {
    ScaleAux = 2 * pow(exp(-X*beta) % Time, gam);
    for(int i=0; i < n; i++) { lambda(i) = RgigDouble(-theta + Event(i), 2, ScaleAux(i));}
  }
  if(mixing == "InvGauss")
  {
    ScaleAux = 2 * pow(exp(-X*beta) % Time, gam) + pow(theta, -2);
    for(int i=0; i < n; i++) { lambda(i) = RgigDouble(Event(i) - 0.5, 1, ScaleAux(i));}
  }

  return lambda;
}

// Sampling of mixing parameters (log-normal mixing only)
arma::mat lambdaUpdate_LN(arma::vec const& lambda0,
                          arma::vec const& prop_var,
                          arma::vec const& Time,
                          arma::vec const& Event,
                          arma::mat const& X,
                          arma::vec const& beta,
                          double const& gam,
                          double const& theta,
                          int const& n)
{
  // CREATING VARIABLES WHERE TO STORE DRAWS
  arma::vec lambda = arma::zeros(n);
  arma::vec ind = arma::zeros(n);

  // PROPOSAL STEP
  arma::vec y_aux = arma::randn(n) % sqrt(prop_var) + lambda0;
  arma::vec y = abs(y_aux);
  arma::vec u = arma::randu(n);

  // ACCEPT/REJECT STEP
  arma::vec log_aux = Event % log(y/lambda0); // -1 cancels out as proposal is in log-scale // (Event-1)*log(y/lambda[l,])
  log_aux -= pow(exp(-X*beta) % Time, gam) % (y - lambda0);
  log_aux -= (1/(2*theta)) * ( pow(log(y),2)- pow(log(lambda0),2) );

  // CREATING OUTPUT VARIABLE & DEBUG
  // DEBUG: Reject very small values to avoid numerical issues
  for (int i=0; i < n; i++)
  {
    if(arma::is_finite(log_aux(i)) & (log(u(i)) < log_aux(i)) & (y_aux(i) > 1e-3))
    {
      ind(i) = 1; lambda(i) = y(i);
    }
    else{ind(i) = 0; lambda(i) = lambda0(i);}
  }

  // OUTPUT
  return join_rows(lambda, ind);
}

// MCMC ALGORITHM
// Argument 'Adapt' added to replace 'MCMC.WEI.NonAdapt' function
// Argument 'FixGam' added to replace 'MCMCR.WEI.gam' function
// Argument 'FixBetaJ' added to replace 'MCMCR.WEI.betaJ.gam' function
// [[Rcpp::export]]
Rcpp::List HiddenRMWreg_MCMC(int N, // Total number of MCMC draws
                             int thin, // Thinning period for MCMC chain
                             int burn, // Burning period for MCMC chain
                             NumericVector Time,
                             NumericVector Event,
                             NumericMatrix X,
                             String mixing,
                             double const& Hyp1Gam,
                             double const& Hyp2Gam,
                             String const& PriorCV,
                             double const& HypTheta,
                             NumericVector beta0,
                             double gam0,
                             double theta0,
                             int Adapt,
                             double ar,
                             int StoreAdapt,
                             int EndAdapt,
                             NumericVector LSbeta0,
                             double LSgam0,
                             double LStheta0,
                             int FixBetaJ, // If FixBeta = 0, all Betas are updated. Else beta[1],...,beta[J] are fixed
                             int FixGam,
                             int FixTheta,
                             int PrintProgress,
                             int lambdaPeriod,
                             int FixLambdaI, // If FixLambda = 0, all lambda's are updated.
                             double RefLambda)
{
  // NUMBER OF REGRESSORS, SAMPLE SIZE AND NUMBER OF STORED DRAWS
  int k = beta0.size(); int n = Time.size(); int Naux = N/thin - burn/thin;

  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat beta = arma::ones(Naux, k);
  arma::vec gam = arma::ones(Naux);
  arma::mat lambda = arma::ones(Naux, n);
  arma::vec theta = arma::ones(Naux);
  arma::mat LSbeta; arma::vec LSgam; arma::vec LStheta;

  // LOG-PROPOSAL VARIANCES
  if(StoreAdapt == 1)
  {
    LSbeta = arma::zeros(Naux, k);
    LSgam = arma::zeros(Naux);
    LStheta = arma::zeros(Naux);
  }

  // TRANSFORMATION ONTO ARMA ELEMENTS
  arma::vec Time_arma = as_arma(Time);
  arma::vec Event_arma = as_arma(Event);
  arma::mat X_arma = as_arma(X);

  // INITIALIZATION OF PARAMETER VALUES FOR MCMC RUN
  arma::mat betaAux = arma::ones(k,2); betaAux.col(0) = as_arma(beta0);
  arma::vec gamAux = arma::ones(2); gamAux(0) = gam0;
  arma::vec thetaAux = arma::ones(2); thetaAux(0) = theta0;
  arma::vec lambdaAux = arma::ones(n);

  // INITIALIZATION OF ADAPTIVE VARIANCES
  arma::vec LSbetaAux = as_arma(LSbeta0);
  double LSgamAux = LSgam0 + 1 - 1; // +-1 to avoid rewriting LSgam0
  double LSthetaAux = LStheta0 + 1 - 1; // +-1 to avoid rewriting LStheta0

  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec betaAccept = arma::zeros(k); arma::vec PbetaAux = arma::zeros(k);
  double gamAccept = 0; double PgamAux = 0;
  double thetaAccept = 0; double PthetaAux = 0;

  // EXTRA PARAMS FOR LOG-NORMAL MIXING ONLY
  arma::vec LSlambdaAux = arma::zeros(n);
  arma::vec lambdaAccept = arma::zeros(n);
  arma::vec PlambdaAux = arma::zeros(n);
  arma::mat lambdaAux_LN = arma::ones(n, 2);

  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i; int j;

  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;

  // START OF MCMC LOOP
  for (i=0; i < N; i++)
  {
    Rcpp::checkUserInterrupt();

    if(i==burn)
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
      Rcpp::Rcout << "End of burn-in period."<< std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
    }

    Ibatch++;

    // UPDATE OF BETA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    for(j=0; j < k; j++)
    {
      if(j > FixBetaJ - 1)
      {
        betaAux.row(j) = betaUpdate(exp(LSbetaAux(j)), j, betaAux.col(0),
                                    Time_arma, Event_arma, X_arma,
                                    gamAux(0), lambdaAux, k, n).t();
      }
      else { betaAux(j,1) = 0; }
    }
    PbetaAux += betaAux.col(1); if(i>=burn) {betaAccept += betaAux.col(1);}

    if(FixGam == 0)
    {
      // UPDATE OF GAM: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
      gamAux = gamUpdate(exp(LSgamAux), gamAux(0),
                         Time_arma, Event_arma, X_arma,
                         betaAux.col(0), thetaAux(0), lambdaAux,
                         PriorCV, HypTheta, Hyp1Gam, Hyp2Gam,
                         mixing, 0.06, FixTheta); // lower_bound
      PgamAux += gamAux(1); if(i>=burn) {gamAccept += gamAux(1);}
    }

    if((FixTheta == 0) & (mixing != "None") & (mixing != "Exponential"))
    {
      thetaAux = thetaUpdate(exp(LSthetaAux), thetaAux(0),
                             gamAux(0), lambdaAux,
                             PriorCV, HypTheta, mixing, n);
      PthetaAux += thetaAux(1); if(i>=burn) {thetaAccept += thetaAux(1);}
    }

    // Updating mixing parameters
    if(i%lambdaPeriod==0)
    {
      if((mixing != "None") & (mixing != "LogNormal"))
      {

        lambdaAux = lambdaUpdate(mixing,
                                 Time_arma, Event_arma, X_arma,
                                 betaAux.col(0), gamAux(0), thetaAux(0), n);
      }
      if(mixing == "LogNormal")
      {
        lambdaAux_LN = lambdaUpdate_LN(lambdaAux, exp(LSlambdaAux),
                                       Time_arma, Event_arma, X_arma,
                                       betaAux.col(0), gamAux(0), thetaAux(0), n);
        lambdaAux = lambdaAux_LN.col(0);
        PlambdaAux += lambdaAux_LN.col(1); if(i>=burn) {lambdaAccept += lambdaAux_LN.col(1);}
      }
    }
    if(FixLambdaI != 0) { lambdaAux(FixLambdaI - 1) = RefLambda; }

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if((i < EndAdapt) & (Adapt == 1))
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PbetaAux /= 50; PbetaAux = -1+2*arma::conv_to<arma::mat>::from(PbetaAux>ar);
        LSbetaAux += PbetaAux*std::min(0.01,1/sqrt(i));
        PbetaAux = arma::zeros(k);

        if(FixGam == 0)
        {
          PgamAux /= 50; PgamAux = -1+2*(PgamAux>ar);
          LSgamAux += PgamAux*std::min(0.01,1/sqrt(i));
          PgamAux = 0;
        }

        if((FixTheta == 0) & (mixing != "None") & (mixing != "Exponential"))
        {
          PthetaAux /= 50; PthetaAux = -1+2*(PthetaAux>ar);
          LSthetaAux += PthetaAux*std::min(0.01,1/sqrt(i));
          PthetaAux = 0;
        }

        if(mixing == "LogNormal")
        {
          PlambdaAux /= (50 / lambdaPeriod); PlambdaAux = -1+2*arma::conv_to<arma::mat>::from(PlambdaAux>ar);
          LSlambdaAux += PlambdaAux*std::min(0.01,1/sqrt(i));
          PlambdaAux = arma::zeros(n);
        }


        Ibatch = 0;
      }
    }

    // STORAGE OF DRAWS
    if((i%thin==0) & (i>=burn))
    {
      beta.row(i/thin - burn/thin) = betaAux.col(0).t();
      gam(i/thin - burn/thin) = gamAux(0);
      theta(i/thin - burn/thin) = thetaAux(0);
      lambda.row(i/thin - burn/thin) = lambdaAux.t();

      if(StoreAdapt == 1)
      {
        LSbeta.row(i/thin - burn/thin) = LSbetaAux.t();
        LSgam(i/thin - burn/thin) = LSgamAux;
        LStheta(i/thin - burn/thin) = LSthetaAux;
      }
    }

    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if((i%(2*thin) == 0) & (PrintProgress == 1))
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
      Rcpp::Rcout << "MCMC iteration " << i << " out of " << N << " has been completed." << std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
      Rcpp::Rcout << "Current draws of some selected parameters are displayed below." << std::endl;
      Rcpp::Rcout << "beta (1st element): " << betaAux(0,0) << std::endl;
      Rcpp::Rcout << "gam: " << gamAux(0) << std::endl;
      Rcpp::Rcout << "theta: " << thetaAux(0) << std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
      Rcpp::Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
      Rcpp::Rcout << "LSbeta (1st element): " << LSbetaAux(0) << std::endl;
      Rcpp::Rcout << "LSgam: " << LSgamAux << std::endl;
      Rcpp::Rcout << "LStheta: " << LSthetaAux << std::endl;
    }
  }

  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "All " << N << " MCMC iterations have been completed." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "Please see below a summary of the overall acceptance rates." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;

  Rcpp::Rcout << "Minimum acceptance rate among beta[j]'s: " << min(betaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among beta[j]'s: " << mean(betaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among beta[j]'s: " << max(betaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Acceptance rate for gam: " << gamAccept/(N-burn) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Acceptance rate for theta: " << thetaAccept/(N-burn) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  if(mixing == "LogNormal")
  {
    Rcpp::Rcout << "Minimum acceptance rate among lambda[i]'s: " << min(lambdaAccept/(N-burn)) << std::endl;
    Rcpp::Rcout << "Average acceptance rate among lambda[i]'s: " << mean(lambdaAccept/(N-burn)) << std::endl;
    Rcpp::Rcout << "Maximum acceptance rate among lambda[i]'s: " << max(lambdaAccept/(N-burn)) << std::endl;
    Rcpp::Rcout << " " << std::endl;
  }
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;

  if(StoreAdapt == 1)
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("beta") = beta,
           Rcpp::Named("gam") = gam,
           Rcpp::Named("theta") = theta,
           Rcpp::Named("lambda") = lambda,
           Rcpp::Named("ls.beta") = LSbeta,
           Rcpp::Named("ls.gam") = LSgam,
           Rcpp::Named("ls.theta") = LStheta));
  }
  else
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("beta") = beta,
           Rcpp::Named("gam") = gam,
           Rcpp::Named("theta") = theta,
           Rcpp::Named("lambda") = lambda));
  }
}

// ################################################################################
// ################################################################################
// # FUNCTIONS REQUIRED TO CALCULATE DIC
// ################################################################################
// ################################################################################

arma::vec logf_RMEIGAM(arma::vec const& Time,
                       arma::vec const& Rate,
                       double const& theta)
{
  arma::vec aux = log(2) + ((theta+1)/2)*log(Rate) + ((theta-1)/2)*log(Time) - R::lgammafn(theta);
  for (int i=0; i < Time.n_elem; i++)
  {
    aux(i) += log(R::bessel_k(2 * sqrt(Rate(i)*Time(i)), -(theta-1), 1));
  }
  return(aux);
}

arma::vec logS_RMEIGAM(arma::vec const& Time,
                       arma::vec const& Rate,
                       double const& theta)
{
  arma::vec aux = log(2) + (theta/2) * log(Rate % Time) - R::lgammafn(theta);
  for (int i=0; i < Time.n_elem; i++)
  {
    aux(i) += log(R::bessel_k(2 * sqrt(Rate(i)*Time(i)), -theta, 1));
  }
  return(aux);
}

double logf_RMEIGAM_d(double const& Time,
                      double const& Rate,
                      double const& theta)
{
  double aux = log(2) + ((theta+1)/2)*log(Rate) + ((theta-1)/2)*log(Time) - R::lgammafn(theta);
  aux += log(R::bessel_k(2 * sqrt(Rate*Time), -(theta-1), 1));

  return(aux);
}

double logS_RMEIGAM_d(double const& Time,
                      double const& Rate,
                      double const& theta)
{
  double aux = log(2) + (theta/2) * log(Rate * Time) - R::lgammafn(theta);
  aux += log(R::bessel_k(2 * sqrt(Rate*Time), -theta, 1));
  return(aux);
}

// LOG-LIKELIHOOD (refer to R code for Mixing = "LogNormal" )
// Debug: function tested to match original R implementation
// [[Rcpp::export]]
double HiddenLogLik(arma::vec const& Time,
                    arma::vec const& Event,
                    arma::mat const& DesignMat,
                    arma::vec const& beta,
                    double const& gam,
                    double const& theta,
                    String const& Mixing,
                    String const& BaseModel)
{
  arma::vec out;
  arma::vec Rate = exp(- gam * DesignMat * beta);

  if(Mixing == "None")
  {
    out = Event % (log(Rate) - Rate % pow(Time, gam)) + (1-Event) % (- Rate % pow(Time,gam));
  }
  if(Mixing == "Exponential")
  {
    out = Event % (log(Rate) - 2*log(Rate % pow(Time,gam) + 1)) + (1-Event) % (-log(Rate % pow(Time,gam) + 1));
  }
  if(Mixing == "Gamma")
  {
    out = Event % (log(Rate) - (theta+1)*log((Rate/theta) % pow(Time,gam)+1)) + (1-Event)%(-theta * log((Rate/theta)%pow(Time,gam)+1));
  }
  if(Mixing == "InvGamma")
  {
    out = Event % logf_RMEIGAM(pow(Time,gam), Rate, theta) + (1-Event) % logS_RMEIGAM(pow(Time,gam), Rate, theta);
  }
  if(Mixing == "InvGauss")
  {
      arma::vec auxIG = 2 * Rate % pow(Time,gam) + pow(theta,-2);
      out = Event%(log(Rate)+1/theta-0.5*log(auxIG)-sqrt(auxIG)) + (1-Event)%(1/theta-sqrt(auxIG));
  }

  // ADJUSTEMENT FOR RMW MODELS
  if(BaseModel == "Weibull")
  {
    out += log(gam) * Event + (gam-1) * Event % log(Time);
  }

  return(sum(out));
}

// LOG-LIKELIHOOD (individual observations)
// Debug: function tested to match original R implementation
double HiddenLogLikInd(double const& Time,
                       double const& Event,
                       arma::vec const& DesignMat,
                       arma::vec const& beta,
                       double const& gam,
                       double const& theta,
                       String const& Mixing,
                       String const& BaseModel,
                       double const& ExtraConst)
{
  double out;
  double Rate = exp(- gam * sum(DesignMat % beta)) * ExtraConst;

  if(Mixing == "None")
  {
    out = Event * (log(Rate) - Rate * pow(Time, gam)) + (1-Event) * (- Rate * pow(Time,gam));
  }
  if(Mixing == "Exponential")
  {
    out = Event * (log(Rate) - 2*log(Rate * pow(Time,gam) + 1)) + (1-Event) * (-log(Rate * pow(Time,gam) + 1));
  }
  if(Mixing == "Gamma")
  {
    out = Event * (log(Rate) - (theta+1)*log((Rate/theta) * pow(Time,gam)+1)) + (1-Event)*(-theta * log((Rate/theta)*pow(Time,gam)+1));
  }
  if(Mixing == "InvGamma")
  {
    out = Event * logf_RMEIGAM_d(pow(Time,gam), Rate, theta) + (1-Event) * logS_RMEIGAM_d(pow(Time,gam), Rate, theta);
  }
  if(Mixing == "InvGauss")
  {
    double auxIG = 2 * Rate * pow(Time,gam) + pow(theta,-2);
    out = Event*(log(Rate)+1/theta-0.5*log(auxIG)-sqrt(auxIG)) + (1-Event)*(1/theta-sqrt(auxIG));
  }

  // ADJUSTEMENT FOR RMW MODELS
  if(BaseModel == "Weibull")
  {
    out += log(gam) * Event + (gam-1) * Event * log(Time);
  }

  return(out);
}

// [[Rcpp::export]]
double HiddenRMWreg_DIC(Rcpp::List const& Chain,
                        arma::vec const& Time,
                        arma::vec const& Event,
                        arma::mat const& DesignMat,
                        String const& Mixing,
                        String const& BaseModel)
{
  // EXTRACTING CHAINS
  arma::mat beta = Rcpp::as<arma::mat>(Chain["beta"]);
  arma::vec gam = Rcpp::as<arma::vec>(Chain["gam"]);
  arma::vec theta = Rcpp::as<arma::vec>(Chain["theta"]);

  // CALCULATING POSTERIOR MEDIANS
  arma::vec beta_hat = median(beta).t();
  double gam_hat = median(gam);
  double theta_hat = median(theta);

  // ELEMENT SIZES
  int N = beta.n_rows;

  // LIKELIHOOD FOR EACH ITERATION
  arma::vec L = arma::zeros(N);
  for (int i=0; i < N; i++)
  {
    L(i) = HiddenLogLik(Time, Event, DesignMat, beta.row(i).t(), gam(i), theta(i), Mixing, BaseModel);
  }
  double pd = -2 * mean(L);
  pd += 2 * HiddenLogLik(Time, Event, DesignMat, beta_hat, gam_hat, theta_hat, Mixing, BaseModel);

  return(-2 * mean(L) + pd);
}

// ################################################################################
// ################################################################################
// # FUNCTIONS REQUIRED TO PSEUDO MARGINAL LIKELIHOODS
// ################################################################################
// ################################################################################

// [[Rcpp::export]]
arma::mat HiddenRMWreg_CaseDeletion(Rcpp::List const& Chain,
                                    arma::vec const& Time,
                                    arma::vec const& Event,
                                    arma::mat const& DesignMat,
                                    String const& Mixing,
                                    String const& BaseModel)
{
  // EXTRACTING CHAINS
  arma::mat beta = Rcpp::as<arma::mat>(Chain["beta"]);
  arma::vec gam = Rcpp::as<arma::vec>(Chain["gam"]);
  arma::vec theta = Rcpp::as<arma::vec>(Chain["theta"]);

  // ELEMENT SIZES
  int n = Time.n_elem; int N = beta.n_rows;
  arma::vec logCPO = arma::zeros(n);
  arma::vec KLaux = arma::zeros(n);

  // Offset to avoid under/overflows
  arma::vec offset = arma::zeros(n);
  arma::vec aux1 = arma::zeros(N); arma::vec aux2 = arma::zeros(N);

  for(int i = 0; i < n; i++)
  {
    offset(i) = HiddenLogLikInd(Time(i), Event(i), DesignMat.row(i).t(),
                         beta.row(0).t(), gam(0), theta(0),
                         Mixing, BaseModel, 1);
    for(int iter = 0; iter < N; iter++)
    {
      aux2(iter) = HiddenLogLikInd(Time(i), Event(i), DesignMat.row(i).t(),
                            beta.row(iter).t(), gam(iter), theta(iter),
                            Mixing, BaseModel, 1);
      aux1(iter) = exp(-aux2(iter) + offset(i));
    }
    logCPO(i) = - log(mean(aux1)) + offset(i);
    KLaux(i) = mean(aux2);
  }

  arma::vec KL = KLaux - logCPO;
  arma::vec CALIBRATION = 0.5 * (1+sqrt(1-exp(-2*KL)));

   return(join_rows(join_rows(logCPO,KL) ,CALIBRATION));
}

// ################################################################################
// ################################################################################
// # FUNCTIONS REQUIRED TO ESTIMATE MARGINAL LIKELIHOODS
// ################################################################################
// ################################################################################

// ACCEPTANCY PROBABILITY FOR GRWMH beta[j]
// USE lambda=rep(1,times=n) FOR WEIBULL MODEL WITH MO MIXTURE
// Debug: function tested to match original R implementation
double AcceptProb_beta_j(arma::vec const& beta0,
                         arma::vec const& beta1,
                         double const& gam,
                         arma::vec const& Time,
                         arma::vec const& Event,
                         arma::mat const& X,
                         arma::vec const& lambda,
                         int const& j)
{
  // AUXILIARY QUANTITIES
  arma::mat Xj = RemoveColMat(X, j); arma::vec bj = RemoveElemVec(beta0, j);
  arma::vec SumOtherEffects = exp( - gam * (Xj * bj) );

  // ACCEPTANCE PROBABILITY
  double log_aux = (beta0(j) - beta1(j)) * gam * sum(X.col(j) % Event);
  log_aux += sum( SumOtherEffects % pow(Time,gam) % lambda % (exp(-gam * beta0(j) * X.col(j)) - exp(-gam * beta1(j) * X.col(j))));
  log_aux = exp(log_aux);

  if(log_aux < 1) { return(log_aux); }
  else { return(1);}
}

// ACCEPTANCY PROBABILITY FOR GRWMH GAMMA
// USE lambda=rep(1,times=n) FOR WEIBULL MODEL WITH MIXTURE
// Debug: function tested to match original R implementation
// Debug: includes alphaRMWLN.gam function as a special case (rather than a separate function)
double AcceptProb_gam(double const& gam0,
                      double const& gam1,
                      arma::vec const& Time,
                      arma::vec const& Event,
                      arma::mat const& X,
                      arma::vec const& beta,
                      double const& theta,
                      arma::vec const& lambda,
                      String const& PriorCV,
                      double const& HypTheta,
                      double const& Hyp1Gam,
                      double const& Hyp2Gam,
                      String const& mixing,
                      double const& lower_bound,
                      int FIX_THETA = 0)
{
  double log_aux;
  if(gam1 >= lower_bound)
  {
    log_aux = sum(Event) * log(gam1/gam0);
    log_aux += (gam1 - gam0) * sum( Event % (log(Time) - X * beta));
    log_aux -= sum(lambda % (pow(exp(-X*beta) % Time, gam1) - pow(exp(-X*beta) % Time, gam0)));
    if(FIX_THETA == 0)
    {
      log_aux += HiddenLogPriorTheta(theta, gam1, HypTheta, PriorCV, mixing);
      log_aux -= HiddenLogPriorTheta(theta, gam0, HypTheta, PriorCV, mixing);
    }
    log_aux += R::dgamma(gam1, Hyp1Gam, Hyp2Gam, 1);
    log_aux -= R::dgamma(gam0, Hyp1Gam, Hyp2Gam, 1);
    if(mixing == "LogNormal")
    {
      log_aux -= Time.n_elem * log(gam1 / gam0);
      Rcpp::Rcout << "Time.n_elem * log(gam1 / gam0)" << Time.n_elem * log(gam1 / gam0) << std::endl;
      log_aux -= (pow(gam1, -2.0) - pow(gam0, -2.0)) * (0.5/theta) * sum(pow(log(lambda), 2.0));
    }
    log_aux = exp(log_aux);
  }
  else { log_aux = 0; }

  if(log_aux < 1) { return(log_aux); }
  else { return(1);}
}

// ACCEPTANCY PROBABILITY FOR GRWMH THETA
// USE lambda=rep(1,times=n) FOR WEIBULL MODEL WITH MIXTURE
double AcceptProb_theta(double const& theta0,
                        double const& theta1,
                        double const& gam,
                        arma::vec const& lambda,
                        String const& PriorCV,
                        double const& HypTheta,
                        String const& mixing,
                        double const& n)
{
  double log_aux;

  if(mixing == "Gamma")
  {
    log_aux = n * (theta1*log(theta1) - R::lgammafn(theta1) - theta0*log(theta0) + R::lgammafn(theta0));
    log_aux += (theta1-theta0) * sum(log(lambda)) - (theta1-theta0) * sum(lambda);
  }
  if(mixing == "InvGamma")
  {
    log_aux = n * (R::lgammafn(theta0) - R::lgammafn(theta1)) + (theta0-theta1) * sum(log(lambda));
  }
  if(mixing == "InvGauss")
  {
    log_aux = 0.5 * (pow(theta0,-2) - pow(theta1,-2)) * sum(lambda);
    log_aux -= n * (pow(theta0,-1) - pow(theta1,-1));
  }
  if(mixing == "LogNormal")
  {
    log_aux = (n/2) * log(theta0/theta1);
    log_aux -= 0.5 * (pow(theta1,-1) - pow(theta0,-1)) * sum(pow(log(lambda),2));
  }
  // PRIOR FACTOR
  log_aux += HiddenLogPriorTheta(theta1, gam, HypTheta, PriorCV, mixing);
  log_aux -= HiddenLogPriorTheta(theta0, gam, HypTheta, PriorCV, mixing);

  log_aux = exp(log_aux);

  // REJECTING INVALID VALUES
  if((mixing == "Gamma") & (theta1 <= 2/gam)) { log_aux = 0; }
  if((mixing == "InvGamma") & (theta1 <= 1)) { log_aux = 0; }

  if(log_aux < 1) { return(log_aux); }
  else { return(1);}
}

// ################################################################################
// ################################################################################
// # FUNCTIONS RELATED TO OUTLIER DETECTION
// ################################################################################
// ################################################################################

//                                 int N, // Total number of MCMC draws
//                                 int burn, // Burning period for MCMC chain

// [[Rcpp::export]]
double HiddenRMWreg_BFoutlierObs(Rcpp::List Chain,
                                 int const& Obs,
                                 arma::vec const& RefLambda,
                                 arma::vec const& Time,
                                 arma::vec const& Event,
                                 arma::mat const& DesignMat,
                                 String const& PriorCV,
                                 double const& HypTheta,
                                 double const& Hyp1Gam,
                                 double const& Hyp2Gam,
                                 String const& Mixing,
                                 String const& BaseModel,
                                 int thin, // Thinning period for MCMC chain
                                 int lambdaPeriod,
                                 double ar)
{
  // EXTRACTING CHAINS
  arma::mat beta = Rcpp::as<arma::mat>(Chain["beta"]);
  arma::vec gam = Rcpp::as<arma::vec>(Chain["gam"]);
  arma::vec theta = Rcpp::as<arma::vec>(Chain["theta"]);

  int N = beta.n_rows;

  double BF = 1;
  arma::vec BF_iter = arma::ones(N);

  if(Mixing == "Exponential")
  {
    for(int iter = 0; iter < N; iter++)
    {
      BF_iter(iter) = HiddenLogLikInd(Time(Obs), Event(Obs), DesignMat.row(Obs).t(),
                                      beta.row(iter).t(), gam(iter), theta(iter),
                                      "None", BaseModel, RefLambda(Obs));
      BF_iter(iter) -= HiddenLogLikInd(Time(Obs), Event(Obs), DesignMat.row(Obs).t(),
                                       beta.row(iter).t(), gam(iter), theta(iter),
                                       Mixing, BaseModel, 1);
      BF_iter(iter) = exp(BF_iter(iter));
    }
    BF = mean(BF_iter);
  }

  return(BF);
}

// [[Rcpp::export]]
arma::vec HiddenRMWreg_BFoutlier(Rcpp::List Chain,
                                 arma::vec const& RefLambda,
                                 arma::vec const& Time,
                                 arma::vec const& Event,
                                 arma::mat const& DesignMat,
                                 String const& PriorCV,
                                 double const& HypTheta,
                                 double const& Hyp1Gam,
                                 double const& Hyp2Gam,
                                 String const& Mixing,
                                 String const& BaseModel,
                                 int thin, // Thinning period for MCMC chain
                                 int lambdaPeriod,
                                 double ar)
{
  int n = Time.n_elem;

  arma::vec BF = arma::ones(n);
  for(int i = 0; i < n; i++)
  {
//    Rcpp::Rcout << "Observation: " << i + 1 << std::endl;
    BF(i) = HiddenRMWreg_BFoutlierObs(Chain, i, RefLambda,
                                      Time, Event, DesignMat,
                                      PriorCV, HypTheta, Hyp1Gam, Hyp2Gam,
                                      Mixing, BaseModel,
                                      thin, lambdaPeriod, ar);
  }

  return(BF);
}
