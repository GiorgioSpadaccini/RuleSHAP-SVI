#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>
#include <cfloat>
#include <limits>
#include <stdexcept>

// [[Rcpp::plugins("cpp17")]]
// [[Rcpp::depends(RcppEigen)]]





/* Data structure for a Chebyshev series over a given interval */
struct cheb_series {
  const double* c;   // coefficients
  int order;         // order of expansion
  double a;          // lower interval point
  double b;          // upper interval point
  int order_sp;      // effective single precision order

  cheb_series(const double* coeffs, int ord, double lower, double upper, int ord_sp)
    : c(coeffs), order(ord), a(lower), b(upper), order_sp(ord_sp) {}
};

/* Chebyshev series data */
static double AE11_data[39] = {
  0.121503239716065790,
  -0.065088778513550150,
  0.004897651357459670,
  -0.000649237843027216,
  0.000093840434587471,
  0.000000420236380882,
  -0.000008113374735904,
  0.000002804247688663,
  0.000000056487164441,
  -0.000000344809174450,
  0.000000058209273578,
  0.000000038711426349,
  -0.000000012453235014,
  -0.000000005118504888,
  0.000000002148771527,
  0.000000000868459898,
  -0.000000000343650105,
  -0.000000000179796603,
  0.000000000047442060,
  0.000000000040423282,
  -0.000000000003543928,
  -0.000000000008853444,
  -0.000000000000960151,
  0.000000000001692921,
  0.000000000000607990,
  -0.000000000000224338,
  -0.000000000000200327,
  -0.000000000000006246,
  0.000000000000045571,
  0.000000000000016383,
  -0.000000000000005561,
  -0.000000000000006074,
  -0.000000000000000862,
  0.000000000000001223,
  0.000000000000000716,
  -0.000000000000000024,
  -0.000000000000000201,
  -0.000000000000000082,
  0.000000000000000017
};
static cheb_series AE11_cs(AE11_data, 38, -1, 1, 20);

static double AE12_data[25] = {
  0.582417495134726740,
  -0.158348850905782750,
  -0.006764275590323141,
  0.005125843950185725,
  0.000435232492169391,
  -0.000143613366305483,
  -0.000041801320556301,
  -0.000002713395758640,
  0.000001151381913647,
  0.000000420650022012,
  0.000000066581901391,
  0.000000000662143777,
  -0.000000002844104870,
  -0.000000000940724197,
  -0.000000000177476602,
  -0.000000000015830222,
  0.000000000002905732,
  0.000000000001769356,
  0.000000000000492735,
  0.000000000000093709,
  0.000000000000010707,
  -0.000000000000000537,
  -0.000000000000000716,
  -0.000000000000000244,
  -0.000000000000000058
};
static cheb_series AE12_cs(AE12_data, 24, -1, 1, 15);

static double E11_data[19] = {
  -16.11346165557149402600,
  7.79407277874268027690,
  -1.95540581886314195070,
  0.37337293866277945612,
  -0.05692503191092901938,
  0.00721107776966009185,
  -0.00078104901449841593,
  0.00007388093356262168,
  -0.00000620286187580820,
  0.00000046816002303176,
  -0.00000003209288853329,
  0.00000000201519974874,
  -0.00000000011673686816,
  0.00000000000627627066,
  -0.00000000000031481541,
  0.00000000000001479904,
  -0.00000000000000065457,
  0.00000000000000002733,
  -0.00000000000000000108
};
static cheb_series E11_cs(E11_data, 18, -1, 1, 13);

static double E12_data[16] = {
  -0.03739021479220279500,
  0.04272398606220957700,
  -0.13031820798497005440,
  0.01441912402469889073,
  -0.00134617078051068022,
  0.00010731029253063780,
  -0.00000742999951611943,
  0.00000045377325690753,
  -0.00000002476417211390,
  0.00000000122076581374,
  -0.00000000005485141480,
  0.00000000000226362142,
  -0.00000000000008635897,
  0.00000000000000306291,
  -0.00000000000000010148,
  0.00000000000000000315
};
static cheb_series E12_cs(E12_data, 15, -1, 1, 10);

static double AE13_data[25] = {
  -0.605773246640603460,
  -0.112535243483660900,
  0.013432266247902779,
  -0.001926845187381145,
  0.000309118337720603,
  -0.000053564132129618,
  0.000009827812880247,
  -0.000001885368984916,
  0.000000374943193568,
  -0.000000076823455870,
  0.000000016143270567,
  -0.000000003466802211,
  0.000000000758754209,
  -0.000000000168864333,
  0.000000000038145706,
  -0.000000000008733026,
  0.000000000002023672,
  -0.000000000000474132,
  0.000000000000112211,
  -0.000000000000026804,
  0.000000000000006457,
  -0.000000000000001568,
  0.000000000000000383,
  -0.000000000000000094,
  0.000000000000000023
};
static cheb_series AE13_cs(AE13_data, 24, -1, 1, 15);

static double AE14_data[26] = {
  -0.18929180007530170,
  -0.08648117855259871,
  0.00722410154374659,
  -0.00080975594575573,
  0.00010999134432661,
  -0.00001717332998937,
  0.00000298562751447,
  -0.00000056596491457,
  0.00000011526808397,
  -0.00000002495030440,
  0.00000000569232420,
  -0.00000000135995766,
  0.00000000033846628,
  -0.00000000008737853,
  0.00000000002331588,
  -0.00000000000641148,
  0.00000000000181224,
  -0.00000000000052538,
  0.00000000000015592,
  -0.00000000000004729,
  0.00000000000001463,
  -0.00000000000000461,
  0.00000000000000148,
  -0.00000000000000048,
  0.00000000000000016,
  -0.00000000000000005
};
static cheb_series AE14_cs(AE14_data, 25, -1, 1, 13);

/* Chebyshev series evaluator */
inline double cheb_eval(const cheb_series& cs, double x)
{
  double d = 0.0, dd = 0.0;
  double y  = (2.0 * x - cs.a - cs.b) / (cs.b - cs.a);
  double y2 = 2.0 * y;

  for (int j = cs.order; j >= 1; --j) {
    double temp = d;
    d = y2 * d - dd + cs.c[j];
    dd = temp;
  }

  return y * d - dd + 0.5 * cs.c[0];
}

/* Exponential integral E1(x) */

double expint_E1(double x, bool scale = false)
{
  using std::exp;
  using std::fabs;
  using std::log;
  using std::isnan;
  using std::numeric_limits;

  if (isnan(x)) return x;

  const double LOG_DBL_MIN = log(numeric_limits<double>::min());
  const double xmaxt = -LOG_DBL_MIN;
  const double xmax  = xmaxt - log(xmaxt);

  if (x < -xmax && !scale) {
    throw std::overflow_error("overflow in expint_E1");
  } else if (x <= -10.0) {
    const double s = 1.0 / x * (scale ? 1.0 : exp(-x));
    const double cheb = cheb_eval(AE11_cs, 20.0 / x + 1.0);
    return s * (1.0 + cheb);
  } else if (x <= -4.0) {
    const double s = 1.0 / x * (scale ? 1.0 : exp(-x));
    const double cheb = cheb_eval(AE12_cs, (40.0 / x + 7.0) / 3.0);
    return s * (1.0 + cheb);
  } else if (x <= -1.0) {
    const double s = (scale ? exp(x) : 1.0);
    const double ln_term = -log(fabs(x));
    const double cheb = cheb_eval(E11_cs, (2.0 * x + 5.0) / 3.0);
    return s * (ln_term + cheb);
  } else if (x == 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  } else if (x <= 1.0) {
    const double s = (scale ? exp(x) : 1.0);
    const double ln_term = -log(fabs(x));
    const double cheb = cheb_eval(E12_cs, x);
    return s * (ln_term - 0.6875 + x + cheb);
  } else if (x <= 4.0) {
    const double s = 1.0 / x * (scale ? 1.0 : exp(-x));
    const double cheb = cheb_eval(AE13_cs, (8.0 / x - 5.0) / 3.0);
    return s * (1.0 + cheb);
  } else if (x <= xmax || scale) {
    const double s = 1.0 / x * (scale ? 1.0 : exp(-x));
    const double cheb = cheb_eval(AE14_cs, 8.0 / x - 1.0);
    const double res = s * (1.0 + cheb);
    if (res == 0.0)
      throw std::underflow_error("underflow in expint_E1");
    return res;
  } else {
    throw std::underflow_error("underflow in expint_E1");
  }
}

































// Forward declaration for pnorm_both
void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);

double pnorm_cpp(double x, double mu = 0.0, double sigma = 1.0,
                 bool lower_tail = true, bool log_p = false) {
  double p, cp;

  // Handle NaNs using R's internal checks
  if (R_IsNA(x) || R_IsNA(mu) || R_IsNA(sigma)) {
    return x + mu + sigma;
  }

  if (!std::isfinite(x) && mu == x) return NA_REAL;

  if (sigma <= 0) {
    if (sigma < 0) {
      Rcpp::warning("NaNs produced: sigma < 0");
      return NA_REAL;
    }
    // Case: sigma = 0
    bool is_lower = (x < mu);
    if (log_p) return is_lower ? R_NegInf : 0.0;
    return is_lower ? 0.0 : 1.0;
  }

  p = (x - mu) / sigma;
  if (!std::isfinite(p)) {
    bool is_lower = (x < mu);
    if (log_p) return is_lower ? R_NegInf : 0.0;
    return is_lower ? 0.0 : 1.0;
  }

  // Call the core worker function
  pnorm_both(p, &p, &cp, (lower_tail ? 0 : 1), (int)log_p);

  return (lower_tail ? p : cp);
}

// Internal worker function maintaining the Cody algorithm
void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p) {
  // Coefficients as defined in Mathlib
  static const double a[5] = {2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249, 0.065682337918207449113};
  static const double b[4] = {47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956};
  static const double c[9] = {0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226, 2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218, 1.0765576773720192317e-8};
  static const double d[8] = {22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755, 18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727};
  static const double p_coef[6] = {0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466, 2.9112874951168792e-5, 0.02307344176494017303};
  static const double q_coef[5] = {1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244, 7.29751555083966205e-5};

  double xden, xnum, temp, del, xsq, y;
  const double eps = DBL_EPSILON * 0.5;
  const double SIXTEN = 16.0;

  if (R_IsNA(x)) { *cum = *ccum = x; return; }

  bool lower = (i_tail != 1);
  bool upper = (i_tail != 0);

  y = std::abs(x);

  if (y <= 0.67448975) {
    if (y > eps) {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (int i = 0; i < 3; ++i) {
        xnum = (xnum + a[i]) * xsq;
        xden = (xden + b[i]) * xsq;
      }
    } else {
      xnum = xden = 0.0;
    }

    temp = x * (xnum + a[3]) / (xden + b[3]);
    if (lower) *cum = 0.5 + temp;
    if (upper) *ccum = 0.5 - temp;
    if (log_p) {
      if (lower) *cum = std::log(*cum);
      if (upper) *ccum = std::log(*ccum);
    }
  }
  else if (y <= 5.65685424949) { // M_SQRT_32 approx
    xnum = c[8] * y;
    xden = y;
    for (int i = 0; i < 7; ++i) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);

    auto do_del = [&](double X) {
      xsq = std::trunc(X * SIXTEN) / SIXTEN;
      del = (X - xsq) * (X + xsq);
      if (log_p) {
        *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + std::log(temp);
        if ((lower && x > 0.) || (upper && x <= 0.))
          *ccum = std::log1p(-std::exp(-xsq * xsq * 0.5) * std::exp(-del * 0.5) * temp);
      } else {
        *cum = std::exp(-xsq * xsq * 0.5) * std::exp(-del * 0.5) * temp;
        *ccum = 1.0 - *cum;
      }
    };

    do_del(y);
    if (x > 0.) { std::swap(*cum, *ccum); }
  }
  else if ((log_p && y < 1e170) || (lower && x > -37.5193 && x < 8.2924) || (upper && x > -8.2924 && x < 37.5193)) {
    xsq = 1.0 / (x * x);
    xnum = p_coef[5] * xsq;
    xden = xsq;
    for (int i = 0; i < 4; ++i) {
      xnum = (xnum + p_coef[i]) * xsq;
      xden = (xden + q_coef[i]) * xsq;
    }
    temp = xsq * (xnum + p_coef[4]) / (xden + q_coef[4]);
    temp = (M_1_SQRT_2PI - temp) / y;

    auto do_del_large = [&](double X) {
      double xsq_inner = std::trunc(X * SIXTEN) / SIXTEN;
      double del_inner = (X - xsq_inner) * (X + xsq_inner);
      if (log_p) {
        *cum = (-xsq_inner * xsq_inner * 0.5) + (-del_inner * 0.5) + std::log(temp);
        if ((lower && x > 0.) || (upper && x <= 0.))
          *ccum = std::log1p(-std::exp(-xsq_inner * xsq_inner * 0.5) * std::exp(-del_inner * 0.5) * temp);
      } else {
        *cum = std::exp(-xsq_inner * xsq_inner * 0.5) * std::exp(-del_inner * 0.5) * temp;
        *ccum = 1.0 - *cum;
      }
    };

    do_del_large(x);
    if (x > 0.) { std::swap(*cum, *ccum); }
  }
  else {
    if (x > 0) {
      *cum = log_p ? 0.0 : 1.0;
      *ccum = log_p ? R_NegInf : 0.0;
    } else {
      *cum = log_p ? R_NegInf : 0.0;
      *ccum = log_p ? 0.0 : 1.0;
    }
  }
}































//' Fit a horseshoe model
//'
//' This function fits a plain horseshoe Bayesian regression for continuous and binary data with
//' a MCMC algorithm based on Gibbs sampling. This function is also compatible with a grouped
//' local shrinkage for linear terms and rule terms separately.
//'
//' @param X a numeric \eqn{n \times p} matrix whose columns are predictors and whose rows are observations.
//' This is used as training data to fit the model.
//' @param y a numeric vector of outcomes used as training data to fit the model. For binary outcome, this should
//' be converted into zeroes and ones.
//' @param p.lin number of columns of X which are linear terms and should thus get their
//' own global shrinkage separately from the last remaining columns of X. By default, all
//' columns are linear terms and therefore all columns share the same shrinkage.
//' @param family A character string denoting the type of (Generalized) linear model
//' that is being fitted
//' @param method.tau.rules character string denoting the prior distribution of the
//' global shrinkage assigned to rules: the choice is between a half-Cauchy prior
//' (\code{"halfCauchy"}),a truncated Cauchy prior (\code{"truncatedCauchy"}) or
//' a fixed value (\code{"fixed"})
//' @param tau.rules a numeric entry denoting which value the global shrinkage for rules should take,
//' in the case where \code{method.tau.rules} is set to \code{"fixed"}. If other values of \code{method.tau.lin}
//' are chosen, this parameter is used as initial value to start the Gibbs sampling from.
//' @param method.tau.lin character string denoting the prior distribution of the
//' global shrinkage assigned to linear terms: the choice is between a half-Cauchy prior
//' (\code{"halfCauchy"}),a truncated Cauchy prior (\code{"truncatedCauchy"}) or
//' a fixed value (\code{"fixed"})
//' @param tau.lin a numeric entry denoting which value the global shrinkage for linear terms should take,
//' in the case where \code{method.tau.lin} is set to \code{"fixed"}. If other values of \code{method.tau.lin}
//' are chosen, this parameter is used as initial value to start the Gibbs sampling from.
//' @param method.tau character string denoting the prior distribution of the
//' overall global shrinkage assigned to all terms: the choice is between a half-Cauchy prior
//' (\code{"halfCauchy"}),a truncated Cauchy prior (\code{"truncatedCauchy"}) or
//' a fixed value (\code{"fixed"})
//' @param tau a numeric entry denoting which value the global shrinkage for linear terms should take,
//' in the case where \code{method.tau} is set to \code{"fixed"}. If other values of \code{method.tau}
//' are chosen, this parameter is used as initial value to start the Gibbs sampling from.
//' @param method.sigma character string denoting which prior distribution the error variance \eqn{\sigma^2}
//' should have: the choice is between Jeffre's prior (\code{"Jeffreys"}) or a fixed
//' value (\code{"fixed"}).
//' @param Sigma2 a numeric entry denoting which value the error variance \eqn{\sigma^2} should take,
//' in the case where \code{method.sigma} is set to \code{"fixed"}. If other values of \code{method.sigma}
//' are chosen, this parameter is ignored.
//' @param prior_lambda a vector of \eqn{p} positive numeric entries representing local shrinkage relaxation
//' (\code{prior_lambda > 1}) or contraction (\code{prior_lambda < 1}). The \eqn{j}-th entry of this vector
//' denotes by how much the local shrinkage should be relaxed for the \eqn{j}-th predictor.
//' @param burn.in Number of burn-in MCMC iterations that are discarded before
//' the draws start to be recorded.
//' @param nmc Number of MCMC iterations after the initial burn-in period
//' @param thin Thinning parameter: if saving all MCMC draws is computationally
//' infeasible, thinning makes the model only record every \code{thin}-th iteration.
//' @param verbose Logical value determining whether updates should be given on the progress
//' @param ping If the model is selected to be verbose, this denotes every how many
//' MCMC draws an update message should be printed
//' @param Beta_init a numeric vector of lenth \eqn{p} representing the initial coefficients to start
//' the Gibbs sampler with.
//' @param Lambda_init a numeric vector of lenth \eqn{p} representing the initial values of
//' local shrinkage to start the Gibbs sampler with.
//'
//' @return BetaSamples A matrix where every column is a MCMC sample of the coefficients.
//' @return BetaHat A vector containing the posterior mean of the coefficients.
//' @return TauLinSamples A vector where every entry is a MCMC sample of the global shrinkage of the linear terms.
//' @return TauLinHat The posterior mean of the global shrinkage of the linear terms.
//' @return TauRulesSamples A vector where every entry is a MCMC sample of the global shrinkage of the rules.
//' @return TauRulesHat The posterior mean of the global shrinkage of the rules.
//' @return LambdaSamples A matrix where every column is a MCMC sample of the local shrinkages, one entry per term.
//' @return LambdaHat A vector containing the posterior mean of the local shrinkages.
//' @export
// [[Rcpp::export]]
Rcpp::List hs_VB (Eigen::VectorXd y,
                  const Eigen::MatrixXd X,
                  const char family,
                  const int p_lin,
                  const Eigen::ArrayXi delta,
                  const Eigen::VectorXd Beta_init,
                  const Eigen::ArrayXd Lambda_init,
                  const Eigen::ArrayXd Omega_scale_init,
                  const Eigen::ArrayXd prior_lambda,
                  const double tau_rules=1,
                  const double tau_lin=1,
                  const double Sigma2=1,
                  bool update_sigma=true,
                  const int batch_size = 2e2,
                  const double lrate = 1e-1,
                  const char threshold_type = 'M',
                  double threshold = 0,
                  const int threshold_update = 10,
                  const int max_iter = 1000,
                  const double eps = 1e-03,
                  const int ping = 100,
                  const bool verbose = 1,
                  const double adapt_p0=0,
                  const double adapt_p1=-4.6e-4,
                  const double step_size=1e-2,
                  const double decay=0.9,
                  const double step_size_sigma=1e-2,
                  const double decay_sigma=0.9,
                  const bool invert_sigma=false,
                  const bool update_tau=true,
                  const bool dvar = true)
{

  ////// Preallocate and initialize

  // Initial quantities and matrix products
  const int n = y.size(), p = X.cols();
  const Eigen::MatrixXd Xt = X.transpose();
  const Eigen::MatrixXd XtX = Xt*X; //Might want to replace this with a version that uses at most N observations, with e.g. N=2e4
  Eigen::VectorXd Xty = Xt*y;
  const Eigen::ArrayXd XtX_diag = XtX.diagonal();
  double norm_y = y.squaredNorm();
  const bool algoBhatt = true; //!!! Just remove this, and remove the option to do other algorithm

  // Thresholding and batching quantities
  double u_i, p_i, threshold_Db=1e10;
  int m_eff = p, S=p;
  Eigen::ArrayXd meff_frac(p), delta_D=Eigen::ArrayXd::Zero(p);
  std::vector<int> batch_i(batch_size), samples_vec(n),
                   active_i(p), passive_i(p);

  std::iota(batch_i.begin(),batch_i.end(),0);
  std::iota(samples_vec.begin(),samples_vec.end(),0);
  std::iota(active_i.begin(),active_i.end(),0);
  
  Eigen::ArrayXd estimated_dvar=Eigen::ArrayXd::Zero(p), D_diag_sorted(p), ranking(p);

  // Shrinkage parameters
  const double shape_tauL2 = ((double)(p_lin+1))/(double)2, shape_tauR2 = ((double)(p-p_lin+1))/(double)2;
  const Eigen::ArrayXd prior_lambda2inv=1/prior_lambda.square();
  double scale_tauL2=std::pow(tau_lin,2)*shape_tauL2,
    scale_tauR2=std::pow(tau_rules,2)*shape_tauR2,
    avg_tauL2inv=shape_tauL2/scale_tauL2, avg_tauR2inv=shape_tauR2/scale_tauR2, avg_tau2inv,
    meandraws_denom, meandraws_tauL, meandraws_tauR, meandraws_tau, log_mass_tau_all,
    log_mass_lambdas, SqrtLogDetSigma;
  Eigen::ArrayXd scale_lambda2=Lambda_init.square(), avg_lambda2inv=1/scale_lambda2,
    denom_lambda(p), expint_lambda(p), avg_tau2inv_all(p), D_diag(p),
    D_diag_old(p), Sigma_diag(p);

  avg_tau2inv_all << Eigen::VectorXd::Constant(p_lin,avg_tauL2inv), Eigen::VectorXd::Constant(p-p_lin,avg_tauR2inv);
  D_diag = avg_tau2inv_all*avg_lambda2inv*prior_lambda2inv;
  D_diag_old = avg_tau2inv_all*avg_lambda2inv*prior_lambda2inv;
  Sigma_diag=1/(XtX_diag+D_diag);
  Eigen::ArrayXd v_vec=Eigen::ArrayXd::Constant(n,(XtX_diag*Sigma_diag).sum()/n), m_vvec=Eigen::ArrayXd::Zero(n), grad_vvec(n);
  Eigen::ArrayXd cpo=Eigen::ArrayXd::Zero(max_iter);

  // Beta, sigma2 paramters
  const double shape_sigma2 = n/2;
  double scale_sigma2=Sigma2*shape_sigma2, avg_sigma2inv=1/Sigma2, log_mass_beta_sigma;
  Eigen::ArrayXd y_hat=X*Beta_init;

  Eigen::ArrayXd mu_beta=Beta_init.array(),avg_beta2_sigma2inv=Sigma_diag+Beta_init.array().square()*avg_sigma2inv;

  // Model fit parameters
  Eigen::LLT<Eigen::MatrixXd> L;
  Eigen::MatrixXd batch_Xt, batch_X, U;
  double elbo_old, elbo = -1e15, elbo_new, first_Delbo=0, logdens, avg_tauL, avg_tauR, avg_tau;
  std::vector<int> lin_index(p_lin), rule_index(p-p_lin);

  std::iota(lin_index.begin(),lin_index.end(),0);
  std::iota(rule_index.begin(),rule_index.end(),p_lin);

  //Draw gammas for estimation of integrals
  const int n_draws = 100000;
  const Eigen::ArrayXd gdraws_y=Rcpp::as<Eigen::VectorXd>(Rcpp::rgamma(n_draws,shape_tauL2,1)).array(),
    gdraws_y_high=Rcpp::as<Eigen::VectorXd>(Rcpp::rgamma(n_draws,shape_tauL2+1,1)).array(),
    gdraws_z=Rcpp::as<Eigen::VectorXd>(Rcpp::rgamma(n_draws,shape_tauR2,1)).array(),
    gdraws_z_high=Rcpp::as<Eigen::VectorXd>(Rcpp::rgamma(n_draws,shape_tauR2+1,1)).array();
  Eigen::ArrayXd ydraws(n_draws), ydraws_high(n_draws), zdraws(n_draws), zdraws_high(n_draws);


  Eigen::MatrixXd VarInvUt=Eigen::MatrixXd::Identity(batch_size,batch_size);
  double shape_omega=1;
  //Eigen::ArrayXd scale_omega=(y_hat.square()+(X.array().square().rowwise()*Sigma_diag.transpose()).rowwise().sum()).sqrt();
  Eigen::ArrayXd scale_omega=Omega_scale_init;
  Eigen::ArrayXd scale_omega_exp=scale_omega.exp();
  Eigen::ArrayXd avg_omega=-(2/(scale_omega_exp+1)-1)/scale_omega*(shape_omega/2);
  Eigen::ArrayXd sqrt_omega=avg_omega(batch_i).sqrt();
  double avg_avg_omega=avg_omega.sum()/n;
  Eigen::ArrayXd first_guess(batch_size), estimated_vars(batch_size);
  double target_sum=1;
  Eigen::ArrayXd logtimes(n);
  std::vector<int> censored_index(n);
  Eigen::VectorXi censored(n);
  int n_cens=2;
  Eigen::ArrayXd ycensvar = Eigen::ArrayXd::Zero(n);

  if(family == 'b'){
    avg_sigma2inv=1;
    update_sigma=false;
    y=y.array()-0.5;
  } else if(family == 's'){
    // Find censored indices
    censored_index.clear();
    for (int j = 0; j < n; ++j) {
      if (!delta(j)){
        censored_index.push_back(j);
      }
    }
    censored = Eigen::Map<Eigen::VectorXi>(censored_index.data(), censored_index.size());
    logtimes = y.array().log();
    y = logtimes + 0;
    n_cens = static_cast<int>(censored.size());
    ycensvar(censored) = Eigen::ArrayXd::Constant(n_cens,1/avg_sigma2inv);
  }
  Eigen::ArrayXd trunc_alpha(n_cens), Z_phi(n_cens); //!!!
  Eigen::VectorXd z(batch_size);

  // Re-adjusting with glm
  Rcpp::Environment stats_pkg = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function glm_func = stats_pkg["glm.fit"];
  Rcpp::Function binomial_func = stats_pkg["binomial"];
  Rcpp::List family_obj = binomial_func(), fit;
  double coeff=1;

  // Adam parameters
  Eigen::ArrayXd grad_lambda2(p), m_lambda2=Eigen::ArrayXd::Zero(p),
    grad_mu(p),m_mu=Eigen::ArrayXd::Zero(p),
    grad_var(p), m_var=Eigen::ArrayXd::Zero(p),
    new_Sigma_diag(p), new_mu=Eigen::ArrayXd::Zero(p),
    grad_omega(n), m_omega=Eigen::ArrayXd::Zero(n),
    grad_ycens(n_cens), m_ycens=Eigen::ArrayXd::Zero(n_cens),
    grad_ycensvar(n_cens), m_ycensvar=Eigen::ArrayXd::Zero(n_cens),
    grad_deltaD(p), m_deltaD(p);
  double b1_t=1, b2_t=1,
    grad_tauL2, m_tauL2=0,
    grad_tauR2, m_tauR2=0,
    grad_sigma2, m_sigma2=0,
    grad_coeff, m_coeff=0,
    grad_cpo1, m_cpo1=0,
    grad_cpo2, m_cpo2=0,
    cpo1=0, cpo1_old=0,
    cpo2=0, cpo2_old=0,
    cpo1_min = 1e10, cpo2_min=1e10;
    
    
    
  //Randomly generate and then set seed (R-c++ compatibility)
  std::uint32_t cpp_seed = static_cast<std::uint32_t>(R::unif_rand() * 4294967296.0);
  std::mt19937 engine(cpp_seed);



  // Start update loop
  for (int i = 0; i < max_iter; i++){

    //Update lambda
    grad_lambda2=scale_lambda2-avg_beta2_sigma2inv*avg_tau2inv_all*prior_lambda2inv/2;
    m_lambda2=decay*m_lambda2+(1-decay)*grad_lambda2;
    scale_lambda2=scale_lambda2-step_size*m_lambda2;
    expint_lambda=scale_lambda2.unaryExpr([](double v){return expint_E1(v,true);});
    denom_lambda=scale_lambda2*expint_lambda;
    avg_lambda2inv=(1-denom_lambda)/denom_lambda;

    if(i==1 || update_tau){
    // Update tau, tauL and tauR jointly
    grad_tauL2=scale_tauL2-(avg_beta2_sigma2inv(lin_index)*
      prior_lambda2inv(lin_index)*avg_lambda2inv(lin_index)).sum()/2;
    m_tauL2=decay*m_tauL2+(1-decay)*grad_tauL2;
    scale_tauL2=scale_tauL2-step_size*m_tauL2;

    grad_tauR2=scale_tauR2-(avg_beta2_sigma2inv(rule_index)*
      prior_lambda2inv(rule_index)*avg_lambda2inv(rule_index)).sum()/2;
    m_tauR2=decay*m_tauR2+(1-decay)*grad_tauR2;
    scale_tauR2=scale_tauR2-step_size*m_tauR2;

    ydraws=gdraws_y/scale_tauL2;
    ydraws_high=gdraws_y_high/scale_tauL2;
    zdraws=gdraws_z/scale_tauR2;
    zdraws_high=gdraws_z_high/scale_tauR2;

    meandraws_denom=(((ydraws-zdraws)+ydraws.sqrt()*(zdraws-1.0)-
      zdraws.sqrt()*(ydraws-1.0))/
        ((ydraws-1.0)*(zdraws-1.0)*(zdraws-ydraws))).sum()/n_draws;

    meandraws_tauL=(((ydraws_high-zdraws)+ydraws_high.sqrt()*(zdraws-1)-
      zdraws.sqrt()*(ydraws_high-1))/
        ((ydraws_high-1)*(zdraws-1)*(zdraws-ydraws_high))).sum()/n_draws;

    meandraws_tauR=(((ydraws-zdraws_high)+ydraws.sqrt()*(zdraws_high-1)-
      zdraws_high.sqrt()*(ydraws-1))/
        ((ydraws-1)*(zdraws_high-1)*(zdraws_high-ydraws))).sum()/n_draws;

    avg_tauL2inv=meandraws_tauL/meandraws_denom*shape_tauL2/scale_tauL2;
    avg_tauR2inv=meandraws_tauR/meandraws_denom*shape_tauR2/scale_tauR2;

    log_mass_tau_all=log(meandraws_denom)-
      shape_tauL2*log(scale_tauL2)-
      shape_tauR2*log(scale_tauR2);
    }



    //Update joint distribution of beta and sigma2
    avg_tau2inv_all << Eigen::VectorXd::Constant(p_lin,avg_tauL2inv), Eigen::VectorXd::Constant(p-p_lin,avg_tauR2inv);
    D_diag=avg_tau2inv_all*avg_lambda2inv*prior_lambda2inv;
    grad_deltaD=delta_D-(D_diag-D_diag_old);
    m_deltaD=decay*m_deltaD+(1-decay)*grad_deltaD;
    delta_D=delta_D-step_size*m_deltaD;
    D_diag_old=(D_diag+0).eval();

    // Determine predictors to update
    if (threshold_type == 'M' || threshold_type == 'b') {
      if(i % threshold_update == 0){
        u_i = R::runif(0, 1);
        p_i = std::exp(adapt_p0 + adapt_p1 * i);

        if (u_i < p_i) {
          if(family=='b'){
            // Make poor assumption of omega_i all equal
            meff_frac = D_diag / (avg_avg_omega*XtX_diag);
          } else{
            meff_frac = D_diag / XtX_diag;
          }

          m_eff = static_cast<int>(std::ceil((1.0 / (meff_frac + 1.0)).sum()));
        }
      
      // m_eff = p - (Sigma_diag*D_diag).sum(); // !!! overwirting definition of m_eff to a hopefully better estimate of degrees of freedom
      }



      // add +10 and clip at p
      m_eff = std::min(p, static_cast<int>(std::ceil(m_eff)));// + 10); !!!
      m_eff = std::max(1,m_eff);

      // take m_eff-th highest shrinkage
      D_diag_sorted=(D_diag+0).eval();
      std::sort(D_diag_sorted.data(), D_diag_sorted.data() + D_diag_sorted.size());

      threshold = D_diag_sorted(m_eff-1);
      
      if(threshold_type=='b'){
      // Now also least recent predictors
      ranking=(delta_D*mu_beta).abs();
      // take m_eff-th smallest value (with minus signs for reverse order
      D_diag_sorted=(0-ranking).eval();
      std::sort(D_diag_sorted.data(), D_diag_sorted.data() + D_diag_sorted.size());
      
      threshold_Db = -D_diag_sorted(m_eff-1);
      }
    }

    active_i.clear();
    passive_i.clear();
    for (int j = 0; j < p; ++j) {
      if (ranking(j) >= threshold_Db || D_diag(j) <= threshold ){
        active_i.push_back(j);
      } else{
        passive_i.push_back(j);
      }
    }
    S = static_cast<int>(active_i.size());

         //Rcpp::Rcout << "a (iter" << i << ") \n";
         

    // Sample the data batch to use to update beta
    batch_i.clear();
    std::sample(samples_vec.begin(), samples_vec.end(), std::back_inserter(batch_i),
                batch_size, engine);

         //Rcpp::Rcout << "b (iter" << i << ") \n";

    if (algoBhatt) {
      // Restrict yourself to batch and define outcome according to family
      if(family=='b'){
        sqrt_omega=avg_omega(batch_i).sqrt();
        batch_X=X(batch_i,active_i).array().colwise() * sqrt_omega;
        batch_Xt=batch_X.transpose();

        z=y(batch_i).array()/sqrt_omega;
      } else{
        batch_Xt=Xt(active_i,batch_i).eval();
        batch_X=X(batch_i,active_i).eval();

        z=y(batch_i);
      }

         //Rcpp::Rcout << "c (iter" << i << ") \n";

      // Compute D^{-1}X^t where beta=solve(X^tX+D,X^ty), but adjust it by n/batch_size
      U = (batch_Xt.array().colwise() * (n/batch_size/D_diag(active_i))).eval();

      // Take LLt cholesky decomposition
      L=(batch_X * U + Eigen::MatrixXd::Identity(batch_size,batch_size)).llt();
      // Compute mu_beta
      VarInvUt=L.solve(U.transpose()); //Not the quickest way but VarInvUt needed anyhow
      new_mu(active_i)=(z.transpose()*VarInvUt).transpose().array();
      // new_mu(active_i)=(U*L.solve(z)).array(); //This is quicker but VarInvUt needed anyhow

         //Rcpp::Rcout << "d (iter" << i << ") \n";

      if(S<p){
        if(family=='b'){
          // Make poor assumption of omega_i all equal
          new_mu(passive_i)=(Xty(passive_i)-XtX(passive_i,active_i)*(avg_avg_omega*mu_beta(active_i)).matrix()).array()/(avg_avg_omega*XtX_diag(passive_i)+D_diag(passive_i)); // Not quick but for now;
        } else{
          new_mu(passive_i)=(Xty(passive_i)-XtX(passive_i,active_i)*mu_beta(active_i).matrix()).array()/(XtX_diag(passive_i)+D_diag(passive_i)); // Not quick but for now
        }
      }
      grad_mu=mu_beta-new_mu;
      m_mu=decay*m_mu+(1-decay)*grad_mu;
      mu_beta=mu_beta-step_size*m_mu;
    }
    else { ///!!! Not needed, with batching this will never happen
      // Take LLt cholesky decomposition
      L=(XtX(active_i,active_i)+D_diag(active_i).matrix().asDiagonal().toDenseMatrix()).llt();

      //Compute mu_beta
      mu_beta(active_i)=mu_beta(active_i)+lrate*L.solve(Xty(active_i)).array();
    }
         //Rcpp::Rcout << "e (iter" << i << ") \n";
    // Use mu_beta to compute variance (adjusted by joint distribution with sigma2)

    //The adjusting factor is precisely the avg of sigma^(-2)
    if(update_sigma){
      if(family == 's'){
        grad_sigma2=scale_sigma2-(norm_y-Xty.dot(mu_beta.matrix())+ycensvar(censored).sum())/2;
      } else{
        grad_sigma2=scale_sigma2-(norm_y-Xty.dot(mu_beta.matrix()))/2;
      }
      m_sigma2=decay*m_sigma2+(1-decay)*grad_sigma2;
      scale_sigma2=scale_sigma2-step_size*m_sigma2;

      avg_sigma2inv=shape_sigma2/scale_sigma2;
    }
    
         //Rcpp::Rcout << "f (iter" << i << ") \n";

    // Compute diagonal of variance of beta (but without sigma^2)
     if(S < p){
       if(family=='b'){
         // Use poor assumption of omega_i all equal
         new_Sigma_diag(passive_i)=
           1/(D_diag(passive_i)+avg_avg_omega*XtX_diag(passive_i));

         estimated_dvar(passive_i)=avg_avg_omega*XtX_diag(passive_i)+D_diag(passive_i)-1/Sigma_diag(passive_i);
       } else{
         new_Sigma_diag(passive_i)=
           1/(D_diag(passive_i)+XtX_diag(passive_i));

         
         estimated_dvar(passive_i)=XtX_diag(passive_i)+D_diag(passive_i)-1/Sigma_diag(passive_i);
       }
       
       estimated_dvar(passive_i)=estimated_dvar(passive_i).unaryExpr([](double v){
         return std::max(v,0.0);
       });

       new_Sigma_diag(passive_i)=new_Sigma_diag(passive_i)+
         new_Sigma_diag(passive_i).square()*estimated_dvar(passive_i);
     }
     
              //Rcpp::Rcout << "g (iter" << i << ") \n";
     
    if(algoBhatt){
      new_Sigma_diag(active_i) = (1/D_diag(active_i)-batch_size/n*
        (VarInvUt.transpose().array()*U.array()).matrix().rowwise().sum().array());
    }
    else{
      Sigma_diag(active_i)=L.solve(Eigen::MatrixXd::Identity(S,S)).diagonal(); // !!! not super efficient, inverting the matrix twice! Should store this so only inverted once
    }

         //Rcpp::Rcout << "h (iter" << i << ") \n";

//    if(invert_sigma){
//    grad_var=1/Sigma_diag-1/new_Sigma_diag;
//    m_var=decay*m_var+(1-decay)*grad_var;
//    Sigma_diag=1/((1/Sigma_diag)-step_size*m_var);
//    }else{
//    grad_var=Sigma_diag-new_Sigma_diag;
//    m_var=decay*m_var+(1-decay)*grad_var;
//   Sigma_diag=Sigma_diag-step_size*m_var;
//    }
    
    if(invert_sigma){
    grad_var=1/Sigma_diag-1/new_Sigma_diag;
    m_var=decay*m_var+(1-decay)*grad_var;
    Sigma_diag=1/((1/Sigma_diag)-step_size*m_var);
    }else{
    grad_var=Sigma_diag-new_Sigma_diag;
    m_var=decay*m_var+(1-decay)*grad_var;
    Sigma_diag=Sigma_diag-step_size*m_var;
    }

    // Update quantities dependent on beta and gamma parameters
    avg_beta2_sigma2inv = avg_sigma2inv*(mu_beta.square())+Sigma_diag;
    y_hat=X*mu_beta.matrix();

         //Rcpp::Rcout << "i (iter" << i << ") \n";

    // Update outcomes for next iteration if GLM
    if(family=='b'){
      // Compute omega for binary case (use old batch_i and VarInvUt)
       first_guess = (VarInvUt.array()*batch_X.array()).rowwise().sum();
       target_sum = std::min(n-1.0,std::max(1.0,0.0+p - (Sigma_diag*D_diag).sum())); //*batch_size/n; !!! // Using fact that this is equal to degrees of freedom
       ///!!! I never use target sum!
       if(S<p){
       estimated_vars = (X(batch_i,passive_i).array().square().rowwise()*Sigma_diag(passive_i).transpose()).rowwise().sum()+(VarInvUt.array()*batch_X.array()).rowwise().sum()/avg_omega(batch_i)*(batch_size/n);
       } else{
       estimated_vars = (VarInvUt.array()*batch_X.array()).rowwise().sum()/avg_omega(batch_i)*(batch_size/n);
       }

       grad_omega(batch_i)=scale_omega(batch_i)-(y_hat(batch_i).square()+
         estimated_vars).sqrt();
       m_omega(batch_i)=decay*m_omega(batch_i)+(1-decay)*grad_omega(batch_i);
       scale_omega(batch_i)=scale_omega(batch_i)-step_size*m_omega(batch_i);
       scale_omega_exp(batch_i)=scale_omega(batch_i).exp();

       avg_omega(batch_i)=-(2/(scale_omega_exp(batch_i)+1)-1)/scale_omega(batch_i)*(shape_omega/2);
       avg_avg_omega=avg_omega.sum()/n;

    } else if(family=='s'){
      // Compute censored outcomes for survival case
      trunc_alpha=(logtimes(censored)-y_hat(censored))*std::sqrt(avg_sigma2inv);
      Z_phi=trunc_alpha.square()/2+std::log(M_PI)/2+std::log(2)/2+
        trunc_alpha.unaryExpr([](double v){
        return pnorm_cpp(v,0,1,false,true);
        });
      grad_ycens = y(censored).array() - (y_hat(censored) + (-Z_phi).exp()/std::sqrt(avg_sigma2inv));
      m_ycens=decay*m_ycens+(1-decay)*grad_ycens;
      y(censored)=y(censored).array()-step_size*m_ycens;
      
      grad_ycensvar=ycensvar(censored)-1/avg_sigma2inv*(1+trunc_alpha*(-Z_phi).exp()-(-2*Z_phi).exp());
      m_ycensvar=decay*m_ycensvar+(1-decay)*grad_ycensvar;
      ycensvar(censored)=ycensvar(censored).array()-step_size*m_ycensvar;
      
      // y changes so update Xty and norm_y (!!! better version: uncensored only needs to be computed once)
      Xty=Xt*y;
      norm_y=y.squaredNorm();
    }
    
    // !!! Need to include this in previous block cause it's repetitive for binary
    if(S<p){
       estimated_vars = (X(batch_i,passive_i).array().square().rowwise()*Sigma_diag(passive_i).transpose()).rowwise().sum()+(VarInvUt.array()*batch_X.array()).rowwise().sum()*(batch_size/n);
       } else{
       estimated_vars = (VarInvUt.array()*batch_X.array()).rowwise().sum()*(batch_size/n);
    }
    
    grad_vvec(batch_i)=v_vec(batch_i) - estimated_vars;
    m_vvec(batch_i)=decay*m_vvec(batch_i)+(1-decay)*grad_vvec(batch_i);
    v_vec(batch_i)=v_vec(batch_i).array()-step_size*m_vvec(batch_i);

         //Rcpp::Rcout << "l (iter" << i << ") \n";

    // Compute ELBO
    if(algoBhatt){
      SqrtLogDetSigma=D_diag(active_i).log().sum()/2-L.matrixL().toDenseMatrix().diagonal().array().log().sum();
      if(S<p){
      SqrtLogDetSigma=SqrtLogDetSigma+Sigma_diag(passive_i).log().sum()/2;
      }
    } else{
      SqrtLogDetSigma=-L.matrixL().toDenseMatrix().diagonal().array().log().sum();
    }

         //Rcpp::Rcout << "m (iter" << i << ") \n";

    log_mass_beta_sigma=n/2*log(avg_sigma2inv)+SqrtLogDetSigma;
    log_mass_lambdas=-expint_lambda.log().sum();

    elbo_new=//likelihood simplifies
      //-(avg_beta2_sigma2inv*avg_lambda2inv*prior_lambda2inv*avg_tau2inv_all).sum()/2+ //beta
      //lambda simplifies
      //tau simplifies
      //tau_lin simplifies
      //tau_rules simplifies
      //sigma simplifies
      //+(D_diag*avg_beta2_sigma2inv).sum()/2+log_mass_beta_sigma+//beta&sigma
      log_mass_beta_sigma+//beta&sigma
      +(scale_lambda2*avg_lambda2inv).sum()+log_mass_lambdas+//lambda
      log_mass_tau_all+ //tau, all of them
      scale_tauL2*avg_tauL2inv+
      scale_tauR2*avg_tauR2inv;
      
     if(family=='s'){
     elbo_new=elbo_new+n_cens*log(avg_sigma2inv)/2+trunc_alpha.unaryExpr([](double v){
        return pnorm_cpp(v,0,1,false,true);
        }).sum(); // already computed above but for now lazy !!!
     }
     
     if(i==0){
     elbo=elbo_new;
     } else {
     elbo=(1-std::pow(lrate,2))*elbo+(1-std::pow(lrate,2))*elbo_new;
     }
     if(i==ping){
       first_Delbo=elbo-elbo_old;
     }
     
     if(i % ping == 0){
       if(std::abs(elbo-elbo_old)<eps*first_Delbo){
         Rcpp::Rcout << "Converged at iteration nr. " << i << "\n";
       }
       elbo_old=elbo;
     }
     

    if(i % 100 == 0 && i < max_iter){
      
      if(family=='b'){
      // Rescale log-odds correctly
      fit = glm_func(Rcpp::Named("x") = y_hat, Rcpp::Named("y") = y.array()+0.5,
      Rcpp::Named("family") = family_obj, Rcpp::Named("intercept") = false);
      // double coeff = decay+(1-decay)*std::abs(Rcpp::as<Eigen::ArrayXd>(fit["coefficients"])(0));
      // double coeff = std::max(0.5,std::abs(Rcpp::as<Eigen::ArrayXd>(fit["coefficients"])(0)));
      
      // grad_coeff=coeff-std::abs(Rcpp::as<Eigen::ArrayXd>(fit["coefficients"])(0));
      // m_coeff=decay*m_coeff+(1-decay)*grad_coeff;
      // coeff=coeff-step_size*m_coeff;
      coeff = Rcpp::as<Eigen::ArrayXd>(fit["coefficients"])(0);
      coeff = decay+std::abs(coeff)*(1-decay);
      
      // 5. Rescale everything
      mu_beta = mu_beta * coeff;
      Sigma_diag = Sigma_diag * (coeff*coeff);
      avg_tauL2inv = avg_tauL2inv / coeff;
      avg_tauR2inv = avg_tauR2inv / coeff;
      avg_lambda2inv = avg_lambda2inv / coeff;
      avg_omega = avg_omega/coeff;
      avg_avg_omega=avg_avg_omega/coeff;
      scale_omega = scale_omega * coeff;
      scale_tauL2 = scale_tauL2 * coeff;
      scale_tauR2 = scale_tauR2 * coeff;
      scale_lambda2 = scale_lambda2 * coeff;
      }
      }

         //Rcpp::Rcout << "n (iter" << i << ") \n";
         
         if(i==1){
         cpo1=((1-v_vec).log().sum()/2-
                                                n*shape_sigma2*std::log(scale_sigma2)+
                                                (shape_sigma2-0.5)*((shape_sigma2-(y.array()-y_hat).square()/(2*(1-v_vec))).log().sum()));
    
         cpo2=(n/batch_size*(1-estimated_vars).log().sum()/2-
                    n*shape_sigma2*std::log(scale_sigma2)+
                    n/batch_size*(shape_sigma2-0.5)*((shape_sigma2-(y(batch_i).array()-y_hat(batch_i)).square()/(2*(1-estimated_vars))).log().sum()));
                    
                    
         } else{
         grad_cpo1=cpo1-((1-v_vec).log().sum()/2-
                                                n*shape_sigma2*std::log(scale_sigma2)+
                                                (shape_sigma2-0.5)*((shape_sigma2-(y.array()-y_hat).square()/(2*(1-v_vec))).log().sum()));
      m_cpo1=0*m_cpo1+(1-0)*grad_cpo1;
      cpo1=cpo1-(5e-4)*m_cpo1;
    
    grad_cpo2=cpo2-(n/batch_size*(1-estimated_vars).log().sum()/2-
                    n*shape_sigma2*std::log(scale_sigma2)+
                    n/batch_size*(shape_sigma2-0.5)*((shape_sigma2-(y(batch_i).array()-y_hat(batch_i)).square()/(2*(1-estimated_vars))).log().sum()));
    m_cpo2=(1-std::pow(1-decay,2))*m_cpo2+(std::pow(1-decay,2))*grad_cpo2;
    cpo2=cpo2-std::pow(step_size,2)*m_cpo2;
         }
    
        cpo(i) = cpo1;

    if(i % ping == 0){
      // Calculate global shrinkages individually
      meandraws_tau=(((zdraws-ydraws)-ydraws.sqrt()*ydraws*(zdraws-1)+
      zdraws.sqrt()*zdraws*(ydraws-1))/
        ((ydraws-1)*(zdraws-1)*(zdraws-ydraws))).sum()/n_draws;

      avg_tau2inv=meandraws_tau/meandraws_denom;
      
      // Note: density is written after change of coordinates!
    
      logdens = -((y.array()-y_hat).square().sum()*avg_sigma2inv+
                  (D_diag*(mu_beta.array().square())).sum()+
                  D_diag.log().sum()+
                  (n+p+1)*std::log(avg_sigma2inv))/2-
                  (1/avg_lambda2inv+1).log().sum()-
                  std::log(1/avg_tau2inv+1)-
                  std::log(avg_tau2inv/avg_tauL2inv+1)-
                  std::log(avg_tau2inv/avg_tauR2inv+1)-
                  //From here on: adjusting for variable change in density
                  1.5*(std::log(avg_sigma2inv)+
                  avg_lambda2inv.log().sum()+
                  std::log(avg_tauL2inv)+
                  std::log(avg_tauR2inv))-
                  std::log(avg_tau2inv)/2;
    
    
    
    if(verbose){
      Rcpp::Rcout << "ELBO after iteration nr. " << i << ": " << elbo << "\n";
      Rcpp::Rcout << "logdensity evaluated at iteration nr. " << i << ": " << logdens << "\n";
      //Rcpp::Rcout << "CPO evaluated at iteration nr. " << i << ": " << cpo(static_cast<int>(std::floor(i/ping))) << "\n";
      //Rcpp::Rcout << "CPO1 evaluated at iteration nr. " << i << ": " << cpo1 << "\n";
      //Rcpp::Rcout << "CPO2 evaluated at iteration nr. " << i << ": " << cpo2 << "\n";
      //Rcpp::Rcout << "CPO1 momentum evaluated at iteration nr. " << i << ": " << m_cpo1 << "\n";
      //Rcpp::Rcout << "CPO2 momentum evaluated at iteration nr. " << i << ": " << m_cpo2 << "\n";
      //Rcpp::Rcout << "CPO1 difference from last ping evaluated at iteration nr. " << i << ": " << std::abs(cpo1_old-cpo1) << "\n";
      //Rcpp::Rcout << "CPO2 difference from last ping evaluated at iteration nr. " << i << ": " << std::abs(cpo2_old-cpo2) << "\n";
      Rcpp::Rcout << "Betas updated at iteration nr. " << i << ": " << S << "\n";
      Rcpp::Rcout << "tauL2inv at iteration nr. " << i << ": " << avg_tauL2inv << "\n";
      Rcpp::Rcout << "tauRinv at iteration nr. " << i << ": " << avg_tauR2inv << "\n";
      if(family == 'b'){
        Rcpp::Rcout << "average omega at iteration nr. " << i << ": " << avg_avg_omega << "\n";
        Rcpp::Rcout << "logodds-rescale at iteration nr. " << i << ": " << coeff << "\n";
      } else{
        Rcpp::Rcout << "sigma2inv at iteration nr. " << i << ": " << avg_sigma2inv << "\n";
      }
      if(family=='s'){
        Rcpp::Rcout << "t_hat before censoring: " << (logtimes(censored).array() > y_hat(censored)).count() << "\n";
      }
    }
    cpo1_old=cpo1;
    cpo2_old=cpo2;
    
    }
    
    
    



  }
  
  
  // !!! This needs to change for binary setting cause then you're not using student t anymore
  Eigen::ArrayXd var_beta=Sigma_diag/avg_sigma2inv*(1+2/(n-2));
  
  Eigen::MatrixXd BetaPrec = (XtX+D_diag.matrix().asDiagonal().toDenseMatrix()).array()*avg_sigma2inv;

  return Rcpp::List::create(Rcpp::Named("BetaHat") = mu_beta,
                            Rcpp::Named("BetaVar") = var_beta,
                            Rcpp::Named("TauL2InvHat") = avg_tauL2inv,
                            Rcpp::Named("TauR2InvHat") = avg_tauR2inv,
                            Rcpp::Named("Sigma2InvHat") = avg_sigma2inv,
                            Rcpp::Named("Lambda2InvHat") = avg_lambda2inv,
                            Rcpp::Named("BetaPrec") = BetaPrec,
                            Rcpp::Named("CPO") = cpo);
}

