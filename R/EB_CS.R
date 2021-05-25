#' Main function used in the paper (Du and Hu, 2020+)
#'
#' @param x a sequence of chi-squared test statistics
#' @param df the degrees of freedom
#' @param qq the quantiles used in spline basis
#' @param method LS: parametric least-squares; PLS: penalized least-squares; g-model: g-modeling
#' @param mixture FALSE: there is no a point mass at zero
#'
#' @return a list: posterior mean, variance, and skewness estimates
#' @examples
#' set.seed(2020)
#' p = 1000
#' k = 7
#' # the prior distribution for lambda
#' alpha = 2
#' beta =  10
#' # lambda
#' lambda = rep(0, p)
#' pi_0 = 0
#' p_0 = floor(p*pi_0)
#' p_1 = p-p_0
#' lambda[(p_0+1):p] = rgamma(p_1, shape = alpha, rate=1/beta)
#' # Generate a Poisson RV
#' J = sapply(1:p, function(x){rpois(1, lambda[x]/2)})
#' X = sapply(1:p, function(x){rchisq(1, k+2*J[x])})
#' qq_set = seq(0.01, 0.99, 0.01)
#' out = EB_CS(X, k, qq=qq_set, method='LS', mixture = FALSE)
#' E = out$E_lambda
#' V = out$V_lambda
#' UP = E+1.645*V^(1/2)
#' LOW = E-1.645*V^(1/2)
#'
#' @export
EB_CS <- function(x,
                  df,
                  qq=c(0.2, 0.4, 0.6, 0.8),
                  method=c('LS', 'PLS', 'g_model'),
                  mixture=FALSE){

  # tweedies' formula for posterior mean
  MF <- function(x, l_1, l_2){
    g_k_2_g_est = 1+2*l_1
    g_k_4_g_est = 4*l_2+(1+2*l_1)^2
    E_J_x_est = x/2*g_k_4_g_est/g_k_2_g_est-(k-4)/2
    E_lambda_est = 2*E_J_x_est*(1+2*l_1)
    return( E_lambda_est )
  }
  # Posterior variance
  VF <- function(x, l_1, l_2, l_3, l_4){
    g_k_2_g_est = 1+2*l_1
    g_k_4_g_est = 4*l_2+(1+2*l_1)^2
    g_k_6_g_est = 8*l_3+12*l_2*(1+2*l_1)+(1+2*l_1)^3
    g_k_8_g_est = 16*l_4+32*l_3*(1+2*l_1)+24*l_2*(1+2*l_1)^2+48*l_2^2+(1+2*l_1)^4

    E_J_x_est = x/2*g_k_4_g_est/g_k_2_g_est-(k-4)/2
    E_J2_x_est = x^2/4*g_k_8_g_est/g_k_4_g_est-(k-6)/2*x*g_k_6_g_est/g_k_4_g_est+(k-6)*(k-4)/4
    V_lambda_est = 16*E_J2_x_est*l_2+(4*E_J2_x_est-4*E_J_x_est^2)*(1+2*l_1)^2
    V_lambda_est = pmax(V_lambda_est, 0)

    return(V_lambda_est)
  }
  # the possible value for Poisson
  J_est <- function(x, l_1){
    g_k_2_g_est = 1+2*l_1
    J_hat = x/2*g_k_2_g_est-(k-2)/2
    return(J_hat)
  }

  # parametric approach: LS
  k = df
  p_out = density_LS(x)
  l_1p = p_out$l_1
  l_2p = p_out$l_2
  l_3p = p_out$l_3
  l_4p = p_out$l_4

  # default method
  method = method[1]
  if(method=='LS'){
    E_lambda_est = MF(x, l_1p, l_2p)
    V_lambda_est = VF(x, l_1p, l_2p, l_3p, l_4p)

    J_hat = J_est(x, l_1p)
    J_hat = pmax(J_hat, 0)
    S_lambda_est <- 2/sqrt(J_hat+0.5)

  }else if(method=='PLS'){
    # Estimate the density derivatives
    D_est = density_PLS(x, qq)
    l_1 = D_est$l_1
    l_2 = D_est$l_2
    E_lambda_est = MF(x, l_1, l_2)
    V_lambda_est = VF(x, l_1, l_2, l_3p, l_4p)

    J_hat = J_est(x, l_1)
    J_hat = pmax(J_hat, 0)
    S_lambda_est <- 2/sqrt(J_hat+0.5)

  }else if(method=='g_model'){
    out_pr = predictive_recursion(x, k)
    pi_0 = out_pr$pi_0
    g_prior = out_pr$g_prior
    lambda_set = out_pr$lambda_set

    den = density_g_model(x, k, pi_0, lambda_set, g_prior)
    g_k = den$g_k
    l_1 = den$l_1
    l_2 = den$l_2
    l_3 = den$l_3
    l_4 = den$l_4
    E_lambda_est = MF(x, l_1, l_2)
    V_lambda_est = VF(x, l_1, l_2, l_3, l_4)

    J_hat = J_est(x, l_1)
    J_hat = pmax(J_hat, 0)
    S_lambda_est <- 2/sqrt(J_hat+0.5)

  }else{ stop('Incorrect method') }

  if(mixture){

    # prior and pi_0
    if(method=='g_model'){
      g_k = g_k
    }else{
      out_pr = predictive_recursion(x, k)
      pi_0 = out_pr$pi_0
      g_prior = out_pr$g_prior
      lambda_set = out_pr$lambda_set
      DD = sapply(1:length(g_prior), function(i){stats::dchisq(x, df=k, ncp=lambda_set[i])})
      g_k = pi_0*stats::dchisq(x, df=k)+(1-pi_0)*DD%*%g_prior
    }
    # the estimated lfdr
    lfdr = pi_0*stats::dchisq(x, k)/g_k


    # posterior mean and variance under the alternative
    E_lambda_a_est = E_lambda_est/(1-lfdr)
    V_lambda_a_est = V_lambda_est/(1-lfdr)-lfdr*E_lambda_a_est^2
    V_lambda_a_est = pmax(V_lambda_a_est, 0)

    E_lambda_est = E_lambda_a_est
    V_lambda_est = V_lambda_a_est

    J_hat = J_hat*(1-lfdr)
    J_hat = pmax(J_hat, 0)
    S_lambda_est <- 2/sqrt(J_hat+0.5)
  }

  # output
  out = list(E_lambda = E_lambda_est, V_lambda = V_lambda_est, S_lambda = S_lambda_est)
  return(out)
}
