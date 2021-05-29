#' Main function used in the paper (Du and Hu, 2020)
#'
#' Give a sequence of chi-squared statistic values, the function computes the posterior mean, variance, and skewness
#' of the noncentrality parameter given the data.
#'
#' @param x a sequence of chi-squared test statistics
#' @param df the degrees of freedom
#' @param qq the quantiles used in spline basis
#' @param method LS: parametric least-squares; PLS: penalized least-squares; g-model: g-modeling
#' @param mixture default is FALSE: there is no point mass at zero.
#'
#' @return a list: posterior mean, variance, and skewness estimates
#'
#' @references Du and Hu (2020), \emph{An Empirical Bayes Method for Chi-Squared Data}, \emph{Journal of American Statistical Association}, forthcoming.
#'
#' @examples
#' p = 1000
#' k = 7
#' # the prior distribution for lambda
#' alpha = 2
#' beta =  10
#' # lambda
#' lambda = rep(0, p)
#' pi_0 = 0.8
#' p_0 = floor(p*pi_0)
#' p_1 = p-p_0
#' lambda[(p_0+1):p] = rgamma(p_1, shape = alpha, rate=1/beta)
#' # Generate a Poisson RV
#' J = sapply(1:p, function(x){rpois(1, lambda[x]/2)})
#' X = sapply(1:p, function(x){rchisq(1, k+2*J[x])})
#' qq_set = seq(0.01, 0.99, 0.01)
#' out = EB_CS(X, k, qq=qq_set, method='LS', mixture = TRUE)
#' E = out$E_lambda
#' V = out$V_lambda
#' S = out$S_lambda
#'
#' @export
EB_CS <- function(x,
                  df,
                  qq=c(0.2, 0.4, 0.6, 0.8),
                  method=c('LS', 'PLS', 'g_model'),
                  mixture=FALSE){
  # Predictive recursion by Newton (2002)
  # used to learn the prior distribution by discretization
  predictive_recursion <- function(x, k){

    p = length(x)
    gamma = (1+1:p)^(-0.67)
    # null and non-null proportions
    pi_0 = sum(x<2*k)/(p*stats::pchisq(2*k, df=k))
    # the grid for theta: integration
    grid_size = 0.2
    theta = seq(0, max(x), grid_size)
    # the prior function under alternative
    pi_1_theta = (1-pi_0)*stats::dgamma(theta, shape = 3, rate = 1/5)

    # repeat the experiment 10 times
    pi_0_final = c()
    pi_theta_final = list()
    n_sim = 10
    for(sim in 1:n_sim){
      # re-sample the data
      z = sample(x, length(x))
      for(i in 1:p){
        m_0 = pi_0*stats::dchisq(z[i], df=k)
        f_1_theta = pi_1_theta*stats::dchisq(z[i], df=k, ncp = theta)
        m_1 = pracma::trapz(theta, f_1_theta)
        # update pi_0 and pi_1_theta
        pi_0 = (1-gamma[i])*pi_0+gamma[i]*m_0/(m_0+m_1)
        pi_1_theta = (1-gamma[i])*pi_1_theta+gamma[i]*f_1_theta/(m_0+m_1)
      }
      # output
      pi_1 = 1-pi_0
      pi_theta = pi_1_theta/pi_1*grid_size
      pi_0_final =c(pi_0_final, pi_0)
      pi_theta_final[[sim]] = pi_theta
    }

    pi_0_final = mean(pi_0_final)
    pi_theta_final = Reduce('+', pi_theta_final)
    pi_theta_final = pi_theta_final/n_sim
    out = list(pi_0=pi_0_final, g_prior = pi_theta_final, lambda_set = theta)
    return(out)
  }

  # the g-modeling method in Efron (2016)
  density_g_model <- function(x, k, pi_0, lambda_set, g_prior){

    # dir: the order of derivative: chi-squared distribution
    chi2pdf_dir <- function(x, k, di){

      # Function: chi-squared density
      chi2pdf_du = function(x, k){
        if(k<=0 & k%%2 == 0){
          out = 0
        }else{
          out = x^(k/2-1)*exp(-x/2)/(2^(k/2)*gamma(k/2))
        }
        return(out)
      }

      # order
      if(di==0){
        g_0 = chi2pdf_du(x, k)
      }else if(di==1){
        g_0 = 1/2*chi2pdf_du(x, k-2)-1/2*chi2pdf_du(x, k)
      }else if(di==2){
        g_0 = 1/4*chi2pdf_du(x, k-4)-2/4*chi2pdf_du(x, k-2)+1/4*chi2pdf_du(x, k)
      }else if(di==3){
        g_0 = 1/8*chi2pdf_du(x, k-6)-3/8*chi2pdf_du(x, k-4)+
          3/8*chi2pdf_du(x, k-2)-1/8*chi2pdf_du(x, k)
      }else if(di == 4){
        g_0 = 1/16*chi2pdf_du(x, k-8)-4/16*chi2pdf_du(x, k-6)+6/16*chi2pdf_du(x, k-4)-
          4/16*chi2pdf_du(x, k-2)+1/16*chi2pdf_du(x, k)
      }
      g_0 = as.vector(g_0)
      return(g_0)
    }

    # the derivative of the non-central chi-squared density
    n_chi2pdf_dir <- function(x, k, lambda, di){
      N = 100
      out = sapply(0:N, function(j){
        chi2pdf_dir(x, k+2*j, di)*exp(-lambda/2)*(lambda/2)^j/factorial(j) })
      if (length(x)==1){
        out = unlist(out)
        ss = sum(out)
      }else{
        ss = apply(out, 1, sum)
      }
      ss = as.vector(ss)
      return(ss)
    }

    # density estimate
    DD = sapply(1:length(g_prior), function(i){n_chi2pdf_dir(x, k, lambda_set[i], 0)})
    g_k = pi_0*chi2pdf_dir(x, k, 0)+(1-pi_0)*DD%*%g_prior
    # first derivative
    DD_1 = sapply(1:length(g_prior), function(i){n_chi2pdf_dir(x, k, lambda_set[i], 1)})
    g_1 = pi_0*chi2pdf_dir(x, k, 1)+(1-pi_0)*DD_1%*%g_prior
    g_1 = g_1/g_k
    # second derivative
    DD_2 = sapply(1:length(g_prior), function(i){n_chi2pdf_dir(x, k, lambda_set[i], 2)})
    g_2 = pi_0*chi2pdf_dir(x, k, 2)+(1-pi_0)*DD_2%*%g_prior
    g_2 = g_2/g_k
    # third derivative
    DD_3 = sapply(1:length(g_prior), function(i){n_chi2pdf_dir(x, k, lambda_set[i], 3)})
    g_3 = pi_0*chi2pdf_dir(x, k, 3)+(1-pi_0)*DD_3%*%g_prior
    g_3 = g_3/g_k
    # 4th derivative
    DD_4 = sapply(1:length(g_prior), function(i){n_chi2pdf_dir(x, k, lambda_set[i], 4)})
    g_4 = pi_0*chi2pdf_dir(x, k, 4)+(1-pi_0)*DD_4%*%g_prior
    g_4 = g_4/g_k

    # from g_k to l_k
    l_1 = g_1
    l_2 = g_2-l_1^2
    l_3 = g_3-(l_2+l_1^2)*l_1-2*l_1*l_2
    l_4 = g_4-g_3*g_1-3*l_1^2*l_2-3*l_2^2-3*l_1*l_3

    den = list()
    den$g_k = g_k
    den$l_1 = l_1
    den$l_2 = l_2
    den$l_3 = l_3
    den$l_4 = l_4
    return(den)
  }

  # Tweedie's formula for posterior mean
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

  # J_hat
  J_est <- function(x, l_1, k){
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
  method = match.arg(method)

  if(method=='LS'){
    E_lambda_est = MF(x, l_1p, l_2p)
    V_lambda_est = VF(x, l_1p, l_2p, l_3p, l_4p)

    J_hat = J_est(x, l_1p, k)
    J_hat = pmax(J_hat, 0)
    S_lambda_est <- 2/sqrt(J_hat+0.5)

  }else if(method=='PLS'){
    # Estimate the density derivatives
    D_est = density_PLS(x, qq)
    l_1 = D_est$l_1
    l_2 = D_est$l_2
    E_lambda_est = MF(x, l_1, l_2)
    V_lambda_est = VF(x, l_1, l_2, l_3p, l_4p)

    J_hat = J_est(x, l_1, k)
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
      pi_0 = pi_0
      g_k = g_k
    }else{
      # other methods also requires g_k
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
