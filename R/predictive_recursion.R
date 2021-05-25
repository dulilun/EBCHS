#' Predictive recursion by Newton (2002)
#'
#' @param x a sequence of chi-squared test statistics
#' @param k degrees of freedom
#'
#' @return  a list: null proportion, prior probability, and lambda-mesh values
#' @examples
#' set.seed(2021)
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
#' lambda[(p_0+1):p] = stats::rgamma(p_1, shape = alpha, rate=1/beta)
#' # Generate a Poisson RV
#' J = sapply(1:p, function(x){rpois(1, lambda[x]/2)})
#' X = sapply(1:p, function(x){rchisq(1, k+2*J[x])})
#' out = predictive_recursion(X, k)
#'
#' @export
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
