#' log-density derivatives--parametric approach
#'
#' Assuming the log density of the chi-squared statistics admits a parametric form, this function
#' estimates up to the fourth order log density derivatives.
#'
#' @param x a sequence of chi-squared test statistics
#' @return a list: the first-to-fourth log density derivatives
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
#' lambda[(p_0+1):p] = stats::rgamma(p_1, shape = alpha, rate=1/beta)
#' # Generate a Poisson RV
#' J = sapply(1:p, function(x){rpois(1, lambda[x]/2)})
#' X = sapply(1:p, function(x){rchisq(1, k+2*J[x])})
#' out = density_LS(X)
#'
#' @export
density_LS <-  function(x){

  ONE = rep(1, length(x))
  ZERO = rep(0, length(x))
  B_1 = cbind(ONE, x)
  B_2 = cbind(ZERO, ONE)
  B_3 = cbind(ZERO, ZERO)
  B_4 = cbind(ZERO, ZERO)

  G = t(B_1)%*%B_1/dim(B_1)[1]
  xx = matrix(rep(x, dim(B_1)[2]), ncol=dim(B_1)[2])
  h = (-1)*apply(B_1+B_2*xx, 2, mean)
  beta = solve(G, h)

  # Output
  # the first to fourth derivative estimation
  d_1 = matrix(rep(1/x, dim(B_1)[2]), ncol = dim(B_1)[2])
  d_2 = matrix(rep(1/x^2, dim(B_1)[2]), ncol = dim(B_1)[2])
  d_3 = matrix(rep(1/x^3, dim(B_1)[2]), ncol = dim(B_1)[2])
  d_4 = matrix(rep(1/x^4, dim(B_1)[2]), ncol = dim(B_1)[2])

  l_1 = (B_1*d_1)%*%beta
  l_2 = (B_2*d_1-B_1*d_2)%*%beta
  l_3 = (B_3*d_1-2*B_2*d_2+2*B_1*d_3)%*%beta
  l_4 = (B_4*d_1-3*B_3*d_2+6*B_2*d_3-6*B_1*d_4)%*%beta

  out = list(l_1 = l_1, l_2= l_2, l_3=l_3, l_4 = l_4)
  return(out)
}
