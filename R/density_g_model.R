#' The l_1 to l_4 derivative from the g-modeling method
#'
#' @param x a sequence of chi-squared test statistics
#' @param k degrees of freedom
#' @param pi_0 the proportion of the null
#' @param lambda_set the set of noncentrality values
#' @param g_prior the prior probability for the noncentrality values
#'
#' @return a list: the margianl density, and its first-to-fourth derivatives
#'
#' @export
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
