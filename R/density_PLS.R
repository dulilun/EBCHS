#' Penalized least-squares method
#'
#' The semiparametric model is employed to estimate the log density derivatives of
#' the chi-squared statistics.
#'
#' @param x a sequence of chi-squared test statistics
#' @param qq the quantiles used for splines
#'
#' @return a list: the first and second density derivatives
#'
#' @examples
#' p = 1000
#' k = 7
#' # the prior distribution for lambda
#' alpha = 2
#' beta =  10
#' # lambda
#' lambda = rep(0, p)
#' pi_0 = 0.5
#' p_0 = floor(p*pi_0)
#' p_1 = p-p_0
#' lambda[(p_0+1):p] = stats::rgamma(p_1, shape = alpha, rate=1/beta)
#' # Generate a Poisson RV
#' J = sapply(1:p, function(x){rpois(1, lambda[x]/2)})
#' X = sapply(1:p, function(x){rchisq(1, k+2*J[x])})
#' qq = c(0.2, 0.4, 0.6, 0.8)
#' out = density_PLS(X, qq)
#'
#' @export
density_PLS <- function(x, qq){

  # Function : construct the basis for natural cubic splines and its derivatice
  basis_natual_spline <- function(x, knots, Boundary.knots, df){

    # at least df>=4
    # natural cubic splines (df=4)
    Aknots <- sort(c(rep(Boundary.knots, df), knots))
    # create QR decomposition of the vectors in the boundary
    const <- splines::splineDesign(Aknots, Boundary.knots, df, derivs = df-2, outer.ok= TRUE)
    qr.const <- qr( t(const) )

    # the basis for natural cubic spline
    basis <- splines::splineDesign(Aknots, x, df, outer.ok= TRUE)
    basis <- as.matrix( (t(qr.qty(qr.const, t(basis))))[, -(1L:2L), drop = FALSE] )
    # 1st derivs
    basis_1 <- splines::splineDesign(Aknots, x, df, derivs = 1, outer.ok= TRUE)
    basis_1 <- as.matrix( (t(qr.qty(qr.const, t(basis_1))))[, -(1L:2L), drop = FALSE] )
    # 2st derivs
    basis_2 <- splines::splineDesign(Aknots, x, df, derivs = 2, outer.ok= TRUE)
    basis_2 <- as.matrix( (t(qr.qty(qr.const, t(basis_2))))[, -(1L:2L), drop = FALSE] )
    # 3st derivs
    basis_3 <- splines::splineDesign(Aknots, x, df, derivs = 3, outer.ok= TRUE)
    basis_3 <- as.matrix( (t(qr.qty(qr.const, t(basis_3))))[, -(1L:2L), drop = FALSE] )
    # 4st derivs
    if ( df>4 ){
      basis_4 <- splines::splineDesign(Aknots, x, df, derivs = 4, outer.ok= TRUE)
      basis_4 <- as.matrix( (t(qr.qty(qr.const, t(basis_4))))[, -(1L:2L), drop = FALSE] )
    }else{
      basis_4 = matrix(0, dim(basis)[1], dim(basis)[2])
    }


    # creat the penalty matrix
    CBB = fda::create.bspline.basis(breaks = sort(c(knots, Boundary.knots)), norder = df)
    Omega = fda::bsplinepen(CBB, Lfdobj = 2)


    Omega_1 = qr.qty(qr.const, Omega)
    Omega_2 = qr.qty(qr.const, t(Omega_1))
    n_q = dim(Omega_2)[1]
    Omega_3 = Omega_2[3:n_q, 3:n_q]

    return ( list(basis = basis, basis_1 = basis_1, basis_2 = basis_2,
                  basis_3 = basis_3, basis_4 = basis_4,
                  Omega_3=Omega_3) )
  }

  # Function : the response vectors
  h_response = function(x, B, B_1, B_2, B_3, B_4, k, l){
    # x^{l}g^{k}(x)/g(x): take k-th derivative
    xx = matrix(rep(x, dim(B)[2]), dim(B)[1], dim(B)[2])
    if ( l==1 ){

      if ( k == 1 ){
        h = apply(B+xx*B_1, 2, mean)
        h = (-1)^(k)*h
      } else if( k==2 ){
        h = apply(2*B_1+xx*B_2, 2, mean)
        h = (-1)^(k)*h
      } else if( k==3 ){
        h = apply(3*B_2+xx*B_3, 2, mean)
        h = (-1)^(k)*h
      } else if( k==4 ){
        h = apply(4*B_3+xx*B_4, 2, mean)
        h = (-1)^(k)*h
      } else{
        stop("Incorrect k and l")
      }

    } else if (l==2){

      if ( k == 1 ){
        h = apply(2*xx*B+xx^2*B_1, 2, mean)
        h = (-1)^(k)*h
      } else if( k==2 ){
        h = apply(2*B+4*xx*B_1+xx^2*B_2, 2, mean)
        h = (-1)^(k)*h
      } else if( k==3 ){
        h = apply(6*B_1+6*xx*B_2+xx^2*B_3, 2, mean)
        h = (-1)^(k)*h
      } else if( k==4 ){
        h = apply(12*B_2+8*xx*B_3+xx^2*B_4, 2, mean)
        h = (-1)^(k)*h
      } else{
        stop("Incorrect k and l")
      }

    } else if (l==3){

      if (k==3){
        h = apply(6*B+18*xx*B_1+9*xx^2*B_2+xx^3*B_3, 2, mean)
        h = (-1)^(l)*h
      } else{
        stop("Incorrect k and l")
      }

    } else if (l==4){

      if (k==4){
        h = apply(24*B+96*xx*B_1+72*xx^2*B_2+16*xx^3*B_3+x^4*B_4, 2, mean)
        h = (-1)^(l)*h
      } else{
        stop("Incorrect k and l")
      }

    } else{
      stop("Incorrect k and l")
    }

    return(h)
  }

  # Function : log-density gradient-smoothing splines method
  density_d_np <-  function(x, y, knots, bb, df, lambda, k, l){
    # inpute x
    # predict y
    # k: the k-th derivative: x^(l)g(x)^(k)/g(x)

    # Step I: Estimation
    # Method I: local approximation
    BB = basis_natual_spline(x, knots, bb, df)
    B = BB[["basis"]]
    B_1  = BB[["basis_1"]]
    B_2  = BB[["basis_2"]]
    B_3  = BB[["basis_3"]]
    B_4  = BB[["basis_4"]]
    Omega = BB[["Omega_3"]]

    # tuning parameter for ridge regression
    # what is the artificial response h?
    h = h_response(x, B, B_1, B_2, B_3, B_4, k, l)
    G = t(B)%*%B/dim(B)[1]
    beta = solve( G+lambda*Omega, h )

    # Step II: prediction
    # method I: local approximations
    BB = basis_natual_spline(y, knots, bb, df)
    B = BB[["basis"]]
    B_1  = BB[["basis_1"]]
    B_2  = BB[["basis_2"]]
    B_3  = BB[["basis_3"]]
    B_4  = BB[["basis_4"]]
    # Output
    yy = matrix(rep(y, dim(B)[2]), ncol = dim(B)[2])
    l_k = ( B*yy^(-l) )%*%beta
    # likelihood function for CV criterion
    h = h_response(y, B, B_1, B_2, B_3, B_4, k, l)
    LIKE_H = t(beta)%*%( t(B)%*%B/dim(B)[1] )%*%beta-2*t(h)%*%beta

    out = list(LIKE_H = LIKE_H, l_k = l_k)
    return(out)
  }

  # Function 5: CV for Splines to choose lambada
  CV = function(x, kn, bb, df, lambda_set, k, l){

    n_lambda = length(lambda_set)
    LIKE_H = rep(0, n_lambda)
    # prune five points in the bounday
    #x = x[x>x[order(x)][5]]

    for (i in 1:n_lambda){

      lambda = lambda_set[i]
      # partition the data-set into 10-fold
      x_permutate = sample(x, length(x), replace = FALSE)
      # 10-fold CV
      LIKE = 0
      for (j in 1:10){

        n_test = floor(length(x)/10)
        b = 1+(j-1)*n_test
        e = j*n_test
        x_test = x_permutate[b:e]
        x_train = setdiff(x, x_test)
        # SS method
        A = density_d_np(x_train, x_test, kn, bb, df, lambda, k, l)
        # Evaluate the likelihood function for the test data
        LIKE = LIKE + A$LIKE_H
      }

      LIKE_H[i] = LIKE

    }

    Index = which.min( LIKE_H )
    lambda_s = lambda_set[ Index ]
    return(lambda_s)
  }


  ###############################################
  # density estimation
  # the nodes and boundary points used in the splines
  kn = stats::quantile(x, qq)
  bb = range(x)

  # semiparametric estimation for the 1st and 2nd orders
  if(length(qq)<=5){
    lambda_1 = 0
    lambda_2 = 0
  }else{
    lambda_set = seq(1, length(qq), 2)
    lambda_1 = CV(x, kn, bb, df=4, lambda_set, k=1, l=1)
    lambda_2 = CV(x, kn, bb, df=4, lambda_set, k=2, l=2)
  }

  A = density_d_np(x, x, kn, bb, df=4, lambda_1, k=1, l=1)
  g_1 = A$l_k

  A = density_d_np(x, x, kn, bb, df=4, lambda_2, k=2, l=2)
  g_2 = A$l_k

  # convert to (l_1 to l_2)
  l_1 = g_1
  l_2 = g_2-l_1^2


  den = list()
  den$l_1 = l_1
  den$l_2 = l_2
  return(den)
}
