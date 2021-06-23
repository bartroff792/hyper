#' AMO intervals
#'
#' Create alpha max optimal acceptance intervals
#' @param N population size
#' @param n sample size
#' @param a alpha
#' @return AMO intervals
#' @examples 
#' AMO_interval <- AMO(100,10,0.05);
#' @export
AMO <- function(N,n,a){
    mid = as.integer(N/2)
    table = matrix(data = 0, nrow =2, ncol = N + 1,dimnames = NULL)
    for ( M in 0: mid){
      x_min =  max(0,(M+n- N))
      x_max = min(n,N)
      center = as.integer((n+1)*(M+1)/(N+2))
      l = center
      u = center
      cdf = dhyper(center,M,N - M,n)
      if(l > x_min){dl = dhyper(l - 1,M,N - M,n)}else{dl = 0}
      if(u < x_max){du = dhyper(u + 1,M,N - M,n)}else{du = 0}
      while(cdf < 1 - a){
        if (du  > dl){
          cdf =  cdf + du
          u = u + 1
          if(u < x_max){du = dhyper(u + 1,M,N - M,n)}else{du = 0}}
        else{cdf =  cdf + dl
        l = l - 1
        if(l > x_min){dl = dhyper(l - 1,M,N - M,n)}else{dl = 0}
        }}
      table[,(M + 1)] = c(l,u)
    }
    for(M in (mid + 1):N){
      table[1,(M + 1)] = n - table[2,(N - M + 1)]
      table[2,(M + 1)] = n - table[1,(N - M + 1)]
    }
    colnames(table) <- c(0:N)
    rownames(table) <- c("accept_l","accept_u")
    return(table)
  }

#' AO intervals
#'
#' Create alpha optimal acceptance intervals
#' @param N population size
#' @param n sample size
#' @param a alpha
#' @return AO intervals
#' @examples 
#' AO_interval <- AO(100,10,0.05);
#' @export
AO <- function(N,n,a){
    # make sure the acceptance region of M, [c,d] is nondecreasing from 0 to N/2
    mid = as.integer(N/2)
    table1 = matrix(data = 0, nrow =2, ncol = N + 1,dimnames = NULL)
    table = matrix(data = 0, nrow =2, ncol = mid + 1,dimnames = NULL)
    table = AMO(N,n,a)[,1:(mid + 1)]
    l_1 = 0
    u_1 = 0
    table1[,1] = c(0,0)
    for (M in 1:mid){
      l = table[1,(M +1)]
      u = table[2,(M +1)]
      if(l < l_1){
        u_1 = u + l_1 - l
      }else{
        l_1 = l
        u_1 = u
      }
      table1[,(M+1)] = c(l_1,u_1)}
    l_1 = table1[1,(mid +1)]
    u_1 = table1[2,(mid +1)]
    for (M in (mid - 1):0){
      l = table1[1,(M +1)]
      u = table1[2,(M +1)]
      if(u > u_1){
        l_1 = l + u_1 - u
      }else{
        l_1 = l
        u_1 = u
      }
      table1[,(M+1)] = c(l_1,u_1)
      table1[,(N - M + 1)] = c(n -u_1, n- l_1)
    }
    if(N/2 == mid){
      C = qhyper(a/2,N/2,N/2,n)
      D = n - C
      table1[, (mid + 1)] = c( C , D)}else{
        table1[,mid+2] = c(n - table1[2, (mid + 1)],n - table1[1, (mid + 1)])
      }
    colnames(table1) <- c(0:N)
    rownames(table1) <- c("accept_l","accept_u")
    return(table1)
  }


#' CI intervals
#'
#' Create confidence intervals by inverting symmetric alpha optimal acceptance intervals
#' @param N population size
#' @param n sample size
#' @param a alpha (confidence level will be 1 - a)
#' @return CI intervals
#' @examples 
#' CI_interval_table <- CI_interval(100,10,0.05);
#' @export
CI_interval <- function(N,n,a){
    t = AO(N,n,a)
    table = NULL
    for (i in 0:n){
      s1 = which(t[1,] <= i)
      s2 = which(t[2,] >= i)
      s = intersect(s1,s2)
      U = max(s) - 1
      L = min(s) - 1
      table = cbind(table, c(L,U))
    }
    colnames(table) <- c(0:n)
    rownames(table) <- c("L","U")
    return(table)}

#' CI intervals
#'
#' Create a confidence interval for x by inverting symmetric alpha optimal acceptance intervals
#' @param N population size
#' @param n sample size
#' @param a alpha (confidence level will be 1 - a)
#' @param x the number of special items 
#' @return the CI interval for x
#' @examples 
#' CI_interval_table <- single_CI_interval(100,10,0.05, 2);
#' @export
single_CI_interval <- function(N,n,a,x){
  t = CI_interval(N,n,a)
  return(t[,x+1])}

#' coverage probabilities
#'
#' Calculate the set of coverage probabilities
#' @param N population size
#' @param n sample size
#' @param CItable the set of confidence intervals
#' @return coverage
#' @export
coverprob <- function(N,n,CItable){
    t = c(rep(0,N+1))
    for (M in 0:N){
      L = max(n - N + M, 0)
      U = min(n,M)
      for (x in L:U){
        if (M >= CItable[1,x+1] && M <= CItable[2,x+1]){
          t[M+1]= t[M+1] + dhyper(x,M,N-M,n, log = FALSE)
        }
      }
    }
    return(t)
  }


#' coverage probability plot
#'
#' plot the set of coverage probabilities
#' @param N population size
#' @param n sample size
#' @param a level
#' @param CItable the set of confidence intervals
#' @return coverage plot
#' @export
cover_plot <- function(N,n,a,CI_table){
    cover = coverprob(N,n,CI_table)
    plot(c(0:N),cover, main = "Coverage Probability", xlab = "M", ylab = "prob",type = "l")
    abline(h = 1 - a,col= "red", lty= 2, lwd= 3)
  }


