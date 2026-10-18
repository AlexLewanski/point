##################################
### SINGLE LOCATION ESTIMATION ###
##################################

#' Find the product of a list of multivariate normal distributions
#'
#' @param mu_list A list of the distributions' means
#' @param Sigma_list A list of the distributions' variance--covariance matrices (in the same order the means in mu_list).
#'
#' @returns A list containing the mean and variance--covariance of the product, which is a new normal distribution.
#' @export
#'
multiply_mvnorms <- function(mu_list, Sigma_list) {

  d <- length(mu_list[[1]])

  precision_sum <- matrix(0, d, d)
  weighted_mean_sum <- rep(0, d)

  for(i in seq_along(mu_list)) {
    P <- solve(Sigma_list[[i]])

    precision_sum <- precision_sum + P
    weighted_mean_sum <- weighted_mean_sum + P %*% mu_list[[i]]
  }

  Sigma <- solve(precision_sum)
  mu <- Sigma %*% weighted_mean_sum

  list(mu = drop(mu),
       Sigma = Sigma)
}


#' Calculate the variance--covariance matrix for a single tree.
#'
#' @param tree A single tree stored as a phylo object.
#'
#' @returns A variance--covariance matrix of tree's tips.
#' @export
#'
calc_vcv_single_tree <- function(tree) {

  tip_order <- tree$tip.label
  vcv_mat <- ape::vcv.phylo(tree)

  return(vcv_mat[match(tip_order, rownames(vcv_mat)),
                 match(tip_order, colnames(vcv_mat))]
  )
}


#' Creates a design matrix.
#'
#' @param mean_vec A vector of means for the treees, which is assumed to organized where the x and y values iterate.
#'
#' @returns A design matrix
#' @export
#'
create_design <- function(mean_vec) {

  cbind(
    rep(c(1, 0), times = length(mean_vec)/2),
    rep(c(0, 1), times = length(mean_vec)/2)
  )

}


#' Weighted least-squares function to calculate ML estimate of location.
#'
#' @param A A design matrix.
#' @param Sigma_list List of variance--covariance matrices.
#' @param mu List of mean values (in same order as Sigma_list).
#'
#' @returns Estimate of coordinates.
#' @export
#'
weighted_least_squares <- function(A, Sigma_list, mu) {
  sigma_inv <- as.matrix(Matrix::bdiag(lapply(Sigma_list, solve)))

  elementone <- solve(t(A)%*%sigma_inv%*%A)
  elementtwo <- t(A)%*%sigma_inv%*%mu

  return(elementone %*% elementtwo)
}


#' Condition a normal distribution.
#'
#' @param mu mean of normal distribution
#' @param V Variance--covariance matrix of normal distribution
#' @param A Matrix encoding equality constraints (double check this).
#'
#' @returns Conditional distribution.
#' @export
#'
condition_mvn <- function(mu, V, A) {

  inv_AVAt <- solve(A %*% V %*% t(A))

  mu_c <- mu - V%*%t(A)%*%inv_AVAt%*%(A%*%mu)
  var_c <- V - V%*%t(A)%*%inv_AVAt%*%(A%*%V)

  list(
    mean = mu_c,
    covariance = var_c
  )
}
