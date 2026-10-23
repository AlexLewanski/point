#' Create a design matrix for multiple groups based on node membership. (This is an internal function)
#'
#' @param max_group The number of groups being considered.
#' @param node_membership I think this is a list of node memberships? Update this description.
#'
#' @returns The design matrix.
#' @export
#'
multigroup_design <- function(max_group, node_membership) {
  design_mat <- matrix(0,
                       nrow = length(node_membership)*2,
                       ncol = max_group*2)

  #based on the group membership
  #design matrix is organized as: x1, y1, x2, y2, x3, y3, x4, y4
  groupind <- unlist(lapply(node_membership, function(x) c((x - 1)*2 + 1, (x - 1)*2 + 2)))

  #design_mat[cbind(1:(length(node_membership)*2), groupind)] <- 1
  design_mat[cbind(seq_len(length(node_membership)*2), groupind)] <- 1

  return(design_mat)
}


#' The maximization step of the soft node assignment em.
#'
#' @param param The coordinates to estimate the likelihood of.
#' @param r_list list of responsibilities for each tree.
#' @param cond_mean_list List of conditional means.
#' @param cond_var_list List of conditional variances.
#' @param group_count The number of groups (the K value).
#'
#' @returns Negative log-likelihood of input coordinates.
#' @export
#'
m_step <- function(param, r_list, cond_mean_list, cond_var_list, group_count) {
  ##m_step for the soft assign EM algorithm where you multiply the likelihoods
  #by the corresponding responsibilities

  if (!is.vector(param) || length(param) != group_count*2) {
    stop("param must be a vector with length equal to two times the group_count (2 corresponds to x and y coordinates per group).")
  }

  Xmat <- matrix(param, ncol = 2, byrow = TRUE)

  log_like_list_r <- list()
  for (TREE in seq_along(cond_mean_list)) {
    log_like_list_r[[TREE]] <- logL_groupmat(Xmat = Xmat,
                                             group_count = group_count,
                                             cond_mean_list[[TREE]],
                                             cond_var_list[[TREE]],
                                             normalize = FALSE)*r_list[[TREE]]
  }

  -sum(unlist(log_like_list_r))
}


#' Function to calculate sum of likelihoods when they are stored as log-likelihoods
#'
#' @param x Vector of log likelihoods
#'
#' @returns Sum of log likelihoods.
#' @export
#'
log_sum_exp <- function(x) {
  if (length(x) == 0L) {
    return(-Inf)
  }

  m <- max(x)

  if (m == Inf | m == -Inf) return(m)

  m + log(sum(exp(x - m)))
}



#' normalize log likelihoods so they sum to one
#'
#' @param loglik_mat matrix of log likelihoods
#'
#' @returns Normalized log likelihood matrix.
#' @export
#'
normalize_loglik <- function(loglik_mat) {
  # Check input
  if (!is.numeric(loglik_mat) || length(loglik_mat) == 0L) {
    stop("'loglik_mat' must be a non-empty numeric object.")
  }

  if (anyNA(loglik_mat) || any(is.nan(loglik_mat))) {
    stop("'loglik_mat' contains NA or NaN values.")
  }

  # Handle infinite values
  if (any(loglik_mat == Inf)) {
    # Assign equal probability to all positive infinities
    out <- (loglik_mat == Inf) / sum(loglik_mat == Inf)
    dim(out) <- dim(loglik_mat)
    return(out)
  }

  if (all(loglik_mat == -Inf)) {
    stop("All log-likelihoods are -Inf; normalization is undefined.")
  }



  #NORMALIZE: normalize the likelihoods so that they sum to 1. Because we are
  #initially calculating log likelihoods, I am using the log-sum-exp approach to
  #avoid underflow issues
  #https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  m <- max(loglik_mat)
  m1 <- m + log(sum(exp(loglik_mat - m)))

  return(exp(loglik_mat - m1))
}


#' Calculate log likelihoods across multiple trees
#'
#' @param Xvec matrix of coordinate to estimate the (log) likelihood of.
#' @param group_count The number of groups (K).
#' @param cond_mean_list List of conditional means.
#' @param cond_var_list List of conditional variance--covariance matrices.
#' @param normalize Boolean indicating whether they should sum to one within each tree.
#'
#' @returns A list of log likelihood matrices.
#' @export
#'
logL_groupmat_multitree <- function(Xvec,
                                    group_count,
                                    cond_mean_list,
                                    cond_var_list,
                                    normalize = FALSE) {

  if (!is.vector(Xvec) || length(Xvec) != group_count*2) {
    stop("Xvec must be a matrix with 2 columns and group_count rows.")
  }

  Xmat <- matrix(Xvec, ncol = 2, byrow = TRUE)
  output_list <- list()
  for (i in seq_along(cond_mean_list)) {

    #normalizing so that the responsibilities all sum to 1. This is simply done
    #by dividing each likelihood by the total sum of likelihoods
    output_list[[i]] <- logL_groupmat(Xmat = Xmat,
                                      group_count = group_count,
                                      cond_mean = cond_mean_list[[i]],
                                      cond_var = cond_var_list[[i]],
                                      normalize = normalize)
  }
  return(output_list)
}


#' Calculate log likelihood of coordinate values.
#'
#' @param Xmat Matrix of coordinates to calcualate the log likelihood of.
#' @param group_count The number of groups (K)
#' @param cond_mean conditional mean.
#' @param cond_var conditional variance--covariance matrix.
#' @param normalize Should the likelihood values be normalized so they sum to one?
#'
#' @returns A matrix of log likelihood values
#' @export
#'
logL_groupmat <- function(Xmat,
                          group_count,
                          cond_mean,
                          cond_var,
                          normalize = FALSE) {
  #recover()
  if (!is.matrix(Xmat) || ncol(Xmat) != 2 || nrow(Xmat) != group_count) {
    stop("Xmat must be a matrix with 2 columns and group_count rows.")
  }

  #QUESTIONS:
  # - the conditional variance/covariance is slighly non-symmetrical. Why is this?
  #   Is this just a precision problem? Currently, I'm dealing with this by
  #   forcing the matrix to be symmetric, but this might be me ignoring a bigger
  #   issue with my conditional variance calculation.


  #matrix to hold loglik values. The matrix is organized with node 1 along the
  #rows and node 2 along the columns and the locations indexed by the column and row
  #indices (e.g., location 1 for node 1 is in row 1; location 1 for node 2 is column 1)
  #For k = 2, the matrix is organized as:
  ##  ________________________    ________________________
  ##| (node1-loc1; node2-loc1)    (node1-loc1; node2-loc2) |
  ##| (node1-loc2; node2-loc1)    (node1-loc2; node2-loc2) |
  ##  ________________________    ________________________

  node_names <- unique(gsub(":[xy]", "", rownames(cond_mean))) #node1, node2
  like_mat <- matrix(data = 0, nrow = group_count, ncol = group_count)
  dimnames(like_mat) <- list(paste0(node_names[1], paste0('loc_', 1:group_count)),
                             paste0(node_names[2], paste0('loc_', 1:group_count)))

  for (NODE1 in seq_len(group_count)) {
    for (NODE2 in seq_len(group_count)) {
      #I'm sacrificing concision here for the sake of readability. We are
      #evaluating the probability at c(x_node1, y_node1, x_node2, y_node2). A one
      #line version to sub into x argument of mvtnorm::dmvnorm is: as.vector(t(X[c(NODE1, NODE2),]))
      NODE1_coords <- Xmat[NODE1,]
      NODE2_coords <- Xmat[NODE2,]
      like_mat[NODE1, NODE2] <- mvtnorm::dmvnorm(x = c(NODE1_coords, NODE2_coords),
                                                 mean = cond_mean,
                                                 sigma = as.matrix(Matrix::forceSymmetric(cond_var)),
                                                 log=TRUE)
    }
  }

  #NORMALIZE: normalize the likelihoods so that they sum to 1. Because we are
  #initially calculating log likelihoods, I am using the log-sum-exp approach to
  #avoid underflow issues
  #https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  if (isTRUE(normalize)) return(normalize_loglik(loglik_mat = like_mat))

  return(like_mat)
}
