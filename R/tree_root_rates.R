###################################################
### ROOT AND BROWNIAN RATE ESTIMATION FUNCTIONS ###
###################################################

#' Calculate estimate of root coordinates
#'
#' @param vcv Variance--covariance matrix of tree
#' @param trait_mat matrix containing tip coordinates.
#'
#' @returns Expected root coordinates.
#' @export
#'
expected_root_state <- function(vcv, trait_mat) {
  #recover()

  one_col_vec <- matrix(rep(1, nrow(trait_mat)), ncol = 1)

  vcv_inv <- solve(vcv)

  return(
    t(solve(t(one_col_vec)%*%vcv_inv%*%one_col_vec)%*%t(one_col_vec)%*%vcv_inv%*%trait_mat)
  )

}

#' Wrapper to calculate expected root coordinates only using the tips that have coordinate information.
#'
#' @param missing_inds Tips that don't have coordinate information.
#' @param tree_covar Variance--covariance matrix of tree tips.
#' @param trait_mat Coordinate information for the tips.
#'
#' @returns Expected root coordinates.
#' @export
#'
calc_root_coords <- function(missing_inds, tree_covar, trait_mat) {
  vcv_subset <- tree_covar[!(rownames(tree_covar) %in% missing_inds),
                           !(colnames(tree_covar) %in% missing_inds)]

  root_state <- expected_root_state(vcv_subset,
                                    trait_mat)

  return(root_state)
}



#' Estimate Brownian rate matrix with (or without) missing individuals
#'
#' @param missing_inds Vector of missing individuals in the tree variance--covariance matrix to omit during rate matrix estimation.
#' @param tree_covar Tree variance--covariance matrix.
#' @param root_vals Coordinates of the tree root (generally estimated)
#' @param trait_mat Tip values
#'
#' @returns Estimated Brownian rate matrix
#' @export
#'
est_brownian_rate_mat <- function(missing_inds = NULL, tree_covar, root_vals, trait_mat) {
  #recover()
  vcv_subset <- tree_covar[!(rownames(tree_covar) %in% missing_inds),
                           !(colnames(tree_covar) %in% missing_inds)]
  one_col_vec <- matrix(rep(1, nrow(trait_mat)),
                        ncol = 1)

  term1 <- trait_mat - one_col_vec%*%t(root_vals)
  return(t(term1)%*%solve(vcv_subset)%*%term1/nrow(trait_mat))
}


#' Wrapper that calculates point estimates of root coordinates and Brownian rate matrix given a tree and tip coordinate and (optionally) individuals to exclude from these estimate.
#'
#' @param tree A tree (stored as a phylo object)
#' @param missing_inds The individuals to exclude from the root and rate matrix calculations
#' @param trait_info Matrix of coordinates for the tips.
#'
#' @returns A list containing estimates of the root coordinates and the Brownian rate matrix.
#' @export
#'
estimate_root_and_rate <- function(tree, missing_inds, trait_info) {
  #recover()
  tree_covar <- ape::vcv.phylo(tree)
  vcv_subset <- tree_covar[!(rownames(tree_covar) %in% missing_inds),
                           !(colnames(tree_covar) %in% missing_inds)]
  trait_reorder <- trait_info[rank(match(rownames(tree_covar), rownames(trait_info)), na.last = NA),]


  root_coords <- calc_root_coords(missing_inds = NULL,
                                  tree_covar = vcv_subset,
                                  trait_mat = trait_reorder)

  brown_rate_mat <- est_brownian_rate_mat(missing_inds = NULL,
                                          tree_covar = vcv_subset,
                                          root_vals = root_coords,
                                          trait_mat = trait_reorder)

  return(
    list(
      root_coords = root_coords,
      brown_rate_mat = brown_rate_mat
    )
  )
}


#' Wrapper that calculates point estimates of root coordinates and Brownian rate matrix across a list of trees.
#'
#' @param tree_list A list of trees (each stored as a phylo object).
#' @param trait A matrix of coordinates for each tip.
#' @param missing_inds Individuals (i.e., tips) to exclude from estimating root coordinates and Brownian rate matrices.
#'
#' @returns A list containing a list of root coordinates, a list of Brownian rate matrices, and the mean Brownian matrix.
#' @export
#'
root_rate_multitree <- function(tree_list, trait, missing_inds) {
  #recover()
  info_list <- lapply(tree_list, function(x, missing, trait) {
    estimate_root_and_rate(tree = x,
                           missing_inds = missing,
                           trait_info = trait)

  }, missing = missing_inds, trait = trait)


  return(
    list(
      root_list = lapply(info_list, function(x) x$root_coords),
      brown_rate_list = lapply(info_list, function(x) x$brown_rate_mat),
      mean_brown_rate = Reduce("+", lapply(info_list, function(x) x$brown_rate_mat))/length(info_list)
    )
  )
}

