#########################################################
### PROCESSING AND DISTRIBUTION CALCULATION FUNCTIONS ###
#########################################################

#' Function to process a tree.
#'
#' @param tree The tree (stored as a phylo object).
#' @param trait Matrix of tip coordinates.
#' @param brown_rate_mat The Brownian rate matrix for lineage movement.
#' @param root_coords The tree's root coordinates.
#' @param missing_inds The tips with missing coordinate information.
#' @param retain The elements to retain in the output.
#'
#' @returns A list of relevant information for the tree.
#' @export
#'
process_tree <- function(tree, trait, brown_rate_mat, root_coords, missing_inds,
                         retain = c('tree_covar', 'trait_reorder', 'root_coords', 'brown_rate_mat', 'cond_distr')) {

  element_list <- list()

  element_list[['tree_covar']] <- ape::vcv.phylo(tree)

  element_list[['trait_reorder']] <- trait[rank(match(rownames(element_list[['tree_covar']]), rownames(trait)), na.last = NA),]

  element_list[['root_coords']] <- root_coords
  element_list[['brown_rate_mat']] <- brown_rate_mat

  element_list[['cond_distr']] <- calc_cond_distr(element_list[['tree_covar']],
                                                  brown_rate_mat = element_list[['brown_rate_mat']],
                                                  root_coords = element_list[['root_coords']],
                                                  trait_mat = element_list[['trait_reorder']],
                                                  missing_inds = missing_inds)

  return(element_list[names(element_list) %in% retain])
}


#' Function to process a list of trees.
#'
#' @param tree_list List of trees (stored as phylo objects).
#' @param trait Matrix of tip coordinates.
#' @param brown_rate_mat The Brownian rate matrix for lineage movement.
#' @param root_coords_list List of the trees root coordinates.
#' @param missing_inds The tips with missing coordinate information.
#' @param progress Boolean indicating whether tree processing progress should be reported.
#'
#' @returns A list processed info for each tree.
#' @export
#'
process_trees <- function(tree_list, trait, brown_rate_mat, root_coords_list, missing_inds, progress = TRUE) {

  if (isTRUE(progress)) prog_bar <- utils::txtProgressBar(min = 0, max = length(tree_list), style = 3, char = "+")

  cond_info_list <- list()
  for (i in 1:length(tree_list)) {
    cond_info_list[[i]] <- process_tree(tree = tree_list[[i]],
                                        trait = trait,
                                        brown_rate_mat = brown_rate_mat,
                                        root_coords = root_coords_list[[i]],
                                        missing_inds = missing_inds,
                                        retain = 'cond_distr')

    if (isTRUE(progress)) utils::setTxtProgressBar(prog_bar, i)
  }

  return(cond_info_list)

}


#' Calculate Kronecker product of tree variance--covariance matrix and Brownian rate matrix
#'
#' @param tree_covar Tree variance--covariance matrix
#' @param rate_mat The Brownian rate matrix for lineage movement.
#' @param missing_inds The tips that are missing (i.e., the non-georeferenced individuals).
#' @param include_inds Boolean indicating whether to include info. on tips that are missing and not missing in the output.
#'
#' @returns A list containing the Kronecker product and (optionally) the tips that are missing and not missing.
#' @export
#'
calc_kronmat <- function(tree_covar, rate_mat, missing_inds = NULL, include_inds = TRUE) {

  output_list <- list()

  output_list[['kron_mat']] <- kronecker(tree_covar, rate_mat, make.dimnames = TRUE)

  if (isTRUE(include_inds)) {
    output_list[['yes_ind']] <- which(!gsub(':[xy]', '', rownames(output_list[['kron_mat']])) %in% missing_inds)
    output_list[['no_ind']] <- which(gsub(':[xy]', '', rownames(output_list[['kron_mat']])) %in% missing_inds)
  }

  return(output_list)
}


#' Calculate the elements necessary to calculate the conditional distribution.
#'
#' @param kron_varcovar_ratemat Kronecker product of tree variance--covariance matrix and Brownian rate matrix.
#' @param no_ind Non-georeferenced individuals.
#' @param yes_ind Georeferenced individuals (the individuals to condition the distribution on).
#'
#' @returns List of elements and info. necessary to calculate the conditional distribution.
#' @export
#'
calc_cond_var_elements <- function(kron_varcovar_ratemat,
                                   no_ind,
                                   yes_ind) {

  covar_11 <- kron_varcovar_ratemat[no_ind,no_ind]
  covar_12 <- kron_varcovar_ratemat[no_ind,yes_ind]
  inv_covar_22 <- solve(kron_varcovar_ratemat[yes_ind,yes_ind])
  covar_21 <- kron_varcovar_ratemat[yes_ind,no_ind]

  #cond_mean <- mean1 + (covar_12 %*% solve(covar_22) %*% (trait_vec -  mean2))
  #cond_var <- covar_11 - covar_12%*%inv_covar_22%*%covar_21

  return(
    list(
      covar_11 = covar_11,
      covar_12 = covar_12,
      inv_covar_22 = inv_covar_22,
      covar_21 = covar_21,
      cond_var = covar_11 - covar_12%*%inv_covar_22%*%covar_21
    )
  )

}

#' Calculate mean of conditional distribution of the non-georeferenced tip locations.
#'
#' @param covar_12 Variance--covariance matrix of non-georeferenced individuals and georeferenced individuals.
#' @param inv_covar22 Inverse of variance--covariance matrix of georeferenced individuals.
#' @param trait_vec Vector of coordinates.
#' @param means Vector of means for the unconditioned distribution of tip locations. This is the root coordinates of the tree.
#' @param no_ind Non-georeferenced individuals.
#' @param yes_ind Georeferenced individuals (the individuals to condition the distribution on).
#'
#' @returns Mean vector of conditional distribution.
#' @export
#'
calc_cond_mean <- function(covar_12, inv_covar22,
                           trait_vec, means,
                           no_ind, yes_ind) {
  mean1 <- rep(means, times = length(no_ind)/2)
  mean2 <- rep(means, times = length(yes_ind)/2)
  return(
    mean1 + (covar_12 %*% inv_covar22 %*% (trait_vec -  mean2))
  )
}





#' Calculate distribution of the non-georeferenced tips conditional on the geo-referenced tips.
#'
#' @param tree_covar Tree variance--covariance matrix.
#' @param brown_rate_mat The Brownian rate matrix for lineage movement.
#' @param root_coords Coordinates of tree's root.
#' @param trait_mat Matrix of tip coordinates.
#' @param missing_inds The tips with missing coordinate info. (i.e., the non-georeferenced individual).
#'
#' @returns list of mean and variance--covariance matrix of the conditional distribution.
#' @export
#'
calc_cond_distr <- function(tree_covar, brown_rate_mat, root_coords, trait_mat, missing_inds) {

  trait_reorder <- trait_mat[rank(match(rownames(tree_covar), rownames(trait_mat)), na.last = NA),]

  element_list <- list()

  kron <- calc_kronmat(tree_covar = tree_covar,
                       rate_mat = brown_rate_mat,
                       missing_inds = missing_inds,
                       include_inds = TRUE)

  element_list[['cond_var_elements']] <- calc_cond_var_elements(
    kron_varcovar_ratemat = kron$kron_mat,
    no_ind = kron$no_ind,
    yes_ind = kron$yes_ind
  )


  #calculate conditional mean
  element_list[['cond_mean']] <- calc_cond_mean(covar_12 = element_list[['cond_var_elements']]$covar_12,
                                                inv_covar22 = element_list[['cond_var_elements']]$inv_covar_22,
                                                trait_vec = c(t(trait_reorder)), #element_list[['data_prep_output']]$vector_trait,
                                                means = root_coords, #element_list[['ls_root']], #currently set at least squares estimate of root coordinates
                                                no_ind = kron$no_ind,
                                                yes_ind = kron$yes_ind)

  #return(element_list[names(element_list) %in% retain])
  return(
    list(
      condmean = element_list[['cond_mean']],
      condvar = element_list[['cond_var_elements']]$cond_var
    )
  )
}
