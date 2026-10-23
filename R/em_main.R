#' Multi-run wrapper for running expectation-maximization using hard node assignments
#'
#' @param times The number of times to run the EM algorithm.
#' @param cond_mean List of conditional means for each tree.
#' @param cond_var List of conditional variance--covariance matrices for each tree.
#' @param k The number of ancestry groups to estimate.
#' @param max_stp The maximum number of steps allowed in the EM algorithm.
#' @param conv_thresh The difference in observed likelihoods between runs to use as the convergence threshold.
#' @param lik_equality_tol Value under which likelihoods are considered equivalent.
#' @param init_method The method to use for initialization of the origin coordinates. The options are "random", "kpp", or "greedy_kpp". The latter two implement k-means++ or greedy k-means++ on the conditional means (default: "random").
#' @param sp_bounds List of two vectors indicating minimum and maximum x and y values.
#' @param retain_all_outputs Boolean indicating whether or not to retain all outputs
#' @param starting_location_seed Optional seed to set when generating initial coordinates.
#' @param n_cand The number of candidate origins to evaluate during each round of choosing origins. This is only used when method = "greedy_kpp".
#'
#' @returns A list of results from the EM runs.
#' @export
#'
em_hardassign_wrapper <- function(times = 1,
                                  cond_mean,
                                  cond_var,
                                  k,
                                  max_stp = 20,
                                  conv_thresh = 1e-5,
                                  lik_equality_tol = 1e-6,
                                  init_method = "random",
                                  sp_bounds = NULL,
                                  retain_all_outputs = TRUE,
                                  starting_location_seed = NULL,
                                  n_cand = 3) {
  #recover()

  if (!is.null(starting_location_seed)) {
    set.seed(starting_location_seed)
  }

  est_list <- list()

  for (est in seq_len(times)) {
    est_list[[est]] <- em_hardassign(cond_mean = cond_mean,
                                     cond_var = cond_var,
                                     k = k,
                                     max_stp = max_stp,
                                     conv_thresh = conv_thresh,
                                     init_method = init_method,
                                     sp_bounds = sp_bounds,
                                     starting_location_seed = NULL,
                                     n_cand = n_cand
    )

    message('finished ', est)

  }


  process_output_list <- process_em_output(output_list = est_list,
                                           lik_tol = lik_equality_tol,
                                           assignment_approach = 'hard')

  if (isTRUE(retain_all_outputs)) process_output_list[['extra_info']][['output_list']] <- est_list
  #if (isTRUE(retain_all_outputs)) process_output_list[['output_list']] <- est_list

  return(process_output_list)
}


#' Multi-run wrapper for running expectation-maximization using soft node assignments
#'
#' @param times The number of times to run the EM algorithm.
#' @param cond_mean List of conditional means for each tree.
#' @param cond_var List of conditional variance--covariance matrices for each tree.
#' @param k The number of ancestry groups to estimate.
#' @param max_stp The maximum number of steps allowed in the EM algorithm.
#' @param conv_thresh The difference in observed likelihoods between runs to use as the convergence threshold.
#' @param lik_equality_tol Value under which likelihoods are considered equivalent.
#' @param init_method The method to use for initialization of the origin coordinates. The options are "random", "kpp", or "greedy_kpp". The latter two implement k-means++ or greedy k-means++ on the conditional means (default: "random").
#' @param sp_bounds List of two vectors indicating minimum and maximum x and y values.
#' @param retain_all_outputs Boolean indicating whether or not to retain all outputs
#' @param starting_location_seed Optional seed to set when generating initial coordinates.
#' @param n_cand The number of candidate origins to evaluate during each round of choosing origins. This is only used when method = "greedy_kpp".
#'
#' @returns A list of results from the EM runs.
#' @export
#'
em_softassign_wrapper <- function(times = 1,
                                  cond_mean,
                                  cond_var,
                                  k,
                                  max_stp = 20,
                                  conv_thresh = 1e-5,
                                  lik_equality_tol = 1e-6,
                                  init_method = "random",
                                  sp_bounds = NULL,
                                  retain_all_outputs = TRUE,
                                  starting_location_seed = NULL,
                                  n_cand = 3) {
  #recover()

  if (!is.null(starting_location_seed)) {
    set.seed(starting_location_seed)
  }


  est_list <- list()

  for (est in seq_len(times)) {
    est_list[[est]] <- em_softassign(cond_mean = cond_mean,
                                     cond_var = cond_var,
                                     k = k,
                                     max_stp = max_stp,
                                     conv_thresh = conv_thresh,
                                     init_method = init_method,
                                     sp_bounds = sp_bounds,
                                     starting_location_seed = NULL,
                                     n_cand = n_cand
    )

    message('finished ', est)

  }


  process_output_list <- process_em_output(output_list = est_list,
                                           lik_tol = lik_equality_tol,
                                           assignment_approach = 'soft')

  if (isTRUE(retain_all_outputs)) process_output_list[['extra_info']][['output_list']] <- est_list
  #if (isTRUE(retain_all_outputs)) process_output_list[['output_list']] <- est_list

  return(process_output_list)
}




#' Expectation-maximization using hard node assignments
#'
#' @param cond_mean List of conditional means for each tree.
#' @param cond_var List of conditional variance--covariance matrices for each tree.
#' @param k The number of ancestry groups to estimate.
#' @param max_stp The maximum number of steps allowed in the EM algorithm.
#' @param conv_thresh The difference in observed likelihoods between runs to use as the convergence threshold.
#' @param init_method The method to use for initialization of the origin coordinates. The options are "random", "kpp", or "greedy_kpp". The latter two implement k-means++ or greedy k-means++ on the conditional means (default: "random").
#' @param sp_bounds List of two vectors indicating minimum and maximum x and y values.
#' @param starting_location_seed Optional seed to set when generating initial coordinates.
#' @param n_cand The number of candidate origins to evaluate during each round of choosing origins. This is only used when method = "greedy_kpp".
#'
#' @returns A list of outputs from the EM run.
#' @export
#'
em_hardassign <- function(cond_mean,
                          cond_var,
                          k,
                          max_stp,
                          conv_thresh,
                          init_method = 'random',
                          sp_bounds = NULL,
                          starting_location_seed = NULL,
                          n_cand = 3) {

  #check to make sure max_stp is at least 3
  stopifnot(max_stp >= 3)

  #recover()

  #=======================
  #=== FUNCTION SET-UP ===
  #=======================

  tree_count <- length(cond_mean)

  #matrix of all pairwise location combos
  #group_pairs <- expand.grid(x1 = 1:K, x2 = 1:K)

  #initialize node membership vec
  group_vec <- rep(1, times = tree_count*2)

  #counter for location initialization
  init_counter <- 1

  #lists to store things
  zvec_list <- list() #node membership vectors
  xvec_list <- list() #location coordinate matrices

  #all_lik_list <- list()

  loglik_list <- list()
  loglik_list[[1]] <- NA


  #===============================
  #=== LOCATION INITIALIZATION ===
  #===============================

  #setting seed before the repeat loop so that a new set of coordinates are generated
  #instead of just fixed coordinates if the set.seed is placed within the loop

  # Set the seed only when one is supplied. Wrappers seed once for the whole
  # batch, so leaving this NULL lets each run draw from the continuing RNG stream.
  if (!is.null(starting_location_seed)) {
    set.seed(starting_location_seed)
  }

  repeat {
    Xold <- as.vector(t(sample_locs(cond_mean_list = cond_mean, k = k, method = init_method, sp_bounds = sp_bounds, n_cand = n_cand)))

    # Xold <- as.vector(t(cbind(stats::runif(k,min=sp_bounds$x[1],max=sp_bounds$x[2]),
    #                           stats::runif(k,min=sp_bounds$y[1],max=sp_bounds$y[2]))))
    if (k == 1) break

    group_info_list <- list()
    for (i in seq_len(tree_count)) {
      lik_info <- logL_groupmat(Xmat = matrix(Xold, ncol = 2, byrow = TRUE),
                                group_count = k,
                                cond_mean[[i]],
                                cond_var[[i]],
                                normalize = FALSE)
      #find the node location pairing that maximizes the likelihood
      group_info_list[[i]] <- which(lik_info == max(lik_info), arr.ind = TRUE)[1,]
    }

    group_vec <- unlist(group_info_list) #vector version of node memberships

    init_counter <- init_counter + 1

    if (length(unique(group_vec)) == k) {
      break # Exit the loop if the condition is no longer met
    }

    if (init_counter >= 500) stop('Coordination initialization failed.')
  }

  message("Initialization attempts: ", init_counter)


  #initial param vals
  #zvec_list[[1]] <- group_vec
  zvec_list[[1]] <- NA
  xvec_list[[1]] <- Xold



  #================
  #=== EM ITERS ===
  #================
  stp <- 2

  while (stp < max_stp) {

    ### EXPECTATION: given location coordinates, what are the node location memberships? ###
    group_info_list <- list()
    # group_like <- list()
    for (i in seq_len(tree_count)) {
      lik_info <- logL_groupmat(Xmat = matrix(Xold, ncol = 2, byrow = TRUE),
                                group_count = k,
                                cond_mean[[i]],
                                cond_var[[i]],
                                normalize = FALSE)
      #find the node location pairing that maximizes the likelihood
      group_info_list[[i]] <- which(lik_info == max(lik_info), arr.ind = TRUE, useNames = TRUE)[1,]

      #overly complicated way of assigning node names to the population labels
      names(group_info_list[[i]]) <- paste0(i, "_", c(gsub("loc.*", "", rownames(lik_info)[group_info_list[[i]][1]]), gsub("loc.*", "", colnames(lik_info)[group_info_list[[i]][2]])))
      #names(group_info_list[[i]]) <- paste0('tree', i, c('_node1', '_node2'))
      # group_like[[i]] <- logL_groupmat(X = matrix(Xold, ncol = 2, byrow = TRUE),
      #                                                        group_count = k,
      #                                                        cond_mean[[i]],
      #                                                        cond_var[[i]],
      #                                                      normalize = FALSE)
      # group_info_list[[i]] <- which(group_like[[i]] == max(group_like[[i]]), arr.ind = TRUE)[1,]

    }
    #all_lik_list[[stp]] <- group_like

    group_vec <- unlist(group_info_list) #vector version of node memberships


    ### MAXIMIZATION: DERIVE NEW GROUP COORDINATES USING THE WEIGHTS ###
    #I'm currently doing this by finding the generalized least squares solution
    #based on the node memberships found during the expectation step
    Xnew <- as.vector(weighted_least_squares(multigroup_design(max_group = k,  group_vec),
                                             Sigma_list = cond_var,
                                             mu = unlist(cond_mean)))


    #list of all log likelihoods using the new parameter values
    like_list <- logL_groupmat_multitree(Xvec = Xnew,
                                         group_count = k,
                                         cond_mean_list = cond_mean,
                                         cond_var_list = cond_var,
                                         normalize = FALSE)

    #add new parameter vals to x_vec and assign Xnew to Xold
    xvec_list[[stp]] <- Xnew
    Xold <- Xnew

    zvec_list[[stp]] <- group_vec

    #the observed data log-likelihood is calculated by
    #(1) summing the likelihoods (not log likelihoods) within trees
    #(2) taking the log of the summed likelihood for each tree
    #(3) summing the log summed likelihood across trees
    #loglik_list[[stp]] <- sum(unlist(lapply(like_list, function(x) {log(sum(exp(x)))}))) #sum(unlist(like_list)) #log likelihood
    loglik_list[[stp]] <- sum(unlist(lapply(like_list, log_sum_exp))) #sum(unlist(like_list)) #log likelihood

    message(paste0('iter ', stp, '; loglik = ', loglik_list[[stp]]))


    ### CHECKS FOR THE NEXT ITER ###
    # proceed to next iter if this was the first iter
    if (stp == 2) {
      stp <- stp + 1
      next
    }

    #if the current lik is less than the previous lik, write a warning
    if (loglik_list[[stp]] < loglik_list[[stp - 1]] ) {
      warning('The likelihood decreased!')

    }

    #if the difference in log liks is below the stopping threshold, break out of loop
    loglik_dif <- loglik_list[[stp]] - loglik_list[[stp - 1]]
    if (abs(loglik_dif) <= conv_thresh) {
      convergence <- TRUE
      break
    }
    #################################

    stp <- stp + 1
  }

  return(
    list(final_x = Xnew,
         final_loglik = loglik_list[[length(loglik_list)]],
         xvec_list = xvec_list,
         iter_count = stp,
         group_membership_list = zvec_list,
         like_vec = unlist(loglik_list))
  )
}


#' Expectation-maximization using soft node assignments
#'
#' @param cond_mean List of conditional means for each tree.
#' @param cond_var List of conditional variance--covariance matrices for each tree.
#' @param k The number of ancestry groups to estimate.
#' @param max_stp The maximum number of steps allowed in the EM algorithm.
#' @param conv_thresh The difference in observed likelihoods between runs to use as the convergence threshold.
#' @param init_method The method to use for initialization of the origin coordinates. The options are "random", "kpp", or "greedy_kpp". The latter two implement k-means++ or greedy k-means++ on the conditional means (default: "random").
#' @param sp_bounds List of two vectors indicating minimum and maximum x and y values.
#' @param starting_location_seed Optional seed to set when generating initial coordinates.
#' @param n_cand The number of candidate origins to evaluate during each round of choosing origins. This is only used when method = "greedy_kpp".
#'
#' @returns A list of outputs from the EM run.
#' @export
#'
em_softassign <- function(cond_mean,
                          cond_var,
                          k,
                          max_stp,
                          conv_thresh,
                          init_method = 'random',
                          sp_bounds = NULL,
                          starting_location_seed = NULL,
                          n_cand = 3) {

  #recover()

  #=======================
  #=== FUNCTION SET-UP ===
  #=======================

  #check to make sure max_stp is at least 3
  stopifnot(max_stp >= 3)

  tree_count <- length(cond_mean) #number of trees

  #lists to store things
  weight_list <- list() #node membership vectors
  weight_list[[1]] <- NA

  xvec_list <- list() #location coordinate matrices
  loglik_list <- list()
  loglik_list[[1]] <- NA

  q_list <- list()

  convergence <- FALSE



  #===============================
  #=== LOCATION INITIALIZATION ===
  #===============================

  #set.seed(3543)
  #currently doing this in the simplest way possible: uniformly sampling across
  #the area defined by sp_bounds

  # Set the seed only when one is supplied. Wrappers seed once for the whole
  # batch, so leaving this NULL lets each run draw from the continuing RNG stream.
  if (!is.null(starting_location_seed)) {
    set.seed(starting_location_seed)
  }

  Xold <- as.vector(t(sample_locs(cond_mean_list = cond_mean, k = k, method = init_method, sp_bounds = sp_bounds, n_cand = n_cand)))

  #Xold <- as.vector(t(cbind(stats::runif(k, min = sp_bounds$x[1], max = sp_bounds$x[2]),
  #                          stats::runif(k, min = sp_bounds$y[1], max = sp_bounds$y[2]))))

  #initial param vals
  xvec_list[[1]] <- Xold



  #================
  #=== EM ITERS ===
  #================
  stp <- 2

  while(stp < max_stp){

    ### EXPECTATION: given location coordinates, what are the node location memberships? ###
    pair_w_list <- logL_groupmat_multitree(Xvec = Xold, #matrix(Xold, ncol = 2, byrow = TRUE),
                                           group_count = k,
                                           cond_mean_list = cond_mean,
                                           cond_var_list = cond_var,
                                           normalize = TRUE)


    ### MAXIMIZATION: DERIVE NEW GROUP COORDINATES USING THE WEIGHTS ###
    #I'm currently doing this by maximizing the product of likelihood calculated
    #across all possible latent variable (location membership) states and the
    #corresponding responsibilities
    #maximizing the expected complete-data log-likelihood
    m_val <- stats::optim(
      par = Xold, #starting parameter values at Xold
      fn = m_step,
      r_list= pair_w_list,
      cond_mean_list = cond_mean,
      cond_var_list = cond_var,
      group_count = k,
      method = "BFGS",
      control = c(maxit = 1000)
    )

    message('convergence: ', m_val$convergence)
    message('message: ', m_val$message)
    Xnew <- m_val$par

    qnew <- m_step(param = Xnew,
                   r_list = pair_w_list,
                   cond_mean_list = cond_mean,
                   cond_var_list = cond_var,
                   group_count = k)

    #The observed-data likelihood does the opposite
    #list of all log likelihoods using the new parameter values
    like_list <- logL_groupmat_multitree(Xvec = Xnew, #matrix(Xnew, ncol = 2, byrow = TRUE),
                                         group_count = k,
                                         cond_mean_list = cond_mean,
                                         cond_var_list = cond_var,
                                         normalize = FALSE)

    #add new parameter vals to x_vec and assign Xnew to Xold
    xvec_list[[stp]] <- Xnew
    q_list[[stp]] <- qnew
    Xold <- Xnew

    weight_list[[stp]] <- pair_w_list #list of pairwise responsibilities

    #the observed data log-likelihood is calculated by
    #(1) summing the likelihoods (not log likelihoods) within trees
    #(2) taking the log of the summed likelihood for each tree
    #(3) summing the log summed likelihood across trees
    #loglik_list[[stp]] <- sum(unlist(lapply(like_list, function(x) {log(sum(exp(x)))}))) #sum(unlist(like_list)) #log likelihood
    loglik_list[[stp]] <- sum(unlist(lapply(like_list, log_sum_exp)))

    message(paste0('iter ', stp, '; loglik = ', loglik_list[[stp]]))


    ### CHECKS FOR THE NEXT ITER ###
    # proceed to next iter if this was the first iter
    if (stp == 2) {
      stp <- stp + 1
      next
    }

    #check to see if the likelihood decreases
    previous_loglik <- loglik_list[[stp - 1]]
    current_loglik <- loglik_list[[stp]]
    loglik_diff <- current_loglik - previous_loglik #this shouldn't be negative

    #currently using a "sliding scale" tolerance that reflects the magnitude of
    #the logliks. The minimum threshold is 1e-10. If a loglik is bigger than
    #1, the lik_tol increases slightly to reflect this magnitude
    lik_tol <- 1e-10 * max(1, abs(previous_loglik), abs(current_loglik))

    #if the current loglik is less than the (negative) tolerance, throw an error
    if (loglik_diff < -lik_tol) {
      stop(sprintf("Log-likelihood decreased at iteration %d: previous = %.12g, current = %.12g, difference = %.12g (tol = %.5g).",
                   stp, previous_loglik, current_loglik, loglik_diff, lik_tol), call. = FALSE)
    }

    #if the difference in log liks is below the stopping threshold, break out of loop
    if (loglik_diff <= conv_thresh) {
      convergence <- TRUE
      break
    }
    #################################

    stp <- stp + 1
  }

  return(
    list(final_x = Xnew, #final location estimate
         final_loglik = loglik_list[[length(loglik_list)]], #final loglik
         convergence = convergence,
         xvec_list = xvec_list, #location estimates from each iter
         group_membership_list = weight_list, #list of weights from each step
         like_vec = unlist(loglik_list), #log liks from each step of EM
         iter_count = stp, #number of iters,
         q_list = q_list
    )
  )

}
