#' Create initial origin coordinates
#'
#' @param cond_mean_list List of conditional means. This is only required for method = "kpp" or "greedy_kpp"
#' @param k The number of origins.
#' @param method The method to use for coordinate initialization. The three options are "random" (random coordinates selected within the sp_bounds), or "kpp" (k-means++ on the conditional means), or "greedy_kpp" (greedy k-means++ on the conditional means).
#' @param sp_bounds The boundaries to choose the initial coordinates. This is only required for method = "random"
#' @param n_cand The number of candidate origins to evaluate during each round of choosing origins. This is only used when method = "greedy_kpp".
#'
#' @returns A matrix (col = 2, row = k) of initial coordinates with x values in the first column and y values in the second.
#' @export
#'
sample_locs <- function(cond_mean_list = NULL, k, method = 'random', sp_bounds = NULL, n_cand = 3) {

  method <- match.arg(arg = method, choices = c('random', 'kpp', 'greedy_kpp'), several.ok = FALSE)

  #check K
  if (length(k) != 1L ||
      !is.numeric(k) ||
      is.na(k) ||
      !is.finite(k) ||
      k <= 0 ||
      k != floor(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  #check random inputs
  if (identical(method, "random")) {
    if (!is.list(sp_bounds) ||
        !all(c("x", "y") %in% names(sp_bounds))) {
      stop("For method = 'random', `sp_bounds` must be a list with x and y bounds.",
           call. = FALSE)
    }

    x_bounds <- sp_bounds[["x"]]
    y_bounds <- sp_bounds[["y"]]

    x_check <- is.numeric(x_bounds) &&
      length(x_bounds) == 2L &&
      all(is.finite(x_bounds)) &&
      x_bounds[1] <= x_bounds[2]

    y_check <- is.numeric(y_bounds) &&
      length(y_bounds) == 2L &&
      all(is.finite(y_bounds)) &&
      y_bounds[1] <= y_bounds[2]

    if (!x_check || !y_check) {
      stop("The x and y bounds must each contain two finite, ordered numbers.",
           call. = FALSE)
    }
  }

  initial_locs <- switch(method,
                         'random' = {
                           cbind(stats::runif(k,min=sp_bounds$x[1],max=sp_bounds$x[2]),
                                 stats::runif(k,min=sp_bounds$y[1],max=sp_bounds$y[2]))
                         },
                         'kpp' = {
                           kpp(cond_mean_list = cond_mean_list, k = k, type = 'standard')
                         },
                         'greedy_kpp' = {
                           kpp(cond_mean_list = cond_mean_list, k = k, type = 'greedy', n_cand = n_cand)
                         })

  return(initial_locs)
}



#' Create initial origin coordinates using k-means++ or greedy k-means++
#'
#' @param cond_mean_list cond_mean_list List of conditional means.
#' @param k The number of origins.
#' @param type The approach to use for coordinate generation. Either the standard k-means++ ('standard'; the default) or the greedy k-means++.
#' @param n_cand The number of candidate origins to evaluate during each round of choosing origins. This is only relevant when type = 'greedy'.
#'
#' @returns A matrix (col = 2, row = k) of initial coordinates with x values in the first column and y values in the second.
#' @export
#'
kpp <- function(cond_mean_list, k, type = 'standard', n_cand = 3) {

  #==========================================
  #== INPUT CHECKS AND INITIAL PROCESSING ===
  #==========================================
  type <- match.arg(type, choices = c('standard', 'greedy'), several.ok = FALSE)

  #checks if type = 'greedy'
  if (type == 'greedy') {
    if (length(n_cand) != 1L ||
        !is.numeric(n_cand) ||
        is.na(n_cand) ||
        !is.finite(n_cand) ||
        n_cand <= 0 ||
        n_cand != floor(n_cand)) {
      stop("`n_cand` must be a single positive integer.", call. = FALSE)
    }

    if (n_cand > 100) {
      stop("`n_cand` maximum is arbitrarily set to 100. Please choose a smaller value.", call. = FALSE)
    }
  }

  if (length(k) != 1L ||
      !is.numeric(k) ||
      is.na(k) ||
      !is.finite(k) ||
      k <= 0 ||
      k != floor(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  cond_mean_df <- do.call(rbind, lapply(cond_mean_list, function(x) matrix(as.vector(x), byrow = TRUE, ncol = 2)))
  n_nodes <- nrow(cond_mean_df)

  if (k > n_nodes) {
    stop("`k` cannot be more than the number of nodes.", call. = FALSE)
  }

  #==========================
  #== LOCATION GENERATION ===
  #==========================
  loc_mat <- matrix(data = NA, nrow = k, ncol = 2)

  selected_indices <- integer(k)

  selected_indices[1] <- sample.int(n_nodes, size = 1)
  loc_mat[1, ] <- cond_mean_df[selected_indices[1], ]

  if (k < 2) return(loc_mat)

  for (loc_ind in 2:k) {

    #2. Calculate the distance of each conditional mean coordinate to each origin coordinate (and take the minimum)
    loc_mat_subset <- loc_mat[!is.na(loc_mat[,1]),, drop = FALSE]
    min_dist_vec <- apply(cond_mean_df, 1, function(x, origin_mat) {min(sqrt(rowSums((sweep(origin_mat, 2, x, "-"))^2)))}, origin_mat = loc_mat_subset)

    #3. Select next coordinate probabilistically with coordinates weighted by squared distance
    weights <- min_dist_vec^2 #weights are the squared distances (don't need to normalize because sample.int does that automatically)
    weights[selected_indices[1:(loc_ind - 1)]] <- 0 #prevent previous conditional means from being chosen (by setting weight to 0)

    #deal with edge case of all weights equaling 0
    if (!any(weights > 0)) {
      weights[] <- 1
      weights[selected_indices[1:(loc_ind - 1)]] <- 0
    }

    if (type == 'greedy') {
      #3b. Extra step when type == 'greedy': select multiple candidates, add them temporarily to the origin matrix, calculate the
      #    potential, and then choose the one with the smallest poential
      candidate_inds <- sample.int(n_nodes, size = n_cand, prob = weights, replace = TRUE) #sample the index to use for the new origin coordinates

      potential_vec <- vector(mode = 'numeric', n_cand)
      for (i in 1:length(candidate_inds)) {
        test_mat <- rbind(loc_mat_subset, cond_mean_df[candidate_inds[i],])
        potential_vec[i] <- sum(apply(cond_mean_df, 1, function(x, origin_mat) {min(sqrt(rowSums((sweep(origin_mat, 2, x, "-"))^2)))},
                                      origin_mat = test_mat)^2)
      }

      #find the candidate the smallest potential
      sampled_index <- candidate_inds[which.min(potential_vec)]
    } else {
      #original
      sampled_index <- sample.int(n_nodes, size = 1, prob = weights) #sample the index to use for the new origin coordinates

    }

    selected_indices[loc_ind] <- sampled_index #add the sampled index

    loc_mat[loc_ind,] <- cond_mean_df[sampled_index,] #add the new origin coordinate to loc_mat

  }

  return(loc_mat)
}



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
