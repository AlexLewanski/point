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
