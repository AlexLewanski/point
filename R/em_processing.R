#' Greedy matching approach to align the coordinates in one matrix (compare_location) to another (ref_location).
#'
#' @param ref_location The matrix used as the reference.
#' @param compare_location The matrix to align to the reference.
#'
#' @returns A list containing the ref_location annd compare_location matrices with info on what coordinates match up and re-ordered to match the ref_location and a dataframe with how to the two matrices match up.
#' @export
#'
ranked_match_up <- function(ref_location, compare_location) {
  #recover()

  colnames(ref_location) <- c('x', 'y')
  colnames(compare_location) <- c('x', 'y')

  dist_mat <- matrix(data = NA,
                     nrow = nrow(ref_location),
                     ncol = nrow(compare_location))

  for (REAL in 1:nrow(ref_location)) {
    for (EST in 1:nrow(compare_location)) {
      dist_mat[REAL,EST] <- stats::dist(rbind(ref_location[REAL,], compare_location[EST,]))
    }
  }

  matchup_df <- data.frame(ref_group = rep(NA, ncol(dist_mat)),
                           compare_group = rep(NA, ncol(dist_mat))
  )

  total_length <- nrow(dist_mat)
  row_ind <- seq_len(total_length)
  col_ind <- seq_len(total_length)

  for (i in seq_len(total_length)) {
    #find the smallest distance
    min_coord <- which(dist_mat == min(dist_mat), arr.ind = TRUE)[1,,drop=FALSE] #just take the first row, necessary for ties

    ref_ind <- row_ind[min_coord[1]]
    compare_ind <- col_ind[min_coord[2]]

    matchup_df[i,] <- c(ref_ind, compare_ind)


    if (i == total_length) break

    #remove groups reference and comparison groups that are already matched
    dist_mat <- dist_mat[-min_coord[1],-min_coord[2],drop=FALSE]
    row_ind <- row_ind[-min_coord[1]]
    col_ind <- col_ind[-min_coord[2]]
  }

  ref_df <- as.data.frame(ref_location)
  ref_df[,'original_group'] <- seq_len(nrow(ref_df))
  ref_df[,'ref_group'] <- seq_len(nrow(ref_df))

  ref_df_reorder <- ref_df[matchup_df$ref_group,]

  compare_df <- as.data.frame(compare_location)
  compare_df[,'original_group'] <- seq_len(nrow(compare_df))
  compare_df_reorder <- compare_df[matchup_df$compare_group,]
  compare_df_reorder[,'ref_group'] <- ref_df_reorder$ref_group

  ref_final <- ref_df_reorder[order(ref_df_reorder$ref_group),]
  compare_final <- compare_df_reorder[order(compare_df_reorder$ref_group),]

  #reset rownames
  rownames(ref_final) <- NULL
  rownames(compare_final) <- NULL

  return(
    list(matchup_df,
         ref_group_reorder = ref_final,
         compare_group_reorder = compare_final
    )
  )
}


#' Aligns up multiple matrices relative to a reference matrix (the first matrix in the input list).
#'
#' @param x_list List of matrices containing the coordinate matrices. The first matrix in the list represents the reference.
#'
#' @returns A list containing the matched up matrices.
#' @export
#'
ranked_match_up_multi <- function(x_list) {

  if (length(x_list) == 1) {
    x_list[[1]] <- as.data.frame(x_list[[1]])
    colnames(x_list[[1]]) <- c('x', 'y')
    x_list[[1]]$original_group <- 1:nrow(x_list[[1]])
    x_list[[1]]$ref_group <- 1:nrow(x_list[[1]])
    return(x_list)
  }

  x_reorder_list <- list()
  reorder_output <- ranked_match_up(x_list[[1]], x_list[[2]])
  x_reorder_list[[1]] <- reorder_output$ref_group_reorder
  x_reorder_list[[2]] <- reorder_output$compare_group_reorder

  if (length(x_list) > 2) {
    for (i in 3:length(x_list)) {
      x_reorder_list[[i]] <- ranked_match_up(x_list[[1]], x_list[[i]])$compare_group_reorder
    }
  }
  return(x_reorder_list)
}



#' Detects whether two matrices are duplicated based on their values. This is the matrix analog of the duplicated function in R.
#'
#' @param matrix_list List of matrices.
#' @param tol Tolerance value to decide whether two matrices are equal.
#'
#' @returns Boolean vector indicating whether the matrices are duplicated.
#' @export
#'
duplicated_matrix <- function(matrix_list, tol = 1e-8) {
  #matrix version of duplicated function to identify duplicated matrices
  equality_vec <- rep(FALSE, length(matrix_list))
  for (i in seq_along(matrix_list)) {
    if (i > 1) {
      for (x in seq_len(i - 1)) {
        #setting check.attributes to FALSE because this will return not equal if attributes (like row names)
        #are different, which we don't care about
        CHECK <- isTRUE(all.equal(matrix_list[[x]], matrix_list[[i]], tolerance = tol, check.attributes = FALSE))
        if (CHECK) {
          equality_vec[i] <- TRUE
          break
        }
      }
    }
  }

  return(equality_vec)
}


#' Convert labels to new labels based on a label map.
#'
#' @param old_label The original labels to convert.
#' @param label_map_original The original label values.
#' @param label_map_new The new labels to map the old values to. This vector is in the same order as label_map_original. So, for example, the value of the first element of label_map_original gets switched to the value of the first element of label_map_new.
#'
#' @returns A vector containing the new, mapped labels.
#' @export
#'
convert_labels <- function(old_label, label_map_original, label_map_new) {

  new_label <- stats::setNames(rep(NA, length(old_label)), nm = names(old_label))

  for (i in seq_along(label_map_original)) {
    original_label <- label_map_original[i]
    updated_label <- label_map_new[i]
    new_label[old_label == original_label] <- updated_label
  }

  return(new_label)
}

#' Convert the labels in the matrices of the membership based on the maps in the matchup list.
#'
#' @param membership_list List of matrices to convert labels (CHECK THAT THIS IS RIGHT).
#' @param matchup_list List of dataframes used to map the labels (CHECK THAT THIS IS RIGHT).
#'
#' @returns A list of converted dataframes(?) with converted labels.
#' @export
#'
convert_labels_multi <- function(membership_list, matchup_list) {

  stopifnot(length(membership_list) == length(matchup_list))

  if (length(membership_list) == 1) {
    return(membership_list)
  }

  membership_list_new <- list()
  membership_list_new[[1]] <- membership_list[[1]]
  for (i in 2:length(membership_list)) {
    membership_list_new[[i]] <- convert_labels(old_label = membership_list[[i]],
                                               label_map_original = matchup_list[[i]]$original_group,
                                               label_map_new = matchup_list[[i]]$ref_group)
  }
  return(membership_list_new)
}


#' Relabel the names of a multiple matrices using a list of map names (THE DOCUMENTATION NEEDS TO IMPROVE.
#'
#' @param multi_mat_list List of maps.
#' @param map_list List of maps.
#' @param order_sequentially Boolean indicating whether to order sequentially.
#'
#' @returns A list of relabeled matrices.
#' @export
#'
relabel_mat_names_list_multi <- function(multi_mat_list, map_list, order_sequentially = FALSE) {
  stopifnot(length(multi_mat_list) == length(map_list))

  output_list <- list()

  for (x in 1:length(multi_mat_list)) {
    output_list[[x]] <- relabel_mat_names_list(mat_list = multi_mat_list[[x]],
                                               map = map_list[[x]])

    if (isTRUE(order_sequentially)) {
      #all
      order_ind <- order(as.numeric(sub(".*loc_", "", rownames(output_list[[x]][[1]]))))
      output_list[[x]] <- lapply(output_list[[x]], function(x, ind) {x[ind,ind]}, ind = order_ind)

    }

  }
  return(output_list)
}


#' Relabel a matrix of names (ADD COMPLETE DESCRIPTION).
#'
#' @param mat_list List of matrices to relabel names.
#' @param map map.
#'
#' @returns List of relabeled matrices.
#' @export
#'
relabel_mat_names_list <- function(mat_list, map) {

  for (i in 1:length(mat_list)) {
    rownames(mat_list[[i]]) <- relabel_loc(name_vec = rownames(mat_list[[i]]), map = map)
    colnames(mat_list[[i]]) <- relabel_loc(name_vec = colnames(mat_list[[i]]), map = map)
  }

  return(mat_list)
}



#' Relabel vector based on a map vector (ADD COMPLETE DESCRIPTION).
#'
#' @param name_vec The name vector to relabel
#' @param map  The map vector
#'
#' @returns The vector of mapped names
#' @export
#'
relabel_loc <- function(name_vec, map) {
  loc <- as.numeric(sub(".*loc_", "", name_vec))
  new_loc <- map$ref_group[match(loc, map$original_group)]

  new_name <- paste0(sub("loc_.*","", name_vec), paste0("loc_", new_loc))
  return(new_name)
}



#' Process expectation-maximization output.
#'
#' @param output_list Output list from (soft or hard) EM function.
#' @param lik_tol Threshold to determine whether likelihoods are equal.
#' @param duplicated_tol Threshold to determine whether coordinates are the same value.
#' @param assignment_approach Whether the EM used soft or hard node assignment.
#'
#' @returns Processed output of EM algorithm.
#' @export
#'
process_em_output <- function(output_list, lik_tol = 1e-6, duplicated_tol = 1e-12, assignment_approach) {
  #recover()
  assignment_approach <- match.arg(assignment_approach, choices = c("soft", "hard"), several.ok = FALSE)

  like_list <- lapply(output_list, function(x) x$like_vec) #likelihood vecs from each run
  top_like_vec <- sapply(like_list, function(x) max(x, na.rm = TRUE)) #top likelihood from each run
  max_indices <- which(abs(top_like_vec - max(top_like_vec)) <= lik_tol) #indices of top likelihoods across runs
  top_runs <- output_list[max_indices] #subset output to top runs

  #extract the top likelihood index from each run
  x_list <- list()
  group_membership_list <- list()
  for (i in 1:length(top_runs)) {
    top_ind <- which.max(top_runs[[i]]$like_vec)[1] #if there are multiple equal, top indices, choose the first one
    x_list[[i]] <- matrix(top_runs[[i]]$xvec_list[[top_ind]], ncol = 2, byrow = 2)
    group_membership_list[[i]] <- top_runs[[i]]$group_membership_list[[top_ind]] #node assignments from the top lik step
  }

  #greedy match up approach of groups based on the matrix
  xlist_ranked_matchup <- ranked_match_up_multi(x_list = x_list)

  convert_labels_list <-  switch(assignment_approach,
                                 'soft' = {
                                   relabel_mat_names_list_multi(multi_mat_list = group_membership_list,
                                                                map_list = xlist_ranked_matchup,
                                                                order_sequentially = TRUE)
                                 },
                                 'hard' = {
                                   convert_labels_multi(membership_list = group_membership_list,
                                                        matchup_list = xlist_ranked_matchup)
                                 },
                                 stop('Not a valid assignment approach.'))


  #convert_labels_list <- convert_labels_multi(membership_list = group_membership_list,
  #                                            matchup_list = xlist_ranked_matchup)

  #remove duplicates
  #compare matrices but only compare the x and y coordinates. The reference column should
  #match up but the original groups column shouldn't be compared.
  #DATE: 10/16/2026 -- Fixing an (seemingly) inconsequential but unnecessary step
  #of passing xlist_ranked_matchup through ranked_match_up_multi again. I am
  #simplifying down to just extracting the x and y coordinates from xlist_ranked_matchup
  #duplicated_indices <- duplicated_matrix(matrix_list = lapply(ranked_match_up_multi(xlist_ranked_matchup), function(x) x[,c('x', 'y')]),
  #                                        tol = duplicated_tol)
  duplicated_indices <- duplicated_matrix(matrix_list = lapply(xlist_ranked_matchup, function(x) x[,c('x', 'y')]),
                                          tol = duplicated_tol)
  xlist_remove_duplicates <- xlist_ranked_matchup[!duplicated_indices]
  group_membership_list_remove_duplicates <- group_membership_list[!duplicated_indices]
  convert_labels_list_remove_duplicates <- convert_labels_list[!duplicated_indices]

  extra_info <- list(maxlik_indices = max_indices,
                     max_lik_group_membership_list_nonreorder = group_membership_list_remove_duplicates,
                     nonduplicated_indices = max_indices[!duplicated_indices],
                     group_membership_list = group_membership_list)

  return(
    list(
      max_lik_xmat_list = xlist_remove_duplicates,
      max_lik_group_membership_list = convert_labels_list_remove_duplicates,
      max_lik = max(top_like_vec),
      extra_info = extra_info
    )
  )
}
