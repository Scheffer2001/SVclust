utils::globalVariables(c("j_loop_var", "."))
#' @title Support Vector Clustering
#'
#' @description
#' This function implements Support Vector Clustering (SVC), a method based on support vector machine theory.
#' The algorithm maps data to a high-dimensional feature space using a kernel function.
#' In this feature space, it finds the smallest hypersphere enclosing the data.
#' Support Vectors (SVs) defining this sphere are identified.
#' Paths between data points are checked for whether they leave the sphere;
#' points connected by paths entirely within the sphere are considered to be in the same cluster.
#'
#' @param data A data frame or matrix with numerical data.
#' @param sigma A numerical value for the sigma parameter of the kernel function.
#'        If NULL (default), sigma is optimized using a grid search based on the `metric` chosen.
#' @param scan_n An integer value for the number of points randomly selected in a `[0,1]` interval.
#'        These points define line segments between pairs of data points and are used to
#'        build an adjacency matrix indicating if two data points are in the same cluster.
#' @param kernel A character string for the type of kernel: "gaussian" (default), "laplacian", or "cauchy".
#' @param n_row An integer value. If `n_row` is specified and less than the total number of rows in `data`,
#'        the dataset is split. SVC is applied to a random sample of `n_row` rows.
#'        The remaining `n - n_row` points are then classified using an SVM trained on the SVC clusters.
#'        Use this for large datasets for faster processing (default is NULL, meaning all data is used for SVC).
#' @param metric A character string for the internal clustering validation metric used for sigma optimization if `sigma` is NULL.
#'        Must be "silhouette" (default), "dunn", or "wss".
#' @param sigma_candidates A numeric vector of sigma values to test during grid search. If NULL (default),
#'        candidates are generated heuristically based on data distances with an expanded search range.
#' @param n_sigma_search An integer, the number of sigma values to generate for the heuristic grid search if
#'        `sigma_candidates` is NULL (default is 5).
#' @param verbose A boolean value. If TRUE (default), progress updates are printed.
#'        OSQP solver messages are always suppressed. If FALSE, these messages are suppressed.
#'
#' @returns A list containing:
#' \itemize{
#'    \item `clusters`: A numeric vector of cluster assignments for the data subset SVC was run on (X_sub).
#'    \item `exec_time`: Total execution time of the function.
#'    \item `SV_clustering`: A binary vector indicating if a point in X_sub is an SVC SV.
#'    \item `SV_indices_svm`: Integer vector of indices (relative to X_sub) of SVM SVs if n_row was used.
#'    \item `final_data`: A data frame, potentially X_sub combined with classified x_out points.
#'    \item `sigma_selected`: The final sigma value used for clustering.
#'    \item `metric_used_for_sigma_search`: The metric used if sigma search was performed.
#'    \item `best_metric_score`: The best score achieved during sigma search.
#'    \item `all_sigma_scores`: Data frame of all sigmas tested and their scores during search.
#' }
#' @importFrom osqp osqp osqpSettings
#' @importFrom clValid dunn
#' @importFrom dplyr sample_n anti_join bind_rows mutate select row_number
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom cluster silhouette
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom e1071 svm
#' @importFrom magrittr %>%
#' @importFrom stats dist runif predict median
#' @export
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("dplyr") && requireNamespace("ggplot2") &&
#'     requireNamespace("gridExtra")) {
#'   library(SVclust)
#'   data(iris)
#'   set.seed(123)
#'   iris_subset_data <- iris[sample(1:nrow(iris), 50), 1:4]
#'   res <- SVC(data = iris_subset_data, scan_n = 5, metric = "silhouette", n_sigma_search = 7)
#'   print(paste("Selected sigma (silhouette):", res$sigma_selected))
#'   print(res$all_sigma_scores)
#'
#'   # Plotting
#'   plot_df <- as.data.frame(iris_subset_data)
#'   plot_df$clusters_assigned <- res$clusters
#'   plot_df$is_SV_clustering <- as.logical(res$SV_clustering)
#'   library(ggplot2)
#'   p1 <- ggplot(plot_df, aes(x = Sepal.Length, y = Sepal.Width,
#'                                  color = factor(clusters_assigned),
#'                                  shape = factor(is_SV_clustering))) +
#'     geom_point(size = 3) +
#'     scale_shape_manual(values=c("FALSE"=19, "TRUE"=8)) +
#'     labs(title = "SVC on Iris Subset", color = "Cluster", shape = "SVC SV") +
#'     theme_minimal()
#'   print(p1)
#'  }
#' }
SVC <- function(data, sigma = NULL, scan_n, kernel = "gaussian", n_row = NULL,
                metric = "silhouette",
                sigma_candidates = NULL,
                n_sigma_search = 5,
                verbose = TRUE) {
  start_total_time <- Sys.time()
  calculate_metric_value <- function(cluster_assignments, data_for_metric, metric_name) {
    if (length(unique(cluster_assignments)) == 0) {
      return(ifelse(metric_name == "wss", Inf, ifelse(metric_name == "silhouette", -1, 0)))
    }
    if (nrow(data_for_metric) == 0) {
      return(ifelse(metric_name == "wss", Inf, ifelse(metric_name == "silhouette", -1, 0)))
    }
    if (length(unique(cluster_assignments)) < 2 && metric_name %in% c("silhouette", "dunn")) {
      return(ifelse(metric_name == "silhouette", -1, 0))
    }
    if (nrow(data_for_metric) <= 1 && metric_name %in% c("silhouette", "dunn")) {
      return(ifelse(metric_name == "silhouette", -1, 0))
    }
    distance_matrix <- NULL
    switch(metric_name,
           "silhouette" = {
             if (length(unique(cluster_assignments)) < 2 || length(unique(cluster_assignments)) > (nrow(data_for_metric) -1) ) return(-1)
             distance_matrix <- stats::dist(data_for_metric)
             sil_obj <- tryCatch(cluster::silhouette(as.integer(cluster_assignments), distance_matrix), error = function(e) NULL)
             if (is.null(sil_obj) || !is.matrix(sil_obj) || ncol(sil_obj) < 3) return(-1)
             return(mean(sil_obj[, 3], na.rm = TRUE))
           },
           "dunn" = {
             if (length(unique(cluster_assignments)) < 2) return(0)
             distance_matrix <- stats::dist(data_for_metric)
             dunn_value <- tryCatch(clValid::dunn(distance_matrix, as.integer(cluster_assignments)), error = function(e) 0)
             if(is.na(dunn_value) || !is.finite(dunn_value)) return(0)
             return(dunn_value)
           },
           "wss" = {
             if (length(unique(cluster_assignments)) == 0) return(Inf)
             total_wss <- 0
             for (k_value in unique(cluster_assignments)) {
               if (is.na(k_value)) next
               current_cluster_points <- data_for_metric[cluster_assignments == k_value, , drop = FALSE]
               if (nrow(current_cluster_points) > 0) {
                 cluster_center <- colMeans(current_cluster_points, na.rm = TRUE)
                 total_wss <- total_wss + sum(rowSums(sweep(current_cluster_points, 2, cluster_center, "-")^2, na.rm = TRUE), na.rm = TRUE)
               }
             }
             return(total_wss)
           },
           stop(paste("Invalid metric. Use 'silhouette', 'dunn', or 'wss'. Provided:", metric_name))
    )
  }
  input_data_df <- as.data.frame(data)
  if (is.vector(data) && ncol(input_data_df) == 1 && is.null(colnames(input_data_df))) {
    colnames(input_data_df) <- "V1"
  }
  data_subset_svc <- NULL
  data_out_of_sample <- NULL
  if (!is.null(n_row) && n_row > 0 && n_row < nrow(input_data_df)) {
    if (verbose) cat("Sampling", n_row, "rows from the data for SVC.\n")
    input_data_df_with_id <- input_data_df %>% dplyr::mutate(..temp_id_for_sampling = dplyr::row_number())
    data_subset_svc_with_id <- dplyr::sample_n(input_data_df_with_id, size = n_row)
    data_out_of_sample_with_id <- dplyr::anti_join(input_data_df_with_id, data_subset_svc_with_id, by = "..temp_id_for_sampling")
    data_subset_svc <- data_subset_svc_with_id %>% dplyr::select(-"..temp_id_for_sampling")
    data_out_of_sample <- data_out_of_sample_with_id %>% dplyr::select(-"..temp_id_for_sampling")
    if (verbose && (nrow(data_out_of_sample) + nrow(data_subset_svc) != nrow(input_data_df))) {
      cat("Warning: Row count mismatch after sampling. SVC subset:", nrow(data_subset_svc),
          "Out-of-sample:", nrow(data_out_of_sample), "Original:", nrow(input_data_df), "\n")
    }
  } else {
    if (verbose && !is.null(n_row) && (n_row <=0 || n_row >= nrow(input_data_df))) {
      cat("n_row is not in a valid range for splitting or is not needed. Using all data for SVC.\n")
    }
    data_subset_svc <- input_data_df
    data_out_of_sample <- NULL
  }
  data_subset_svc_numeric <- as.matrix(data_subset_svc)
  selected_kernel_function <- switch(kernel,
                                     "gaussian" = function(sig_param, r1, r2) exp(-sig_param * sum((r1 - r2)^2)),
                                     "cauchy" = function(sig_param, r1, r2) 1 / (1 + sig_param * sum((r1 - r2)^2)),
                                     "laplacian" = function(sig_param, r1, r2) exp(-sig_param * sqrt(sum((r1 - r2)^2))),
                                     stop("Invalid kernel. Use 'gaussian', 'laplacian', or 'cauchy'.")
  )
  run_core_svc_algorithm <- function(current_numeric_data, current_sigma_value, kernel_function, num_scan_points_param, num_connectivity_cores) {
    num_points_current_data <- nrow(current_numeric_data)
    if (num_points_current_data == 0) {
      warning("No data points to cluster in core SVC algorithm.")
      return(list(clusters = integer(0), lagrange_multipliers = numeric(0), is_support_vector = logical(0), radius_value = NA, adjacency_matrix = matrix(0,0,0)))
    }
    if (num_points_current_data == 1) {
      return(list(clusters = 1L, lagrange_multipliers = 1, is_support_vector = FALSE, radius_value = 0, adjacency_matrix = matrix(0,1,1)))
    }
    kernel_matrix <- outer(1:num_points_current_data, 1:num_points_current_data, Vectorize(function(i, j) kernel_function(current_sigma_value, current_numeric_data[i, ], current_numeric_data[j, ])))
    P_matrix_osqp <- kernel_matrix
    q_vector_osqp <- diag(kernel_matrix)
    A_matrix_constraints <- rbind(
      rep(1, num_points_current_data),
      diag(num_points_current_data)
    )
    lower_bounds_constraints <- c(1, rep(0, num_points_current_data))
    upper_bounds_constraints <- c(1, rep(1, num_points_current_data))
    osqp_model <- osqp::osqp(P = P_matrix_osqp, q = q_vector_osqp,
                             A = A_matrix_constraints,
                             l = lower_bounds_constraints, u = upper_bounds_constraints,
                             pars = osqp::osqpSettings(verbose = FALSE, eps_abs=1e-05, eps_rel=1e-05, max_iter=10000L))
    osqp_solution <- osqp_model$Solve()
    if(osqp_solution$info$status_val != 1 && osqp_solution$info$status_val != 2 && verbose) {
      cat("Warning: OSQP solver did not converge optimally for sigma =", current_sigma_value, ". Status:", osqp_solution$info$status, "\n")
    }
    lagrange_multipliers_betas <- osqp_solution$x
    lagrange_multipliers_betas[lagrange_multipliers_betas < 1e-5] <- 0
    lagrange_multipliers_betas[lagrange_multipliers_betas > (1 - 1e-5)] <- 1
    is_sv_on_boundary <- (lagrange_multipliers_betas > 1e-5) & (lagrange_multipliers_betas < (1 - 1e-5))
    is_sv_all_types <- lagrange_multipliers_betas > 1e-5
    term3_sphere_equation <- lagrange_multipliers_betas %*% kernel_matrix %*% lagrange_multipliers_betas
    radius_squared_function <- function(point_vector) {
      k_point_point <- kernel_function(current_sigma_value, point_vector, point_vector)
      k_data_point_vector <- sapply(1:num_points_current_data, function(i) kernel_function(current_sigma_value, current_numeric_data[i, ], point_vector))
      term2_sphere_equation <- -2 * sum(lagrange_multipliers_betas * k_data_point_vector)
      return(as.numeric(k_point_point + term2_sphere_equation + term3_sphere_equation))
    }
    radius_sq_values_at_data_points <- apply(current_numeric_data, 1, radius_squared_function)
    radius_sq_threshold_value <- NA
    if (sum(is_sv_on_boundary) > 0) {
      radius_sq_threshold_value <- mean(radius_sq_values_at_data_points[is_sv_on_boundary], na.rm = TRUE)
    } else if (sum(is_sv_all_types) > 0) {
      radius_sq_threshold_value <- mean(radius_sq_values_at_data_points[is_sv_all_types], na.rm = TRUE)
      if (verbose) cat("Note: No SVs strictly on boundary (0 < beta < 1) for sigma =", current_sigma_value, ". Using all SVs (beta > 0) for R threshold.\n")
    } else {
      if (verbose) cat("Warning: No Support Vectors found for sigma =", current_sigma_value, ". Clustering may be trivial.\n")
      radius_sq_threshold_value <- Inf
    }
    if(is.na(radius_sq_threshold_value) || !is.finite(radius_sq_threshold_value)){
      if(verbose) cat("Warning: R_sq_threshold is NA/Inf for sigma =", current_sigma_value, ". Defaulting to R_sq_threshold that likely forms one cluster.\n")
      radius_sq_threshold_value <- Inf
    }
    point_adjacency_matrix <- matrix(0, num_points_current_data, num_points_current_data)
    if (num_points_current_data > 1) {
      data_point_pairs <- utils::combn(num_points_current_data, 2)
      actual_num_cores <- min(num_connectivity_cores, parallel::detectCores() -1)
      actual_num_cores <- max(1, actual_num_cores)
      parallel_cluster_obj <- NULL
      t_values_for_scan <- stats::runif(num_scan_points_param)
      if (actual_num_cores > 1 && ncol(data_point_pairs) > 0) {
        parallel_cluster_obj <- parallel::makeCluster(actual_num_cores)
        doParallel::registerDoParallel(parallel_cluster_obj)
      } else {
        foreach::registerDoSEQ()
      }
      connectivity_check_results <- foreach::foreach(j_loop_var = 1:ncol(data_point_pairs), .combine = 'rbind', .packages=c("dplyr")) %dopar% {
        idx1 <- data_point_pairs[1, j_loop_var]
        idx2 <- data_point_pairs[2, j_loop_var]
        point_vector1 <- current_numeric_data[idx1, ]
        point_vector2 <- current_numeric_data[idx2, ]
        points_on_segment <- t(sapply(t_values_for_scan, function(t_val) t_val * point_vector1 + (1 - t_val) * point_vector2))
        if(!is.matrix(points_on_segment) && length(point_vector1) > 0) points_on_segment <- matrix(points_on_segment, ncol=length(point_vector1))
        max_radius_sq_on_segment <- -Inf
        if(nrow(points_on_segment) > 0){
          max_radius_sq_on_segment <- max(apply(points_on_segment, 1, radius_squared_function), na.rm = TRUE)
        }
        if (is.finite(max_radius_sq_on_segment) && max_radius_sq_on_segment < radius_sq_threshold_value) {
          return(c(idx1, idx2))
        } else {
          return(NULL)
        }
      }
      if(!is.null(parallel_cluster_obj)) parallel::stopCluster(parallel_cluster_obj)
      if (!is.null(connectivity_check_results) && nrow(connectivity_check_results) > 0) {
        for (i_result_row in 1:nrow(connectivity_check_results)) {
          point_adjacency_matrix[connectivity_check_results[i_result_row, 1], connectivity_check_results[i_result_row, 2]] <- 1
        }
      }
    }
    point_adjacency_matrix <- point_adjacency_matrix + t(point_adjacency_matrix)
    graph_object <- igraph::graph_from_adjacency_matrix(point_adjacency_matrix, mode = "undirected")
    connected_components <- igraph::components(graph_object)
    output_cluster_assignments <- connected_components$membership
    if(!is.null(rownames(current_numeric_data)) && length(output_cluster_assignments) == nrow(current_numeric_data)) {
      names(output_cluster_assignments) <- rownames(current_numeric_data)
    }
    return(list(clusters = output_cluster_assignments,
                lagrange_multipliers = lagrange_multipliers_betas,
                is_support_vector = is_sv_all_types,
                radius_value = sqrt(radius_sq_threshold_value),
                adjacency_matrix = point_adjacency_matrix))
  }
  best_sigma_parameters <- list(sigma = NULL, score = NULL, clusters = NULL, lagrange_multipliers = NULL, is_support_vector = NULL, radius_value = NULL)
  if (is.null(sigma)) {
    if (verbose) cat("Sigma not provided. Starting grid search for optimal sigma using '", metric, "' metric.\n")
    if (is.null(sigma_candidates)) {
      if (nrow(data_subset_svc_numeric) > 1 && ncol(data_subset_svc_numeric) > 0) {
        pairwise_distances <- stats::dist(data_subset_svc_numeric)
        median_distance <- stats::median(pairwise_distances[pairwise_distances > 0], na.rm = TRUE)
        base_sigma_heuristic_value <- 1.0
        if (!is.na(median_distance) && median_distance > 1e-9) {
          temp_heuristic <- 1.0
          if (kernel == "gaussian") {
            temp_heuristic <- 1 / (2*median_distance^2)
          } else if (kernel == "laplacian") {
            temp_heuristic <- 1 / median_distance
          } else if (kernel == "cauchy") {
            temp_heuristic <- 1 / (median_distance^2)
          }
          if (is.finite(temp_heuristic) && temp_heuristic > 1e-4 && temp_heuristic < 1e4) {
            base_sigma_heuristic_value <- temp_heuristic
          } else if (verbose) {
            cat("Note: Heuristic from median_dist (",temp_heuristic,") was extreme. Using default base_sigma=1.0 instead.\n")
          }
        } else {
          if(verbose) {
            cat("Warning: Median distance is zero, NA, or very small. Using default base_sigma = 1.0 for heuristic.\n")
          }
        }
        if(!is.finite(base_sigma_heuristic_value) || base_sigma_heuristic_value <= 1e-9) {
          if(verbose) cat("Warning: Calculated base_sigma_heuristic was problematic (", base_sigma_heuristic_value, "). Resetting to 1.0.\n")
          base_sigma_heuristic_value = 1.0
        }
        sigma_candidates <- base_sigma_heuristic_value * (10^seq(-2.0, 2.0, length.out = n_sigma_search))
      } else {
        sigma_candidates <- 10^seq(-2, 1, length.out = n_sigma_search)
      }
      sigma_candidates <- pmax(1e-6, sigma_candidates)
      if (verbose) cat("Generated sigma candidates:", paste(signif(sigma_candidates,3), collapse=", "), "\n")
    } else {
      if (verbose) cat("Using provided sigma candidates:", paste(signif(sigma_candidates,3), collapse=", "), "\n")
    }
    best_metric_score_value <- ifelse(metric %in% c("wss"), Inf, -Inf)
    num_cores_for_sigma_search_conn <- max(1, floor((parallel::detectCores() -1) / 2) )
    all_sigma_search_scores_df <- data.frame()
    for (sigma_test_val in sigma_candidates) {
      if (verbose) cat("  Testing sigma =", signif(sigma_test_val,4), "...")
      core_svc_results_gs <- run_core_svc_algorithm(data_subset_svc_numeric, sigma_test_val, selected_kernel_function, scan_n, num_cores_for_sigma_search_conn)
      current_metric_score <- calculate_metric_value(core_svc_results_gs$clusters, data_subset_svc, metric_name = metric)
      all_sigma_search_scores_df <- rbind(all_sigma_search_scores_df, data.frame(sigma=sigma_test_val, score=current_metric_score, num_clusters=length(unique(core_svc_results_gs$clusters))))
      if (verbose) cat(" Score:", signif(current_metric_score,4), ", Clusters:", length(unique(core_svc_results_gs$clusters)), "\n")
      should_update_best_score <- FALSE
      if (metric == "wss") {
        if (is.finite(current_metric_score) && current_metric_score < best_metric_score_value) should_update_best_score <- TRUE
      } else {
        if (is.finite(current_metric_score) && current_metric_score > best_metric_score_value) should_update_best_score <- TRUE
      }
      if ((metric == "wss" && !is.finite(best_metric_score_value) && is.finite(current_metric_score)) ||
          (metric != "wss" && !is.finite(best_metric_score_value) && is.finite(current_metric_score))) {
        should_update_best_score <- TRUE
      }
      if (should_update_best_score) {
        best_metric_score_value <- current_metric_score
        best_sigma_parameters$sigma <- sigma_test_val
        best_sigma_parameters$clusters <- core_svc_results_gs$clusters
        best_sigma_parameters$lagrange_multipliers <- core_svc_results_gs$lagrange_multipliers
        best_sigma_parameters$is_support_vector <- core_svc_results_gs$is_support_vector
        best_sigma_parameters$radius_value <- core_svc_results_gs$radius_value
      }
    }
    if (is.null(best_sigma_parameters$sigma)) {
      if(nrow(all_sigma_search_scores_df) > 0 && any(is.finite(all_sigma_search_scores_df$score))) {
        stop("Sigma grid search failed. While some scores might be finite, no optimal sigma was selected. Try different sigma_candidates or check data/metric behavior.")
      } else {
        stop("Sigma grid search failed. No valid (finite) scores were generated. Check metric calculation or data.")
      }
    }
    sigma <- best_sigma_parameters$sigma
    if (verbose) {
      cat("Grid search completed. Best sigma:", signif(sigma,4), "with", metric, "score:", signif(best_metric_score_value,4), "\n")
      cat("Grid search summary:\n")
      print(all_sigma_search_scores_df)
    }
    final_cluster_assignments <- best_sigma_parameters$clusters
    final_lagrange_multipliers <- best_sigma_parameters$lagrange_multipliers
    final_is_support_vector_flags <- best_sigma_parameters$is_support_vector
  } else {
    if (verbose) cat("Using provided sigma =", sigma, "\n")
    num_cores_for_main_run_conn <- max(1, parallel::detectCores() - 1)
    main_core_svc_run_results <- run_core_svc_algorithm(data_subset_svc_numeric, sigma, selected_kernel_function, scan_n, num_cores_for_main_run_conn)
    final_cluster_assignments <- main_core_svc_run_results$clusters
    final_lagrange_multipliers <- main_core_svc_run_results$lagrange_multipliers
    final_is_support_vector_flags <- main_core_svc_run_results$is_support_vector
    best_metric_score_value <- NA
    all_sigma_search_scores_df <- NULL
  }
  if(length(final_cluster_assignments) != nrow(data_subset_svc)){
    warning("Mismatch between length of final cluster assignments and rows in SVC data subset. This should not happen.")
  }
  final_output_dataframe <- NULL
  svm_support_vector_indices <- NULL
  if (!is.null(data_out_of_sample) && nrow(data_out_of_sample) > 0) {
    if (verbose) cat("Classifying out-of-sample points using SVM.\n")
    svm_training_data <- data_subset_svc
    if (length(final_cluster_assignments) == nrow(svm_training_data)) {
      svm_training_data$cluster_label <- factor(final_cluster_assignments)
    } else {
      warning("Length of final cluster assignments does not match rows in SVM training data. Cannot proceed with SVM.")
      final_output_dataframe <- data_subset_svc
      if(length(final_cluster_assignments) == nrow(final_output_dataframe)) {
        final_output_dataframe$cluster_label <- factor(final_cluster_assignments)
      } else {
        final_output_dataframe$cluster_label <- NA
      }
      svm_support_vector_indices <- NULL
    }
    if (!is.null(svm_training_data$cluster_label) && length(unique(svm_training_data$cluster_label)) > 1 && nrow(svm_training_data) > (ncol(svm_training_data) -1) ) { # ncol needs to account for added cluster_label
      svm_training_features <- svm_training_data[, !(names(svm_training_data) %in% "cluster_label"), drop=FALSE]
      svm_classification_model <- tryCatch({
        e1071::svm(x = svm_training_features, y = svm_training_data$cluster_label, probability = FALSE)
      }, error = function(e) {
        if (verbose) cat("Error training SVM for out-of-sample classification:", e$message, "\n")
        return(NULL)
      })

      if (!is.null(svm_classification_model)) {
        predicted_labels_out_of_sample <- predict(svm_classification_model, newdata = data_out_of_sample)
        classified_out_of_sample_data <- data_out_of_sample
        classified_out_of_sample_data$cluster_label <- predicted_labels_out_of_sample
        data_subset_svc_with_svc_labels <- data_subset_svc
        data_subset_svc_with_svc_labels$cluster_label <- factor(final_cluster_assignments)
        final_output_dataframe <- dplyr::bind_rows(data_subset_svc_with_svc_labels, classified_out_of_sample_data)
        svm_support_vector_indices <- svm_classification_model$index
      } else {
        if (verbose) cat("SVM training failed. Out-of-sample points will not be assigned to clusters by SVM.\n")
        data_subset_svc_with_svc_labels <- data_subset_svc
        data_subset_svc_with_svc_labels$cluster_label <- factor(final_cluster_assignments)
        unclassified_out_of_sample_data <- data_out_of_sample
        unclassified_out_of_sample_data$cluster_label <- NA
        final_output_dataframe <- dplyr::bind_rows(data_subset_svc_with_svc_labels, unclassified_out_of_sample_data)
        svm_support_vector_indices <- NULL
      }
    } else {
      if (verbose && (!is.null(svm_training_data$cluster_label) || length(final_cluster_assignments) > 0) ) cat("Not enough classes or data to train SVM for out-of-sample points.\n")
      else if (verbose) cat("SVM training data preparation issue for out-of-sample points.\n")

      data_subset_svc_with_svc_labels <- data_subset_svc
      if(length(final_cluster_assignments) == nrow(data_subset_svc_with_svc_labels)){
        data_subset_svc_with_svc_labels$cluster_label <- factor(final_cluster_assignments)
      } else {
        data_subset_svc_with_svc_labels$cluster_label <- NA
      }
      unclassified_out_of_sample_data <- data_out_of_sample
      unclassified_out_of_sample_data$cluster_label <- NA
      final_output_dataframe <- dplyr::bind_rows(data_subset_svc_with_svc_labels, unclassified_out_of_sample_data)
      svm_support_vector_indices <- NULL
    }
  } else {
    final_output_dataframe <- data_subset_svc
    if (length(final_cluster_assignments) == nrow(final_output_dataframe)) {
      final_output_dataframe$cluster_label <- factor(final_cluster_assignments)
    } else {
      final_output_dataframe$cluster_label <- NA
      warning("Final cluster assignments length mismatch even with no out-of-sample data.")
    }
    svm_support_vector_indices <- NULL
  }
  end_total_time <- Sys.time()
  total_execution_time <- end_total_time - start_total_time
  if (verbose) cat("Total execution time:", format(total_execution_time), "\n")
  return(list(
    clusters = final_cluster_assignments,
    exec_time = total_execution_time,
    SV_clustering = ifelse(final_is_support_vector_flags, 1, 0),
    SV_indices_svm = svm_support_vector_indices,
    final_data = final_output_dataframe,
    sigma_selected = sigma,
    metric_used_for_sigma_search = if(is.null(best_sigma_parameters$sigma) && is.null(all_sigma_search_scores_df)) NA else metric,
    best_metric_score = if(is.null(best_sigma_parameters$sigma) && is.null(all_sigma_search_scores_df)) NA else best_metric_score_value,
    all_sigma_scores = all_sigma_search_scores_df
  ))
}
