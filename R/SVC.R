utils::globalVariables(c("j_loop_var", "j", "."))
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
#' @param gamma A numerical value or "auto". If "auto" (default), gamma is estimated
#' using the inverse of median squared distance of our data. If a numeric value is provided, it is used directly.
#' @param scan_n An integer value for the number of points randomly selected in a `[0,1]` interval.
#'        These points define line segments between pairs of data points and are used to
#'        build an adjacency matrix indicating if two data points are in the same cluster.
#' @param kernel A character string for the type of kernel: "gaussian" (default), "laplacian", or "cauchy".
#' @param n_obs An integer value. If `n_obs` is specified and less than the total number of rows in `data`,
#'        the dataset is split. SVC is applied to a random sample of `n_obs` rows.
#'        The remaining `n - n_obs` points are then classified using an SVM trained on the SVC clusters.
#'        Use this for large datasets for faster processing (default is NULL, meaning all data is used for SVC).
#' @param verbose A boolean value. If TRUE (default), progress updates are printed.
#'        OSQP solver messages are always suppressed. If FALSE, these messages are suppressed.
#' @param n_cores An integer specifying the number of cores to use for parallel processing.
#' If NULL (default), it automatically detects available cores.
#'
#' @returns A list containing:
#' \itemize{
#'    \item `clusters`: A numeric vector of cluster assignments for the data subset SVC was run on (X_sub).
#'    \item `exec_time`: Total execution time of the function.
#'    \item `SV_clustering`: A binary vector indicating if a point in X_sub is an SVC SV.
#'    \item `SV_indices_svm`: Integer vector of indices (relative to X_sub) of SVM SVs if n_obs was used.
#'    \item `final_data`: A data frame, potentially X_sub combined with classified x_out points.
#'    \item `gamma_selected`: The final gamma value used for clustering.
#' }
#' @importFrom osqp osqp osqpSettings
#' @importFrom dplyr sample_n anti_join bind_rows mutate select row_number
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom e1071 svm
#' @importFrom magrittr %>%
#' @importFrom stats dist runif predict median na.omit quantile
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
#'   res <- SVC(data = iris_subset_data, scan_n = 5, n_cores = 2)
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
SVC <- function(data, gamma = NULL,
                scan_n,
                kernel = "gaussian",
                n_obs = NULL,
                verbose = TRUE,
                n_cores = NULL) {

  start_total_time <- Sys.time()

  if (is.null(n_cores)) {
    if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") == "TRUE") {
      actual_num_cores <- 2
    } else {
      actual_num_cores <- max(1, parallel::detectCores() - 1)
    }
  } else {
    actual_num_cores <- n_cores
  }

  input_data_df <- as.data.frame(data)
  if (is.vector(data) && ncol(input_data_df) == 1 && is.null(colnames(input_data_df))) {
    colnames(input_data_df) <- "V1"
  }

  data_subset_svc <- NULL
  data_out_of_sample <- NULL
  if (!is.null(n_obs) && n_obs > 0 && n_obs < nrow(input_data_df)) {
    if (verbose) cat("Sampling", n_obs, "rows from the data for SVC.\n")
    input_data_df_with_id <- input_data_df %>% dplyr::mutate(..temp_id_for_sampling = dplyr::row_number())
    data_subset_svc_with_id <- dplyr::sample_n(input_data_df_with_id, size = n_obs)
    data_out_of_sample_with_id <- dplyr::anti_join(input_data_df_with_id, data_subset_svc_with_id, by = "..temp_id_for_sampling")
    data_subset_svc <- data_subset_svc_with_id %>% dplyr::select(-"..temp_id_for_sampling")
    data_out_of_sample <- data_out_of_sample_with_id %>% dplyr::select(-"..temp_id_for_sampling")
  } else {
    data_subset_svc <- input_data_df
    data_out_of_sample <- NULL
  }

  data_subset_svc_numeric <- as.matrix(data_subset_svc)

  if (is.null(gamma)) {
    if (verbose) cat("gamma not provided. Estimating optimal gamma \n")
    sig_estimates <- as.numeric(dist(data_subset_svc_numeric)^2)
    gamma <- 1/quantile(sig_estimates[sig_estimates > 0], 0.5)
    if (verbose) cat("Estimated gamma:", signif(gamma, 4), "\n")
  }

  selected_kernel_function <- switch(kernel,
                                     "gaussian" = function(sig_param, r1, r2) exp(-sig_param * sum((r1 - r2)^2)),
                                     "cauchy"   = function(sig_param, r1, r2) 1 / (1 + sig_param * sum((r1 - r2)^2)),
                                     "laplacian" = function(sig_param, r1, r2) exp(-sig_param * sqrt(sum((r1 - r2)^2))),
                                     stop("Invalid kernel."))

  run_core_svc_algorithm <- function(current_numeric_data, current_gamma_value, kernel_function, num_scan_points_param, n_cores_to_use) {
    num_points <- nrow(current_numeric_data)
    if (num_points <= 1) return(list(clusters = 1L, is_support_vector = FALSE))

    kernel_matrix <- outer(1:num_points, 1:num_points, Vectorize(function(i, j)
      kernel_function(current_gamma_value, current_numeric_data[i, ], current_numeric_data[j, ])))

    osqp_model <- osqp::osqp(P = kernel_matrix, q = diag(kernel_matrix),
                             A = rbind(rep(1, num_points), diag(num_points)),
                             l = c(1, rep(0, num_points)), u = c(1, rep(1, num_points)),
                             pars = osqp::osqpSettings(verbose = FALSE, eps_abs=1e-05, eps_rel=1e-05))

    betas <- osqp_model@Solve()$x
    betas[betas < 1e-5] <- 0
    is_sv <- betas > 1e-5
    is_sv_boundary <- (betas > 1e-5) & (betas < (1 - 1e-5))

    term3 <- as.numeric(betas %*% kernel_matrix %*% betas)
    radius_sq_func <- function(p) {
      k_pp <- kernel_function(current_gamma_value, p, p)
      k_dp <- sapply(1:num_points, function(i) kernel_function(current_gamma_value, current_numeric_data[i, ], p))
      return(k_pp - 2 * sum(betas * k_dp) + term3)
    }

    r_sq_vals <- apply(current_numeric_data, 1, radius_sq_func)
    r_sq_limit <- if(any(is_sv_boundary)) mean(r_sq_vals[is_sv_boundary]) else mean(r_sq_vals[is_sv])

    pairs <- utils::combn(num_points, 2)
    adj <- matrix(0, num_points, num_points)
    t_vals <- stats::runif(num_scan_points_param)

    if (n_cores_to_use > 1) {
      cl <- parallel::makeCluster(n_cores_to_use)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
    } else {
      foreach::registerDoSEQ()
    }

    results <- foreach::foreach(j = 1:ncol(pairs), .combine = 'rbind') %dopar% {
      p1 <- current_numeric_data[pairs[1, j], ]; p2 <- current_numeric_data[pairs[2, j], ]
      seg <- t(sapply(t_vals, function(t) t * p1 + (1 - t) * p2))
      if(max(apply(seg, 1, radius_sq_func)) < r_sq_limit) return(c(pairs[1, j], pairs[2, j]))
      return(NULL)
    }

    if (!is.null(results)) for(i in 1:nrow(results)) adj[results[i,1], results[i,2]] <- 1
    adj <- adj + t(adj)
    list(clusters = igraph::components(igraph::graph_from_adjacency_matrix(adj, "undirected"))$membership,
         is_sv = is_sv)
  }

  main_res <- run_core_svc_algorithm(data_subset_svc_numeric, gamma, selected_kernel_function, scan_n, actual_num_cores)

  final_df <- data_subset_svc
  final_df$cluster_label <- factor(main_res$clusters)
  svm_idx <- NULL

  if (!is.null(data_out_of_sample)) {
    if (verbose) cat("Classifying out-of-sample points using SVM.\n")
    model_svm <- e1071::svm(x = data_subset_svc, y = final_df$cluster_label)
    preds <- stats::predict(model_svm, newdata = data_out_of_sample)
    out_df <- data_out_of_sample; out_df$cluster_label <- preds
    final_df <- dplyr::bind_rows(final_df, out_df)
    svm_idx <- model_svm$index
  }

  total_time <- Sys.time() - start_total_time
  if (verbose) cat("Total execution time:", format(total_time), "\n")

  return(list(clusters = main_res$clusters,
              exec_time = total_time,
              SV_clustering = as.numeric(main_res$is_sv),
              SV_indices_svm = svm_idx,
              final_data = final_df,
              gamma_selected = gamma))
}
