#' @title Support Vector Clustering
#'
#' @description
#' Implements Support Vector Clustering (SVC) as proposed by Ben-Hur, Horn,
#' Siegelmann and Vapnik (2001). The data are mapped to a high-dimensional
#' feature space by a kernel function, where the smallest hypersphere enclosing
#' the data is found by solving the dual quadratic program.
#' Support Vectors (SVs) lie on the surface of this sphere.
#' Paths between data points are checked for whether they leave the sphere;
#' points connected by paths entirely within the sphere are considered to be in the same cluster.
#'
#' @param data A numeric data frame, matrix or vector, with at least two rows and
#'        no missing or infinite values.
#' @param gamma The kernel width parameter: a single positive number, or `NULL`
#'        (default) or `"auto"` to estimate it with the median heuristic described in Details.
#' @param scan_n A positive integer: the number of points randomly selected in the `[0,1]` interval.
#'        They define the positions checked on the line segment between each pair of data points
#'        (the same positions are used for every pair) to build the adjacency matrix indicating
#'        if two data points are in the same cluster.
#' @param kernel A character string for the type of kernel: "gaussian" (default), "laplacian", or "cauchy".
#' @param n_obs An integer value of at least 2. If `n_obs` is specified and less than the total number of rows
#'        in `data`, SVC is applied to a random sample of `n_obs` rows, and the remaining rows are
#'        classified with an SVM trained on the SVC clusters (if SVC finds a single cluster, all of them
#'        are assigned to it). Use this for large datasets for faster processing
#'        (default is NULL, meaning all data is used for SVC).
#' @param verbose A boolean value. If TRUE (default), progress messages are printed.
#'        OSQP solver messages are always suppressed.
#' @param n_cores An integer specifying the number of cores to use for parallel processing.
#'        If NULL (default), all detected cores but one are used (two when the environment variable
#'        `_R_CHECK_LIMIT_CORES_` is set to `TRUE`, as during CRAN checks).
#'
#' @details
#' **Automatic gamma.** With `gamma = NULL` or `gamma = "auto"`, gamma is estimated on the rows
#' used for SVC as the inverse of the median pairwise squared distance (Gaussian and Cauchy
#' kernels) or of the median pairwise distance (Laplacian kernel, whose exponent uses the
#' distance itself), so that the estimate does not depend on the units of the data. For the
#' Gaussian kernel this is the middle value returned by `kernlab::sigest()` with
#' `scaled = FALSE`, computed over all pairs instead of a random sample.
#'
#' This heuristic sets the kernel width at the overall scale of the data, where SVC usually
#' finds a single cluster. Treat it as a starting point: Ben-Hur et al. (2001, Section 4.2)
#' start from a small gamma, where there is a single cluster, and increase it to observe the
#' formation of an increasing number of clusters.
#'
#' **Reproducibility.** The positions checked on the segments and, when `n_obs` is used, the
#' rows used for SVC are drawn with R's random number generator. Call [set.seed()] before
#' `SVC()` for reproducible results.
#'
#' **Parallel backend.** When more than one core is used, a PSOCK cluster is registered as the
#' `foreach` backend during the call. On exit, the cluster is stopped and the sequential backend
#' is registered.
#'
#' @returns A list containing:
#' \itemize{
#'    \item `clusters`: A numeric vector of cluster assignments for the rows SVC was run on,
#'          in the order given by `subset_index`.
#'    \item `exec_time`: Total execution time of the function.
#'    \item `SV_clustering`: A binary vector indicating if each row SVC was run on is an SVC SV,
#'          in the order given by `subset_index`.
#'    \item `SV_indices_svm`: Integer vector of indices (relative to the rows SVC was run on) of
#'          the SVM SVs, or `NULL` if `n_obs` was not used or SVC found a single cluster.
#'    \item `final_data`: `data` as a data frame, in its original row order, with a factor column
#'          `cluster_label` holding the cluster of every row (SVC clusters for the rows SVC was
#'          run on, SVM predictions for the other rows).
#'    \item `gamma_selected`: The gamma value used for clustering.
#'    \item `subset_index`: Integer vector with the rows of `data` SVC was run on;
#'          `seq_len(nrow(data))` when `n_obs` is not used.
#' }
#'
#' @references
#' Ben-Hur, A., Horn, D., Siegelmann, H. T., and Vapnik, V. (2001). Support vector clustering.
#' *Journal of Machine Learning Research*, 2, 125-137.
#' <https://www.jmlr.org/papers/v2/horn01a.html>
#'
#' @importFrom osqp osqp osqpSettings
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom e1071 svm
#' @importFrom stats dist median predict runif
#' @export
#'
#' @examples
#' # Two well-separated groups of points
#' set.seed(1)
#' x <- rbind(matrix(rnorm(40, mean = 0, sd = 0.3), ncol = 2),
#'            matrix(rnorm(40, mean = 3, sd = 0.3), ncol = 2))
#' res <- SVC(x, gamma = 1, scan_n = 5, n_cores = 1, verbose = FALSE)
#' table(res$clusters)
#' plot(x, col = res$clusters, pch = ifelse(res$SV_clustering == 1, 8, 19),
#'      main = "SVC clusters (stars: support vectors)")
#'
#' # Run SVC on a random subset of 25 rows and classify the other 15 with an SVM
#' res_sub <- SVC(x, gamma = 1, scan_n = 5, n_obs = 25, n_cores = 1, verbose = FALSE)
#' table(res_sub$final_data$cluster_label)
SVC <- function(data, gamma = NULL,
                scan_n,
                kernel = "gaussian",
                n_obs = NULL,
                verbose = TRUE,
                n_cores = NULL) {

  start_total_time <- Sys.time()

  is_count <- function(x, min = 1) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= min && x == round(x)
  }

  kernel <- match.arg(kernel, c("gaussian", "laplacian", "cauchy"))
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }
  if (missing(scan_n) || !is_count(scan_n)) {
    stop("`scan_n` must be a single positive integer.")
  }
  if (!is.null(n_obs) && !is_count(n_obs, min = 2)) {
    stop("`n_obs` must be NULL or a single integer greater than or equal to 2.")
  }
  if (is.null(gamma) || identical(gamma, "auto")) {
    gamma <- NULL
  } else if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0) {
    stop("`gamma` must be NULL, \"auto\" or a single positive number.")
  }

  if (is.null(n_cores)) {
    if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") == "TRUE") {
      actual_num_cores <- 2
    } else {
      detected_cores <- parallel::detectCores()
      actual_num_cores <- if (is.na(detected_cores)) 1 else max(1, detected_cores - 1)
    }
  } else if (is_count(n_cores)) {
    actual_num_cores <- n_cores
  } else {
    stop("`n_cores` must be NULL or a single positive integer.")
  }

  input_data_df <- as.data.frame(data)
  if (nrow(input_data_df) < 2L || ncol(input_data_df) < 1L) {
    stop("`data` must have at least two rows and one column.")
  }
  if (!all(vapply(input_data_df, is.numeric, logical(1)))) {
    stop("All columns of `data` must be numeric.")
  }
  if (!all(is.finite(as.matrix(input_data_df)))) {
    stop("`data` must not contain missing or infinite values.")
  }

  n_total <- nrow(input_data_df)
  if (!is.null(n_obs) && n_obs < n_total) {
    if (verbose) cat("Sampling", n_obs, "rows from the data for SVC.\n")
    subset_index <- sample.int(n_total, n_obs)
  } else {
    subset_index <- seq_len(n_total)
  }
  data_subset_svc <- input_data_df[subset_index, , drop = FALSE]
  data_subset_svc_numeric <- as.matrix(data_subset_svc)

  if (is.null(gamma)) {
    if (verbose) cat("gamma not provided. Estimating gamma with the median heuristic.\n")
    pairwise_dist <- as.numeric(dist(data_subset_svc_numeric))
    pairwise_dist <- pairwise_dist[pairwise_dist > 0]
    if (length(pairwise_dist) == 0L) {
      stop("Cannot estimate `gamma`: all rows used for SVC are identical.")
    }
    gamma <- if (kernel == "laplacian") 1 / median(pairwise_dist) else 1 / median(pairwise_dist^2)
    if (verbose) cat("Estimated gamma:", signif(gamma, 4), "\n")
  }

  selected_kernel_function <- switch(kernel,
                                     "gaussian" = function(sig_param, r1, r2) exp(-sig_param * sum((r1 - r2)^2)),
                                     "cauchy"   = function(sig_param, r1, r2) 1 / (1 + sig_param * sum((r1 - r2)^2)),
                                     "laplacian" = function(sig_param, r1, r2) exp(-sig_param * sqrt(sum((r1 - r2)^2))))

  run_core_svc_algorithm <- function(current_numeric_data, current_gamma_value, kernel_function, num_scan_points_param, n_cores_to_use) {
    num_points <- nrow(current_numeric_data)

    kernel_matrix <- outer(1:num_points, 1:num_points, Vectorize(function(i, j)
      kernel_function(current_gamma_value, current_numeric_data[i, ], current_numeric_data[j, ])))

    # Tight tolerances plus polishing (refinement of the ADMM solution on its active set).
    # With looser settings, small gamma values (nearly constant kernel matrix) leave spurious
    # nonzero betas above the 1e-5 cut-off: extra support vectors off the sphere surface and
    # spurious singleton clusters.
    osqp_model <- osqp::osqp(P = kernel_matrix, q = diag(kernel_matrix),
                             A = rbind(rep(1, num_points), diag(num_points)),
                             l = c(1, rep(0, num_points)), u = c(1, rep(1, num_points)),
                             pars = osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-08, eps_rel = 1e-08,
                                                       max_iter = 20000L, polishing = TRUE))

    solution <- osqp_model@Solve()
    if (solution$info$status_val != 1L) {
      warning("The QP solver did not converge (OSQP status: ", solution$info$status,
              "); support vectors and clusters may be inaccurate.", call. = FALSE)
    } else if (solution$info$status_polish != 1L) {
      warning("The QP solution could not be polished; support vectors and clusters may be ",
              "inaccurate. Consider a different `gamma`.", call. = FALSE)
    }
    betas <- solution$x
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
      on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      })
      doParallel::registerDoParallel(cl)
    } else {
      foreach::registerDoSEQ()
    }

    results <- foreach::foreach(j = 1:ncol(pairs), .combine = 'rbind') %dopar% {
      p1 <- current_numeric_data[pairs[1, j], ]; p2 <- current_numeric_data[pairs[2, j], ]
      # outer() keeps one row per scan point even for one-dimensional data
      seg <- outer(t_vals, p1) + outer(1 - t_vals, p2)
      if(max(apply(seg, 1, radius_sq_func)) < r_sq_limit) return(c(pairs[1, j], pairs[2, j]))
      return(NULL)
    }

    # A single connected pair comes back from foreach as a vector, not a matrix
    if (!is.null(results)) adj[matrix(results, ncol = 2)] <- 1
    adj <- adj + t(adj)
    list(clusters = igraph::components(igraph::graph_from_adjacency_matrix(adj, "undirected"))$membership,
         is_sv = is_sv)
  }

  main_res <- run_core_svc_algorithm(data_subset_svc_numeric, gamma, selected_kernel_function, scan_n, actual_num_cores)

  svc_labels <- factor(main_res$clusters)
  final_df <- input_data_df
  final_df$cluster_label <- factor(rep(NA, n_total), levels = levels(svc_labels))
  final_df$cluster_label[subset_index] <- svc_labels
  svm_idx <- NULL

  if (length(subset_index) < n_total) {
    out_index <- setdiff(seq_len(n_total), subset_index)
    if (nlevels(svc_labels) > 1L) {
      if (verbose) cat("Classifying out-of-sample points using SVM.\n")
      model_svm <- e1071::svm(x = data_subset_svc, y = svc_labels)
      final_df$cluster_label[out_index] <- stats::predict(model_svm, newdata = input_data_df[out_index, , drop = FALSE])
      svm_idx <- model_svm$index
    } else {
      if (verbose) cat("SVC found a single cluster; assigning it to the out-of-sample points.\n")
      final_df$cluster_label[out_index] <- levels(svc_labels)
    }
  }

  total_time <- Sys.time() - start_total_time
  if (verbose) cat("Total execution time:", format(total_time), "\n")

  return(list(clusters = main_res$clusters,
              exec_time = total_time,
              SV_clustering = as.numeric(main_res$is_sv),
              SV_indices_svm = svm_idx,
              final_data = final_df,
              gamma_selected = gamma,
              subset_index = subset_index))
}
