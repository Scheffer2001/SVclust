two_blobs <- function(n_per = 15, seed = 1) {
  set.seed(seed)
  rbind(matrix(rnorm(2 * n_per, 0, 0.3), ncol = 2),
        matrix(rnorm(2 * n_per, 3, 0.3), ncol = 2))
}

# TRUE when two labelings describe the same partition (up to relabeling)
same_partition <- function(a, b) {
  tb <- table(a, b)
  all(rowSums(tb > 0) == 1) && all(colSums(tb > 0) == 1)
}

test_that("SVC separates two well-separated groups", {
  x <- two_blobs()
  set.seed(2)
  res <- SVC(x, gamma = 1, scan_n = 5, n_cores = 1, verbose = FALSE)
  expect_true(same_partition(res$clusters, rep(1:2, each = 15)))
  expect_length(res$SV_clustering, 30)
  expect_true(all(res$SV_clustering %in% c(0, 1)))
  expect_identical(res$subset_index, 1:30)
  expect_null(res$SV_indices_svm)
  expect_true(same_partition(res$final_data$cluster_label, res$clusters))
})

test_that("gamma = 'auto' is the same as gamma = NULL", {
  x <- two_blobs()
  set.seed(3); a <- SVC(x, gamma = "auto", scan_n = 5, n_cores = 1, verbose = FALSE)
  set.seed(3); b <- SVC(x, gamma = NULL, scan_n = 5, n_cores = 1, verbose = FALSE)
  expect_identical(a$gamma_selected, b$gamma_selected)
  expect_identical(a$clusters, b$clusters)
  expect_null(names(a$gamma_selected))
})

test_that("automatic gamma is the median heuristic and does not depend on units", {
  x <- two_blobs()
  d <- as.numeric(dist(x))
  expect_equal(SVC(x, scan_n = 1, n_cores = 1, verbose = FALSE)$gamma_selected, 1 / median(d^2))
  expect_equal(SVC(x, scan_n = 1, kernel = "laplacian", n_cores = 1, verbose = FALSE)$gamma_selected,
               1 / median(d))
  for (k in c("gaussian", "laplacian", "cauchy")) {
    set.seed(4); a <- SVC(x, scan_n = 5, kernel = k, n_cores = 1, verbose = FALSE)
    set.seed(4); b <- SVC(100 * x, scan_n = 5, kernel = k, n_cores = 1, verbose = FALSE)
    expect_true(same_partition(a$clusters, b$clusters), info = k)
  }
})

test_that("small gamma gives a single cluster without spurious support vectors", {
  # Loose solver settings used to leave spurious nonzero betas at small gamma: support
  # vectors off the sphere surface and many singleton clusters.
  x <- two_blobs(n_per = 30)
  set.seed(11)
  res <- SVC(x, scan_n = 5, n_cores = 1, verbose = FALSE)
  expect_length(unique(res$clusters), 1)
  expect_lt(sum(res$SV_clustering), 10)
})

test_that("one-dimensional input gives the same clusters as the same points in 2-D", {
  set.seed(5)
  x1 <- c(rnorm(10, 0, 0.1), rnorm(10, 3, 0.1))
  set.seed(6); a <- SVC(x1, gamma = 1, scan_n = 5, n_cores = 1, verbose = FALSE)
  set.seed(6); b <- SVC(cbind(x1, 0), gamma = 1, scan_n = 5, n_cores = 1, verbose = FALSE)
  expect_true(same_partition(a$clusters, b$clusters))
  expect_true(same_partition(a$clusters, rep(1:2, each = 10)))
})

test_that("n_obs keeps the original row order and classifies the other rows", {
  x <- as.data.frame(two_blobs())
  rownames(x) <- paste0("obs", seq_len(nrow(x)))
  set.seed(7)
  res <- SVC(x, gamma = 2, scan_n = 5, n_obs = 20, n_cores = 1, verbose = FALSE)
  expect_length(res$subset_index, 20)
  expect_length(res$clusters, 20)
  expect_identical(rownames(res$final_data), rownames(x))
  expect_equal(res$final_data[, 1:2], x)
  expect_false(anyNA(res$final_data$cluster_label))
  expect_identical(as.integer(res$final_data$cluster_label[res$subset_index]),
                   as.integer(factor(res$clusters)))
  expect_true(same_partition(res$final_data$cluster_label, rep(1:2, each = 15)))
  expect_true(all(res$SV_indices_svm %in% seq_along(res$subset_index)))
})

test_that("n_obs assigns every row when SVC finds a single cluster", {
  x <- two_blobs()
  set.seed(8)
  res <- SVC(x, gamma = 1e-4, scan_n = 5, n_obs = 20, n_cores = 1, verbose = FALSE)
  expect_length(unique(res$clusters), 1)
  expect_null(res$SV_indices_svm)
  expect_false(anyNA(res$final_data$cluster_label))
  expect_length(unique(res$final_data$cluster_label), 1)
})

test_that("two points (a single pair) work", {
  res <- SVC(matrix(c(0, 0, 0.1, 0.1), nrow = 2, byrow = TRUE), gamma = 1, scan_n = 5,
             n_cores = 1, verbose = FALSE)
  expect_length(res$clusters, 2)
})

test_that("all kernels run", {
  x <- two_blobs()
  for (k in c("gaussian", "laplacian", "cauchy")) {
    expect_no_error(SVC(x, gamma = 1, scan_n = 3, kernel = k, n_cores = 1, verbose = FALSE))
  }
})

test_that("verbose = TRUE prints progress", {
  expect_output(SVC(two_blobs(), scan_n = 3, n_cores = 1), "Estimated gamma")
})

test_that("invalid inputs give informative errors", {
  x <- two_blobs()
  expect_error(SVC(x), "scan_n")
  expect_error(SVC(x, scan_n = 0), "scan_n")
  expect_error(SVC(x, scan_n = 2.5), "scan_n")
  expect_error(SVC(x, scan_n = 5, gamma = -1), "gamma")
  expect_error(SVC(x, scan_n = 5, gamma = "median"), "gamma")
  expect_error(SVC(x, scan_n = 5, kernel = "rbf"))
  expect_error(SVC(x, scan_n = 5, n_obs = 1), "n_obs")
  expect_error(SVC(x, scan_n = 5, n_cores = 0), "n_cores")
  expect_error(SVC(x, scan_n = 5, verbose = NA), "verbose")
  expect_error(SVC(x[1, , drop = FALSE], scan_n = 5), "two rows")
  expect_error(SVC(data.frame(a = 1:3, b = letters[1:3]), scan_n = 5), "numeric")
  expect_error(SVC(data.frame(a = c(1, NA, 3)), scan_n = 5), "missing")
  expect_error(SVC(matrix(1, 3, 2), scan_n = 5, verbose = FALSE), "identical")
})

test_that("parallel and sequential runs agree and leave the sequential backend registered", {
  skip_on_cran()
  x <- two_blobs()
  set.seed(9); a <- SVC(x, gamma = 1, scan_n = 5, n_cores = 1, verbose = FALSE)
  set.seed(9); b <- SVC(x, gamma = 1, scan_n = 5, n_cores = 2, verbose = FALSE)
  expect_identical(a$clusters, b$clusters)
  expect_identical(foreach::getDoParName(), "doSEQ")
})
