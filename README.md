# SVclust: Support Vector Clustering in R

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**SVclust** is an R package for non-linear cluster analysis using the **Support Vector Clustering (SVC)** algorithm of Ben-Hur et al. (2001). By mapping data into a high-dimensional feature space, SVC can identify clusters of arbitrary shapes that traditional methods like K-means often fail to capture.

This implementation leverages the **OSQP** solver for efficient optimization and provides a heuristic starting value for the kernel width.

## Features

* **Non-linear Boundaries:** Effectively clusters data with complex geometries (e.g., concentric circles).
* **Automatic Starting Value for $\gamma$:** A median-distance heuristic (the middle value of `kernlab::sigest()`) estimates the kernel width, as a starting point for exploring larger values.
* **Scalability:** Hybrid approach using SVC on data subsets and SVM for out-of-sample classification.
* **Performance:** High-speed Quadratic Programming (QP) solving via the `osqp` package.

## Theory

Support Vector Clustering maps data points into a high-dimensional feature space using a non-linear kernel. The goal is to find the smallest enclosing sphere that contains most of the data points.

### Primal Problem

The optimization seeks to minimize the radius $R$ while allowing for some points to lie outside the sphere through slack variables $\xi_i$:

$$\min_{R, \mathbf{a}, \xi} \; R^2 + C \sum_{i=1}^n \xi_i$$

Subject to:

$$\|\Phi(\mathbf{x}_i) - \mathbf{a}\|^2 \le R^2 + \xi_i, \quad \xi_i \ge 0$$

where $\mathbf{a}$ is the center of the sphere and $C$ is the penalty constant.

### Dual Problem

Using the Lagrangian formulation and the kernel trick, we solve the following dual problem to find the multipliers $\beta_i$:

$$\max_{\beta} \sum_{i=1}^n \beta_i K(\mathbf{x}_i, \mathbf{x}_i) - \sum_{i,j=1}^n \beta_i \beta_j K(\mathbf{x}_i, \mathbf{x}_j)$$

Subject to:

$$0 \le \beta_i \le C, \quad \sum_{i=1}^n \beta_i = 1$$

The available kernels are the Gaussian $K(\mathbf{x}, \mathbf{x}') = e^{-\gamma \|\mathbf{x} - \mathbf{x}'\|^2}$, the Laplacian $e^{-\gamma \|\mathbf{x} - \mathbf{x}'\|}$ and the Cauchy $1 / (1 + \gamma \|\mathbf{x} - \mathbf{x}'\|^2)$.

### Distance and Connectivity

For any point $\mathbf{x}$, its distance to the center of the sphere in the feature space is calculated by:

$$R^2(\mathbf{x}) = K(\mathbf{x}, \mathbf{x}) - 2 \sum_{i=1}^n \beta_i K(\mathbf{x}_i, \mathbf{x}) + \sum_{i,j=1}^n \beta_i \beta_j K(\mathbf{x}_i, \mathbf{x}_j)$$

Two points $\mathbf{x}$ and $\mathbf{x}'$ belong to the same cluster if the segment connecting them remains within the sphere:

$$\|\Phi(t\mathbf{x} + (1-t)\mathbf{x}') - \mathbf{a}\|^2 \le R^2, \quad \forall t \in [0,1]$$

## Installation

You can install the development version of **SVclust** directly from GitHub:

```r
# If you don't have devtools installed:
# install.packages("devtools")

devtools::install_github("Scheffer2001/SVclust")
```

## Quick Start

```r
library(SVclust)

# Create a sample non-linear dataset: two concentric circles
set.seed(42)
n <- 100
theta <- seq(0, 2*pi, length.out = n)
c1 <- cbind(cos(theta), sin(theta)) + matrix(runif(n*2, -0.1, 0.1), ncol=2)
c2 <- cbind(2*cos(theta), 2*sin(theta)) + matrix(runif(n*2, -0.1, 0.1), ncol=2)
df <- as.data.frame(rbind(c1, c2))

results <- SVC(data = df, gamma = 2, scan_n = 10)

table(results$clusters)
#>
#>   1   2
#> 100 100
```

## Choosing gamma

`gamma` sets the scale at which the data is probed: small values give smooth boundaries and few clusters, large values tighter boundaries and more clusters. With `gamma = NULL` (the default), `SVC()` estimates it with the median heuristic, which places the kernel width at the overall scale of the data and usually yields a single cluster; on the circles above it gives `gamma = 0.26` and one cluster.

Following Ben-Hur et al. (2001, Section 4.2), use the automatic value as a starting point and increase `gamma` to observe the formation of clusters, preferring solutions with few support vectors (`mean(results$SV_clustering)`) that remain stable over a range of values. On the circles, the two circles appear at `gamma = 2` and remain at `gamma = 4`; at `gamma = 8` the outer circle starts to split.

The positions checked on each segment are drawn at random, so call `set.seed()` before `SVC()` for reproducible results.

## References

1. Ben-Hur, A., Horn, D., Siegelmann, H. T., & Vapnik, V. (2001). Support vector clustering. *Journal of Machine Learning Research*, 2(Dec), 125-137. <https://www.jmlr.org/papers/v2/horn01a.html>
2. Schölkopf, B., Williamson, R. C., Smola, A. J., Shawe-Taylor, J., & Platt, J. C. (1999). Support vector method for novelty detection. *Advances in Neural Information Processing Systems*, 12.
3. Stellato, B., Banjac, G., Goulart, P., Bemporad, A., & Boyd, S. (2020). OSQP: An Operator Splitting Solver for Quadratic Programs. *Mathematical Programming Computation*, 12(4), 637-672.
