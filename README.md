\# SVclust: Support Vector Clustering in R



\[!\[License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



\*\*SVclust\*\* is an R package for non-linear cluster analysis using the \*\*Support Vector Clustering (SVC)\*\* algorithm. By mapping data into a high-dimensional feature space, SVC can identify clusters of arbitrary shapes that traditional methods like K-means often fail to capture.



This implementation leverages the \*\*OSQP\*\* solver for efficient optimization and provides automated hyperparameter tuning.



\## Features



\* \*\*Non-linear Boundaries:\*\* Effectively clusters data with complex geometries (e.g., concentric circles).

\* \*\*Automated Tuning:\*\* Grid search for the optimal kernel parameter ($\\sigma$) using Silhouette, Dunn, or WSS metrics.

\* \*\*Scalability:\*\* Hybrid approach using SVC on data subsets and SVM for out-of-sample classification.

\* \*\*Performance:\*\* High-speed Quadratic Programming (QP) solving via the `osqp` package.



\## Theory



Support Vector Clustering maps data points into a high-dimensional feature space using a non-linear kernel. The goal is to find the smallest enclosing sphere that contains most of the data points.



\### Primal Problem

The optimization seeks to minimize the radius $R$ while allowing for some points to lie outside the sphere through slack variables $\\xi\_i$:



$$R^2 + C \\sum\_{i=1}^n \\xi\_i$$



Subject to:

$$\\|\\Phi(\\mathbf{x}\_i) - \\mathbf{a}\\|^2 \\le R^2 + \\xi\_i, \\quad \\xi\_i \\ge 0$$



where $\\mathbf{a}$ is the center of the sphere and $C$ is the penalty constant.



\### Dual Problem

Using the Lagrangian formulation and the kernel trick, we solve the following dual problem to find the multipliers $\\beta\_i$:



$$\\max\_{\\beta} \\sum\_{i=1}^n \\beta\_i K(\\mathbf{x}\_i, \\mathbf{x}\_i) - \\sum\_{i,j=1}^n \\beta\_i \\beta\_j K(\\mathbf{x}\_i, \\mathbf{x}\_j)$$



Subject to:

$$0 \\le \\beta\_i \\le C, \\quad \\sum\_{i=1}^n \\beta\_i = 1$$







\### Distance and Connectivity

For any point $\\mathbf{x}$, its distance to the center of the sphere in the feature space is calculated by:



$$R^2(\\mathbf{x}) = K(\\mathbf{x}, \\mathbf{x}) - 2 \\sum\_{i=1}^n \\beta\_i K(\\mathbf{x}\_i, \\mathbf{x}) + \\sum\_{i,j=1}^n \\beta\_i \\beta\_j K(\\mathbf{x}\_i, \\mathbf{x}\_j)$$



Two points $\\mathbf{x}$ and $\\mathbf{x}'$ belong to the same cluster if the segment connecting them remains within the sphere:



$$\\|\\Phi(t\\mathbf{x} + (1-t)\\mathbf{x}') - \\mathbf{a}\\|^2 \\le R^2, \\quad \\forall t \\in \[0,1]$$







\## Installation



You can install the development version of \*\*SVclust\*\* directly from GitHub:



```r

\# If you don't have devtools installed:

\# install.packages("devtools")



devtools::install_github("Scheffer2001/SVclust")

```



\## Quick Start



```r

library(SVclust)

\# Create a sample non-linear dataset

set.seed(42)

n <- 100

theta <- seq(0, 2\*pi, length.out = n)

c1 <- cbind(cos(theta), sin(theta)) + matrix(runif(n\*2, -0.1, 0.1), ncol=2)

c2 <- cbind(2\*cos(theta), 2\*sin(theta)) + matrix(runif(n\*2, -0.1, 0.1), ncol=2)

df <- as.data.frame(rbind(c1, c2))



results <- SVC(data = df, scan\_n = 10)



print(results$gamma\_selected)

table(results$clusters)

```



\## References



1\. Ben-Hur, A., Horn, D., Siegelmann, H. T., \& Vapnik, V. (2001). Support vector clustering. \*Journal of Machine Learning Research\*, 2(Dec), 125-137.

2\. Schölkopf, B., Williamson, R. C., Smola, A. J., Shawe-Taylor, J., \& Platt, J. C. (1999). Support vector method for novelty detection. \*Advances in Neural Information Processing Systems\*, 12.

3\. Stellato, B., Banjac, G., Goulart, P., Bemporad, A., \& Boyd, S. (2020). OSQP: An Operator Splitting Solver for Quadratic Programs. \*Mathematical Programming Computation\*, 12(4), 637-672.

