# try the check for fold = 1 
# Tests runArpackSVD().
# library(testthat); library(BiocSingular); source("setup.R"); source("test-arpack-svd.R")


library(testthat)
library(RSpectra)  

set.seed(80000)
test_that("Arpack SVD works on input matrices", {
  y <- matrix(rnorm(50000), ncol=200)
  set.seed(100)
  out <- runArpackSVD(y, k=5, nu=5, nv=5)
  set.seed(100)
  ref <- RSpectra::svds(as.matrix(y), k=5, nu=5, nv=5)
  expect_equal_svd(out, ref[c("d", "u", "v")])
  
  # Handles truncation.
  set.seed(100)
  out <- runArpackSVD(y, k=10, nv=5, nu=3)
  set.seed(100)
  ref <- RSpectra::svds(as.matrix(y), k=10, nv=5, nu=3)
  ref$v <- ref$v[, 1:5]
  expect_equal_svd(out, ref[c("d", "u", "v")])
})

set.seed(80001)
test_that("Random SVD works on thin matrices", {
  y <- matrix(rnorm(10000), ncol=10)
  set.seed(200)
  out <- runArpackSVD(y, k=3)
  set.seed(200)
  ref <- RSpectra::svds(as.matrix(y), k=3, nv=3, nu=3)
  expect_equal_svd(out, ref)
  
  # Handles truncation.
  set.seed(200)
  out <- runArpackSVD(y, k=5, nv=3, nu=2)
  set.seed(200)
  ref <- RSpectra::svds(as.matrix(y), k=5, nu=2, nv=3)
  ref$v <- ref$v[,1:3]
  expect_equal_svd(out, ref)
  
  set.seed(200)
  out <- runArpackSVD(y, k=2, nv=2, nu=2)
  set.seed(200)
  ref <- RSpectra::svds(as.matrix(y), k=2, nu=2, nv=2)
  ref$d <- ref$d[1:2]
  expect_equal_svd(out, ref)
})

set.seed(80002)
test_that("Arpack SVD works on fat matrices", {
  y <- matrix(rnorm(10000), nrow=10)
  
  # Test for k = 4
  set.seed(300)
  out <- runArpackSVD(y, k=4)
  set.seed(300)
  ref <- RSpectra::svds(y, k=4, nu=4, nv=4)  
  expect_equal_svd(out, ref)
  
  # Handles truncation
  set.seed(300)
  out <- runArpackSVD(y, k=5, nv=4, nu=2)
  set.seed(300)
  ref <- RSpectra::svds(y, k=5, nv=4, nu=2)  
  expect_equal_svd(out, ref)
  
  # Test for k = 6
  set.seed(300)
  out <- runArpackSVD(y, k=6, nu=2, nv=6)
  set.seed(300)
  ref <- RSpectra::svds(y, k=6, nu=2, nv=6)  
  expect_equal_svd(out, ref)
})


set.seed(80004)
test_that("Arpack SVD works with centering and scaling", {
  y <- matrix(rnorm(10000), ncol=50)
  center <- runif(ncol(y))
  scale <- runif(ncol(y))
  
  set.seed(100)
  out <- runArpackSVD(y, k=5, center=center, scale=scale)
  set.seed(100)
  ry <- scale(y, center=center, scale=scale)
  ref <- RSpectra::svds(ry, k=5, nv=5, nu=5)
  expect_equal_svd(out, ref)
  
  # Works with the cross-product.
  y <- matrix(rnorm(10000), ncol=10)
  center <- runif(ncol(y))
  scale <- runif(ncol(y))
  
  ry <- scale(y, center=center, scale=scale)
  set.seed(200)
  ref <- RSpectra::svds(ry, k=6, nu=6, nv=6)
  set.seed(200)
  out <- runArpackSVD(y, k=6, center=center, scale=scale)
  expect_equal_svd(out, ref)
  
  # Works with the deferred operations. 
  set.seed(100)
  ref <- runArpackSVD(ry, k=7)
  set.seed(100)
  out <- runArpackSVD(y, k=7, center=center, scale=scale, deferred=TRUE)
  expect_equal_svd(out, ref)
})

set.seed(800041)
test_that("Arpack SVD handles named inputs", {
  y <- matrix(rnorm(10000), ncol=50)
  rownames(y) <- sprintf("THING_%i", seq_len(nrow(y)))
  colnames(y) <- sprintf("STUFF_%i", seq_len(ncol(y)))
  
  out <- runArpackSVD(y, k=3)
  expect_identical(rownames(out$u), rownames(y))
  expect_identical(rownames(out$v), colnames(y))
})

set.seed(80005)
test_that("Arpack SVD fails gracefully with silly inputs", {
  y <- matrix(rnorm(10000), ncol=50)
  expect_error(runArpackSVD(y, k=-1), "non-negative")
  expect_error(runArpackSVD(y, nu=-1), "non-negative")
  expect_error(runArpackSVD(y, nv=-1), "non-negative")
  
  expect_warning(runArpackSVD(y, k=1e6), "requested than available")
  expect_warning(runArpackSVD(y, nu=1e6), "requested than available")
  expect_warning(runArpackSVD(y, nv=1e6), "requested than available")
})

set.seed(8000001)
test_that("Arpack SVD not handles zeroes", {
  y <- matrix(rnorm(50000), ncol=200)
  expect_error(runArpackSVD(y, k=0, nv=0, nu=0), "try svd()")
  expect_error(runArpackSVD(y, k=0), "try svd()")
})


