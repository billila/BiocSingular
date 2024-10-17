#' @export
#' @importFrom BiocParallel bpstart bpstop SerialParam
#' @importFrom utils head
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom RSpectra svds
runArpackSVD <- function(x, k=min(dim(x)), nu=k, nv=k, center=FALSE, scale=FALSE, deferred=FALSE, fold=Inf, BPPARAM=SerialParam())
{
  if (!is(BPPARAM, "SerialParam"))
    stop("Parallel Computation is not supported")
  
  if (!missing(fold))
    stop("Fold argument is not used")
  
  checked <- check_numbers(x, k=k, nu=nu, nv=nv)
  k <- checked$k
  nv <- checked$nv
  nu <- checked$nu
  
  # Setting up the parallelization environment.
  old <- getAutoBPPARAM()
  setAutoBPPARAM(BPPARAM)
  on.exit(setAutoBPPARAM(old))
  
  if (!.bpisup2(BPPARAM)) {
    bpstart(BPPARAM)
    on.exit(bpstop(BPPARAM), add=TRUE)
  }
  
  x <- standardize_matrix(x, center=center, scale=scale, deferred=deferred, BPPARAM=BPPARAM)
  if (use_crossprod(x, fold)) {
    res <- svd_via_crossprod(x, k=k, nu=nu, nv=nv, FUN=RSpectra::svds)
  } else {
    res <- RSpectra::svds(as.matrix(x), k=k, nu=nu, nv=nv)
    res$d <- head(res$d, k)
    res <- standardize_output_SVD(res, x)
  }
  
  res
}
