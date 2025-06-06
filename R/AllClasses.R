#' @export
#' @import methods
setClass("BiocSingularParam", contains="VIRTUAL", slots=c(deferred="logical", fold="numeric"))

#' @export
setClass("ExactParam", contains="BiocSingularParam")

#' @export
setClass("IrlbaParam", contains="BiocSingularParam", slots=c(extra.work="integer", args="list"))

#' @export
setClass("RandomParam", contains="BiocSingularParam", slots=c(args="list"))

#' @export
setClass("FastAutoParam", contains="BiocSingularParam")

#' @export
setClass("LowRankMatrixSeed", slots=c(rotation="ANY", components="ANY"))

#' @export
setClass("ArpackParam", contains="BiocSingularParam", slots=c(args="list"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("LowRankMatrix",
    contains="DelayedMatrix",
    representation(seed="LowRankMatrixSeed")
)
