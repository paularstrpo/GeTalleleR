#' Fit V(pr) for a single layer to the model
#' @param x
#' @param models
#' @param bin_edges
#' @return The V(pr) for a single layer for the model.
fit_vpr <- function(x,models,bin_edges){
    Nx <- as.numeric(hist(x, breaks=bin_edges, plot=FALSE)$counts)
    lx <- sum(Nx)
    csx <- cumsum(Nx/lx)
    csx_mat <- rep(csx, times=51)
    emd_v <- sum(abs(csx_mat - models), 2)
    vpr <- (min(emd_v)+49)/100
    # return(min(emd_v))
    return(vpr)
}
