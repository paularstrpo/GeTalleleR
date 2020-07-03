#' Fit V(pr) for a single layer to the model
#' @param x
#' @param models
#' @param bin_edges
#' @return The V(pr) for a single layer for the model.
fit_vpr <- function(x,models,bin_edges){
    Nx <- as.numeric(hist(x,
               breaks=100, #bin_edges,
               plot=FALSE)$counts)
    lx <- sum(Nx)
    csx <- cumsum(Nx/lx)
    csx_mat <- matrix(rep(csx,each=51),ncol = 51)
    emd_v <- abs(csx_mat - models[1:nrow(csx_mat),]) # coerce dimensions to work until the hist counts bug is fixed
    vpr <- (min(emd_v)+49)/100
    return(vpr)
}
