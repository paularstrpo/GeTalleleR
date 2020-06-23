#' Generate V(pr) model distributions
#' @param total_reads a numeric vector with total reads per site for one layer.
#' @param bins_edges A vector of bin edges outputted from `farey_bins()`.
#' @param rm.hom Logical. Remove homozygous sites before making the cumulative density function? Default is TRUE.
#' @param nboots A number. How many bootstrap iterations to perform? Default is 10000.
#' @return A list with the bootstrapped V(pr) model distributions and empirical cumulative density functions for one layer.

gen_ideal_hist <- function(total_reads, bin_edges, nboots=10000, rm.hom=TRUE){

    # bootstrap binomial distributions for Vpr
    vaf_model <- sapply(50:100, function(k){
        nb_reads <- as.numeric(total_reads)

        p <- (k/100)*rep(1, nboots)
        new_count <- sample(nb_reads,nboots, replace=TRUE)
        y <- rbinom(new_count, size=nboots, prob=p)
        return(y/new_count)
    })

    # remove homozygous sites.
    if (rm.hom){
        idx <- vaf_model < 0.1 | vaf_model > 0.9
        vaf_model[idx] <- NA
    }

    # Compute ECDF manually to ensure the same order of arguments in computing EMD
    vaf_cdf <- apply(vaf_model, MARGIN = 2, function(z){
        # TODO: can't seem to get farey sequence bins to work here...
        # instead for now specifying number of bins to make.
        # For some reason this computes a different length result for each z...
        return(hist(z, length(bin_edges), plot=FALSE)$density) 
    })

    # return the bootstrapped Vpr distribution and the cumulative density function (cdf) for this layer
    # Only need first column from vaf_model - remove rest for memory reasons
    return(list(cdf=vaf_cdf, dist=vaf_model[, 1]))
}
