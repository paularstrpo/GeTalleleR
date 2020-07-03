#' Get Windows and V(pr) for a chromosome
#'
#' @param pos Integer or Numeric vector with genomic positions (within the same chromosome) to segment.
#' @param vaf Numeric vector of VAF (variant allele frequencies) for a given layer. Must correspond to sites in pos.
#' @param seg_start Where to start segmentation? Default is 1.
#' @param seg_end Where to end segmentation? Default is the maximum value of `pos`.
#' @param sensitivity Parameter for segmentation via findchangepoints
#' @param models_cdf cumulative density for a VAF model generated by gen_ideal_hist
#' @return a DF or list with the information: chromosome | layer_name | seg_start |  seg_end | genomic_length | num_sites | merged_window_idx | merged_window_start | merged_window_end | V(pr)
calc_window_vpr <- function(bin_edges){

    # compute V(pr)
    # fitted_vprs <- zeros(1,nrow(m))

    folded_m <- abs(m[, layer]-0.5)+0.5
    fitted_vprs <- fit_vpr(x = folded_m, models = models_mat[, layer], bin_edges = bin_edges)

    # change to return a df of layer + positons in each segments
    # [sorted_fitted_vprs, idx_sorted] <- order(fitted_vprs)

    m <- m[, idx_sorted]
    d_fitted_vprs <- diff(sorted_fitted_vprs)
    idx_same_vprs <- which(d_fitted_vprs<=0.011)

    # cuts data into segments using the change points
    # if : segment length is less than threshold, and/or variance of V(pr) is too low,
    # then: add it to the previous segment.
    while (!is.na(idx_same_vprs)){
        idx_rm <- zeros(1,numel(m))

        for (i in idx_same_vprs(}:-1:1)){
            m(i).win_data <- [m(i).win_data; m(i+1).win_data]
            m(i).window <- [m(i).window m(i+1).window]
            m(i).win_bp_length <- m(i).win_bp_length+m(i+1).win_bp_length
            idx_rm(i) <- i+1
        }
        # m(idx_rm(idx_rm>0)) <- []
        for (i in 1:nrow(m)){
            folded_m <- abs(m(i).win_data(:,layer+1)-0.5)+0.5 # change to use layer map convention
            fitted_vprs(i) <- fit_vpr(folded_m,models_mat(layer).models,bin_edges)
        }

        [sorted_fitted_vprs, idx_sorted] <- order(fitted_vprs)
        m <- m[idx_sorted, ]
        d_fitted_vprs <- diff(sorted_fitted_vprs))
        idx_same_vprs <- which(d_fitted_vprs<=0.011)

# TODO: MOVE STATISTICS TO ITS OWN FUNCTION
# fitting vpr in merged windows and analysis of the other layers
# change this to be able to analyze/compare to a user specified number of layers

# compare all pairs of layers
all_l_cmb <- combn(nlayers, 2)

# Fit V(pr) after segmentation and compare to other layers - code
# This loop iterates over all the segments that were found to check:
#  - how many samples
#  - length of segment in terms of a) num points and b) genomic length
for (i in 1:numel(m)){
    res$win_dp_length(i) <- numel(folded_m) # keep this
    res$win_bp_length(i) <- m(i).win_bp_length # keep this
    # fit vprs and check if different than model with vpr=0.5
    # fits V(pr) again to available layers & run KS test to check whether distributions are different
    for (k in layer_map){ #
        folded_m <- abs(m(i).win_data(:,k+1)-0.5)+0.5 # change to reflect layer map & keep

        # fits V(pr) for a single layer to the corresponding CDF models
        res$fitted_vprs(k,i) <- fit_vpr(folded_m,models_mat(k).models,bin_edges) # keep this

        # KS test of whether your layer is different than the null, which we expect to be 0.5 throughout
        [!,res$pv_ks_05(k,i)] <- kstest2(folded_m,abs(model_05(k).model-0.5)+0.5)# move to own fun
    }

    # move to own function
    for (k in 1:6){
        s1 <- abs(m(i).win_data(:,all_l_cmb(k,1)+1)-0.5)+0.5
        s2 <- abs(m(i).win_data(:,all_l_cmb(k,2)+1)-0.5)+0.5
        # KS test of the originally segmented layer versus all the other layers -d ifferent or same dist?
        [!,res$pv_ks_ll(k,i)] <- kstest2(s1,s2)
    }
}


    }