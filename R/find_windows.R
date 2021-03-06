#' Get windows for V(pr) calculation in a chromosome
#'
#' @param pos Integer or Numeric vector with genomic positions (within the same chromosome) to segment.
#' @param vaf Numeric vector of VAF (variant allele frequencies) for a given layer. Must correspond to sites in pos.
#' @param seg_start Where to start segmentation? Default is 1.
#' @param seg_end Where to end segmentation? Default (0) uses the maximum value of `pos`.
#' @param sensitivity Penalty value for segmentation via `changepoint::cpt.meanvar`
#' @return a DF or list with information about windows and the vafs for SNPs within each window. (TODO: flesh this doc out more)

find_windows <- function(pos, vaf, seg_start=1, seg_end=0, sensitivity=0.2){
    # The authors make no representations about the suitability of this software
    # for any purpose. It is provided "as is" without express or implied warranty.

    # data on the chromosome containing the segment of interest
    data <- data.frame(pos=pos, vaf=vaf) # data is a matrix with base-pairs and all layers
    rm(pos, vaf)

    # if seg_end is not specified then use the largest positional value in the dataset.
    if (seg_end==0){
        seg_end <- max(as.numeric(data$pos))
    }

    # Grab data in the segment we are looking at
    seg_filt <- (data$pos >= seg_start) & (data$pos <= seg_end)

    if (length(data$pos[seg_filt])==0){
        return(print('Error: specified segment has zero sites.'))
    }

    # processeing of segments with no data
    data_sgmnt <- data[seg_filt,] # trim data to fit specified segment
    ds_size <- nrow(data_sgmnt) # length of segment in terms of sites

    # GENERATE WINDOWS
    # note: could instead set 'MinDistance', Lmin, to 11 ponints
    # Uses the changepoint package:
    iptsT <- changepoint::cpt.meanvar(abs(data_sgmnt$pos-0.5)+0.5, penalty='Manual', pen.value=sensitivity)
    # grab indeces for window edges
    idx_all_edges <- c(seg_start, iptsT@cpts[1:iptsT@ncpts.max], ds_size)

    # merge windows with less than a threshold number of points (or a threshold Vpr)
    # always keep the start point and merge with the second window instead.
    # TODO: refine this!
    idx_gap <- diff(idx_all_edges<10)
    if (idx_gap[1]==1){
        idx_gap[1] <- 0
        idx_gap[2] <- 1
    }

    # grab indeces for window starts & ends
    idx_win_start <- idx_all_edges[1:length(idx_all_edges)-1]
    idx_win_end <- idx_all_edges[2:length(idx_all_edges)]

    # shift all but the last window ends by one to avoid double-counting the breakpoint
    idx_win_end[1:(length(idx_win_end)-1)] <- idx_win_end[1:(length(idx_win_end)-1)] + 1

    # save window information!
    m <- data.frame(window_start=data$pos[idx_win_start],
                    window_end=data$pos[idx_win_end])

    # genomic length of window
    m$window_bp_length <- m$window_end - m$window_start

    # length in terms of number of SNPs residing in the window
    m$window_site_length <- sapply(1:nrow(m), function(x){
        length(data_sgmnt$pos[idx_win_start[x]:idx_win_end[x]])
        })

    # Append layer vafs
    # TODO: silence the warnings for this!
    mlist <- lapply(1:nrow(m), function(x){
        meta <- m[x,]
        # grab vafs that fall into a particular segment and
        # bind them to the window bound information
        vafs <-  data_sgmnt[data_sgmnt$pos >= meta$window_start &
                                data_sgmnt$pos <= meta$window_end,]
        return(cbind(meta, vafs))
    })

    # return results in list format
    res <- list(segment_inputs=data_sgmnt,
                total_length=ds_size,
                total_bp_length=seg_end-seg_start,
                num_windows=nrow(m),
                window_metadata=m,
                result=mlist)

    return(res)
}
