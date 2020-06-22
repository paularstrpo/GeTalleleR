#' Create bins with a farey sequence
#' @param n A number.
#' @return A list with the bin edges and bin centres for the resulting sequence.
farey_bins <- function(n){
    # original matlab code:
    # bin_centres=farey_sequence(n);

    # bin_edges = bin_centres(:)';

    # binwidth = diff(bin_edges);
    # bin_edges = [0, bin_edges(1:end-1)+binwidth/2, 1];
    # bin_edges = full(real(bin_edges));

    # % Shift bins so the interval is ( ] instead of [ ) for

    # bin_edges = bin_edges + eps(bin_edges);
    # bin_edges(1) = 0;
    # bin_edges(end) = 1;

    # Use the farey() function from the elliptic package
    bin_centres <- elliptic::farey(n)
    bin_edges <- as.numeric(t(bin_centres))
    binwidth <- as.numeric(diff(bin_edges))

    bin_edges <- c(0, bin_edges[-1]+binwidth/2, 1)
    bin_edges <- as.matrix(as.numeric(bin_edges))

    # Shift bins so the interval is ( ] instead of [ ) 
    bin_edges <- bin_edges + 0L
    bin_edges[c(1, length(bin_edges))] <- c(0, 1)

    return(list(edges=bin_edges, centres=bin_centres))
}