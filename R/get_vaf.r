#' Calculate VAF
#' 
#' @param ref A numeric vector with number of reference reads for a set of loci
#' @param alt A numeric vector with numer of variant reads for a set of loci
#' @return The Variant Allele Frequency of a locus given reference and altered (variant) reads for a given locus
#' @examples
#' alt <- sample(10:1000, 5)
#' ref <- sample(10:1000, 5)
#' vaf <- get_vaf(ref, alt)
get_vaf <- function(ref, alt){
    return(alt / (ref + alt))
}

