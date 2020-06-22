gen_ideal_hist <- function(matrix, bin_edges, layer_map){

# Demand that the input matrix be a data frame with the following format:
# CHR | POS | REF | ALT | Ttr | Ntr | Nex | Tex

# Layer map must be a named character map column names to layer type.
# Example:
#                              Ttr                              Ntr
# "X013_Ttr_BRCA_TCGA.BH.A0B8.01A" "X013_Ntr_BRCA_TCGA.BH.A0B8.11A"
#                              Nex                              Tex
# "X013_Nex_BRCA_TCGA.BH.A0B8.11A" "X013_Tex_BRCA_TCGA.BH.A0B8.01A"

# extract chromosome
chr <- matrix[, 1]
model <- list()

for (layer in layer_map){

    layer <- matrix[, layer]

    # Is this grabbing the column with tumor / normal vafs (or reads)? confused about the original input format
    # l1 <- layer+6
    # l2 <- l1+4

    idx_c <- which(chr<=22)
    # total reads?
    # counts_vector <- matrix[idx_c, l1] + matrix[idx_c,l2]

    # TODO: maybe vectorize this with *apply and use boot package to bootstrap?
    # (https://stats.idre.ucla.edu/r/faq/how-can-i-generate-bootstrap-statistics-in-r/)
    # This way, user can manually specify number of bootstrap iterations..
    # model_samples <- list()

    # bootstrap binomial distributions for Vpr
    model_samples <- sapply(50:100, function(k){
        nb_reads <- as.numeric(counts_vector)

        p <- (k/100)*rep(1, 5000)
        new_count <- sample(nb_reads,5000, replace=TRUE)
        y <- rbinom(new_count, size=5000, prob=p)
        surr1 <- y/new_count

        p <- (1-k/100)*rep(1, 5000)
        new_count <- sample(nb_reads,5000, replace=TRUE)
        y <- rbinom(new_count, size=5000, prob=p)
        surr2 <- y/new_count

        return(c(surr1, surr2))
    })

    # mdl_cdf <- matrix(0, nrow=51, ncol=length(bin_edges))

    # remove homozygous sites from Nex
    if (names(layer) == 'Nex'){
        idx <- model_samples < 0.1 | model_samples > 0.9
        model_samples[idx] <- NA
    }

    # wait, is this just an ECDF for the models?
    mdl_cdf <- apply(model_samples, MARGIN = 2, function(z){
        f_sample <- abs(z-0.5)+0.5
        Ni <- sparse(histc(f_sample, bin_edges))
        li <- sum(Ni)
        cumsum(Ni)/li
    })

    # return the bootstrapped Vpr distribution and the cumulative density function (cdf) for this layer
    model[[names(layer)]] <- list(cdf=mdl_cdf, dist=model_samples)
}
# return results!
return(model)
}
