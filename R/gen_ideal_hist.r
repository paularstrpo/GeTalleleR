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

for (layer in layer_map){
    
    layer <- matrix[, layer]

    # is this grabbing the column with tumor / normal? confused about the original input format
    
    l1 <- layer+6
    l2 <- l1+4


    idx_c <- which(chr<=22)
    counts_vector <- matrix[idx_c, l1] + matrix[idx_c,l2]

    model_samples <- struct('smpl', cell(1,51))

    for (k in 50:100){
        nb_reads <- numel(counts_vector)

        p <- k/100
        p <- p*rep(1, 5000)
        new_count_idx <- randi(nb_reads,1,5000)
        new_count <- counts_vector(new_count_idx)
        y <- rbinom(new_count(:),p(:))
        surr1 <- y(:)./new_count(:)

        p <- 1-k/100
        p <- p*ones(1,5000)
        new_count_idx <- randi(nb_reads,1,5000)
        new_count <- counts_vector(new_count_idx)
        y <- rbinom(new_count(:),p(:))
        surr2 <- y(:)./new_count(:)

        model_samples(k-49).smpl <- [surr1; surr2]
    }

    mdl_cdf <- zeros(51,numel(bin_edges))
    for (z in 1:51){

        if (layer==1){
           idx <- model_samples(z).smpl<0.1 | model_samples(z).smpl>0.9
           model_samples(z).smpl(idx) <- []
        }

        f_sample <- abs(model_samples(z).smpl-0.5)+0.5
        Ni <- sparse(histc(f_sample,bin_edges))
        li <- sum(Ni)
        mdl_cdf(z,:) <- cumsum(Ni)/li
    }

    models_cdf(layer).models <- mdl_cdf
    model_05(layer).model <- model_samples(1).smpl
    all_models(layer).model <- model_samples
}
