# find_win_and_vpr <- function(data,chr,start_pos,end_pos,layer,sensitivity,model_05,models_mat,bin_edges){
# # The authors make no representations about the suitability of this software for any purpose.
# # It is provided "as is" without express or implied warranty.
# #
# # input:
# #   data - matrix with chromosome, base pair positions and with VAF values
# #          for analysis from RNA2DNAlign
# #
# #   chr,start_pos,end_pos,layer - coordinates of the segment in the data:
# #                                   chr=chosen chromosome,
# #                                   start_pos=first base pair,
# #                                   end_pos  =last base pair,
# #                                   layer    =signal for windows generation
# #                                             (convention: Nex=1, Ntr=2, Tex=3, Ttr=4)
# #
# #
# # output:
# # results - structure with fields:
#
# #
# # If you have any questions please contact: p.m.slowinski@exeter.ac.uk
#
# # START
# # here we find indices that correspond to the chromosome that contains the segment
# results.chr <- chr
# results.layer <- layer
#
# idx_chr <- which(data[, 1]==chr)
# if (sum(diff(idx_chr)>1)>0){
#     print(paste0('warning: something is wrong with the chromosome: ', chr))
#     # quick check if a chromosome was in consecutive rows
# }
#
# # data on the chromosome containing the segment of interest
# data_chr <- data[idx_chr, 2:6]# data_chr is a matrix with base-pairs and all layers
#
# # grab data in the segment
# idx_sgmnt_start <- which(data_chr[,1] >= start_pos)# ,1,'first');
# idx_sgmnt_end <- which(data_chr[,1] <= end_pos)# ,1,'last');
#
# # processeing of segments with no data
# if (is.na(idx_sgmnt_start) || is.na(idx_sgmnt_length(data_chr)) || data_chr[idx_sgmnt_length(data_chr),1]<data_chr(idx_sgmnt_start,1)){
#     results.chr <- NA
#     results.layer <- NA
#     results.vpr_data <- NA
#     #disp('no data in segment')
#     return
# }
#
# data_sgmnt <- data_chr(idx_sgmnt_start:idx_sgmnt_},:)
# ds_size <- size(data_sgmnt,1)
#
# results.segment <- data_sgmnt
# results.total_length <- ds_size
# results.total_bp_length <- }_pos-start_pos
#
# ########## GENERATION OF WINDOWS ##########
# smpl <- data_sgmnt(:,layer+1)
#
# results.layer <- layer
# iptsT <- findchangepts(abs(smpl-0.5)+0.5,'MinThreshold',sensitivity)#could instead set 'MinDistance', Lmin, to 11 ponints
#
# idx_all_edges <- [1, t(iptsT), ds_size]
# idx_gap <- diff(idx_all_edges)<10){
#
# # always keep the start point and merge with the second window instead
# if (idx_gap(1)==1){
#    idx_gap(1) <- 0
#    idx_gap(2) <- 1
# }
#
# idx_all_edges(idx_gap) <- []
# idx_win_start <- idx_all_edges(1:}-1)
# idx_win_} <- [idx_all_edges(2:}-1)-1 ds_size]
#
# nb_of_windows <- numel(idx_win_start)
# ########## PUT DATA FROM WINDOWS INTO A STRUCTURE ##########
# m <- struct([])
# for (i in 1:nb_of_windows){
#     idx_window <- idx_win_start(i):idx_win_}(i)
#     m(i).win_data <- data_sgmnt(idx_window,:)
#     m(i).window <- idx_window
#     m(i).win_bp_length <- m(i).win_data(},1)-m(i).win_data(1,1)
# }
#
# ########## COMPUTE VPRS ##########
# #pre-fitting
# fitted_vprs <- zeros(1,numel(m))
#
# for (i in 1:numel(m)){
#     folded_m <- abs(m(i).win_data(:,layer+1)-0.5)+0.5
#     fitted_vprs(i) <- fit_vpr(folded_m,models_mat(layer).models,bin_edges)
# }
#
# [sorted_fitted_vprs, idx_sorted] <- sort(fitted_vprs)
# m <- m(idx_sorted)
# d_fitted_vprs <- diff(sorted_fitted_vprs)){
# idx_same_vprs=find(d_fitted_vprs<=0.011)
#
# while (!is.na(idx_same_vprs)){
#     idx_rm <- zeros(1,numel(m))
#
#     for (i in idx_same_vprs(}:-1:1)){
#         m(i).win_data <- [m(i).win_data; m(i+1).win_data]
#         m(i).window <- [m(i).window m(i+1).window]
#         m(i).win_bp_length <- m(i).win_bp_length+m(i+1).win_bp_length
#         idx_rm(i) <- i+1
#     }
#     m(idx_rm(idx_rm>0)) <- []
#
#     fitted_vprs <- zeros(1,numel(m))
#
#     for (i in 1:numel(m)){
#         folded_m <- abs(m(i).win_data(:,layer+1)-0.5)+0.5
#         fitted_vprs(i) <- fit_vpr(folded_m,models_mat(layer).models,bin_edges)
#     }
#
#     [sorted_fitted_vprs, idx_sorted] <- sort(fitted_vprs)
#      m <- m(idx_sorted)
#     d_fitted_vprs <- diff(sorted_fitted_vprs)){
#     idx_same_vprs=find(d_fitted_vprs<=0.011)
# }
#
# # fitting vpr in merged windows and analysis of the other layers
# all_l_cmb <- nchoosek(1:4,2)
#
# for (i in 1:numel(m)){
#     results.win_dp_length(i) <- numel(folded_m)
#     results.win_bp_length(i) <- m(i).win_bp_length
#     #fit vprs and check if different than model with vpr=0.5
#     for (k in 1:4){
#         folded_m <- abs(m(i).win_data(:,k+1)-0.5)+0.5
#         results.fitted_vprs(k,i) <- fit_vpr(folded_m,models_mat(k).models,bin_edges)
#         [!,results.pv_ks_05(k,i)] <- kstest2(folded_m,abs(model_05(k).model-0.5)+0.5)
#     }
#     for (k in 1:6){
#         s1 <- abs(m(i).win_data(:,all_l_cmb(k,1)+1)-0.5)+0.5
#         s2 <- abs(m(i).win_data(:,all_l_cmb(k,2)+1)-0.5)+0.5
#         [!,results.pv_ks_ll(k,i)] <- kstest2(s1,s2)
#     }
# }
#
# results.vpr_data <- m
# }
#
# ######### SOME EXTRA FUNCTIONS #########
# fit_vpr <- function(x1,models,bin_edges){
#
# Nx <- sparse(histc(x1,bin_edges))
# lx <- sum(Nx)
# csx <- t(cumsum(Nx/lx))
# csx_mat <- repmat(csx,51,1)
#
# emd_v <- sum(abs(csx_mat-models),2)
#
# [!,i_vpr] <- min(emd_v)
# vpr <- (i_vpr+49)/100
#
# }
