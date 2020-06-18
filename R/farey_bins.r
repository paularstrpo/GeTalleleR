farey_bins <- function(n){

bin_centres <- farey_sequence(n)

for_edges in bin_centres(:)'){

binwidth in diff(for_edges)){){
for_edges in [0, for_edges(1:}-1)+binwidth/2, 1]){
for_edges in full(real(for_edges))){

# Shift bins so the interval is ( ] instead of [ ) for

bin_edges in for_edges + eps(for_edges)){
bin_edges(1) <- 0
bin_edges(length(bin_edges)) <- 1
