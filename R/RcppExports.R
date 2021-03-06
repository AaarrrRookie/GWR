# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

cross_mat <- function(Match, ii, p, Sv, n, nc) {
    .Call('GWR_cross_mat', PACKAGE = 'GWR', Match, ii, p, Sv, n, nc)
}

cross_mat_vec <- function(Match, Sv, Y, n, nc) {
    .Call('GWR_cross_mat_vec', PACKAGE = 'GWR', Match, Sv, Y, n, nc)
}

cross_vec <- function(Match, Sv, n, nc) {
    .Call('GWR_cross_vec', PACKAGE = 'GWR', Match, Sv, n, nc)
}

cross_vec_mat <- function(Match, Sv, Y, n, nc) {
    .Call('GWR_cross_vec_mat', PACKAGE = 'GWR', Match, Sv, Y, n, nc)
}

GWR_core1_Deck <- function(Match, D, Y, Z) {
    .Call('GWR_GWR_core1_Deck', PACKAGE = 'GWR', Match, D, Y, Z)
}

GWR_core2 <- function(Match, D, i, p, n, nc) {
    .Call('GWR_GWR_core2', PACKAGE = 'GWR', Match, D, i, p, n, nc)
}

GWR_core_time <- function(Match, D, Y) {
    .Call('GWR_GWR_core_time', PACKAGE = 'GWR', Match, D, Y)
}

MoranI <- function(Match, plz, Dplz, D, Y, TSS) {
    .Call('GWR_MoranI', PACKAGE = 'GWR', Match, plz, Dplz, D, Y, TSS)
}

MoranI_MC <- function(Match, plz, Dplz, D, Y, TSS) {
    .Call('GWR_MoranI_MC', PACKAGE = 'GWR', Match, plz, Dplz, D, Y, TSS)
}

