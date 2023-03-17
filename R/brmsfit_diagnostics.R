#' check brmsfit diagnostics
#' @param brmsfit a brmsfit object
#' @return logical vector. First element is whether divergences passed; second
#' is whether rhat passed
#' @export
check_brmsfit_diagnostics <- function(brmsfit) {
  out <- rep(T, 2)
  out[1] <- rstan::get_num_divergent(brmsfit$fit) == 0
  out[2] <- max(brms::rhat(brmsfit)) < 1.05
  out
}
