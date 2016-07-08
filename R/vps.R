
getWaveletFilters <- function(filter.number = 1, family = "DaubExPhase"){
  h <- filter.select(filter.number = filter.number,
                     family=family)$H
  len <- length(h)
  # Definition of g[n]:
  # g[n] = (-1)^(1-n)h[1-n]
  index <- 1 + - (len - 1):0
  g <- rev(h) * (-1)^(1 - index)
  list(h = h, g = g)
}

# expectedVps -------------------------------------------------------------
#' @export
expectedVps <- function(H, sigma2, Ts, filter.number, family, nlevels) {
  filters = getWaveletFilters(filter.number,family)
  filter_len = length(filters$g)

  # Normal convolution is
  # convolve(filters$g, rev(filters$g), type = "o")
  # but we have to convolve
  # convolve(filters$g, rev(rev(filters$g)), type = "o") and thus...
  g_acf = convolve(filters$g, filters$g, type = "o")
  h_acf = convolve(filters$h, filters$h, type = "o")

  vps = numeric(nlevels)
  max_D_index = (2 ^ nlevels - 1)*(length(filters$g) - 1) + 1
  D_index = (-max_D_index):max_D_index
  D = abs(D_index) ^ (2 * H)

  for (i in seq_len(nlevels)) {
    use_D = D[D_index %in% (-filter_len + 1):(filter_len - 1)]
    vps[[i]] = -sigma2 * Ts ^ (2 * H) /  2 * (use_D %*% g_acf)
    # convolution ranging from 2 * (-filter_len + 1) to 2 * (filter_len -1)
    # so we can downsanmple the signal using samples from 1, 3, 5, 7
    D =
      convolve(D, rev(h_acf),
               type = "o")#[seq(1, length.out = (2 * filter_len - 1), by = 2)]
    last_sample_index = (length(D) - 1) / 2
    D_index = (-last_sample_index):last_sample_index
    if (D_index[[1]] %% 2 == 0) {
      # if first index is pair, start there
      # End at zero to build a symmetric D signal, as supposed to be
      use_index = seq(1, which(D_index == 0), by = 2)
      first_D_index = D_index[[1]]
    } else {
      # End at zero to build a symmetric D signal, as supposed to be
      use_index = seq(2, which(D_index == 0), by = 2)
      first_D_index = D_index[[2]]
    }
    # Reflect the signal to get a symmetric D, as supposed to be
    D = c(D[use_index], rev(head(D[use_index], -1)))
    # D_index[[1]] / 2 is due to the downsampling
    D_index = seq(first_D_index / 2, length.out = length(D), by = 1)
  }
  vps
}
