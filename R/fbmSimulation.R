#' @export 
fbmSim = function(n, H) {
  .Call('fracdet_simulateFbmCpp', PACKAGE = 'fracdet', n, H)
}

#' @export
#' @importFrom FGN SimulateFGN
fgnSim = function(n, H) {
  FGN::SimulateFGN(n,H)
}
