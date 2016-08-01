# fracdet constructor -----------------------------------------------------

#' @export
Fracdet = function(x, fbmPars) {
  UseMethod("Fracdet", x)
}

#' @export
Fracdet.wd = function(x, fbmPars) {
  attr(x, "WaveletVar") = WaveletVar.wd(x)
  # obtain a matrix where the first column represents the estimate of the fbmPars
  # and the second column the Std. Error
  parsMatrix = extractFbmPars(fbmPars)
  attr(x, "fbmPars") = parsMatrix
  class(x) = c("Fracdet", class(x))
  x
}

extractFbmPars = function(fbmPars) {
  UseMethod("extractFbmPars", fbmPars)
}

extractFbmPars.matrix = function(fbmPars) {
  if (!all(dim(fbmPars) == c(2, 2))) {
    stop("Incorrect fbmPars dimensions (expected 2x2)")
  }
  rownames(fbmPars) = c("H", "sigma2")
  colnames(fbmPars) = c("Estimate", "Std. Error")
  fbmPars
}

extractFbmPars.nls = function(fbmPars) {
  extractFbmParsFromFit(fbmPars)
}

extractFbmPars.lm = function(fbmPars) {
  extractFbmParsFromFit(fbmPars)
}


extractFbmParsFromFit = function(fbmPars) {
  pars = summary(fit)$coefficients[,1:2]
  stopifnot(rownames(pars) == c("H", "sigma2"))
  pars
}


# Basic functionality -----------------------------------------------------
#' @export
getWaveletVar = function(x) {
  UseMethod("getWaveletVar", x)
}

#' @export
getWaveletVar.Fracdet = function(x) {
  attr(x, "WaveletVar")
}

#' @export
getFbmPars = function(x) {
  UseMethod("getFbmPars", x)
}

#' @export
getFbmPars.Fracdet = function(x) {
  attr(x, "fbmPars")
}


#' @export
as.wd = function(x) {
  UseMethod("as.wd", x)
}

#' @export
as.wd.Fracdet = function(x){
  class(x) = "wd"
  x
}

# methods of class wd -----------------------------------------------------
# Required due to an error in the "wd" class, that checks the object class by
# direct comparation instead of using inherits

accessC.Fracdet = function(x, ...) { accessC.wd(as.wd(x), ...) }
accessD.Fracdet = function(x, ...) { accessD.wd(as.wd(x), ...) }
convert.Fracdet = function(x, ...) { convert.wd(as.wd(x), ...) }
draw.Fracdet = function(x, ...) { draw.wd(as.wd(x), ...) }
image.Fracdet = function(x, ...) { image.wd(as.wd(x), ...) }
IsEarly.Fracdet = function(x) { IsEarly.wd(as.wd(x)) }
LocalSpec.Fracdet = function(x, ...) { LocalSpec.wd(as.wd(x), ...) }
modernise.Fracdet = function(x, ...) { modernise.wd(as.wd(x), ...) }
nullevels.Fracdet = function(x, ...) { nullevels.wd(as.wd(x), ...) }
plot.Fracdet = function(x, ...) { plot.wd(as.wd(x), ...) }
print.Fracdet = function(x, ...) { print.wd(as.wd(x), ...) }
putC.Fracdet = function(x, ...) { putC.wd(as.wd(x), ...) }
putD.Fracdet = function(x, ...) { putD.wd(as.wd(x), ...) }
summary.Fracdet = function(x, ...) { summary.wd(as.wd(x), ...) }
threshold.Fracdet = function(x, ...) { threshold.wd(as.wd(x), ...) }
WaveletVar.Fracdet = function(x, ...) { WaveletVar.wd(as.wd(x), ...) }
wr.Fracdet = function(x, ...) { wr.wd(as.wd(x), ...) }


# estimateDetSignal -----------------------------------------------------------
