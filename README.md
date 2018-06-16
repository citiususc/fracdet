# fracdet: estimation of deterministic and fractal components in time series

## What's `fracdet`?
`fracdet` is an R package implementing a simple, yet realistic, fractal-deterministic model
that may be able to capture key features of complex time series. The package also provides 
functionality that, under the assumptions of this model, is able to characterize the fractal
component and to provide an estimation of the deterministic components present in a given time 
series. All the theoretical background is detailed in the paper:

> García, C.A., Otero, A., Félix, P., Presedo, J. & Márquez D.G., (2016). **Simultaneous estimation of deterministic and fractal components in non-stationary time series**, Physica D: Nonlinear Phenomena
Volumes 374–375, 1 July 2018, Pages 45-57. [link](https://www.sciencedirect.com/science/article/pii/S0167278916304018)

## Installation
`fracdet` is not currently available on [CRAN](http://cran.r-project.org/), but it may be installed directly from github using [devtools](https://github.com/hadley/devtools).

```r
library("devtools")
install_version("FGN", "2.0-12")
install_github("citiususc/fracdet")
```

## Quickstart
For a quickstart guide, see this [vignette](https://citiususc.github.io/fracdet/pages/vignette.html).
