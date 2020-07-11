\name{wgaim-internal}
\alias{constructCM}
\alias{getQTL}
\alias{gwasreml.default}
\alias{findQTLPeaks}
\alias{theme_scatter}
\title{Internal wgaim functions}
\description{
These are internal \code{gwasreml} functions that are not intended to be
directly called by the user.
}
\usage{
constructCM(genoData, scale.method = "diag")
getQTL(object, genObj)
findQTLPeaks(wald, qtl.exc, eff.exc)
theme_scatter(base_size = 11, base_family = "")
\method{gwasreml}{default}(baseModel, \ldots)
}
\author{Julian Taylor}
\keyword{internal}