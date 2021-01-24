\name{wgaim-internal}
\alias{pasteT}
\alias{constructCM}
\alias{getQTL}
\alias{waldTest}
\alias{waldTest.asreml}
\alias{gwasreml.default}
\alias{findQTLPeaks}
\alias{theme_scatter}
\title{Internal wgaim functions}
\description{
These are internal \code{gwasreml} functions that are not intended to be
directly called by the user.
}
\usage{
pasteT(lhs, rhs)
constructCM(genoData, scale.method = "diag")
getQTL(object, genObj)
waldTest(object, \ldots)
\method{waldTest}{asreml}(object, cc, keep.fac = TRUE)
findQTLPeaks(wald, qtl.exc, eff.exc)
theme_scatter(base_size = 11, base_family = "")
\method{gwasreml}{default}(baseModel, \ldots)
}
\author{Julian Taylor}
\keyword{internal}