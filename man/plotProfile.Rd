\name{plotProfile}
\alias{plotProfile}
\title{
Graphically displays Wald or LOD statistics profiles of marker effects
}
\description{
Graphically displays Wald or LOD statistics profiles of marker effects
from an object of class \code{"gwasreml"}.   
}
\usage{
plotProfile(object, genObj, chr = names(object$QTL$geno), by.trait =
                 FALSE, LOD = TRUE, thresh = TRUE, use.dist = TRUE,
                 chr.lines = TRUE, annotate.peaks = TRUE)
}
\arguments{
  \item{object}{An object of class \code{"gwasreml"} obtained from
    running \code{gwasreml.asreml()} (see \code{?gwasreml.asreml} for
    more information). 
}
\item{genObj}{An object inheriting the class \code{"cross"} used in the
  call to \code{gwasreml.asreml()} (see \code{?gwasreml.asreml} for
    more information). 
}
\item{chr}{Subset of chromosomes to display.
}
\item{by.trait}{The type of profile to display. If \code{by.trait = FALSE} then
  the appropriate Wald statistics are displayed for the overall effect
  across \code{Traits}. If \code{by.trait = TRUE} then Wald statistics
  are displayed for each level of \code{Trait}. The latter is only possible if
  \code{main.effects = FALSE} was set in the original \code{gwasreml.asreml()} call.
}
  \item{LOD}{If \code{TRUE} then convert all Wald statistics to LOD
    scores. The default is \code{FALSE}.
}
  \item{thresh}{If \code{TRUE} then plot calculated significance
    thresholds. The default is \code{TRUE}.
}
  \item{use.dist}{If \code{TRUE} then use actual distances map
    distances. Setting to \code{FALSE} will just place markers
    equidistant from each other. The default = \code{TRUE}.
}
  \item{chr.lines}{If \code{TRUE} then plot vertical lines between
    chromosomes to enhance visual separation. The default is \code{TRUE}.
}
\item{annotate.peaks}{If \code{TRUE} then annotate putative QTL peaks on
  profiles. The default is \code{TRUE}.
}
}
\value{
\code{ggplot()} functionality is used to display Wald or LOD statistics
in the appropriate faceted panels with putative QTL annotated on the peaks.
}
\author{
Julian Taylor
}
\seealso{
\code{\link{gwasreml.asreml}}, \code{\link{summary.gwasreml}}
}
\examples{
\dontrun{
# convert "cross" object to "interval" object

mapi <- cross2int(map, impute = "Martinez", id = "Genotype")

# baseline multi-environment model

hwt.fa2 <- asreml(hectlitwt_500 ~ 1, random = ~ fa(Expt, 2):Genotype + at(Expt):Range
                  + at(Expt):Row, residual = ~ dsum(~ ar1(Row):ar1(Range) | Expt),
                  G.param = hwt.fa1$G.param, R.param = hwt.fa1$R.param,
                  data = heath, maxit = 30)

# Run scanning function across chromosome 2B

hwt.qtl <- gwasreml(hwt.fa2, mapi, merge.by = "Genotype", chr = "2B",
           Trait = "Expt", main.effects = TRUE, n.fa = 0, na.action =
           na.method(x = "include"))

# Plot main eeffct and interaction profiles with QTL peaks annotated.

plotProfile(hwt.qtl, mapi, LOD = TRUE, thresh = TRUE, annotate.peaks = TRUE)

}
}
\keyword{hplot}
