\name{summary.gwasreml}
\alias{summary.gwasreml}
\title{Summary and print methods for the class \code{"gwasreml"}.
}
\description{
Print a QTL summary from the \code{"gwasreml"} object in presentable
format.
}
\usage{
\method{summary}{gwasreml}(object, genObj, LOD = TRUE, \ldots)
}
\arguments{
  \item{object}{an object of class \code{"gwasreml"} (see Details)}
  \item{genObj}{a genetic data structure inheriting the class \code{"cross"}}
  \item{LOD}{logical value. If TRUE LOD scores for QTL are calculated, defaults to \code{TRUE}}
  \item{\ldots}{further arguments passed to or from other methods}
}
\details{
  It is important that the \code{genObj} is not missing in
  \code{summary.wgaim()} as it contains vital summary information about
  each of the QTL detected.
}
\value{A summary of the QTL component of the \code{"gwasreml"} object is
  printed to the screen. For each QTL detected, if an \code{"interval"}
  analysis was performed then \code{summary.gwasreml()} prints the
  chromosome, name and distance of each flanking marker, size,
  p-value and LOD score if desired. If a \code{"marker"} analysis was
  performed then the chromosome, name and distance of the associated
  marker, size, p-value and LOD score are printed. 
}
\author{Julian Taylor
}
\seealso{\code{\link{gwasreml.asreml}}}
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

# summarise

summary(hwt.qtl, mapi)
}
}
\keyword{regression}
\keyword{methods}
