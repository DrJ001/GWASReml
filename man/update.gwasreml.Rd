\name{update.gwasreml}
\alias{update.gwasreml}
\title{
Update function for \code{"gwasreml"} objects.
}
\description{
Update function for \code{"gwasreml"} objects that respects previous
results from a \code{gwasreml()} run.
}
\usage{
\method{update}{gwasreml}(object, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"gwasreml"} that inherits the
    \code{"asreml"} class. 
}
\item{\dots}{Any arguments to \code{gwasreml.asreml} (see
  \code{?gwasreml.asreml} for more information). This can also be used
  to pass through additional \code{asreml()} or \code{asreml.options()}
  arguments (see either \code{?asreml} or \code{?asreml.options} for
  more information).
}
}
\details{
This is a simple updating function to be used in conjunction with
\code{gwasreml.asreml()} to assist in determining significant putative main
effect QTL and putative \code{Trait} by QTL interactions. There are three
very useful arguments from \code{gwasreml.asreml()} that can be used for
this task. The \code{gwasreml} object can be updated by adding another
chromosome using the \code{chr} argument. The function checks
results of previous runs and will only scan chromosomes that have
not been previously scanned. Results obtained from these new chromosomes
are appended to results from previously run chromosomes. If
at least one \code{qtl.window} and \code{effects.window} are altered in
\code{update.gwasreml()} then the putative QTL interaction and main
effects are recalculated using the new rules and the final model is
re-fitted with the amended putative QTL. 
}
\value{
Updated version of the final model with any necessary changes.
}
\author{
Julian Taylor
}
\seealso{
\code{\link{summary.gwasreml}}, \code{\link{gwasreml.asreml}}  
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

# Update with additional chromsomes 

hwt.qtl2 <- update(hwt.qtl, chr = c("1D","7A"))

# Update with additional chromsomes and alter QTL window

hwt.qtl2 <- update(hwt.qtl, chr = c("1D","7A"), qtl.window = 30)
}
}
\keyword{regression}
