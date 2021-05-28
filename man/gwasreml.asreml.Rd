\name{gwasreml.asreml}
\alias{gwasreml.asreml}
\alias{gwasreml}
\title{Genome wide association analysis of single and multi-section
  experiments using ASReml-R V4
}
\description{
A comprehensive marker scanning function for genome wide association
analysis of single and multi-section (environment, trait, treatment) experiments using ASReml-R V4
}
\usage{
\method{gwasreml}{asreml}(baseModel, genObj, merge.by = NULL,
       Trait = NULL, covariate = NULL, type = "main+int",
       fix.lines = TRUE, chr = names(nmar(genObj)), gen.type =
       "interval", n.fa = 0, breakout = -1, thresh = "b-corr",
       TypeI = 0.05, qtl.window = 40, effects.window = 20,
       trace = TRUE, \ldots)
}
\arguments{
\item{baseModel}{a linear mixed model object of class \code{"asreml"} usually representing
    a base \code{asreml()} model to be extended.
}
\item{genObj}{a list object containing the genotypic data in the form of
  a \pkg{qtl} \code{"cross"} object (see Details).
}
 \item{merge.by}{a character string or name of the column(s) in
   \code{baseModel} phenotypic data \code{genObj} to merge the
   phenotypic and genotypic data sets.
 }
\item{Trait}{a character string column name in the phenotypic
       data used in \code{baseModel} to identify to the multiple
       sections of the data (see Details).
     }
\item{covariate}{a character string column name in the phenotypic
       data used in \code{baseModel} to identify the covariate
       to be used in the marker scan (see Details).
}
\item{type}{a character string namng the type of fixed effects to
  be fitted when individually scanning markers/intervals. Types are \code{"main"}
  for main effects only, \code{"int"} for interaction effects only and
  \code{"main+int"} for additive main effects and their interaction with
  \code{Trait} (see Details).
}      
\item{fix.lines}{a logical value. If \code{TRUE} then lines existing in
  the phenotype data that do not exist in \code{genObj} are fixed and placed in
  the \code{fixed} component of the \code{asreml()} models (see
  Details). It is recommended to set this to TRUE. Defaults to \code{TRUE}.
}
\item{chr}{a character string vector of chromosome names from \code{genObj}. All markers
  within these named chromosomes will be analysed.
 }
\item{gen.type}{a character string determining the type of genetic data to
  be used in the analysis. Possibilities are \code{"interval"} and
  \code{"markers"}. The default is \code{"interval"}.
}
\item{n.fa}{number of factors for the Factor Analytic whole genome random marker by
  \code{Trait} interaction model. Default = 0 (see Details).
}
\item{breakout}{a numerical integer equivalent to the iteration
 where the algorithm breaks out. The default is
  -1 which ensures the algorithm scans all markers within the specified
  chromosomes. (see Details)
}
\item{thresh}{the thresholding method for determining marker
  significance in Wald/LOD profiles (see Details).
}   
\item{qtl.window}{For each significant QTL, the genetic distance (in
  map units) to the left and right side excluded from further
  analysis. The default = 20 (see Details).
}
\item{effects.window}{For each significant interaction QTL, the genetic
  distance (in map units) to the left and right side to exclude main
  effect QTL. The default = 20 (see Details). 
}
\item{TypeI}{a numerical value determining the familywise alpha level of significance
    for detecting a QTL. The default is 0.05.
  }
\item{trace}{a automatic tracing facility. If \code{trace = TRUE} then
    all \code{asreml} output is piped to the screen during the analysis.
  If \code{trace = "file.txt"}, then output from all asreml models is
  piped to "\code{file.txt}". Both trace mechanisms will display a
  message if a QTL is detected.
}
\item{\dots}{
Any extra arguments to be passed to the \code{asreml}
  calls. within the function (see \code{?asreml} and \code{?asreml.options} for more
  information).
}
}
\details{
In the initial call to \code{gwasreml.asreml()}, the marker or interval
information is collated from \code{genObj}. If \code{gen.type =
  "interval"} then interval marker data is collated from the
\code{"interval.data"} components of the chromosomes in
\code{genObj}. If \code{gen.type = "markers"} then markers are
collated from the \code{"imputed.data"} components of the chromosomes in
\code{genObj}.

The \code{Trait} argument is used as a placeholder for potentially
different types of sectioning in the data. For example, \code{Trait} may
be the column factor describing different experiments or
environments. It may also be the column factor partitioning multiple
measurements of the same experimental unit. If \code{Trait} is set to
\code{NULL} then the function assumes there is no sectioning of the
data and proceeds with a single environment genome wide analysis. If
this is true and \code{n.fa} is greater than zero the it will
be set to \code{n.fa = 0} by the function. Additionally, if \code{type} is
not set to \code{"main"} then an error is returned for the user to
adjust the argument accordingly. 

If \code{Trait} is non-NULL then \code{type} can be flexibly set to the
users analysis requirements. For each of the different \code{type}
arguments the algorithm is specified below. 
\itemize{
  \item \code{type = "main"} then the algorithm individually fits
  marker/interval main effects terms only and initially assesses them for
  significance. Marker/intervals are then retained for the final model 
  according to the \code{qtl.window} rules outlined below.
  \item \code{"int"} then the algorithm follows Malosetti et al.
  (2013) and initially fits marker/interval by \code{Trait}
  interaction effects and assesses them for significance using the
  \code{qtl.window} rules outlined below. If significant
  then then the model is refitted with both main and interaction
  terms and assessed again. If the interaction is non-significant it is
  removed from the final model and if an interaction is significant then the
  corresponding main effect is removed from the final model to ensure
  marker/interval effects for all levels of \code{Trait} are available.
  \item \code{"main+int"} then the algorithm fits additive main effects
  and their interaction with \code{Trait} and assesses both terms
  simultaneously using the \code{qtl.window} and \code{effects.window}
  rules outlined below. Retained terms from this process are fitted in a
  final model. 
}

If \code{covariate} is non-NULL then the \code{covariate} will replace
\code{Trait} in the marker scan and \code{type="main+int"} is required
to be set. This ensures the algorithm defined above is followed with the
exception that, if an interaction term is significant the final model
will contain both main and interaction effects.   

For all \code{type} arguments, Wald/LOD statistic profiles can be
graphically displayed using \code{plotProfile()}. Currently there is
only one familywise QTL thresholding method, \code{thresh = 
"b-corr"}, a Bonferroni corrected threshold of Li & Ji (2005). Note if
\code{type = "main+int"} then the threshold for the interaction effects
may differ from the main effects if the effects have a differing number of 
degrees of freedom (i.e. the number of levels in \code{Trait} is
greater than two).

The function has in-built automatic detection of putative interaction and main
effect QTL peaks in each chromosome. This automation occurs
in several ways. Firstly, peaks are required to be above the adjusted familywise
Bonferroni corrected Wald/LOD values calculated using \code{thresh =
  "b-corr"}. Within a chromosome, Wald/LOD
statistic peaks are identified and ordered from highest to lowest.
A \code{qtl.window} is then placed to the left and right of the highest
peak and other peaks falling within this window are removed. If a peak
falls outside this window and remains part of the higher peak being assessed
then the peak outside the window is removed. For \code{type =
  "main+int"} a final assessment is
conducted by placing a secondary window around the
interaction peak and removing any significant main effect terms that fall
within \code{effects.window}. The final set of putative QTL are then
placed as fixed effects and a final model is refitted. (Note: these two
window arguments have low numerical defaults to match a \code{genObj}
with centimorgan map units. These would need to be dramatically
increased if map units, for example, were base-pairs from a physical map.)

Users can break out of the algorithm using the \code{breakout}
argument. If a numerical value greater than zero is given then the
scanning algorithm will breakout after the specified number of
iterations. This may be useful for quick investigation of a particular
markers within a chromosome. This feature will be enhanced in the future
versions of the function.

It is recommended to set \code{fix.lines = TRUE} to ensure additive and
non-additive genetic variances are estimated from lines in the
\code{merge.by} component of the phenotypic data that have genetic
marker data in \code{genObj}. Lines in the phenotype
\code{merge.by} factor not existing in \code{genObj} will be placed
as a fixed factor (called \code{Gomit}) in the \code{asreml} model.
Note, if there are others factors in the model that have some
potential confounding with \code{Gomit} then \code{asreml} will
indicate this with a simple message \code{'Terms with zero df listed
in attribute 'zerodf' of the wald table'} at the end of its iterative
maximisation. This confounding will have no effect on the outcome and
can be safely ignored. If \code{fix.lines = FALSE} is set then all
available lines in the \code{merge.by} component of the phenotypic data
will be used to estimate the non-additive genetic variance.
In this instance, users also need to be aware that \code{asreml} will
output a large number of warnings due to an inherent mismatch in the levels
of the lines contained in the phenotype data compared to the lines in
\code{genObj}.

The \code{n.fa} argument provides flexibility in specifying the Factor
Analytic model for the whole genome random marker by \code{Trait}
interaction effects. For \code{n.fa = 0} a \code{diag(Trait}) prefix term is
used and separate additive genetic variances are assumed for all levels specified
in \code{Trait}. For \code{n.fa = n} then a \code{fa(Trait, n)} prefix
term is used then a Factor Analytic model of order n is used as a
parsimonious approximation to the interaction effects. To ensure model
stability the Factor Analytic models are built sequentially. For
example, if \code{n.fa = 2} is set then the models for \code{n.fa = 0}
and \code{n.fa = 1} are fitted in sequence before the model with
\code{n.fa = 2}.

It is recommended that \code{trace = "file.txt"} be used to pipe the
sometimes invasive tracing of \code{asreml} licensing and fitting
numerics for each model to a file. Errors, warnings and messages will
still appear on screen during this process. Note some warnings that
appear may be passed through from an asreml call and are outputted upon
exit. These may be ignored as they are handled during the execution of
the function.
}
\value{An object of class \code{"gwasreml"} which also inherits the class
  \code{"asreml"} by default. The object returned is actually an \code{asreml}
  object (see \code{asreml.object}) with the addition of components from
  the QTL detection listed below.

\item{QTL}{A list of components from the significant QTL detected
    including a character vector of the significant QTL. The complete
    set of Wald statistics are also returned for graphically displaying
    using \code{plotProfile()}.
  }
 }
\references{
Li, J. and Ji, L. (2005) Adjusting multiple testing in multilocus
analyses using the eigenvalues of a correlation matrix.
\emph{Heredity}, \bold{95}, 221-227.
}
\author{
Julian Taylor
}
\seealso{
\code{\link{summary.gwasreml}}, \code{\link{update.gwasreml}}   
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
           Trait = "Expt", type = "main+int", n.fa = 0, na.action =
           na.method(x = "include"))

# summarise

summary(hwt.qtl, mapi)
}
}
\keyword{regression}
