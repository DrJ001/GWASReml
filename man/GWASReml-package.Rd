\name{GW-ASReml-package}
\alias{GW-ASReml-package}
\docType{package}
\title{Genome wide multi-section analysis using ASReml-R V4}
\description{
This package provides a flexible computational implementation of
multi-trait QTL analysis using new technological advances in ASReml-R
V4 functionality. Here, "multi-trait" is universal term encapsulating
three possible different types, multi-environment, multi-variate or multi-treatment.  
}
\details{
\tabular{ll}{
Package: \tab GW-ASReml\cr
Type: \tab Package\cr
Version: \tab 1.x\cr
Date: \tab 2021-01-22\cr
License: \tab GPL 2\cr
}
Welcome to version 1.x of \pkg{GW-ASReml}! The documentation given in this help file is
only brief and users should consult individual help files for the main functions of the
package. The package requires users have an installed and fully licensed
version of ASReml-R V4 (see
\url{https://www.vsni.co.uk/software/asreml-r} for more information).

It should be noted the package contains no functions for preparing your
genetic marker data for analysis. This can be done externally in many
ways but the genetic data object \code{genObj} to be passed into
\code{\link{gwasreml.asreml}} should be a \pkg{qtl} \code{"cross"}
object. Within each of the chromosomes of \code{genObj} there at least needs to be an
\code{imputed.data} or \code{interval.data} element containing complete
(no missing values) set of markers. It may be helpful to consult the
\code{cross2int()} function in the \pkg{wgaim} package.

Multi-trait QTL analysis is conducted using the function
\code{\link{gwasreml.asreml}} which, as its first argument, requires a
multi-section \code{asreml} base model. For a set of genetic
markers in \code{genObj}, the function builds on the supplied base model
in two fundamental ways. Given as a set of markers being scanned in a
chromosome, an additive relationship matrix is generated from all other
markers in the remaining chromosomes and incorporated into the variance
of a multi-trait term in the model. With this term in place, each
marker in the chromosome being tested are placed in the fixed model with
a main marker effect term and a trait by marker interaction term. Once
all markers have been scanned \code{gwasreml.asreml()} automatically
finds putative interaction and main effect QTL. A degree of flexibility
is provided for this automation by altering the two arguments
\code{qtl.window} and \code{effects.window} in
\code{gwasreml.asreml()}. To provide further flexibility, updating of the
returned object can occur using \code{\link{update.gwasreml}} where
the arguments of \code{gwasreml.asreml()} can be altered and the
analysis re-run.

A summary method is available for the returned \code{"gwasreml"}
object and provide users with a detailed report on the putative QTL
found by the algorithm including their size, the closest linked markers,
significance (including LOD score if desired). The effect profiles
returned in the \code{"gwasreml"} object and can be plotted using the
function \code{\link{plotProfile}}. This function allows
subsetting by chromosome and annotates putative QTL on the profiles where
appropriate. 
}
\author{
Author: Julian Taylor
Maintainer: Julian Taylor <julian.taylor@adelaide.edu.au>
}
\keyword{package}
\seealso{
\code{\link{qtl-package}} 
}
