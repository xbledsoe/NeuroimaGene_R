% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NIDP_names.R
\name{listNIDPs}
\alias{listNIDPs}
\title{List NIDPs}
\usage{
listNIDPs(modality = NA, atlas = NA, filename = NA, verbose = FALSE)
}
\arguments{
\item{modality}{Neuroimaging modality. Defaults to NA; see README for additional options}

\item{atlas}{Neuroimaging parcellation atlas for NIDP query. Defaults to NA; see README for additional options}

\item{filename}{optional filename for writing data to a table}

\item{verbose}{print runtime messages to R console. Default to FALSE}
}
\value{
a list of NIDP names satisfying the required criteria
}
\description{
Supplies a list of all NIDPs by name for any given modality or atlas
}
\examples{
dk_names <- listNIDPs(modality = 'T1', atlas = 'Desikan')


}
\keyword{reference}
