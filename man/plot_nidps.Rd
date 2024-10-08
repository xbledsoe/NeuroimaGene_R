% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_nidps.R
\name{plot_nidps}
\alias{plot_nidps}
\title{NIDP summary plot of NeuroimaGene object}
\usage{
plot_nidps(ng_obj, maxNidps = 30, title = NA, verbose = FALSE)
}
\arguments{
\item{ng_obj}{NeuroimaGene Object}

\item{maxNidps}{maximum number of NIDPs to visualize. default=30}

\item{title}{optional title tag for the plot}

\item{verbose}{print runtime messages to R console. Default to FALSE}
}
\value{
a ggplot class object detailing mean effect size magnitude per NIDP, colored by brain region
}
\description{
Generate overview plot of the neuroimagene object according to nidps
}
\examples{
gene_list <- c('TRIM35', 'PROSER3', 'EXOSC6', 'PICK1', 'UPK1A', 'ESPNL', 'ZIC4')
ng <- neuroimaGene(gene_list, atlas = NA, mtc = 'BH', vignette = TRUE)
plot_nidps(ng)

}
\keyword{neuroimaging}
