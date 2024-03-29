% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_scanbase.R
\name{global_from_concat}
\alias{global_from_concat}
\title{((scan->nlin)->global) -> (nlin->global)}
\usage{
global_from_concat(scan_to_global, nlin_to_global, clobber = TRUE, dry = FALSE)
}
\arguments{
\item{scan_to_global}{path to an xfm file composed of a nonlinear scan->nlin and a nonlinear
nlin->global. The xfm file should have four component transforms linear->grid->linear->grid.
The final two are extracted and written out.}

\item{nlin_to_global}{A filename for the resultant transform, to be written to the directory
of \code{scan_to_global}}

\item{dry}{whether to do a dry run, meaning not run the command.}
}
\value{
the path to the output file invisibly.
}
\description{
Deconstruct a scan to global transform into an nlin to global transform
The output file is created in the same directory as the input file.
}
