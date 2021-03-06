% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/target_msms_purify.R
\name{target_msms_purify}
\alias{target_msms_purify}
\title{target_msms_purify}
\usage{
target_msms_purify(file.name = NULL, abs.min.count = 100,
  rela.min.count = 1)
}
\arguments{
\item{file.name}{The csv file name of spectra.}

\item{abs.min.count}{The minimnum value of absolute intensity. Default: 100.}

\item{rela.min.count}{The minimnum value of relative intensity (1000 as benchmark). Default: 1.}
}
\description{
Purify msms spectra of targets
}
\author{
Zhiwei Zhou
\email{zhouzw@sioc.ac.cn}
}
