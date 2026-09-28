\name{print.metaweave}
\alias{print.assemblage}
\alias{print.distribution_collection}
\alias{print.ecological_network}
\alias{print.inference_model}
\alias{print.network_collection}
\title{Print a MetaWeave object}
\description{Prints a compact summary of an assemblage, distribution collection,
ecological network, inference model, or network collection. The object is
returned invisibly. These methods are normally called through \code{print()}.}

\usage{
\method{print}{assemblage}(x, ...)
\method{print}{distribution_collection}(x, ...)
\method{print}{ecological_network}(x, ...)
\method{print}{inference_model}(x, ...)
\method{print}{network_collection}(x, ...)
}
\arguments{
\item{x}{A metaweave object of the corresponding class.}
\item{...}{Additional arguments, currently unused.}
}
\value{The input object, returned invisibly after printing a compact summary.}
