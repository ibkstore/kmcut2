% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractcolumns.r
\name{extractcolumns}
\alias{extractcolumns}
\title{Extract a sub-set of columns (such as a sub-set of samples) from a data table. All rows will be preserved.}
\usage{
extractcolumns(fnamein, fids, fnameout, wdir)
}
\arguments{
\item{fnamein}{character vector that specifies the name of tab-delimited text file with the input data table.}

\item{fids}{character vector that specifies the name of text file with column ids (such as sample ids).
The file must contain one column id per line, without any trailing spaces or any other additional symbols.}

\item{fnameout}{character vector that specifies the name of output file where the new data table will be saved.}

\item{wdir}{character vector that specifies the name of the working directory for the input/output files.}
}
\value{
no return value
}
\description{
Extract a sub-set of columns (such as a sub-set of samples) from a data table. All rows will be preserved.
}
\examples{

Basic usage:
extractcolumns(fnamein="table.txt", fids="columnids.txt", fnameout="newtable.txt", wdir="c:/test")

Example with built-in data files:

library(data.table)
library(kmcut)

# Load example gene expression data table for 2 genes
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
# Load a list that contains column (sample) ids
idlist = system.file("extdata", "columnids.txt", package="kmcut")
# Run the function
extractcolumns(fnamein=fdat, fids=idlist, fnameout="example_samples_subset.txt", wdir="c:/test")

This will create a tab-delimited text file in directory "c:/test":
"example_samples_subset.txt"
}
