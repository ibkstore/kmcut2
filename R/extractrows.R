#' Extract a sub-set of rows (such as a group of gene ids) from a data table. All columns will be preserved.
#' @param fnamein character vector that specifies the name of tab-delimited text file with the input data table.
#' @param fids character vector that specifies the name of text file with row ids (such as gene ids).
#' The file must contain one row id per line, without any trailing spaces or any other additional symbols.
#' @param fnameout character vector that specifies the name of output file where the new data table will be saved.
#' @param wdir character vector that specifies the name of the working directory for the input/output files.
#' @return no return value
#'
#' @export
#' @examples
#'
#' Basic usage:
#' extractrows(fnamein="table.txt", fids="geneids.txt", fnameout="newtable.txt", wdir="c:/test")
#'
#' Example with built-in data files:
#'
#' library(data.table)
#' library(kmcut)
#'
#' # Load example gene expression data table for 2 genes
#' fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
#' # Load a list that contains one gene id (MYCN)
#' idlist = system.file("extdata", "rowids.txt", package="kmcut")
#' # Run the function
#' extractrows(fnamein=fdat, fids=idlist, fnameout="example_genes_subset.txt", wdir="c:/test")
#'
#' This will create a tab-delimited text file with one row "MYCN" in directory "c:/test":
#' "example_genes_subset.txt"

extractrows<-function(
  # The name of tab-delimited text file with the input data table
  fnamein,
  # The name of text file with row ids (such as gene ids)
  fids,
  # The name of output file where the new data table will be saved
  fnameout,
  # Working directory for the input/output files
  wdir
)
# begin function
{
setwd(wdir)

dat = read.delim(fnamein, header = TRUE, stringsAsFactors = FALSE)
dat[,1] = make.names(dat[,1], unique = TRUE)
length(rownames(dat))

# The gene list
ids = unique(as.vector(unlist(read.delim(fids, header = FALSE))))

length(ids)

dat = dat[dat[,1] %in% ids,]

if(length(rownames(dat)) == 0)
{
  warning("ERROR: no matching row ids found in the data table...")
}

write.table(dat, file = fnameout, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
# end function

