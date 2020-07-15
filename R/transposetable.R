#' Transpose a data table: convert rows to columns and columns to rows.
#' Row names will become column names, and column names will become row names.
#' @param fnamein character vector that specifies the name of tab-delimited text file with the input data table.
#' @param fnameout character vector that specifies the name of output file where the transposed data table will be saved.
#' @param wdir character vector that specifies the name of the working directory for the input/output files.
#' @return no return value
#'
#' @export
#' @examples
#'
#' Basic usage:
#' transposetable(fnamein="table.txt", fnameout="transposedtable.txt", wdir="c:/test")
#'
#' Example with built-in data files:
#'
#' library(data.table)
#' library(kmcut)
#'
#' # Load example gene expression data table for 2 genes
#' fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
#'
#' transposetable(fnamein=fdat, fnameout="example_genes_295_transposed.txt", wdir="c:/test")
#'
#' This will create a tab-delimited text file with the transposed table in directory "c:/test":
#' "example_genes_295_transposed.txt"

transposetable<-function(
  # The name of tab-delimited text file with the input data table
  fnamein,
  # The name of output file where the transposed data table will be saved
  fnameout,
  # Working directory for the input/output files
  wdir
)
# begin function
{
setwd(wdir)

dat = read.delim(fnamein, header = TRUE, row.names = 1, stringsAsFactors = FALSE)

length(colnames(dat))
length(rownames(dat))

df = as.data.frame(t(dat))
colnames(df) = row.names(dat)

# Convert to table and convert row names into a column
setDT(df, keep.rownames=TRUE)
colnames(df)[1] = "tracking_id"

length(colnames(df))
length(rownames(df))

write.table(df, file = fnameout, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
# end function
