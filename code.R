#<- Direct installation from github.com ->
# Install package "devtools"
install.packages("devtools")
# Install package "kmcut"
library(devtools)
install_github("ibkstore/kmcut2")

# <- Download kmcut2.tar.gz and install from a local folder ->
install.packages("C:/test/kmcut2.tar.gz", repos = NULL, type = "source")
install.packages(c("survival", "stringr", "data.table","pracma"))

# View the example survival data file
library(kmcut)
sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
writeLines(readLines(con=sdat))

# View the example expression data file
library(kmcut)
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
writeLines(readLines(con=fdat))

# Example of how to use "kmoptpermcut" with the data files included in the package
library(stats)
library(survival)
library(stringr)
library(data.table)
library(tools)
library(pracma)
library(kmcut)
# Load example gene expression data and survival data for 2 genes and 295 samples
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
# Run the permutation test for 100 iterations for each gene, with fixed seed=1234
kmoptpermcut(fname=fdat, sfname=sdat, wdir="c:/test", n_iter=100, seed=1234)

# Example of how to use "kmoptscut" with the data files included in the package
library(stats)
library(survival)
library(stringr)
library(data.table)
library(tools)
library(pracma)
library(kmcut)
# Load example gene expression data and survival data for 2 genes and 295 samples
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
# Run the function
kmoptscut(fname=fdat, sfname=sdat, wdir="c:/test")

# Example of how to use "kmqcut" with the data files included in the package
library(stats)
library(survival)
library(stringr)
library(data.table)
library(tools)
library(pracma)
library(kmcut)
# Load example gene expression data and survival data for 2 genes and 295 samples
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
# Use the 50th quantile (the median) to stratify the samples
kmqcut(fname=fdat, sfname=sdat, wdir="c:/test", quant=50)

# Example of how to use "kmucut" with the data files included in the package
library(stats)
library(survival)
library(stringr)
library(data.table)
library(tools)
library(pracma)
library(kmcut)
# Load example gene expression data and survival data for 2 genes and 295 samples
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
# Use the cutoff=5 to stratify the samples and remove features that have less than 90% unique # values (this eliminates MYH2 gene)
kmucut(fname=fdat, sfname=sdat, wdir="c:/test", cutoff=5, min_uval=90)

# Example of how to use "kmvalcut" with the data files included in the package 
library(stats)
library(survival)
library(stringr)
library(data.table)
library(tools)
library(pracma)
library(kmcut)
# Load example gene expression data and survival data for 2 genes and 295 samples
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
# "example_genes_295_KM_quant_50.txt" is a file with cutoffs created by 'kmqcut' and must exist in directory "c:/test"
kmvalcut(input1="example_genes_295_KM_quant_50.txt", input2=fdat, sfname=sdat, wdir="c:/test")

# Example of how to use "ucoxbatch" with the data files included in the package
library(stats)
library(survival)
library(stringr)
library(data.table)
library(tools)
library(pracma)
library(kmcut)
# Load example gene expression data and survival data for 2 genes and 295 samples
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
# Run the function
ucoxbatch(fname = fdat, sfname = sdat, wdir = "c:/test")

# Example of how to use "ucoxpred" with the data files included in the package
library(stats)
library(survival)
library(stringr)
library(data.table)
library(tools)
library(pracma)
library(kmcut)
# Load example gene expression data and survival data for 2 genes and 295 samples
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
# Use the example dataset as both a training dataset and as a test dataset (that is, perform a 
# resubstitution). Note that in a real world application training and test datasets will be different from one another.
ucoxpred(fname1=fdat, sfname1=sdat, fname2=fdat, sfname2=sdat, wdir="c:/test")

# Example of how to use the function 'extractrows' with the data files included in the package
library(data.table)
library(kmcut)
# Load example gene expression data table for 2 genes (2 rows)
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
# Load a list that contains one gene id (MYCN)
idlist = system.file("extdata", "rowids.txt", package="kmcut")
# Run the function
extractrows(fnamein=fdat, fids=idlist, fnameout="example_genes_subset.txt", wdir="c:/test")

# Example of how to use the function "extractcolumns" with the data files included in the package
library(data.table)
library(kmcut)
# Load example gene expression data table for 2 genes
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
# Load a list that contains column (sample) ids
idlist = system.file("extdata", "columnids.txt", package="kmcut")
# Run the function
extractcolumns(fnamein=fdat, fids=idlist, fnameout="example_samples_subset.txt", wdir="c:/test")

# Example of how to use the function "transposetable" with the data files included in the package
library(data.table)
library(kmcut)
# Load example gene expression data table for 2 genes. In this file, genes are rows and samples are columns.
fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
# Run the function
transposetable(fnamein=fdat, fnameout="example_genes_295_transposed.txt", wdir="c:/test")
