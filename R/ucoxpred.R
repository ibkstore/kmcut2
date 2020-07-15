#' For each feature in the input file, this function fits a univariate Cox regression model on a training dataset and then uses the model
#' to predict risk scores for a test dataset.
#'
#' @param fname1 character vector that specifies the name of the training file with feature(s) for each sample. The file must be tab-delimited,
#' where features are in rows and samples are in columns. First column must contain feature names. Column names must contain sample ids.
#' @param sfname1 character vector that specifies the name of the file with right-censored survival time data for training dataset. The file must be tab-delimited,
#' where samples are in rows. First column must be named sample_id and contain sample ids that match those in fname1. The file must
#' contain columns called stime and scens, with survival time and censoring variable (0 or 1), respectively.
#' @param fname2 character vector that specifies the name of the test file with feature(s) for each sample. The file must be tab-delimited,
#' where features are in rows and samples are in columns. First column must contain feature names. Column names must contain sample ids.
#' @param sfname2 character vector that specifies the name of the file with right-censored survival time data for test dataset. The file must be tab-delimited,
#' where samples are in rows. First column must be named sample_id and contain sample ids that match those in fname2. The file must
#' contain columns called stime and scens, with survival time and censoring variable (0 or 1), respectively.
#' @param wdir character vector that specifies the name of the working directory for the input/output files.
#' Three output files are automatically created by adding:\cr'_cox_train_sum.txt' to 'fname1',
#' '_train_score.txt' to 'fname1';\cr'_test_score.txt' to 'fname2';
#' @param min_uval numeric value that specifies the minimal percentage of unique values per feature (default is 50).
#' Features that have less than 'min_uval' percent unique values are excluded from the analysis.
#' @param psort logical value whether to sort the output table by p-values in increasing order (default is FALSE).
#' @param verbose logical value whether to print progress (default is TRUE).
#' @return no return value
#'
#' @export
#' @examples
#'
#' Basic usage:
#' ucoxpred(fname1="train.txt", sfname1="survivaltrain.txt", fname2="test.txt", sfname2="survivaltest.txt", wdir="c:/test")
#'
#' Example with built-in data files:
#'
#' library(stats)
#' library(survival)
#' library(stringr)
#' library(data.table)
#' library(tools)
#' library(pracma)
#' library(kmcut)
#'
#' # Load example gene expression data and survival data for 2 genes and 295 samples
#' fdat = system.file("extdata", "example_genes_295.txt", package="kmcut")
#' sdat = system.file("extdata", "survival_data_295.txt", package="kmcut")
#'
#' # Use the example dataset as a training dataset and as a test dataset (perform resubstitution).
#' # NOTE: In a typical real world application training and test datasets will be different from one another.
#' ucoxpred(fname1=fdat, sfname1=sdat, fname2=fdat, sfname2=sdat, wdir="c:/test")
#'
#' This will create three output files in directory "c:/test":
#' 1) Tab-delimited text file with Cox summary for the training data: "example_genes_295_cox_train_sum.txt"
#' 2) Tab-delimited text file with risk scores for training data: "example_genes_295_train_score.txt"
#' 3) Tab-delimited text file with risk scores for test data: "example_genes_295_test_score.txt"

ucoxpred<-function(
  # The training file with feature(s) for each sample (samples are in columns, features are in rows)
  fname1,
  # The file with survival time data for training dataset
  sfname1,
  # The test file with feature(s) for each sample (samples are in columns, features are in rows)
  fname2,
  # The file with survival time data for test dataset
  sfname2,
  # Working directory for the input/output files
  wdir,
  # Min percentage of unique values in ]0, 100] for each feature
  min_uval = 50,
  # Option to sort the output table by p-values in increasing order (TRUE by default)
  psort = FALSE,
  # Print progress (TRUE by default)
  verbose = TRUE
)
# begin function
{

if(min_uval <= 0 || min_uval > 100)
{
  stop("min_uval must be in ]0, 100]")
}

setwd(wdir)

# Read training data
train_edat = read.delim(fname1, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
# Read test data
test_edat = read.delim(fname2, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
# Name of the output file with COXPH summary for the training data
cox_out_train = basename(file_path_sans_ext(fname1))
cox_out_train = sprintf("%s_cox_train_sum.txt",cox_out_train)
# Name of the output file with the risk scores for training data
risk_out_train = basename(file_path_sans_ext(fname1))
risk_out_train = sprintf("%s_train_score.txt",risk_out_train)
# Name of the output file with the risk scores for test data
risk_out_test = basename(file_path_sans_ext(fname2))
risk_out_test = sprintf("%s_test_score.txt",risk_out_test)

# The training survival time data
train_sdat_init = read.delim(sfname1, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
l = length(rownames(train_sdat_init))
l
# Must have "status" and "time" variables
train_sdat = data.frame(train_sdat_init$sample_id, train_sdat_init$scens, train_sdat_init$stime)
colnames(train_sdat) = c("sample_id","status","time")

# The test survival time data
test_sdat_init = read.delim(sfname2, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
l = length(rownames(test_sdat_init))
l
# Must have "status" and "time" variables
test_sdat = data.frame(test_sdat_init$sample_id, test_sdat_init$scens, test_sdat_init$stime)
colnames(test_sdat) = c("sample_id","status","time")

# Remove genes that have less than U% of unique values in the training dataset
n_genes1 = length(rownames(train_edat))
N = length(colnames(train_edat))
rv = logical(length = n_genes1)
for(i in 1:n_genes1)
{
  un = 100 * length(unique(unlist(train_edat[i,]))) / N
  if(un < min_uval)
  {
    rv[i] = FALSE
  }else
  {
    rv[i] = TRUE
  }
}
train_edat = train_edat[rv,]
length(rownames(train_edat))

# Remove genes that have less than U% of unique values in the test dataset
n_genes1 = length(rownames(test_edat))
N = length(colnames(test_edat))
rv = logical(length = n_genes1)
for(i in 1:n_genes1)
{
  un = 100 * length(unique(unlist(test_edat[i,]))) / N
  if(un < min_uval)
  {
    rv[i] = FALSE
  }else
  {
    rv[i] = TRUE
  }
}
test_edat = test_edat[rv,]
length(rownames(test_edat))

# Merge train expression data with survival data
train_edat = as.data.frame(t(train_edat))

length(colnames(train_edat))

# Convert row names into first column
rnames = rownames(train_edat)
train_edat = cbind(rnames,train_edat)
colnames(train_edat)[1] = "sample_id"

train_edat = merge(train_edat, train_sdat, by.x = 1, by.y = 1, all = FALSE)
rownames(train_edat) = unlist(train_edat["sample_id"])
train_edat["sample_id"] = NULL
colnames(train_edat) = make.names(colnames(train_edat), unique = TRUE)

n_train_samples = length(rownames(train_edat))

# Check for missing and non-numeric elements
row.has.na = apply( train_edat, 1, function(x){any(is.na(x) || is.nan(x) || is.infinite(x))} )
s = sum(row.has.na)
s

# Recalculate the number of features in the table
n_genes_train = length(colnames(train_edat)) - 2
n_genes_train

# Merge test expression data with survival data
test_edat = as.data.frame(t(test_edat))

length(colnames(test_edat))

# Convert row names into first column
rnames = rownames(test_edat)
test_edat = cbind(rnames,test_edat)
colnames(test_edat)[1] = "sample_id"

test_edat = merge(test_edat, test_sdat, by.x = 1, by.y = 1, all = FALSE)
rownames(test_edat) = unlist(test_edat["sample_id"])
test_edat["sample_id"] = NULL
colnames(test_edat) = make.names(colnames(test_edat), unique = TRUE)

n_test_samples = length(rownames(test_edat))

# Check for missing and non-numeric elements
row.has.na = apply( test_edat, 1, function(x){any(is.na(x) || is.nan(x) || is.infinite(x))} )
s = sum(row.has.na)
s

# Recalculate the number of features in the table
n_genes_test = length(colnames(test_edat)) - 2
n_genes_test

# Table with coxph summary for the training data
cox_train = matrix(data = NA, nrow = n_genes_train, ncol = 4)
colnames(cox_train) = c("CC","HR","P","FDR_P")
rownames(cox_train) = colnames(train_edat)[1:n_genes_train]

# Table with risk scores for training data
results_train = matrix(data = NA, nrow = n_genes_train, ncol = n_train_samples)
colnames(results_train) = rownames(train_edat)
rownames(results_train) = colnames(train_edat)[1:n_genes_train]

# Table with risk scores for test data
results_test = matrix(data = NA, nrow = n_genes_test, ncol = n_test_samples)
colnames(results_test) = rownames(test_edat)
rownames(results_test) = colnames(test_edat)[1:n_genes_test]

# Calculate the risk score for each feature in the training and test datasets
for(i in 1:n_genes_train)
{
  if(verbose == TRUE)
  {
    print(paste("Processing ",i," of ",n_genes_train))
  }

  res = coxph(as.formula(paste('Surv(time, status)~', colnames(train_edat)[i])), data=train_edat, model = TRUE)
  r = summary(res)
  cox_train[i,"CC"] = unlist(r$concordance["C"])
  cox_train[i,"HR"] = r$coefficients[2]
  df = as.data.frame(r["logtest"])
  cox_train[i,"P"] = df["pvalue",1]

  pred_resubst = predict(res, newdata = train_edat, type="risk")
  results_train[i,] = pred_resubst

  crow = which(colnames(train_edat) == colnames(test_edat)[i])
  if(length(crow) != 1) next

  pred_test = predict(res, newdata = test_edat, type="risk")
  results_test[crow,] = pred_test
}

cox_train[,"FDR_P"]  = p.adjust(cox_train[,"P"], method = "fdr")

# Convert to table and convert row names into a column
df = as.data.frame(cox_train)
setDT(df, keep.rownames=TRUE)
colnames(df)[1] = "tracking_id"
write.table(df, file = cox_out_train, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

df = as.data.frame(results_train)
setDT(df, keep.rownames=TRUE)
colnames(df)[1] = "tracking_id"
write.table(df, file = risk_out_train, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

df = as.data.frame(results_test)
setDT(df, keep.rownames=TRUE)
colnames(df)[1] = "tracking_id"
write.table(df, file = risk_out_test, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

}
# end function
