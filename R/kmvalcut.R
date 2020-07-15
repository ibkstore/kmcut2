#' Creates Kaplan-Meier survival curves for each feature from a validation data set by using a file with previously determined
#' stratification thresholds (one threshold per feature), and calculates the log-rank test p-value.
#'
#' @param input1 character vector that specifies the name of tab-delimited file with the table that contains one or more feature
#' and a stratification threshold for each feature (this table is produced by kmoptscut, kmoptpermcut, kmqcut or kmucut).
#' The file with previously determined stratification thresholds must have first two columns named as 'tracking_id' and 'CUTOFF'.
#' The 'tracking_id' column contains feature names, the 'CUTOFF' column contains stratification threshold for each feature.
#' @param input2 character vector that specifies the name of the file with feature(s) for each sample. The file must be tab-delimited,
#' where features are in rows and samples are in columns. First column must contain feature names. Column names must contain sample ids.
#' Feature names must exactly match the ones in 'input1' file.
#' @param sfname character vector that specifies the name of the file with right-censored survival time data. The file must be tab-delimited,
#' where samples are in rows. First column must be named 'sample_id' and contain sample ids that match those in 'fname'. The file must contain columns called 'stime' and 'scens',
#' with survival time and censoring variable (0 or 1), respectively.
#' @param wdir character vector that specifies the name of the working directory for the input/output files.
#' Output file names are automatically created by adding\cr'_KM_val' and corresponding extension to 'input2'.
#' @param min_uval numeric value that specifies the minimal percentage of unique values per feature (default is 50).
#' Features that have less than 'min_uval' percent unique values are excluded from the analysis.
#' @param psort logical value whether to sort the output table by p-values in increasing order (default is FALSE).
#' @param wlabels logical value whether to write a CSV file with low/high (below/above the cutoff) group sample labels (default is TRUE).
#' @param wpdf logical value whether to write a PDF file with survival curves and plots (default is TRUE).
#' @return no return value
#'
#' @export
#' @examples
#'
#' Basic usage:
#' kmvalcut(input1="cutoffs.txt",input2="validation.txt",sfname="survival.txt",wdir="c:/test")
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
#' # "example_genes_295_KM_quant_50.txt" is a file with cutoffs created by 'kmqcut'
#' # and must exist in directory "c:/test"
#' kmvalcut(input1="example_genes_295_KM_quant_50.txt",input2=fdat,sfname=sdat,wdir="c:/test")
#'
#' This will create three output files in directory "c:/test":
#' 1) PDF file with plots: "example_genes_295_KM_val.pdf"
#' 2) Tab-delimited text file with the results: "example_genes_295_KM_val.txt"
#' 3) CSV file with low/high sample labels: "example_genes_295_KM_val_labels.csv"

kmvalcut<-function(
   # File with the table that contains feature(s) and cutoff(s) (table is produced by kmoptpermcut, kmoptcut, kmqcut or kmucut)
  input1,
  # The file with feature(s) for each sample (samples are in columns, features are in rows)
  input2,
  # The file with survival time data
  sfname,
  # Working directory with the input/output files
  wdir,
  # Min percentage of unique values in ]0, 100%] for each feature
  min_uval = 50,
  # Option to sort the output table by p-values in increasing order (TRUE by default)
  psort = FALSE,
  # Write a CSV file with low/high sample labels (TRUE by default)
  wlabels = TRUE,
  # Write a PDF file with survival curves (TRUE by default)
  wpdf = TRUE
)
# begin function
{

if(min_uval <= 0 || min_uval > 100)
{
  stop("min_uval must be in ]0, 100]")
}

setwd(wdir)

# Name of the output PDF file
pdf_file = basename(file_path_sans_ext(input2))
pdf_file = sprintf("%s_KM_val.pdf",pdf_file)

# Name of the output TXT file
txt_file = basename(file_path_sans_ext(input2))
txt_file = sprintf("%s_KM_val.txt",txt_file)

# Name of the output CSV file with low/high sample labels
csv_file = basename(file_path_sans_ext(input2))
csv_file = sprintf("%s_KM_val_labels.csv",csv_file)

sdat = read.delim(sfname, header = TRUE, stringsAsFactors = FALSE)

# The input table with cutoffs
cdat = read.delim(input1, header = TRUE)

# The input data table with features to test
edat = read.delim(input2, header = TRUE, row.names = 1)

ids = intersect(sdat$sample_id,colnames(edat))

sdat = sdat[sdat$sample_id %in% ids,]

edat = edat[, colnames(edat) %in% ids]

sdat = sdat[order(sdat$sample_id), ]
edat = edat[,order(colnames(edat))]

########################### Remove features that have less than min_uval% of unique values
N = length(colnames(edat))
n_genes1 = length(rownames(edat))
rv = logical(length = n_genes1)

for(i in 1:n_genes1)
{
  un = 100 * length(unique(unlist(edat[i,]))) / N
  if(un < min_uval)
  {
    rv[i] = FALSE
  }else
  {
    rv[i] = TRUE
  }
}

edat = edat[rv,]
###########################

# Check for missing and non-numeric elements
row.has.na = apply( edat, 1, function(x){any(is.na(x) || is.nan(x) || is.infinite(x))} )
s = sum(row.has.na)
if(s > 0)
{
  stop("The input data table has missing or non-numeric elements")
}

# Convert expression data table into a matrix
edat = as.matrix(edat)

# The number of genes in the table with cutoffs
n_genes = length(rownames(cdat))

# The number of samples in the table
n_samples = length(colnames(edat))

results = matrix(data = NA, nrow = n_genes, ncol = 6)

sample_labels = matrix(data = NA, nrow = n_samples, ncol = n_genes)
colnames(sample_labels) = cdat[,1]
rownames(sample_labels) = colnames(edat)

colnames(results) = c("CUTOFF","CHI_SQ","LOW_N","HIGH_N","P","FDR_P")
rownames(results) = cdat[,1]

if(wpdf == TRUE)
{
  pdf(pdf_file)
}

for(i in 1:n_genes)
{
  if(is.na(cdat[i,2])) next

  cur_low = NA
  cur_high = NA
  cur_p = NA
  cur_chi = NA
  cur_cutoff = NA

  crow = which(rownames(edat) == cdat[i,1])
  if(length(crow) != 1) next

  cur_cutoff = cdat[i,2]

    # Create the Kaplan-Meier curvs for two groups using the cutoff
    labels = as.numeric(edat[crow,] > cur_cutoff) + 1
    sample_labels[,i] = labels
    fact = factor(labels)
    cur_low = sum(fact == 1)
    cur_high = sum(fact == 2)
    if(cur_low > 0 && cur_high > 0)
    {
      sur = survdiff(Surv(time = sdat$stime, event = sdat$scens, type = 'right') ~ fact, rho = 0)
      cur_p = pchisq(sur$chisq, 1, lower.tail=FALSE)
      cur_chi = sur$chisq

      if(!is.na(cur_p) && wpdf == TRUE)
      {
        title = rownames(edat)[crow]
        title = paste(title,"; Cutoff = ", sep = "")
        title = paste(title,sprintf("%G",cur_cutoff), sep = "")
        title = paste(title,"; P-value = ", sep = "")
        title = paste(title,sprintf("%G",cur_p), sep = "")

        sfit = survfit(Surv(time=sdat$stime, sdat$scens==1) ~ strata(fact))
        plot(sfit, lty=c(1, 2), xlab="Time", ylab="Survival Probability", mark.time=TRUE)
        title(main = title, sub = basename(input2), cex.sub = 0.6, cex.main = 0.6)
        ll = paste("Low, n=", cur_low, sep = "")
        lh = paste("High, n=", cur_high, sep = "")
        legend(x = "topright", y = NULL, c(ll, lh), lty=c(1, 2))
      }
    }

  results[i,"CUTOFF"] = cur_cutoff
  results[i,"CHI_SQ"] = cur_chi
  results[i,"P"] = cur_p
  results[i,"LOW_N"] = cur_low
  results[i,"HIGH_N"] = cur_high
}

results[,"FDR_P"]  = p.adjust(results[,"P"], method = "fdr")

if(wpdf == TRUE)
{
  dev.off()
}

if(length(rownames(results)) > 1 && psort == TRUE)
{
  results = results[order(results[,"P"]),]
}

# Convert to table and convert row names into a column
df = as.data.frame(results)
setDT(df, keep.rownames=TRUE)
colnames(df)[1] = "tracking_id"
write.table(df, file = txt_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

if(wlabels == TRUE)
{
  df = as.data.frame(sample_labels)
  setDT(df, keep.rownames=TRUE)
  colnames(df)[1] = "sample_id"
  write.table(df, file = csv_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
}

}
# end function
