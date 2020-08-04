#' For each feature finds a cutoff that optimally stratifies samples into 2 groups, plots Kaplan-Meier survival curves
#' and observed vs. expected optimization plot. Then, uses permutation test to estimate the statistical significance of the optimal cutoff.
#'
#' @param fname character vector that specifies the name of the file with feature(s) for each sample. The file must be tab-delimited,
#' where features are in rows and samples are in columns. First column must contain feature names. Column names must contain sample ids.
#' @param sfname character vector that specifies the name of the file with right-censored survival time data. The file must be tab-delimited,
#' where samples are in rows. First column must contain sample ids that match those in 'fname'. The file must contain columns called 'stime' and 'scens',
#' with survival time and censoring variable (0 or 1), respectively.
#' @param wdir character vector that specifies the name of the working directory for the input/output files.
#' Output file names are automatically created by adding\cr"KMoptp_minf_.2f_iter_d"
#' and corresponding extension to 'fname'.
#' @param seed an integer value that specifies the seed for random number generator
#' @param min_fraction numeric value that specifies the minimal fraction of samples in the smaller group (default is 0.1).
#' @param min_up_down numeric value that specifies the minimal number of up/down points on either side of the peak for
#' pracma::findpeaks function (default is 1).
#' @param n_iter numeric value that specifies the number of iterations for the permutation test.
#' The default is n_iter=100 for fast calculations. Recommended is n_iter=10000 (slow, especially for a large number of samples/features).
#' @param peak_tolerance numeric value that specifies the maximal difference between in heigth between top peaks.
#' The peak within 'peak_tolerance' closest to the median value is selected.
#' @param min_uval numeric value that specifies the minimal percentage of unique values per feature (default is 50).
#' Features that have less than 'min_uval' percent unique values are excluded from the analysis.
#' @param psort logical value whether to sort the output table by p-values in increasing order (default is FALSE).
#' @param wlabels logical value whether to write a CSV file with low/high (below/above the cutoff) group sample labels (default is TRUE).
#' @param wpdf logical value whether to write a PDF file with survival curves and plots (default is TRUE).
#' @param verbose logical value whether to print progress (default is TRUE).
#' @return no return value
#'
#' @export
#' @examples
#'
#' Basic usage:
#' kmoptpermcut(fname="table.txt",sfname="survival.txt",wdir="c:/test",n_iter=1000,seed=1234)
#'
#' Example with data files included in the package:
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
#' # Search for optimal cutoffs and run the permutation tests
#' kmoptpermcut(fname=fdat, sfname=sdat, wdir="c:/test", n_iter=100, seed=1234)
#'
#' This will create three output files in directory "c:/test":
#' 1) PDF file with plots: "example_genes_295_KMoptp_minf_0.10_iter_100.pdf"
#' 2) Tab-delimited text file with the results: "example_genes_295_KMoptp_minf_0.10_iter_100.txt"
#' 3) CSV file with low/high sample labels: "example_genes_295_KMoptp_minf_0.10_iter_100_labels.csv"

kmoptpermcut<-function(
  # The full path to the file with feature(s) for each sample (samples are in columns, features are in rows)
  fname,
  # The full path to the file with survival time data
  sfname,
  # Working directory where the output files will be created
  wdir,
  # Random number generator seed (the default is NULL, which means set seed automatically using system time)
  seed = NULL,
  # If the fraction of samples in the smaller group is below this value, skip this partitioning
  min_fraction = 0.1,
  # Min number of up/down points on either side of the peak
  min_up_down = 1,
  # Number of random permutations
  n_iter = 100,
  # Peak tolerance
  peak_tolerance = 0.1,
  # Option to sort the output table by p-values in increasing order (TRUE by default)
  psort = FALSE,
  # Min percentage of unique values in ]0, 100] for each feature
  min_uval = 50,
  # Write a CSV file with low/high sample labels (TRUE by default)
  wlabels = TRUE,
  # Write a PDF file with survival curves (TRUE by default)
  wpdf = TRUE,
  # Print progress (TRUE by default)
  verbose = TRUE
)
# begin function
{

if(min_fraction < 0 || min_fraction >= 0.5)
{
  stop("min_fraction must be in [0, 0.5[")
}
if(min_up_down <= 0)
{
  stop("min_up_down must be greater than 0")
}
if(n_iter <= 0)
{
  stop("n_iter must be greater than 0")
}
if(peak_tolerance <= 0)
{
  stop("peak_tolerance must be greater than 0")
}
if(min_uval <= 0 || min_uval > 100)
{
    stop("min_uval must be in ]0, 100]")
}
if(is.null(seed)==TRUE)
{
  seed = ceil(as.numeric(Sys.time()))
}
set.seed(seed)
print(paste("Seed=",seed,sep=""))

setwd(wdir)

# Name of the output PDF file
pdf_file = basename(file_path_sans_ext(fname))
pdf_file = sprintf("%s_KMoptp_minf_%.2f_iter_%d.pdf",pdf_file,min_fraction,n_iter)

# Name of the output TXT file
txt_file = basename(file_path_sans_ext(fname))
txt_file = sprintf("%s_KMoptp_minf_%.2f_iter_%d.txt",txt_file,min_fraction,n_iter)

# Name of the output CSV file with low/high sample labels
csv_file = basename(file_path_sans_ext(fname))
csv_file = sprintf("%s_KMoptp_minf_%.2f_iter_%d_labels.csv",csv_file,min_fraction,n_iter)

# The survival time data
sdat = read.delim(sfname, header = TRUE, stringsAsFactors = FALSE)

# The input data table
edat = read.delim(fname, header = TRUE, row.names = 1)

# The number of genes in the table
n_genes1 = length(rownames(edat))
n_genes1

# The number of samples in the table
N = length(colnames(edat))

# Remove genes that have less than MIN_UVAL% of unique values
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
length(rownames(edat))

ids = intersect(sdat$sample_id,colnames(edat))

sdat = sdat[sdat$sample_id %in% ids,]

edat = edat[, colnames(edat) %in% ids]

sdat = sdat[order(sdat$sample_id), ]
edat = edat[,order(colnames(edat))]

# Check for missing and non-numeric elements
row.has.na = apply( edat, 1, function(x){any(is.na(x) || is.nan(x) || is.infinite(x))} )
s = sum(row.has.na)
if(s > 0)
{
  stop("The input data table has missing or non-numeric elements")
}

# Convert expression data table into a matrix
edat = as.matrix(edat)

# The number of features in the table
n_genes = length(rownames(edat))

# The number of samples in the table
n_samples = length(colnames(edat))

sample_labels = matrix(data = NA, nrow = n_samples, ncol = n_genes)
colnames(sample_labels) = rownames(edat)
rownames(sample_labels) = colnames(edat)

results = matrix(data = NA, nrow = n_genes, ncol = 7)

colnames(results) = c("CUTOFF","CHI_SQ","LOW_N","HIGH_N","R","P","FDR_P")
rownames(results) = rownames(edat)

# Go over each feature and find the partitioning with the largest Chi-square

if(wpdf == TRUE)
{
  pdf(pdf_file)
}

for(i in 1:n_genes)
{
  if(verbose == TRUE)
  {
    print(paste("Processing ",i," of ",n_genes))
  }

  cutoffs = numeric()
  chi_values = numeric()

  ucutoffs = as.numeric(unique(edat[i,]))
  ucutoffs = sort(ucutoffs, decreasing = FALSE)

  n_points = length(ucutoffs)

  if(n_points <= 1)
  {
    next
  }

  for(j in 1:n_points)
  {
    cutoff = ucutoffs[j]
    labels = as.numeric(as.vector(unlist(edat[i,])) > cutoff) + 1
    fact = factor(labels)
    low = sum(fact == 1)
    high = sum(fact == 2)

    if(low < min_fraction*n_samples || high < min_fraction*n_samples) next

    sur = survdiff(Surv(time = sdat$stime, event = sdat$scens, type = 'right') ~ fact, rho = 0)

    chi_values = c(chi_values,sur$chisq)
    cutoffs = c(cutoffs,cutoff)
  }

  # Peak finder
  x = cutoffs
  y = chi_values
  peaks = findpeaks(y, nups = min_up_down, ndowns = min_up_down, npeaks=5, sortstr=TRUE,
                    minpeakdistance = 1, zero = "+")

  if(is.null(peaks)) next

  # If there are peaks with nearly equal height, select the one that results in the most balanced (closest to the median) split
  tolerance = peak_tolerance
  best_peaks_i = peaks[abs(peaks[,1] - peaks[1,1]) < tolerance, 2]
  mid_i = ceil(length(y)/2)
  if(length(best_peaks_i) > 1)
  {
    vd = abs(best_peaks_i - mid_i)
    best_i = best_peaks_i[vd == min(vd)][1]
  }else
  {
    best_i = best_peaks_i[1]
  }

  # Correlation between the observed and the ideal plots
  # The ideal plot is defined by two lines with intercepts ans slopes alpha1,alpha2 and beta1,beta2
  max_chi = y[best_i]
  max_i = length(y)
  pks = c( ceil(max_i / 4), ceil(max_i / 2), ceil(3 * max_i / 4) )
  vd = abs(pks - best_i)
  pksi = pks[vd == min(vd)]
  if(length(pksi) == 1)
  {
    id_peak = pksi
  }else
  {
    id_peak = ceil(max_i / 2)
  }
  beta1 = max_chi / (id_peak - 1)
  alpha1 = -beta1
  beta2 = max_chi / (id_peak - max_i)
  alpha2 = -beta2 * max_i
  ideal_y = fit_ideal_linear_plot(alpha1,beta1,alpha2,beta2,id_peak,max_i)

  r = cor(y, ideal_y, method = "spearman", use ="na.or.complete")
  best_p = pchisq(max_chi, 1, lower.tail=FALSE)

  # Create the Kaplan-Meier curvs for two groups using the best cutoff
  best_cutoff = x[best_i]
  labels = as.numeric(edat[i,] > best_cutoff) + 1
  fact = factor(labels)
  best_low = sum(fact == 1)
  best_high = sum(fact == 2)
  if(wlabels == TRUE)
  {
    sample_labels[,i] = labels
  }

  results[i,"CUTOFF"] = best_cutoff
  results[i,"CHI_SQ"] = max_chi
  results[i,"LOW_N"] = best_low
  results[i,"HIGH_N"] = best_high
  results[i,"R"] = r

  # Skip the plotting if low or high group size is at the min limit
  # if(min(c(best_low,best_high)) < ceiling(min_fraction*n_samples))
  # {
  #   next
  # }

  # Estimate the joint random distribution of chisq and r
  rv = find_surv_chi_r_random(mcrf.stime = as.vector(unlist(sdat$stime)), mcrf.scens = as.vector(unlist(sdat$scens)),
                              mcrf.featdata = as.vector(unlist(edat[i,])), mcrf.ideal_y = ideal_y,
                              mcrf.niter = n_iter, mcrf.minfrac = min_fraction, mcrf.min_up_dn = min_up_down,
                              mcrf.tolerance = peak_tolerance, mcrf.verbose = verbose)
  rv = na.omit(rv)
  if(dim(rv)[1] > 0)
  {
     joint_p = sum( as.numeric(rv[,1] >= max_chi) + as.numeric(rv[,2] >= r) == 2) / n_iter
  }else
  {
    joint_p = 0
  }

  results[i,"P"] = joint_p

if(wpdf == TRUE)
{
  plot(x, y, type="b", col="navy", xlab="Cutoff", ylab="Chi-sq", ylim=c(0, max_chi+1))
  grid()
  title(main = paste(rownames(edat)[i],"; R=", round(r, digits = 3), sep = ""), sub = basename(fname),
        cex.sub = 0.6, cex.main = 0.6)
  points(x[best_i], y[best_i], pch=20, col="red", cex = 2)
  points(x, ideal_y, pch = 24, col = "green",type="b")

  title = rownames(edat)[i]
  title = paste(title,"; Cutoff=", sep = "")
  title = paste(title,sprintf("%G",best_cutoff), sep = "")
  title = paste(title,"; P=",sprintf("%G",joint_p), sep = "")

  sfit = survfit(Surv(time=sdat$stime, event = sdat$scens, type = 'right') ~ strata(fact))
  plot(sfit, lty=c(1, 2), xlab="Time", ylab="Survival Probability", mark.time=TRUE)
  title(main = title, sub = basename(fname), cex.sub = 0.6, cex.main = 0.6)
  ll = paste("Low, n=", best_low, sep = "")
  lh = paste("High, n=", best_high, sep = "")
  legend(x = "topright", y = NULL, c(ll, lh), lty=c(1, 2))
}

}

if(wpdf == TRUE)
{
  dev.off()
}

results[,"FDR_P"]  = p.adjust(results[,"P"], method = "fdr")

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
