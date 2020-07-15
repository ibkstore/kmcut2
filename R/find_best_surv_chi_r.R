# This function returns the best Chi-square from K-M minimization based on pracma findpeaks()
# and Spearman correlation between the observed and ideal linear plot Chi-sq vs cutoffs
# Returns a vector (best_cutoff, chi_sq, corr)
# Variables 'stime' and 'scens' are vectors with the survival times and censor variable, respectively
# Variable 'featdata' contains values of the feature for each sample (such as gene expression)
# 'stime', 'scens', and 'featdata' must be sorted in the same order
# 'ideal_y' - the ideal plot for the same number of points
# 'minfrac' - the value of the smallest allowed group size for K-M partitioning (in [0, 1])
# 'min_up_dn' - the number of up/down steps before/after peak for findpeaks

find_best_surv_chi_r <- function(mc.stime, mc.scens, mc.featdata, mc.ideal_y, mc.minfrac = 0.1, mc.min_up_dn = 1, mc.tolerance = 0.1)
{

  if(length(mc.stime) != length(mc.scens) || length(mc.stime) != length(mc.featdata) || length(mc.scens) != length(mc.featdata))
  {
    stop("Unequal vector lengths in function call")
  }
  if(mc.minfrac < 0)
  {
    stop("minfrac must be greater than or equal to 0")
  }

  mc.cutoffs = numeric()
  mc.chi_values = numeric()

  mc.min_n_s = mc.minfrac*length(mc.featdata)
  # print(mc.min_n_s)

  mc.uc = unique(mc.featdata)
  mc.n = length(mc.uc)
  mc.uc = sort(mc.uc, decreasing = FALSE)
  # print(mc.uc)

  if(mc.n <= 1)
  {
    warning("Only one unique cutoff point found")
    return(c(NA,NA,NA))
  }

  for(mc.j in 1:mc.n)
  {
    mc.cutoff = mc.uc[mc.j]
    mc.labels = as.numeric(mc.featdata > mc.cutoff) + 1
    mc.fact = factor(mc.labels)
    mc.low = sum(mc.fact == 1)
    mc.high = sum(mc.fact == 2)

    if(mc.low < mc.min_n_s || mc.high < mc.min_n_s) next

    mc.sur = survdiff(Surv(time=mc.stime, event=mc.scens, type = 'right') ~ mc.fact, rho = 0)

    mc.chi_values = c(mc.chi_values,mc.sur$chisq)
    mc.cutoffs = c(mc.cutoffs,mc.cutoff)
  }

  # Peak finder
  mc.peaks = findpeaks(mc.chi_values, nups = mc.min_up_dn, ndowns = mc.min_up_dn, npeaks=5, sortstr=TRUE,
                       minpeakdistance = 1, zero = "+")

  if(is.null(mc.peaks))
  {
    mc.best_cutoff = NA
    mc.max_chi = NA
  }else
  {

    # If there are peaks with nearly equal height, select the one that results in the most balanced (closest to the median) split
    # mc.tolerance = .Machine$double.eps ^ 0.5
    mc.best_peaks_i = mc.peaks[abs(mc.peaks[,1] - mc.peaks[1,1]) < mc.tolerance, 2]
    mc.mid_i = ceil(length(mc.chi_values)/2)
    if(length(mc.best_peaks_i) > 1)
    {
      mc.vd = abs(mc.best_peaks_i - mc.mid_i)
      mc.best_i = mc.best_peaks_i[mc.vd == min(mc.vd)][1]
    }else
    {
      mc.best_i = mc.best_peaks_i[1]
    }

    mc.max_chi = mc.chi_values[mc.best_i]
    mc.best_cutoff= mc.cutoffs[mc.best_i]
  }

  # Correlation between the observed and the ideal plots
  # The ideal plot is defined by two lines with intercepts and slopes alpha1,alpha2 and beta1,beta2
  mc.r = cor(mc.chi_values, mc.ideal_y, method = "spearman", use ="na.or.complete")

  return(c(mc.best_cutoff,mc.max_chi,mc.r))
}
