# Returns a matrix of best random Chi-square values and correlations between the observed and ideal curves
# First column - chi-sq, Second column - corr, each row - one iteration of random sampling
# This function shuffles the feature values
# Variables 'stime' and 'scens' are vectors with the survival times and censor variable, respectively
# Variable 'featdata' contains values of the feature for each sample (such as gene expression)
# 'stime', 'scens', and 'featdata' must be sorted in the same order
# 'ideal_y' - the ideal plot for the same number of points
# 'niter' - number of iterations for random sampling
# 'minfrac' - the value of the smallest allowed group size for K-M partitioning (in [0, 1])
# 'min_up_dn' - the number of up/down steps before/after peak for findpeaks

find_surv_chi_r_random <- function(mcrf.stime, mcrf.scens, mcrf.featdata, mcrf.ideal_y, mcrf.niter = 100,
                                   mcrf.minfrac = 0.1, mcrf.min_up_dn = 1, mcrf.tolerance = 0.1, mcrf.verbose = TRUE)
{
  
  if(length(mcrf.stime) != length(mcrf.scens) || length(mcrf.stime) != length(mcrf.featdata) || length(mcrf.scens) != length(mcrf.featdata))
  {
    stop("Unequal vector lengths in function call")
  }
  if(mcrf.niter <= 0)
  {
    stop("The number of iterations must be greater than 0")
  }
  if(mcrf.minfrac < 0)
  {
    stop("minfrac must be greater than or equal to 0")
  }
  
  mcrf.results = matrix(data = NA, nrow = mcrf.niter, ncol = 2)
  colnames(mcrf.results) = c("chisq","corr")
  
  # Shufffle the features
  mcrf.rd = mcrf.featdata
  mcrf.rt = mcrf.stime
  mcrf.rc = mcrf.scens
  
  for(mcrf.i in 1:mcrf.niter)
  {
    if(mcrf.verbose == TRUE && mcrf.i %% 100 == 0) 
    {
      print(mcrf.i)
    }
    mcrf.rd = sample(mcrf.rd, replace = FALSE, prob = NULL)
    mcrf.rt = sample(mcrf.rt, replace = FALSE, prob = NULL)
    mcrf.rc = sample(mcrf.rc, replace = FALSE, prob = NULL)
    mcrf.r = find_best_surv_chi_r(mc.stime = mcrf.rt, mc.scens = mcrf.rc, mc.featdata = mcrf.rd, mc.ideal_y = mcrf.ideal_y,
                                  mc.minfrac = mcrf.minfrac, mc.min_up_dn = mcrf.min_up_dn, mc.tolerance = mcrf.tolerance)
    mcrf.results[mcrf.i,1] = mcrf.r[2] # Chi-sq
    mcrf.results[mcrf.i,2] = mcrf.r[3] # Correlation
  }
  
  return(mcrf.results)
}
