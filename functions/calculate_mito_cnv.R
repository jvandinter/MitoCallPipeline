# Function: outputs DF with mean, std error, median, MAD per sample per patient
# Jip van Dinter
# 25-02-2020

calculate_mito_cnv <- function(wd, mito_reference_df) {
  
  setwd(wd)
  
  nfiles = list.files(pattern = "WGSMetrics", recursive = T)
  mtfiles = list.files(pattern = ".tsv", recursive = T)
  
  # Read data
  ncov = lapply(nfiles, read.table, header = T, fill = T, stringsAsFactors = F)
  ncov = lapply(ncov, function(x) x[nrow(x)-(nrow(x)-1):nrow(x),])
  mtcov = lapply(mtfiles, read.table, header = T, stringsAsFactors = F)
  
  # Extract mean, median, SD, MAD
  #Nuclear
  nextr = lapply(ncov, function(x) {
      med = x[4]
      mad = x[5]
      mean = as.numeric(x[2])
      sdmean = x[3]
      df = data.frame("n_mean" = mean,"n_sd" = sdmean,'n_median' = med,'n_MAD' = mad)
      return(df)
    })
    
  # Add all samples to one DF
  ntable = bind_rows(nextr, .id = "df")
  nl = sapply(nfiles, function(x) {
    a = str_split_fixed(x,
                        pattern = "_dedup",
                        n = 2)
    return(a[1])
  })
  ntable$df = nl
  
  # MT
  mtextr = lapply(mtcov, function(x) {
    med = median(x[,4])
    mt_mean = mean(x[,4])
    sdmean = sd(x[,4])
    mad = mad(x[,4])
    df = data.frame("mt_mean" = mt_mean, "mt_sd" = sdmean, "mt_median" = med, "mt_MAD" = mad)
    return(df)
  })
  
  mttable = bind_rows(mtextr, .id = "df")
  sl = sapply(mtfiles, function(x) {
    a = str_split_fixed(x,
                        pattern = "/",
                        n = 2)
    return(a[1])
  })
  mttable$df = sl
  
  # Merge tables and compute CNV
  
  df <- left_join(mttable, ntable)
  colnames(df) = c("sample","mt_mean","mt_sd", "mt_median","mt_MAD","n_mean","n_sd",'n_median','n_MAD')
  df <- mutate(df, cnv_mean = mt_mean / n_mean) 
  df <- mutate(df, cnv_median = mt_median / n_median)
  
  # Add additional informational columns
  ## Add patient identifier to all its clones
  df <- left_join(df, mito_reference_df, by = "sample")

  return(df)
}
