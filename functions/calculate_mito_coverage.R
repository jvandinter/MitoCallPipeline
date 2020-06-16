# Calculate Mitochondrial coverage
# Jip van Dinter
# 25-02-2020

# supply list of of patients and the dir that contains the tabular coverage data
calculate_mt_coverage <- function(wd, mito_reference_df) {
  setwd(wd)
  f = list.files(pattern = ".tsv", all.files = T, recursive = T)
  cf = lapply(f, read.table, header = T)
  dat = bind_rows(cf, .id = "df")
  
  # Add sampatient_liste ID to relevant rows
  sl = paste(sapply(f,function(x) {
    a = str_split_fixed(x, pattern = "/",n = 2)
    return(a[1])
  }))
  
  dat$df = rep(sl, sapply(cf,nrow))
  dat$target = NULL
  colnames(dat) = c("sample",'chrom','position','coverage')
  
  # Add patient identifier to all its clones
  rdat <- left_join(dat, mito_reference_df, by = "sample")
  
  return(rdat)
}
