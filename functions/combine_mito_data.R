combine_mito_data <- function(shared_mutation_output, 
                              vaf_output, 
                              low_vaf_filter, 
                              high_vaf_filter) {
  
  # Join informational DFs
  cdf = left_join(shared_mutation_output, vaf_output)
  df = filter(cdf, 
                    cdf$VAF >= low_vaf_filter & cdf$VAF <= high_vaf_filter)
  
  mut_pos <- sapply(df$variant, function(x) {
    a = str_split_fixed(x,
                        pattern = "_",
                        n = 2)
    b = str_split_fixed(a[1],
                        pattern = ":",
                        n = 2)
    return(as.integer(b[2]))
  })
  
  ref = sapply(df$variant, function(x) {
    a = str_split_fixed(x,
                        pattern = "_",
                        n = 2)
    b = str_split_fixed(a[2],
                        pattern = "/",
                        n = 2)
    return((b[1]))
  })
  
  alt = sapply(df$variant, function(x) {
    a = str_split_fixed(x,
                        pattern = "_",
                        n = 2)
    b = str_split_fixed(a[2],
                        pattern = "/",
                        n = 2)
    return((b[2]))
  })
  
  df$mut_pos = mut_pos
  df$ref = ref
  df$alt = alt
  
  return(df)
}
