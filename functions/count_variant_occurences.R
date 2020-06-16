count_variant_occurences <- function(combined_mito_df, 
                                     mito_reference_df) {
  oc = as.data.frame(table(unlist(combined_mito_df$sample)))
  colnames(oc) <- c('sample','freq')
  
  # Add patient data
  oc <- left_join(mito_reference_df, oc, by = 'sample')
  oc[is.na(oc)] <- 0
  return(oc)
}
