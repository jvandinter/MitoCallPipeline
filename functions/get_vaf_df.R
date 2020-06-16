# Function: Get VAF in df
# Jip van Dinter & Freek Manders
# 25-02-2020

get_vaf_df <- function(vcf_list) {
  
  # Adds point mutation context
  vaf = lapply(vcf_list, function(x) { 
    ad = geno(x)$AD[,1]
    vaf = sapply(ad, 
                 function(x) x[[2]] / sum(x))
    vaf[is.na(vaf)] = 0 #For sites with 0 reads
    vaf = sapply(vaf, unlist)
    return(vaf)
    })
  
  # Create df of all mutations
  dfvs = lapply(vaf, function(x) {
    a = as.data.frame(x)
    return(a)
  })
  
  dfv = do.call(rbind,dfvs)
  
  # Add rownames as new column
  dfv$sample = rownames(dfv)
  
  #split column in samples and chromosome locations
  dfv = separate(dfv, sample, sep = "\\.", into = c("sample","variant"), extra = "merge")
  colnames(dfv) = c("VAF","sample","variant")
  return(dfv)
}
