# Filter for unique mutations in list of VCFs
# Jip van Dinter
# 25-02-2020

genotype_df <- function(vcf_list) {
  
  gtl = lapply(vcf_list, 
               function(x) {
                 gt = geno(x)$GT %>% as.data.frame
                 gt <- rownames_to_column(gt, var = "variant")
                 return(gt)
               })
  
  dfms = lapply(gtl, function(x) {
    a = as.data.frame(x)
    colnames(a) = c("variant","genotype")
    return(a)
  })
  return(dfms)
}

split_patient_df <- function(geno_df, mito_reference_df) {
  
  dfm = do.call(rbind,geno_df)
  
  # Add rownames as new column
  dfm$sample = rownames(dfm)
  dfm$sample = sapply(dfm$sample, function(x) {
    a = str_split_fixed(x, pattern = "\\.", n = 2)
    return(a[1])
  })
  # Add patient data
  dfm <- left_join(dfm, mito_reference_df, by = 'sample')
  
  # split df in df per patient
  padf = split(dfm,dfm$patient)
  dfum = do.call(rbind,padf)
  return(dfum)
}

count_vars <- function(patient_df, mito_reference_df) {
  # Count number of occurences of variant in patient
  spm = lapply(patient_df, function(x) {
    xx = unique(x)
    oc = as.data.frame(table(unlist(xx$variant)))
    colnames(oc) <- c('variant','freq')
    # Connect number of occurences to variant
    a = left_join(xx, oc, by = "variant")
    return(a)
  })
  
  # Add them back together
  dfum = do.call(rbind,spm)
  oc = as.data.frame(table(unlist(mito_reference_df$patient)))
  colnames(oc) <- c('patient','max_freq')
  dfum = left_join(dfum, oc, by = 'patient')
  dfum = filter(dfum, dfum$genotype %in% c('0/1','0|1'))
  return(dfum)
}

create_variant_df <- function(vcf_list, mito_reference_df) {

  geno_df = genotype_df(vcf_list = vcf_list)
  
  padf = split_patient_df(geno_df = geno_df, 
                         mito_reference_df = mito_reference_df)
  
  #df = count_vars(patient_df = pdf, 
  #                mito_reference_df = mito_reference_df)
  
  return(padf)
}

filter_variants <- function(variant_df, vcf_all_df, sample_reference_df) {
  
  # Add frequencies to variant table
  
  vdf <- left_join(variant_df, vcf_all_df)
  oc <- as.data.frame(table(unlist(sample_reference_df$patient)))
  colnames(oc) <- c('patient','max_freq')
  vdf <- left_join(vdf, oc)
  
  vdf$single_sample <- ifelse(vdf$max_freq == 1, 'yes','no')
  # Filter for correct genotype
  vdf <- filter(vdf, genotype %in% c('0/1','0|1')) %>%
    # Filter for occurence in patient
    filter(single_sample == 'yes' | patient_freq < max_freq) %>%
    # Filter for global occurence
    filter(overall_freq <= patient_freq)
  vdf$single_sample <- NULL
  return(vdf)
}
