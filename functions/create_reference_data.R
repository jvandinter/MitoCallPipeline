create_reference_data <- function(wd, vcf_filter, patient_list, patient_age_df, trisomy_list, target_table) {
  
  setwd(wd)
  
  f <- list.files(pattern = vcf_filter, 
                  recursive = T)
  
  ff <- grep(pattern = "merged", 
             x = f, 
             val = T, 
             invert = T)
  
  sample_list = paste(sapply(ff, 
                              function(x) {
                                a = str_split_fixed(x, 
                                                    pattern = "/", 
                                                    n = 2)
                                return(a[1])
                              }))
  
  df <- data.frame(sample_list)
  
  colnames(df) = 'sample'
  df$patient = ''

  df <- df %>% mutate(patient = case_when(
    grepl('(^AC)([^0-9]{2})', sample) ~ patient_list[1],
    grepl(patient_list[2], sample) ~ patient_list[2],
    grepl(patient_list[3], sample) ~ patient_list[3],
    grepl(patient_list[4], sample) ~ patient_list[4],
    grepl(patient_list[5], sample) ~ patient_list[5],
    grepl(patient_list[6], sample) ~ patient_list[6],
    grepl(patient_list[7], sample) ~ patient_list[7],
    grepl(patient_list[8], sample) ~ patient_list[8],
    grepl(patient_list[9], sample) ~ patient_list[9],
    grepl(patient_list[10], sample) ~ patient_list[10],
    grepl(patient_list[11], sample) ~ patient_list[11],
    grepl(patient_list[12], sample) ~ patient_list[12],
    grepl(patient_list[13], sample) ~ patient_list[13],
    grepl(patient_list[14], sample) ~ patient_list[14]))
  
  df <- left_join(df, patient_age_df, by = "patient")
  # Label healthy and diseased states
  df$health = ifelse(grepl("AML",df$sample),"AML", 'healthy')
  df$ploid = ifelse(df$patient %in% trisomy_list, 'triploid',"diploid")
  df$state = ifelse(df$health == 'AML', "AML", 
                    ifelse(df$ploid == 'triploid', 'trisomy', 'healthy'))
  # label mesenchymal bulk samples
  df$bulk = ifelse(grepl('*.(BULK)',df$sample) & df$state == 'healthy','BULK',"clone")

  
  return(df)
}
  