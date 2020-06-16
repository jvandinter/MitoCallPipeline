check_ubiquitous_mut <- function(df, vcfl) {
  ndf <- mutate(df, ubiq = df$max_freq - df$freq)
  ndf$dual <- ifelse(ndf$max_freq == 2, '2S','MS')
  smdf <- subset(ndf, ubiq <= 1 & dual == 'MS')
  ndf$dual <- NULL
  smdf$pass <- apply(smdf, MARGIN = 1, function(x) {
    p <- as.character(x[4])
    v <- as.character(x[1])
    ps <- subset(vcf_all, patient == p & variant == v)
    a <- ifelse(length(ps) > as.numeric(x[10]), 0, 1)
    return(a)
    })
  fdf <- subset(smdf, pass == 1)
  output <- anti_join(ndf, fdf, by = 'variant')
  return(output)
}

check_mutation_frequency <- function(vcf_all_list) {
  return(vcf_all_list)
}
#subset table on patient


#look for freq compared to number of patients

