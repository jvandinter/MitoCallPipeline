# Load VCFs and filter variants that PASSED previous filtering

filter_passed_variants <- function(ref_genome, wd, vcf_filter) {
  
  library(ref_genome, character.only = TRUE)
  setwd(wd)
  
  # Load VCFs to analyse
  # Merged VCFs of same sample from Mito pipeline WDL
  
  f <- list.files(pattern = vcf_filter, 
                  recursive = T)
  
  ff <- grep(pattern = "merged", 
             x = f, 
             val = T, 
             invert = T)
  
  # Filter only for PASS variants
  
  vcfl_pass <- lapply(ff, 
                      function(x) {
                        v = readVcf(x)
                        passed = grep("PASS", 
                                      filt(v), 
                                      fixed = T)
                        return(v[passed,])
                      })
  
  # Add sample names
  
  names(vcfl_pass) <- paste(sapply(ff, 
                                   function(x) {
                                     a = str_split_fixed(x, 
                                                         pattern = "/", 
                                                         n = 2)
                                     return(a[1])
                                   }))
  return(vcfl_pass)
}

# Return all called variants in a DF

grab_all_variants<- function(ref_genome, wd, vcf_filter, mito_reference_df) {
  
  library(ref_genome, character.only = TRUE)
  setwd(wd)
  
  # Load VCFs to analyse

  f = list.files(pattern = vcf_filter, 
                  recursive = T)
  
  ff = grep(pattern = "merged", 
             x = f, 
             val = T, 
             invert = T)
  
  vcfs = lapply(ff, function(x) {
    v = readVcf(x)
    return(v)
    })
  
  names(vcfs) = paste(sapply(ff, 
                                   function(x) {
                                     a = str_split_fixed(x, 
                                                         pattern = "/", 
                                                         n = 2)
                                     return(a[1])
                                   }))
  
  gtl = lapply(vcfs, 
               function(x) {
                 flt = fixed(x)$FILTER
                 gt = geno(x)$GT %>% 
                   as.data.frame %>% 
                   rownames_to_column(var = "variant")
                 gt["filter"] <- flt
                 return(gt)
               })

  dfms = lapply(gtl, function(x) {
    a = as.data.frame(x)
    colnames(a) = c("variant","genotype",'filter')
    return(a)
  })

  dfm = do.call(rbind,dfms)
  
  # Add rownames as new column
  dfm$sample = rownames(dfm)
  dfm$sample = sapply(dfm$sample, function(x) {
    a = str_split_fixed(x, pattern = "\\.", n = 2)
    return(a[1])
  })
  # Add patient data
  dfm <- left_join(dfm, mito_reference_df, by = 'sample')
  dfm <- filter(dfm, complete.cases(dfm))
  dfm$varsamp <- paste(dfm$variant, dfm$patient, sep = '_')
  
  # Add frequency of variant per patient
  oc = as.data.frame(table(unlist(dfm$varsamp)))
  colnames(oc) <- c('varsamp','patient_freq')
  dfm <- left_join(dfm, oc, by = 'varsamp')
  dfm$varsamp <- NULL
  
  # Add frequency of variant overall
  oc = as.data.frame(table(unlist(dfm$variant)))
  colnames(oc) <- c('variant','overall_freq')
  dfm <- left_join(dfm, oc, by = 'variant')

  return(dfm)
}

# Load variant annotations in a tibble

get_exon_table = function(vcf) {    #Get gene annotation
  ann = info(vcf)$ANN %>%
    as.list()
  ann[elementNROWS(ann) == 0] = ""
  ann = purrr::map_chr(ann, str_c, collapse = ";")    #Filter for coding muts
  coding_f = str_detect(ann, "synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant")
  if (!sum(coding_f)){
    return(NULL)
  }
  ann_exon = ann[coding_f]    #Get info for gene table
  ann_l = str_split(ann_exon, pattern = "\\|")
  effect_type = purrr::map_chr(ann_l, 2) #Missense, nonnsense ect.
  effect = purrr::map_chr(ann_l, 3)
  gene = purrr::map_chr(ann_l, 4)
  vcf = vcf[coding_f]
  gt = geno(vcf)$GT
  tb = gt %>%
    as.data.frame() %>%
    rownames_to_column( var = 'variant') %>%
    dplyr::mutate(sample = colnames(gt[1]),
                  gene = gene,
                  effect = effect,
                  effect_type = effect_type)
  tb$sample = colnames(tb)[2]
  colnames(tb)[2] = 'genotype'
    return(tb)
}

grab_variant_annotations <- function(vcf_filter) {
  
  f <- list.files(pattern = vcf_filter, 
                  recursive = T)
  
  ff = grep(pattern = "merged", 
            x = f, 
            val = T, 
            invert = T)
  
  vcfs = lapply(ff, function(x) {
    v = readVcf(x)
    return(v)
  })
  
  names(vcfs) = paste(sapply(ff, 
                             function(x) {
                               a = str_split_fixed(x, 
                                                   pattern = "/", 
                                                   n = 2)
                               return(a[1])
                             }))
  
  annotated_variants = lapply(vcfs, get_exon_table)
  
  dfav = do.call(rbind,annotated_variants)
  
  return(dfav)
}

single_type_occurence <- function(grange, ref_genome, subset_variant_df) {
  df <- data_frame()
  vcf <- grange
  types <- mut_type(vcf)
  CpG = c("ACG", "CCG", "TCG", "GCG")
  column_names = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G",
                   "C>T at CpG", "C>T other")
  CT_context = 0
  CT_at_CpG = 0
  CT_at_other = 0
  CT_muts = which(types == "C>T")
  
  if (length(CT_muts) > 0) {
    CT_context = type_context(vcf[CT_muts], ref_genome)[[2]]
    CT_at_CpG = sum(!(is.na(BiocGenerics::match(CT_context,CpG))))
    CT_at_other = length(CT_muts) - CT_at_CpG
  }
  
  # Construct a table and handle missing mutation types.
  full_table = table(factor(types, levels = column_names))
  full_table["C>T at CpG"] = CT_at_CpG
  full_table["C>T other"] = CT_at_other
  df = BiocGenerics::rbind(df, full_table)
  colnames(df) = names(full_table)
  #row.names(df) = subset_variant_df$sample
  return(df)
}
