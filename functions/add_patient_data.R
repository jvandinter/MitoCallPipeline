add_patient_data <- function(df, patient_list, age_list) {

  # Initialise columns
  a = df

  # Add data to columns based on sample column
  a = a %>% mutate(patient = case_when(
    grepl('(^AC)([^0-9]{2})', sample) ~ patient_list[1],
    grepl(patient_list[2], sample) ~ patient_list[2],
    grepl(patient_list[3], sample) ~ patient_list[3],
    grepl(patient_list[4], sample) ~ patient_list[4],
    grepl(patient_list[5], sample) ~ patient_list[5],
    grepl(patient_list[6], sample) ~ patient_list[6],
    grepl(patient_list[7], sample) ~ patient_list[7],
    grepl(patient_list[8], sample) ~ patient_list[8],
    grepl(patient_list[9], sample) ~ patient_list[9]))

  a = a %>% mutate(age = case_when(
    grepl('(^AC)([^0-9]{2})', sample) ~ age_list[1],
    grepl(patient_list[2], sample) ~ age_list[2],
    grepl(patient_list[3], sample) ~ age_list[3],
    grepl(patient_list[4], sample) ~ age_list[4],
    grepl(patient_list[5], sample) ~ age_list[5],
    grepl(patient_list[6], sample) ~ age_list[6],
    grepl(patient_list[7], sample) ~ age_list[7],
    grepl(patient_list[8], sample) ~ age_list[8],
    grepl(patient_list[9], sample) ~ age_list[9]))
  return(a)
}
