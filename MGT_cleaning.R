library(strex)
library(tidyverse)


#read in df of all patients and remove uninformative columns or columns with PDFs
full_df <- read_csv("~/data/patient_data/Atrium/MGT/NMTRC009MolecularGui_DATA_LABELS_2022-06-08_1522.csv") %>% 
  discard(function(x) all(is.na(x)) | any(grep("pdf", x))) 

#split df into a list for each patient
full_list <- split(full_df, full_df$`Subject ID`)
full_list <- full_list[grep("MGT", names(full_list))]


#make new cleaned df
patient_cleaned <- data.frame(matrix(nrow = length(full_list), dimnames = list(names(full_list), "MGT_ID")))
patient_cleaned$MGT_ID <- rownames(patient_cleaned)

#######columns we have to track
interim_therapy_columns <- colnames(full_df)[grep("Drug Name|Unit|Route|Frequency|Start Date|Stop Date", colnames(full_df))]
interim_therapy_name_columns <- colnames(full_df)[grep("Drug Name", colnames(full_df))]

medication_columns <- colnames(full_df)[grep("Medication|Reason For Administration|Unit|Route|Frequency|Start Date|Stop Date", colnames(full_df))]
medication_name_columns <- colnames(full_df)[grep("Medication", colnames(full_df))]
start_date_columns<- colnames(full_df)[grep("Start Date", colnames(full_df))]
end_date_columns<- colnames(full_df)[grep("Stop Date", colnames(full_df))]

metastasis_columns <- colnames(full_df)[grep("Metastatic Site", colnames(full_df))]

specific_response_columns <- colnames(full_df)[grep("Select most appropriate response for this time point", colnames(full_df))]
response_columns <- colnames(full_df)[grep("Indicate time|Specify time|Percent Tumor|Select most appropriate response for this time point", colnames(full_df))]

off_therapy_columns <- colnames(full_df)[grep("Event Name|Date first therapy began|Date last therapy ended|Reason for going off  therapy", colnames(full_df))]

########start making new columns
for (i in 1: length(full_list)){

patient_cleaned[names(full_list[i]), 'initial_diagnosis_type'] <- full_list[[i]] %>% 
  filter(`Event Name` == "Screening", is.na(`Repeat Instrument`)) %>% 
  select(`Select the most appropriate description of the patient's disease at the time of enrollment onto this study:`)

patient_cleaned[names(full_list[i]), 'treatment_start'] <- full_list[[i]] %>% 
  filter(`Event Name` == "Screening", is.na(`Repeat Instrument`)) %>% 
  select(`Date Treatment Started`)

patient_cleaned[names(full_list[i]), 'initial_primary_tumor_location'] <- full_list[[i]] %>% 
  filter(`Event Name` == "Screening", is.na(`Repeat Instrument`)) %>% 
  select(`Primary Tumor Site`)


#underscore separate if multiple metastasis
patient_cleaned[names(full_list[i]), 'initial_metastasis_location'] <- full_list[[i]] %>% 
  filter(`Event Name` == "Screening", is.na(`Repeat Instrument`)) %>% 
  select(all_of(metastasis_columns)) %>% 
  unite('merge', sep = "_", na.rm = TRUE)



patient_cleaned[names(full_list[i]), 'intial_stage'] <- full_list[[i]] %>% 
  filter(`Event Name` == "Screening", is.na(`Repeat Instrument`)) %>% 
  select(`Stage`)


patient_cleaned[names(full_list[i]), 'initial_histology'] <- full_list[[i]] %>% 
  filter(`Event Name` == "Screening", is.na(`Repeat Instrument`)) %>% 
  select(`Histology`)


patient_cleaned[names(full_list[i]), 'n-myc_amplified'] <- full_list[[i]] %>% 
  filter(`Event Name` == "Screening", is.na(`Repeat Instrument`)) %>% 
  select(`N-myc amplification`)

#each drug is underscore separated name, unit, route, Frequency, Start Date, Stop Date and then a plus sign
#start end end dated have "start-" or "end-" before date name, respectively, since some are missing one or the other
  temp <- full_list[[i]] %>% 
  filter(`Event Name` == "Interim Therapy (1 cycle)") %>% 
  select(all_of(interim_therapy_columns)) %>%
  discard(~all(is.na(.))) %>%
  mutate(across(any_of(interim_therapy_name_columns), ~ gsub("(.*)", "\\+\\1", .))) %>%
  mutate(across(any_of(start_date_columns), ~ gsub("(.*)", "start-\\1", .))) %>%
  mutate(across(any_of(end_date_columns), ~ gsub("(.*)", "end-\\1", .)))

if(length(temp) == 0){
  patient_cleaned[names(full_list[i]), 'post_biopsy_pre_treatment_chemotherapy'] <- NA
} else {patient_cleaned[names(full_list[i]), 'post_biopsy_pre_treatment_chemotherapy'] <- temp %>%
  unite('merge', sep = "_") %>%
  gsub("^\\+", "", .)
}

#disease status before treatment
  temp <- full_list[[i]] %>% 
  filter(`Event Name` == "Interim Therapy (1 cycle)") %>% 
  select(`Which statement best describes the patient's disease status at the end of interim therapy prior to beginning the study treatment plan?`)
if(nrow(temp) == 0){
  patient_cleaned[names(full_list[i]), 'pre_treatment_status'] <- NA
} else {patient_cleaned[names(full_list[i]), 'pre_treatment_status'] <- temp
  
}


#each medication is underscore separated name, reason, unit, route, Frequency, Start Date, Stop Date and then a plus sign
#start end end dated have "start-" or "end-" before date name, respectively, since some are missing one or the other
  temp <- full_list[[i]] %>% 
  filter(`Event Name` == "All Visits") %>% 
  select(all_of(medication_columns)) %>%
  discard(~all(is.na(.))) %>%
  mutate(across(any_of(medication_name_columns), ~ gsub("(.*)", "\\+\\1", .))) %>%
  mutate(across(any_of(start_date_columns), ~ gsub("(.*)", "start-\\1", .))) %>%
  mutate(across(any_of(end_date_columns), ~ gsub("(.*)", "end-\\1", .)))
if(length(temp) == 0){
  patient_cleaned[names(full_list[i]), 'Medications'] <- NA
} else {patient_cleaned[names(full_list[i]), 'Medications'] <- temp %>%
  unite('merge', sep = "_") %>%
  gsub("^\\+", "", .)
}


#Response column
#Select most appropriate response for this time point refers to ct/mri response, MIBG PET-CT response, and bone marrow response
#resulting column has the time point(cycle), % change of tumor ("unknown%" if empty), and the 3 response columns with
#ctmri-/mibgpet-/bm- in front ("unknown" if empty). "_" separates all entries for 1 time point and "+" separates time points 
patient_cleaned[names(full_list[i]), 'response'] <- paste(
  full_list[[i]] %>% 
  select(all_of(response_columns)) %>%
  filter(!if_all(everything(), ~is.na(.))) %>%
  mutate(`Percent Tumor Change from Enrollment` = replace_na(`Percent Tumor Change from Enrollment`, "unknown")) %>%
  mutate(`Percent Tumor Change from Enrollment` =  gsub("(.*)", "\\1%", `Percent Tumor Change from Enrollment`)) %>%
  mutate(across(all_of(specific_response_columns), ~replace_na(.x, "unknown"))) %>%
  mutate(across(specific_response_columns[1], ~ gsub("(.*)", "ctmri-\\1", .))) %>%
  mutate(across(specific_response_columns[2], ~ gsub("(.*)", "mibgpet-\\1", .))) %>%
  mutate(across(specific_response_columns[3], ~ gsub("(.*)", "bm-\\1", .))) %>%
  mutate(`Indicate time point` = gsub("Other \\(specify\\)", NA, `Indicate time point`)) %>%
  mutate(`Indicate time point` = coalesce(`Indicate time point`,`Specify time point for response evaluation`)) %>%
  select(-`Specify time point for response evaluation`) %>%
  unite("merge", sep = "_") %>%
  pull(1),
  collapse = "+"
  )

##Off therapy column, event name, has date first therapy began, date last therapy ended, reason for going off therapy. seperated by "_"
patient_cleaned[names(full_list[i]), 'off_therapy'] <- paste(
  full_list[[i]] %>% 
  filter(grepl("Off", `Event Name`)) %>%
  select(all_of(off_therapy_columns)) %>%
  mutate("Reason for going off  therapy" = replace_na(`Reason for going off  therapy`, "unknown")) %>%
  unite("merge", sep = "_") %>%
  pull(1),
  collapse = "+"
  )



###PFS column with event name, pfs in days of most recent previous relapse therapy, then pfs on this therapy in days, then ratio of the latter over the former
patient_cleaned[names(full_list[i]), 'pfs'] <- paste(full_list[[i]] %>% 
  select(`Event Name`, grep("PFS", colnames(.))) %>%
  filter(!if_all(grep("PFS", colnames(.)), ~is.na(.))) %>%
  unite("merge", sep = "_") %>%
  pull(1),
collapse = "+"
)




i <- i+1
}

patient_cleaned <- patient_cleaned %>% mutate_all(na_if,"")

