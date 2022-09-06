library(redcapAPI)
library(tidyverse)

options(redcap_api_url = "https://rc2e.atriumhealth.org/redcap/api/")

rcon <- redcapConnection(token = "9EE4A3BD18ED2725C5353772E4CC0DAA", project = "NMTRC 009 - Molecular Guided Therapy Year 3")
test <- exportRecords(rcon)

exportFiles(rcon, record = "MGT9-107-08", field = "upload_prev_tx_hx", event = "screening_arm_1", dir = "~/Desktop/test")


for (patient in patients){
  skip <- FALSE
  dir.create(paste("/Volumes/delve/Patients/MGT", patient, "scans", sep = "/"))
  event_names <- test %>% filter(subject_initials == patient) %>% filter(primary_source_doc != "<NA>") %>% pull(redcap_event_name)
  for (event_name in event_names){
 tryCatch(exportFiles(rcon, record = patient, field = "primary_source_doc", event = event_name, 
                      dir = paste("/Volumes/delve/Patients/MGT", patient, "scans", sep = "/")),
  error = function(e){skip <- TRUE})
  
  if(skip){next}
  }
}

for (patient in temp$subject_initials){
  Sys.sleep(3)
  skip <- FALSE
  #dir.create(paste("/Volumes/delve/Patients/MGT", patient, "scans", sep = "/"))
  tryCatch(exportFiles(rcon, record = patient, field = "primary_source_doc", event = "additional_biopsy2_arm_1",
                       dir = paste("/Volumes/delve/Patients/MGT", patient, "scans", sep = "/")),
           error = function(e){skip <- TRUE})

    if(skip){next}
  }

for(i in 1:nrow(events)){
  temp <- test %>% select(subject_initials, redcap_event_name, pe_upload) %>% 
    filter(pe_upload != "<NA>", subject_initials %in% patients, redcap_event_name == events$redcap_event_name[i])
  for (patient in temp$subject_initials){
    Sys.sleep(3)
    skip <- FALSE
    #dir.create(paste("/Volumes/delve/Patients/MGT", patient, "systems_review", sep = "/"))
    tryCatch(exportFiles(rcon, record = patient, field = "pe_upload", event = events$redcap_event_name[i],
                         dir = paste("/Volumes/delve/Patients/MGT", patient, "systems_review", sep = "/")),
             error = function(e){skip <- TRUE})
    
    if(skip){next}
  }
}




