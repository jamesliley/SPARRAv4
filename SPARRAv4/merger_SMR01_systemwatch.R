## Merge SMR01 with systemwatch. TO BE COMMENTED OR REPLACED WITH SIMON'S CODE.

code_folder="./"
#source(paste0(code_folder,"utilities.R"))

library(dplyr)

## Not-in operator
`%notin%` = negate(`%in%`)



######## Functions

find_all_codes=function(.data) {
  condition_columns=c("main_condition",paste0("other_condition_",1:5))
  
  all_unique=c()
  for (col in condition_columns) {
    new_unique=unique(.data[[col]])
    new_unique=new_unique[new_unique != ""]
    all_unique=unique(c(all_unique,new_unique))
  }
  
  return(list(all_unique))
}


add_all_codes_SMR01=function(SMR01){
  SMR01 = SMR01 %>% 
    group_by(id, cis_marker) %>%
    mutate(all_stay_code_list = find_all_codes(cur_data())) %>%
    ungroup()
  return(SMR01)
}

bespoke_mapply_function = function(item_1,item_2,item_3,item_4,item_5,item_6) {
  w=unique(c(item_1,item_2,item_3,item_4,item_5,item_6))
  w=w[w !=""]
  return(list(w))
}


condition_to_list=function(input_df) {
  output_df=input_df %>%
    mutate(
      code_list=mapply(
        bespoke_mapply_function,
        main_condition,
        other_condition_1,
        other_condition_2,
        other_condition_3,
        other_condition_4,
        other_condition_5
      )
    )
  return(output_df)
}



##### This function added by JL - SystemWatch seems to be different between
#####  Azure and NSH
condition_to_list_systemwatch=function(input_df) {
  output_df=input_df %>%
    mutate(
      code_list=mapply(
        bespoke_mapply_function,
        main_condition,
        diagnosis_2,
        diagnosis_3,
        diagnosis_4,
        diagnosis_5,
        diagnosis_6
      )
    )
  return(output_df)
}




######## Check codes

check_codes = function(.data,group_set,suffixcheck=TRUE) {
  if (!suffixcheck) {
    return(
      ifelse(
        .data$main_condition %in% group_set |
          .data$other_condition_1 %in% group_set |
          .data$other_condition_2 %in% group_set |
          .data$other_condition_3 %in% group_set |
          .data$other_condition_4 %in% group_set |
          .data$other_condition_5 %in% group_set, 1, 0
      )
    )
  } else {
    return(
      ifelse(
        gsub("([A-Za-z0-9]+).*","\\1",.data$main_condition) %in% group_set |
          gsub("([A-Za-z0-9]+).*","\\1",.data$other_condition_1) %in% group_set |
          gsub("([A-Za-z0-9]+).*","\\1",.data$other_condition_2) %in% group_set |
          gsub("([A-Za-z0-9]+).*","\\1",.data$other_condition_3) %in% group_set |
          gsub("([A-Za-z0-9]+).*","\\1",.data$other_condition_4) %in% group_set |
          gsub("([A-Za-z0-9]+).*","\\1",.data$other_condition_5) %in% group_set, 1, 0
      )
    )
  }
}



######## Prepare SMR01

prepare_smr01=function(SMR01) {
  print("Preparing SMR01")
  
  print("Adding code list column")
  SMR01 = condition_to_list(SMR01)
  
  print("Adding column with all codes in stay")
  SMR01 = add_all_codes_SMR01(SMR01)
  
  ## alcohol, substance, falls, self harm
  falls_codes=get_falls_codes()
  selfharm_codes=get_selfharm_codes()
  alcohol_drug_groupings=get_icd10_grouping_drug_alcohol_selfharm()
  print("Adding alcohol, drug, self-harm flags")
  SMR01 = SMR01 %>% 
    mutate(
      alcohol_admin=check_codes(.,alcohol_drug_groupings$alcohol),
      drug_admin=check_codes(.,alcohol_drug_groupings$drug),
      falls_admin=check_codes(.,falls_codes),
      selfharm_admin=check_codes(.,selfharm_codes)
    )
  
  print("Merging CIS marker")
  SMR01 = SMR01 %>%
    group_by(id,cis_marker) %>%
    mutate(
      length_of_stay = sum(length_of_stay),
      alcohol_admin = ifelse(
        sum(alcohol_admin) > 0, 1, 0
      ),
      drug_admin=ifelse(
        sum(drug_admin) > 0, 1, 0
      ),
      falls_admin=ifelse(
        sum(falls_admin) > 0, 1, 0
      ),
      selfharm_admin=ifelse(
        sum(selfharm_admin) > 0, 1, 0
      ),
      time_discharge = max(time_discharge),
      temp_column = abs(admission_type - 18)
    ) %>%
    arrange(time, time_discharge, -temp_column) %>%
    filter(row_number() ==1) %>%
    ungroup() %>%
    mutate(
      length_of_stay = as.numeric(
        difftime(time_discharge,time, units=c("days"))
      )
    ) %>%
    select(-temp_column)
  
  print("Adding location code to SMR01")
  SMR01 = SMR01 %>% mutate(location_code= location)
  
  return(SMR01)
}



########## Prepare systemwatch

prepare_systemwatch=function(
  SystemWatch,
  time_cutoff,
  use_diagnosis_codes=FALSE
) {
  
  print("preparing systemwatch")

  
  
  ################## ADDED BY JL ##########################
  print("Changing SystemWatch column names")
  SystemWatch = SystemWatch %>% 
    rename(
      other_condition_1=diagnosis_2,
      other_condition_2=diagnosis_3,
      other_condition_3=diagnosis_4,
      other_condition_4=diagnosis_5,
      other_condition_5=diagnosis_6
    )
  
  print("Changing SystemWatch column type")
  SystemWatch$management_of_patient = as.character(SystemWatch$management_of_patient)
  SystemWatch$admission_transfer_from = as.character(SystemWatch$admission_transfer_from)
  ################## ADDED BY JL ##########################
  
  print("Adding code list column")
  SystemWatch = condition_to_list(SystemWatch)
  SystemWatch$all_stay_code_list = SystemWatch$code_list
  
  print("Replacing NA with time cutoff")
  SystemWatch = SystemWatch %>%
    tidyr::replace_na(list(time_discharge = time_cutoff))
  
  print("computing length of stay")
  SystemWatch$length_of_stay = as.numeric(
    difftime(SystemWatch$time_discharge, SystemWatch$time,units=c("days"))
  )
  SystemWatch$inpatient_daycase_identifier = "I"
  
  if (use_diagnosis_codes) {
    print("Using SW diagnosis codes to catch flags")
    
    falls_codes=get_falls_codes()
    selfharm_codes=get_selfharm_codes()
    alcohol_drug_groupings=get_icd10_grouping_drug_alcohol_selfharm()
    print("Adding alcohol, drug, self-harm flags")
    SystemWatch = SystemWatch %>% 
      mutate(
        alcohol_admin=check_codes(.,alcohol_drug_groupings$alcohol),
        drug_admin=check_codes(.,alcohol_drug_groupings$drug),
        falls_admin=check_codes(.,falls_codes),
        selfharm_admin=check_codes(.,selfharm_codes)
      )
  } else {
    print("Ignoring systemwatch diagnosis codes")
    SystemWatch = SystemWatch %>%
      mutate(
        alcohol_admin=0,
        drug_admin=0,
        falls_admin=0,
        selfharm_admin=0
      )
  }
    
 return(SystemWatch)  
  
}



##### Prepare SMR01E

prepare_smr01e = function(SMR01E) {
  
  
  print("Preparing SMR01E for merge")
  
  print("Adding code list column")
  SMR01E = condition_to_list(SMR01E)
  SMR01E$all_stay_code_list = SMR01E$code_list
  
  ################## ADDED BY JL
  print("Changing SMR01E column types")
  SMR01E$management_of_patient=as.character(SMR01E$management_of_patient)
  
  
  print("Using SMR01E diagnosis codes to add flags")
  falls_codes=get_falls_codes()
  selfharm_codes=get_selfharm_codes()
  alcohol_drug_groupings=get_icd10_grouping_drug_alcohol_selfharm()
  print("Adding alcohol, drug, self-harm flags")
  SMR01E = SMR01E %>% 
    mutate(
      alcohol_admin=check_codes(.,alcohol_drug_groupings$alcohol),
      drug_admin=check_codes(.,alcohol_drug_groupings$drug),
      falls_admin=check_codes(.,falls_codes),
      selfharm_admin=check_codes(.,selfharm_codes)
    ) %>% mutate(location_code= location)
  
  return(SMR01E)
}



######### Merge tables

merge_tables=function(SMR01,SystemWatch) {
  
  print("Performing anti-join for merge")
  sw_rows_to_add = anti_join(
    SystemWatch,
    SMR01,
    by=c("id","time","location_code")
  )
  print("Adding rows found via antijoin")
  SMR01 = SMR01 %>% bind_rows(sw_rows_to_add)
  
  return(SMR01)
}


######## Emergency indicators


emergency_indicator=function(.data) {
  
  emergency_admission_types = c(20,21,22,30,31,32,33,34,35,36,39)
  return(
    ifelse(
      .data$admission_type %in% emergency_admission_types |
        (.data$admission_type ==38 & 
           !stringr::str_sub(.data$admission_transfer_from,1,1) %in% c("4","5")),1,0
    )
  )
}

elective_indicator = function(.data) {
  elective_admission_types = c(10:12,19)
  return(
    ifelse(
      (.data$admission_type %in% elective_admission_types) &
        (.data$inpatient_daycase_identifier =="I"),1,0
    )
  )
}

final_smr01_indicators =function(SMR01) {
  
  print("Adding emergency, elective, dc, alcohol, emergency flags to merged")
  output = SMR01 %>%
    mutate(
      emergency_admin = emergency_indicator(.),
      elective_admin = elective_indicator(.),
      dc_admin = ifelse(inpatient_daycase_identifier =="D",1,0)
    ) %>%
    mutate(
      emergency_alcohol_substance = ifelse(
        emergency_admin ==1 & (alcohol_admin == 1|drug_admin==1),1,0
      ),
      emergency_alcohol = ifelse(
        emergency_admin ==1 & alcohol_admin ==1,1,0
      ),
      emergency_selfharm = ifelse(
        emergency_admin ==1 & selfharm_admin == 1,1,0
      )
    )
  
  return(output)
}





#### Main merge


smr01_systemwatch_merge = function(
  list_of_data_tables,
  time_cutoff,
  use_sw_diagnosis_codes = TRUE,
  include_smr01e= TRUE
) {
  
  if ("SMR01" %notin% names(list_of_data_tables)) stop("No SMR01")
  if ("SystemWatch" %notin% names(list_of_data_tables)) stop("No SystemWatch")
  
  n_smr01_rows_orig = nrow(list_of_data_tables$SMR01)
  
  list_of_data_tables$SMR01 = prepare_smr01(list_of_data_tables$SMR01)
  n_smr01_rows=nrow(list_of_data_tables$SMR01)
  
  
  
  if (include_smr01e) {
    if (!"SMR01E" %in% names(list_of_data_tables)) {
      print("No SMR01E, continuing without")
    } else {
      list_of_data_tables$SMR01E = prepare_smr01e(list_of_data_tables$SMR01E)
      list_of_data_tables$SMR01 = list_of_data_tables$SMR01 %>%
        bind_rows(list_of_data_tables$SMR01E)
    }
  }
  
  
  list_of_data_tables$SystemWatch = prepare_systemwatch(
    list_of_data_tables$SystemWatch,
    time_cutoff,
    use_diagnosis_codes=use_sw_diagnosis_codes
  )
  n_sw_rows = nrow(list_of_data_tables$SystemWatch)
  
  
  
  
  merged_smr01 = merge_tables(
    list_of_data_tables$SMR01,
    list_of_data_tables$SystemWatch
  )
  n_merged=nrow(merged_smr01)
  
  
  output=final_smr01_indicators(merged_smr01)
  print("Done")
  print(paste0("Orig SMR01 rows: ",n_smr01_rows_orig,", after cis merging: ",
               n_smr01_rows,", orig systemwatch: ", n_sw_rows,", after merging: ",
               n_merged))

  if(include_smr01e) print(paste0("N SMR01E rows: ",nrow(list_of_data_tables$SMR01E)))  
  print(paste0("Final rows: ",nrow(output)))
  
  return(output)
}



