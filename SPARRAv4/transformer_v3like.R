## New v3 transformer
## Simon Rogers, Oct 21


#' transformer_sparrav3approx_features
#' This transformer approximately replicates the features used by SPARRAv3. An exact replication of 
#' v3 features is given in the function v3_feature_creation().
#' 
#' The changes from the exact v3 specificiations are:
#'  1. A universal 3-year lookback on all tables except LTC, and an arbitrary lookback (30y) there
#'  2. Identical BNF section counts across cohorts
#'  3. Removed binarised versions of bnf_vitamins or bnf_bandages
#'  4. epilepsy_indicated and num_ltc identical across cohorts
#'  5. Removal of constant columns blood_indicated, endocrine_indicated, congenital_indicated and digestive_indicated
#'  6. Added pis_antiepileptics, pis_parkinsonism, pis_respiratory_corticosteroids, pis_minerals
#'  7. Added pis_bandages and pis_vitamins as multi-BNF categories
#'  8. Removed redundant pis_substance and corrected spelling of pis_anticoagulants
#'  9. Added 'indicated' covariates for other LTC-listed covariates, and removed BNF/LTC dual features epilepsy_indicated etc
#'  10. Used a single unified BNF grouping, adding multiple pis_ features.
#'  
#' This function was adapted from that written by Simon Rogers, Oct 21.
#' 
#' 
#' Usage example
#' transformer_sparrav3(patients, episodes, list_of_data_tables, time_cutoff)
#' 
transformer_v3like = function(patients,episodes,list_of_data_tables,
                              time_cutoffs,
                              emergency_stay_threshold=26,
                              elective_stay_threshold=19,
                              other_stay_threshold=6,
                              force_one_year_lookback=FALSE) {
  
	# Cutoff data and mutate to unixtime
	episodes = episodes %>% mutate(unixtime = as.numeric(time))
	for (i in 1:length(list_of_data_tables)) {
		list_of_data_tables[[i]] = list_of_data_tables[[i]] %>% mutate(unixtime = as.numeric(time))
	}
	
	# Lookback in days
	one_year_lookback = 1 * 366 * 3600 * 24
	three_year_lookback = (365+365+366)* 3600*24 # 3 * 366  * 3600 * 24
	ltc_lookback = 30 * 366 * 3600 * 24 # 30 years
	
	# Max lookback
	if (force_one_year_lookback) {
	  max_lookback_vector= c(
	    one_year_lookback,
	    one_year_lookback,
	    one_year_lookback,
	    one_year_lookback,
	    ltc_lookback,
	    one_year_lookback,
	    one_year_lookback,
	    one_year_lookback,
	    ltc_lookback,
	    one_year_lookback)
	} else {
	  max_lookback_vector= c(
	    three_year_lookback,
	    three_year_lookback,
	    three_year_lookback,
	    three_year_lookback,
	    ltc_lookback,
	    three_year_lookback,
	    three_year_lookback,
	    three_year_lookback,
	    ltc_lookback,
	    three_year_lookback)
	}
	
	lookback_table=data.frame(
		source_table=factor(c(
			"SMR01",
			"SMR00",
			"SMR01E",
			"SystemWatch",
			"SPARRALTC",
			"AE2",
			"SMR04",
			"PIS",
			"deaths",
			"SMR01M" # Merged SMR01, SMR01E and SystemWatch
			),
			levels=as.character(levels(episodes$source_table))),
		max_lookback = max_lookback_vector
	) %>% mutate(max_lookback = as.numeric(max_lookback))

	print(lookback_table)
	episodes = episodes %>% left_join(lookback_table, by =c("source_table"))
	
	for (t in time_cutoffs) {
		
		# Filter episodes and data tables
		filtered_episodes = episodes %>%
			filter(
				unixtime < t,
				unixtime >= t - max_lookback
			) %>%
			select(-unixtime)
		
		# PIS cutoff is one month earlier
		pis_cutoff_posix = as.POSIXct(t, origin = "1970-01-1") %m-% months(1)
		pis_cutoff = as.numeric(pis_cutoff_posix)
		
		filtered_list_of_data_tables = list_of_data_tables
		for (table_name in names(filtered_list_of_data_tables)) {
			max_l = lookback_table %>% filter(source_table == table_name)
			tc = t
			if (table_name == "PIS") {
				tc = pis_cutoff
			}
			print(paste0(table_name," ",max_l$max_lookback))
			print(table_name)
			filtered_list_of_data_tables[[table_name]] = filtered_list_of_data_tables[[table_name]] %>%
				filter(
					unixtime < tc,
					unixtime >= t - max_l$max_lookback
				)
		}

		filtered_patients = patients %>% filter(id %in% filtered_episodes$id)

		# Output data structure
		out = NULL
		
		
		
		# SMR01 + Systemwatch merge
		# This merges data and adds flags for
		# emergency, alcohol, substance abuse, daycase, elective
		print(paste0("Merging SMR01"))
		SMR01_merged=filtered_list_of_data_tables$SMR01M #merge_smr01(filtered_list_of_data_tables,time_cutoff=t)
		print(paste0("Completed merge"))
		
		
		
		# Emergency / elective admissions and bed days (needs time cutoff)
		print(paste0("Emergency / elective admissions and bed days"))
		tmp = admissions_and_bed_days(SMR01_merged,filtered_list_of_data_tables$SMR01,
		                              filtered_list_of_data_tables$SMR01E,t,three_year_lookback,
		                              emergency_stay_threshold=emergency_stay_threshold,
		                              elective_stay_threshold=elective_stay_threshold)
		out = join_transformer_outputs(out,tmp)
		

		# Number of LTCs resulting in admission
		print(paste0("Number of LTCs resulting in admissions"))
		tmp=icd10_admission(filtered_list_of_data_tables$SMR01M)
		out = join_transformer_outputs(out,tmp)
		
		
		# Features that indicate a LTC based on either an LTC or a PIS entry
		print(paste0("Indicated features"))
		tmp = indicated_features(filtered_list_of_data_tables)
		out=join_transformer_outputs(out,tmp)
		
		# LTC features
		# - counts of LTCs are different per cohort
		# - this makes a general_ltc_count (for LTC, FE) and a yed_count (for YED): TODO (U16?)
		print(paste0("LTC_features"))
		tmp = ltc_features(filtered_list_of_data_tables$SPARRALTC)
		out = join_transformer_outputs(out, tmp)
		
		# General outpatient attendances
		print(paste0("General outpatient"))
		tmp = general_outpatient(filtered_list_of_data_tables$SMR00)
		out= join_transformer_outputs(out,tmp)

		# Psych outpatient
		print(paste0("Psych outpatient"))
		tmp = psych_outpatient(filtered_list_of_data_tables$SMR00)
		out = join_transformer_outputs(out,tmp)
						
		# Psych admissions
		print(paste0("Psych admissions"))
		tmp = psych_admissions(filtered_list_of_data_tables$SMR04)
		out = join_transformer_outputs(out,tmp)
		
		
		# AE2 attendances
		# Note - had to request additional attendance code column from Fraser
		print(paste0("AE2"))
		tmp = ae2_attendances(filtered_list_of_data_tables$AE2)
		out = join_transformer_outputs(out,tmp)
		
		# PIS section counts. These are the individual PIS features.
		# e.g. pis_diabetes
		print(paste0("PIS section counts"))
		for (pis_feature in names(bnf_groupings)) {
			print(paste0(pis_feature))
			tmp = bnf_count(
				filtered_list_of_data_tables$PIS,
				toupper(bnf_groupings[[pis_feature]]),
				pis_feature
			)
			if (dim(filtered_list_of_data_tables$PIS)[1]>2e6) gc()
			out=join_transformer_outputs(out,tmp)
		}
		
		# In the LTC cohort there are two pis features computed differently
		# These convert the number to binary and then add (see p45 of specs)
		#print(paste0("PIS binary OR features"))
		#for (pis_feature in names(bnf_binary)) {
		#		print(paste0(pis_feature))
		#	tmp = bnf_binary_count(
		#		filtered_list_of_data_tables$PIS,
		#		bnf_binary[[pis_feature]],
		#		pis_feature
		#	)
		#	out=join_transformer_outputs(out,tmp)
		#}

		# Total BNF count
		print(paste0("NUM BNF TOTAL"))
		tmp = bnf_total(filtered_list_of_data_tables$PIS)
		out = join_transformer_outputs(out,tmp)
		
				
		# Number of BNF sections
		# Notes: these are based on the tables starting on p31 of the function specs.
		# Initial interpretation was that the ones listed as grouped, eg. skin, would
		# be binary (any one of the group present vs all of the group present). But it
		# seems that they are counted separately.
		print(paste0("NUM BNF SECTIONS"))
		tmp = bnf_num_sections(filtered_list_of_data_tables$PIS)
		out = join_transformer_outputs(out,tmp)

		
		# Individual LTCs for U16 cohort
		tmp = individual_LTCs(filtered_list_of_data_tables$SPARRALTC)
		out = join_transformer_outputs(out,tmp)
		
		# Age, sex and SIMD (added by JL for NSH)
		print(paste0("Age, sex and SIMD"))
		tmp=transformer_patient_info(patients,episodes,list_of_data_tables,time_cutoff=t)
		tmp$matrix$gender=as.integer(tmp$matrix$gender=="Male")
		tmp$matrix =tmp$matrix %>% rename(sexM=gender)
		names(tmp$missing_default)=gsub("gender","sexM",names(tmp$missing_default))
		out = join_transformer_outputs(out,tmp)
	}
	
	out$matrix[is.na(out$matrix)] = 0
	return(out$matrix)
}
	
	

#**********************************************************#  
# Utility methods --------------------------------
#**********************************************************#  


#**********************************************************#  
# Feature joiner --------------------------------
#**********************************************************#  

join_transformer_outputs = function(out,new_inp) {
	if (is.null(out)) {
		out = new_inp
	} else {
		out$matrix = out$matrix %>%
			full_join(new_inp$matrix,by="id")
		out$missing_default = c(out$missing_default,new_inp$missing_default)
	}
	return(out)
}



#**********************************************************#  
# BNF total count ------------------------------------------
#**********************************************************#  

bnf_total = function(PIS) {
  tmp = list()
  # Extract relevant columns, binarise number sold, keep just one of each section
  pis_sub = PIS %>%
    select(id, bnf_section,num_items) %>%
    group_by(id) %>%
    summarise(bnf_total_count=sum(num_items)) %>% #filter(row_number() == 1) %>%
    ungroup()
  
  
  # make the features
  tmp=list()
  tmp$matrix = pis_sub %>%
    transmute(
      id=id,
      num_bnf_total = bnf_total_count
    )
  tmp$missing_default = list(
    num_bnf_total = 0
  )
  return(tmp)
}




#**********************************************************#  
# BNF section count --------------------------------
#**********************************************************#  

bnf_num_sections = function(PIS) {
	tmp = list()
	# Extract relevant columns, binarise number sold, keep just one of each section
	pis_sub = PIS %>%
		select(id, bnf_section,num_items) %>%
		mutate(num_items = ifelse(num_items > 0 , 1 ,0)) %>%
		group_by(id, bnf_section) %>%
		filter(row_number() == 1) %>%
		ungroup()
	
	# Add the grouped ones (see e.g. combined nutrition, p33 of spec)
	count_data =NULL
	for (group_name in names(bnf_section_ref$groups)) {
		print(paste0(group_name))
		temp = pis_sub %>%
			filter(bnf_section %in% bnf_section_ref$groups[[group_name]]) %>%
			group_by(id) %>%
			summarise(!!group_name := n_distinct(bnf_section)) # count everything in the group
		if (is.null(count_data)) {
			count_data=temp
		} else {
			count_data = count_data %>%
				full_join(temp, by= c("id"))
		}
 	}

	# Add the individual (all, we remove the excluded ones later)
	temp0 = pis_sub %>%
		filter(bnf_section %in% bnf_section_ref$individual) #%>%
	  #	group_by(id) #%>%
	#summarise(all_sections = n_distinct(bnf_section)) # This is very slow on the NSH: rewritten below
	temp=as_tibble(data.frame(
	  id=unique(temp0$id),
	  all_sections=as.integer(
	    table(temp0$id[which(!duplicated(
	      paste0(temp0$id,as.character(temp0$bnf_section))))]
	      ))))
	
	count_data = count_data %>%
		full_join(temp, by = c("id"))

	
	# Count the sum of excluded ones
	for (excluded_name in names(bnf_section_ref$excluded)) {
		print(paste0(excluded_name))
		temp = pis_sub %>%
			filter(bnf_section %in% bnf_section_ref$excluded[[excluded_name]]) %>%
			group_by(id) %>%
			summarise(!!excluded_name := n_distinct(bnf_section))
		count_data = count_data %>%
			full_join(temp, by =c("id"))
	}
	# Do the total - i.e. all sections in the table
	temp0 = pis_sub #%>%
	#	group_by(id) #%>%
	#	summarise(num_bnf_total = n_distinct(bnf_section)) # This is very slow on the NSH: rewritten below
	ptab=table(temp0$id[which(!duplicated(paste0(temp0$id,as.character(temp0$bnf_section))))]) # Frequency table for ID-BNF_section pairs
	temp=as_tibble(data.frame(
	  id=as.integer(names(ptab)),
	  num_bnf_sections=as.integer(ptab)))
	
	
	
	count_data = count_data %>%
		full_join(temp, by=c("id"))
	
	# Fill in NAs with 0
	count_data[is.na(count_data)] = 0
	
	# make the features
	tmp$matrix = count_data %>%
		transmute(
			id=id,
			num_bnf_sections = num_bnf_sections
		)
	tmp$missing_default = list(
		num_bnf_all = 0,
		num_bnf_sections = 0
	)
	return(tmp)
}


#**********************************************************#  
# Standard count BNF features ---------------------------
#**********************************************************#  

bnf_count = function(PIS, section_list, feature_name) {
	tmp = list()
	tmp$matrix = PIS %>%
		select(id, bnf_section, num_items) %>%
		filter(bnf_section %in% section_list) %>%
		group_by(id) %>%
		summarise(!!feature_name := sum(num_items))
	tmp$missing_default = list()
	tmp$missing_default[feature_name] = 0
	return(tmp)
}


#**********************************************************#  
# AE2 attendances --------------------------------
#**********************************************************#  

# Check if a particular code is in attendance record
check_codes_ae2 = function(.data,group_set) {
  return(
    ifelse(
      .data$DIAGNOSIS_1 %in% group_set |
        .data$DIAGNOSIS_2 %in% group_set |
        .data$DIAGNOSIS_3 %in% group_set, 1, 0
    )
  )
}



ae2_attendances = function(AE2) {
  
  tmp = list()
  tmp$matrix = AE2 %>% 
    mutate(
      alcohol_drug_attendance=check_codes_ae2(.,2),
      psych_attendance=check_codes_ae2(.,16)
    ) %>%
    group_by(id) %>%
    summarise(num_ae2_attendances = n(),
              num_alcohol_drug_attendances = sum(alcohol_drug_attendance),
              num_psych_attendances = sum(psych_attendance))
  tmp$missing_default=list(
    num_ae2_attendances = 0,
    num_alcohol_drug_attendances = 0,
    num_psych_attendances = 0
  )
  return(tmp)
}




#**********************************************************#  
# Psych admissions --------------------------------
#**********************************************************#  

psych_admissions = function(SMR04) {
	tmp = list()
	tmp$matrix = SMR04 %>%
		group_by(id) %>%
		summarise(
			num_psych_admissions = n_distinct(CIS_MARKER)
		)
	tmp$missing_default = list(
		num_psych_admissions = 0
	)
	return(tmp)
}



#**********************************************************#  
# Admissions and bed days --------------------------------
#**********************************************************#  

admissions_and_bed_days = function(SMR01_merged,SMR01,SMR01E,time_cutoff,
                                   max_lookback=round(3*365.25),
                                   emergency_stay_threshold=26,
                                   elective_stay_threshold=19,
                                   other_stay_threshold=6) {
	tmp = list()
	
	# Construct length of stay variable, truncating appropriately
	xlos=round((SMR01_merged$time_discharge-SMR01_merged$time)/(24*3600))
	trm=SMR01_merged$time>=time_cutoff
	trunc=SMR01_merged$time_discharge>=time_cutoff
	em=which(trunc==TRUE & trm==FALSE & SMR01_merged$emergency_admin==1)
	el=which(trunc==TRUE & trm==FALSE & SMR01_merged$elective_admin==1)
	eo=which(trunc==TRUE & trm==FALSE & SMR01_merged$emergency_admin==0 & SMR01_merged$elective_admin==0)
	for (i in 1:length(em)){
	  if ((i %% 100)==0) print(paste0("Emergency length of stay check ",i," of ",length(em)))
	  idx=SMR01_merged$id[em[i]]; 
	  tx=SMR01_merged$time[em[i]];
	  xcis=SMR01_merged$cis_marker[em[i]]
	  
	  if (!is.finite(xcis)) {
	    los2=(time_cutoff-as.numeric(SMR01_merged$time[em[i]]))/(24*3600)
	    xlos[em[i]]=pmin(los2,emergency_stay_threshold)
	  } else {
	    # SMR01, SystemWatch and SMR01E tables corresponding to this admission
	    x01=SMR01[which(SMR01$id==idx),]; x01=x01[which(x01$CIS_MARKER==xcis & x01$time>=tx & x01$time<= time_cutoff),]
	    x01e=SMR01E[which(SMR01E$id==idx),]; x01e=x01e[which(x01e$CIS_MARKER==xcis & x01e$time>tx & x01e$time<= time_cutoff),]
	    
	    # All discharge dates corresponding to this admission
	    dtx=c(x01$time_discharge, x01e$time_discharge)
	    dtx=dtx[which(dtx < time_cutoff)]
	    
	    # Count all known full episodes (admission+discharge date) towards length of stay,
	    #  plus a truncated estimate for final stay before time cutoff
	    if (length(dtx)>0) {
	      xlos[em[i]]=(as.numeric(max(dtx))-as.numeric(tx))/(24*3600) + 
	        pmin((as.numeric(time_cutoff)-as.numeric(max(dtx)))/(24*3600),emergency_stay_threshold)
	    } else {
	      los2=(time_cutoff-as.numeric(SMR01_merged$time[em[i]]))/(24*3600)
	      xlos[em[i]]=pmin(los2,emergency_stay_threshold)
	    }  
	  }
	}
	for (i in 1:length(el)){
	  if ((i %% 100)==0) print(paste0("Elective length of stay check ",i," of ",length(el)))
	  idx=SMR01_merged$id[el[i]]; 
	  tx=SMR01_merged$time[el[i]];
	  xcis=SMR01_merged$cis_marker[el[i]]
	  
	  if (!is.finite(xcis)) {
	    los2=(time_cutoff-as.numeric(SMR01_merged$time[el[i]]))/(24*3600)
	    xlos[el[i]]=pmin(los2,elective_stay_threshold) 
	  } else {
	    # SMR01, SystemWatch and SMR01E tables corresponding to this admission
	    x01=SMR01[which(SMR01$id==idx),]; x01=x01[which(x01$CIS_MARKER==xcis & x01$time>=tx & x01$time<= time_cutoff),]
	    x01e=SMR01E[which(SMR01E$id==idx),]; x01e=x01e[which(x01e$CIS_MARKER==xcis & x01e$time>tx & x01e$time<= time_cutoff),]
	    
	    # All discharge dates corresponding to this admission
	    dtx=c(x01$time_discharge, x01e$time_discharge)
	    dtx=dtx[which(dtx < time_cutoff)]
	    
	    # Count all known full episodes (admission+discharge date) towards length of stay,
	    #  plus a truncated estimate for final stay before time cutoff
	    if (length(dtx)>0) {
	      xlos[el[i]]=(as.numeric(max(dtx))-as.numeric(tx))/(24*3600) + 
	        pmin((as.numeric(time_cutoff)-as.numeric(max(dtx)))/(24*3600),elective_stay_threshold)
	    } else {
	      los2=(time_cutoff-as.numeric(SMR01_merged$time[el[i]]))/(24*3600)
	      xlos[el[i]]=pmin(los2,elective_stay_threshold)
	    }  
	  }
	}
	for (i in 1:length(eo)){
	  if ((i %% 100)==0) print(paste0("Other length of stay check ",i," of ",length(eo)))
	  idx=SMR01_merged$id[eo[i]]; 
	  tx=SMR01_merged$time[eo[i]];
	  xcis=SMR01_merged$cis_marker[eo[i]]
	  
	  if (!is.finite(xcis)) {
	    los2=(time_cutoff-as.numeric(SMR01_merged$time[eo[i]]))/(24*3600)
	    xlos[eo[i]]=pmin(los2,other_stay_threshold) 
	  } else {
	    # SMR01, SystemWatch and SMR01E tables corresponding to this admission
	    x01=SMR01[which(SMR01$id==idx),]; x01=x01[which(x01$CIS_MARKER==xcis & x01$time>=tx & x01$time<= time_cutoff),]
	    x01e=SMR01E[which(SMR01E$id==idx),]; x01e=x01e[which(x01e$CIS_MARKER==xcis & x01e$time>tx & x01e$time<= time_cutoff),]
	    
	    # All discharge dates corresponding to this admission
	    dtx=c(x01$time_discharge, x01e$time_discharge)
	    dtx=dtx[which(dtx < time_cutoff)]
	    
	    # Count all known full episodes (admission+discharge date) towards length of stay,
	    #  plus a truncated estimate for final stay before time cutoff
	    if (length(dtx)>0) {
	      xlos[eo[i]]=(as.numeric(max(dtx))-as.numeric(tx))/(24*3600) + 
	        pmin((as.numeric(time_cutoff)-as.numeric(max(dtx)))/(24*3600),other_stay_threshold)
	    } else {
	      los2=(time_cutoff-as.numeric(SMR01_merged$time[eo[i]]))/(24*3600)
	      xlos[eo[i]]=pmin(los2,other_stay_threshold)
	    }  
	  }
	}
	
	xlos[which(trm==TRUE)]=0
	SMR01_merged$length_of_stay=as.numeric(xlos)
	
	# Remove records beginning before earliest time threshold
	w=which(SMR01_merged$unixtime >= as.numeric(time_cutoff) - as.numeric(max_lookback))
	SMR01_merged=SMR01_merged[w,]
	
	
	
	tmp$matrix = SMR01_merged %>%
		group_by(id) %>%
		summarise(
			num_emergency_admissions = sum(emergency_admin),
			emergency_bed_days = sum(emergency_admin * length_of_stay),
			num_elective_admissions = sum(elective_admin),
			elective_bed_days = sum(elective_admin * length_of_stay),
			num_other_admissions = sum((1-emergency_admin)*(1-elective_admin)),
			other_bed_days = sum((1-elective_admin)*(1-emergency_admin)*length_of_stay),
		  num_alcohol_substance_admissions = sum(emergency_alcohol_substance),
		  num_dc_admissions = sum(dc_admin),
			num_alcohol_admissions = sum(emergency_alcohol),
		  num_emergency_selfharm=sum(emergency_selfharm)
		)
	tmp$missing_default = list(
		num_emergency_admissions = 0,
		emergency_bed_days = 0,
		num_elective_admissions = 0,
		elective_bed_days = 0,
		num_other_admissions = 0,
		other_bed_days = 0,
		alcohol_substance = 0,
		num_dc_admissions = 0,
		num_alcohol_admissions = 0,
		num_elective_dc_admissions = 0,
	  num_emergency_selfharm = 0
	)
	
	return(tmp)
}


#**********************************************************#  
# ICD10 groups resulting in admission ---------------------
#**********************************************************#  

#' Inputs
#' @param patients A tibble of patients with related basic information
#' @param episodes A tibble of individual records with related basic information
#' @param list_of_data_tables A list of the individual data tables with all detailed columns from the raw data
#' 
#' Output
#' @return A named list with two elements: \itemize{
#'    \item{matrix}{a PATIENTS (rows) x FEATURES (columns) matrix, one patient per row, the two columns represent 
#'                        "number of alcohol/drug related admissions" and "number of different LTCs which have resulted in admission" }
#'    \item{missing_default}{a FEATURES (length) list of default value (0 here) per feature to set for all patients that have that particular feature as missing}
#' }
#' 
#' Adapted from a function by Sam Oduro
#' 
#' @export
icd10_admission <- function(SMR01M){

  ltc_groupings <- get_icd10_grouping_ltc()

  # Add "ltc_" prefix to names and "_admission" postfix
  names(ltc_groupings) <- sapply(names(ltc_groupings), function(x) str_c(str_c("ltc_", x), "_admission"))
  
  data_inpatient <- SMR01M  %>% 
    select(id,time,code_list,emergency_admin,elective_admin)

  # regex extracts first word; this ensures ICD10 codes such as 'C180 D' count as 'C180'.
  codelist=data_inpatient$code_list
  codelist=lapply(codelist,function(x) gsub("([A-Za-z0-9]+).*","\\1",x))
  
  # Transform
  for (ltc_name in names(ltc_groupings)) {
    lname=unlist(lapply(codelist,
                 function(x) ifelse(any(x %in% ltc_groupings[[ltc_name]]),1,0)))
    data_inpatient <- data_inpatient %>% mutate(!!ltc_name:= lname)
  }

  # Num LTCs resulting in any admission
  di1=data_inpatient %>% group_by(id) %>%
    summarise(across(starts_with("ltc"),~sum(.))) %>%
    mutate(across(starts_with("ltc"), function(x) x>0)) %>%
    transmute(id=id,
              numLTCs_resulting_in_admin = rowSums(select(., starts_with("ltc_")))
    )

  # Num LTCs resulting in emergency admission
  di2=data_inpatient; 
  for (i in 6:dim(di2)[2]) di2[,i]=di2[,i]*di2$emergency_admin
  di2=di2 %>% group_by(id) %>%
    summarise(across(starts_with("ltc"),~sum(.))) %>%
    mutate(across(starts_with("ltc"), function(x) x>0)) %>%
    transmute(id=id,
              numLTCs_resulting_in_emergency_admin = rowSums(select(., starts_with("ltc_")))
    )
  
  # Num LTCs resulting in elective admission
  di3=data_inpatient; 
  for (i in 6:dim(di3)[2]) di3[,i]=di3[,i]*di3$elective_admin
  di3=di3 %>% group_by(id) %>%
    summarise(across(starts_with("ltc"),~sum(.))) %>%
    mutate(across(starts_with("ltc"), function(x) x>0)) %>%
    transmute(id=id,
              numLTCs_resulting_in_elective_admin = rowSums(select(., starts_with("ltc_")))
    )

  # Num LTCs resulting in other admission
  di4=data_inpatient; 
  for (i in 6:dim(di4)[2]) di4[,i]=di4[,i]*(1-di4$emergency_admin)*(1-di4$elective_admin)
  di4=di4 %>% group_by(id) %>%
    summarise(across(starts_with("ltc"),~sum(.))) %>%
    mutate(across(starts_with("ltc"), function(x) x>0)) %>%
    transmute(id=id,
              numLTCs_resulting_in_other_admin = rowSums(select(., starts_with("ltc_")))
    )
  
  
    
  dix = di1 %>% full_join(di2,by="id") %>% full_join(di3,by="id") %>% full_join(di4,by="id")
  
  out=list(matrix=dix,
           missing_default=list(
             num_LTCs_resulting_in_admin=0,
             num_LTCs_resulting_in_emergency_admin=0,
             num_LTCs_resulting_in_elective_admin=0,
             num_LTCs_resulting_in_other_admin=0))
  
  # Output a list containing the final matrix and the defaults (0 for each category)
  return(out)
}





#**********************************************************#  
# LTC/prescription indicated features -----------
#**********************************************************#  

indicated_features = function(list_of_data_tables) {
	
  # MS
	indicated_table = indicated(
		list_of_data_tables$SPARRALTC,
		list_of_data_tables$PIS,
		"FIRST_MULTIPLE_SCLEROSIS_EPISODE",
		"num_bnf_1002"
	) %>% transmute(id = id, MS_indicated = indicated)
	
	# Parkinsons
	indicated_table = indicated_table %>%
		full_join(
			indicated(
				list_of_data_tables$SPARRALTC,
				list_of_data_tables$PIS,
				"FIRST_PARKINSON_DISEASE_EPISODE",
				"num_bnf_0409"
			) %>% transmute(id = id, parkinsons_indicated = indicated),
			by = c("id")
		)

	# Epilepsy
	indicated_table = indicated_table %>%
		full_join(
			indicated(
				list_of_data_tables$SPARRALTC,
				list_of_data_tables$PIS,
				"FIRST_EPILEPSY_EPISODE",
				"num_bnf_0408"
			) %>% transmute(id = id, epilepsy_indicated = indicated),
			by = c("id")
		)

		
	# Dementia
	indicated_table = indicated_table %>%
		full_join(
			indicated(
				list_of_data_tables$SPARRALTC,
				list_of_data_tables$PIS,
				"FIRST_DEMENTIA_EPISODE",
				"num_bnf_0411"
			) %>% transmute(id = id, dementia_indicated = indicated),
			by = c("id")
		)

	# Diabetes
	indicated_table = indicated_table %>%
	  full_join(
	    indicated(
	      list_of_data_tables$SPARRALTC,
	      list_of_data_tables$PIS,
	      "FIRST_DIABETES_EPISODE",
	      "num_bnf_0601"
	    ) %>% transmute(id = id, diabetes_indicated = indicated),
	    by = c("id")
	  )


	indicated_table[is.na(indicated_table)] = 0
	
	tmp = list()
	tmp$matrix = indicated_table
	tmp$missing_default = list(
		MS_indcated = 0,
		parkinsons_indicated = 0,
	  epilepsy_indicated = 0,
	  epilepsy_indicated_yed = 0,
	  diabetes_indicated = 0
	)
	return(tmp)	
}
	
# Method that adds a 1 to an ID if a particular ltc_term is in the ltc_data
# OR a particular prescription is present at least once
indicated = function(ltc_data,pis_data,ltc_term,pis_term) {
	# Check for the LTC
	indicated_ltc = ltc_data %>% filter(LTC_TYPE == ltc_term) %>%
		transmute(id = id,ltc = 1)
	# Check the PIC
	if (!is.null(pis_term)) {
	indicated_pis = pis_data %>% filter(bnf_section == pis_term) %>%
		group_by(id) %>%
		filter(row_number() == 1) %>% # note - just take the first if there is more than one
		ungroup() %>%
		transmute(id=id,pis=1)
	
	# Combine the two
	indicated = indicated_ltc %>%
		full_join(indicated_pis, by = c("id")) %>%
		mutate(
			ltc = replace_na(ltc,0),
			pis = replace_na(pis,0)
		) %>%
		transmute(id = id, indicated = ltc + pis)
	} else {
	  indicated=indicated_ltc %>%
	    mutate(
	      ltc = replace_na(ltc,0),
	    ) %>%
	    transmute(id = id, indicated = ltc)
	}
	return(indicated)
}




#**********************************************************#  
# General outpatient features -----------------------------
#**********************************************************#  

# Note: REFERRAL_TYPE < 3

general_outpatient = function(SMR00) {
	codes_to_include = c(
		"A1", "A2", "A3", "A6", "A7", "A8", "A81", "A82", 
		"A9", "AB", "AD", "H2", "AF", "CA", "AG", "AH", "AM",
		"AP", "AQ", "AR", "C1", "C11", "C12", "C13", "C3",
		"C31", "C4", "C41", "C42", "C5", "C51", "C6", "C7", "C8",
		"C9", "CB", "D3", "D4", "D6", "E12", "F2", "G1", "G1A",
		"G2", "G21", "G22", "G3", "G4", "G5", "G6", "H1", "J3",
		"J4", "J5"
	)
	num_op_appointments = SMR00 %>%
		filter(
			SPECIALTY %in% codes_to_include,
			REFERRAL_TYPE < 3
		) %>%
		group_by(id) %>%
		summarise(num_outpatient_appointment_general = n())
	num_op_appointments_followup = SMR00 %>%
	  filter(
	    SPECIALTY %in% codes_to_include,
	    REFERRAL_TYPE == 3
	  ) %>%
	  group_by(id) %>%
	  summarise(num_outpatient_appointment_followup_general = n())
	
	tmp= list()
	tmp$matrix = num_op_appointments
	tmp$matrix = full_join(tmp$matrix,num_op_appointments_followup,by="id")
	tmp$missing_default = list(num_outpatient_appointment_general = 0,
	                           num_outpatient_appointment_followup_general = 0)
	
	return(tmp)
}



#**********************************************************#  
# Psychiatric outpatient features -------------------------
#**********************************************************#  

# Note: REFERRAL_TYPE < 3

psych_outpatient = function(SMR00) {
  codes_to_include = c(
    "G1", "G1A", "G2", "G21", "G22", "G3", "G4", "G5", "G6"
  )
  num_psych_op_appointments = SMR00 %>%
    filter(
      SPECIALTY %in% codes_to_include,
      REFERRAL_TYPE < 3
    ) %>%
    group_by(id) %>%
    summarise(num_outpatient_appointment_psych = n())
  num_psych_op_appointments_followup = SMR00 %>%
    filter(
      SPECIALTY %in% codes_to_include,
      REFERRAL_TYPE == 3
    ) %>%
    group_by(id) %>%
    summarise(num_outpatient_appointment_followup_psych = n())
  
  tmp=list()
  tmp$matrix = num_psych_op_appointments
  tmp$matrix = full_join(tmp$matrix,num_psych_op_appointments_followup,by="id")
  tmp$missing_default = list(num_outpatient_appointment_psych = 0,
                             num_outpatient_appointment_followup_psych = 0)
  

  return(tmp)
}







#**********************************************************#  
# LTC features --------------------------------
#**********************************************************#  

general_ltc_names=c("FIRST_ARTHRITIS_EPISODE", 
	"FIRST_ASTHMA_EPISODE",
	"FIRST_ATRIAL_FIBRILLATION_EPISODE", 
	"FIRST_COPD_EPISODE", 
	"FIRST_CANCER_EPISODE", 
	"FIRST_CEREBROVASCULAR_DISEASE_EPISODE", 
	"FIRST_CHRONIC_LIVER_DISEASE_EPISODE",
	"FIRST_DEMENTIA_EPISODE", 
	"FIRST_DIABETES_EPISODE",
	"FIRST_EPILEPSY_EPISODE", 
	"FIRST_HEART_DISEASE_EPISODE", 
	"FIRST_HEART_FAILURE_EPISODE", 
	"FIRST_MULTIPLE_SCLEROSIS_EPISODE", 
	"FIRST_PARKINSON_DISEASE_EPISODE", 
	"FIRST_RENAL_FAILURE_EPISODE",
	"FIRST_CONGENITAL_PROBLEMS_EPISODE",
	"FIRST_DISEASES_OF_THE_BLOOD_AND_BLOOD_FORMING_ORGANS_EPISODE",
	"FIRST_OTHER_ENDOCRINE_METABOLIC_DISEASES_EPISODE",
	"FIRST_OTHER_DISEASES_OF_DIGESTIVE_SYSTEM_EPISODE"
)

ltc_features = function(ltc_table) {
	# ltc_table should be the sparse-style table
	ltc_sub_table =ltc_table %>%
		select(id, LTC_TYPE, time)
	
	wide_table = ltc_sub_table %>%
		pivot_wider(names_from = "LTC_TYPE",values_from="time")
	
	general_ltc_names_sub = c()
	# filter names as some might not be present
	for (name in general_ltc_names) {
		if (name %in% names(wide_table)) {
			general_ltc_names_sub = c(general_ltc_names_sub,name)
		}
	}

	wide_table$num_general_ltc = rowSums(!is.na(wide_table[general_ltc_names_sub]))

	wide_table = wide_table %>%
		select(id,num_general_ltc)
	
	tmp = list()
	tmp$matrix = wide_table
	tmp$missing_default = list(
		num_general_ltc = 0
	)
	return(tmp)
}



#**********************************************************#  
# Individual LTCs --------------------------------
#**********************************************************#  




individual_LTCs = function(ltc_table) {
  name_map = list(
    "ARTHRITIS"="FIRST_ARTHRITIS_EPISODE", 
    "ASTHMA"="FIRST_ASTHMA_EPISODE",
    "ATRIAL_FIBRILLATION"="FIRST_ATRIAL_FIBRILLATION_EPISODE", 
    "COPD"="FIRST_COPD_EPISODE", 
    "CANCER"="FIRST_CANCER_EPISODE", 
    "CEREBROVASCULAR_DISEASE"="FIRST_CEREBROVASCULAR_DISEASE_EPISODE", 
    "CHRONIC_LIVER_DISEASE"="FIRST_CHRONIC_LIVER_DISEASE_EPISODE",
    "DEMENTIA"="FIRST_DEMENTIA_EPISODE", 
    "DIABETES"="FIRST_DIABETES_EPISODE",
    "EPILEPSY"="FIRST_EPILEPSY_EPISODE", 
    "HEART_DISEASE"="FIRST_HEART_DISEASE_EPISODE", 
    "HEART_FAILURE"="FIRST_HEART_FAILURE_EPISODE", 
    "MULTIPLE_SCLEROSIS"="FIRST_MULTIPLE_SCLEROSIS_EPISODE", 
    "PARKINSON_DISEASE"="FIRST_PARKINSON_DISEASE_EPISODE", 
    "RENAL_FAILURE"="FIRST_RENAL_FAILURE_EPISODE",
    "CONGENITAL_PROBLEMS"="FIRST_CONGENITAL_PROBLEMS_EPISODE",
    "DIS_BLOOD"="FIRST_DISEASES_OF_THE_BLOOD_AND_BLOOD_FORMING_ORGANS_EPISODE",
    "ENDOCRINE_MET"="FIRST_OTHER_ENDOCRINE_METABOLIC_DISEASES_EPISODE",
    "OTHER_DIGESTIVE"="FIRST_OTHER_DISEASES_OF_DIGESTIVE_SYSTEM_EPISODE"
  )
  tmp = list()
  tmp$matrix = NULL
  tmp$missing_default=list()
  for (feature_name in names(name_map)) {
    ltc_name = name_map[[feature_name]]
    print(paste0(ltc_name," to ",feature_name))
    tt = ltc_table %>%
      filter(LTC_TYPE==ltc_name) %>%
      group_by(id) %>%
      summarise(!!feature_name := 1)
    if (is.null(tmp$matrix)) {
      tmp$matrix = tt
    } else {
      tmp$matrix = full_join(tmp$matrix,tt,by="id")
    }
    tmp$missing_default[[feature_name]]=0
  }
  tmp$matrix[is.na(tmp$matrix)]=0
  return(tmp)
}



#**********************************************************#  
# PIS definitions --------------------------------
#**********************************************************#  

bnf_groupings=list()
bnf_groupings$pis_gastro_int=c(
  "num_bnf_0101","num_bnf_0102","num_bnf_0103","num_bnf_0105",
  "num_bnf_0106","num_bnf_0107","num_bnf_0109"
)
bnf_groupings$pis_respiratory = c("num_bnf_0301", "num_bnf_0302", 
                                  "num_bnf_0303", "num_bnf_0304", "num_bnf_0305", "num_bnf_0306",
                                  "num_bnf_0307", "num_bnf_0308", "num_bnf_0309", "num_bnf_0310")
bnf_groupings$pis_cns = c("num_bnf_0401", "num_bnf_0402", "num_bnf_0403",
                          "num_bnf_0404", "num_bnf_0405", "num_bnf_0406", "num_bnf_0407",
                          "num_bnf_0408", "num_bnf_0409", "num_bnf_0410", "num_bnf_0411")
bnf_groupings$pis_infections = c("num_bnf_0501", "num_bnf_0502", "num_bnf_0503", "num_bnf_0504",
                                 "num_bnf_0505")
bnf_groupings$pis_endocrine = c( "num_bnf_0601", "num_bnf_0602", "num_bnf_0603",
                                 "num_bnf_0604", "num_bnf_0605", "num_bnf_0606", "num_bnf_0607")
bnf_groupings$pis_incontinence = c( "num_bnf_2201", "num_bnf_2202",
                                    "num_bnf_2205", "num_bnf_2210", "num_bnf_2215", "num_bnf_2220",
                                    "num_bnf_2230", "num_bnf_2240", "num_bnf_2250", "num_bnf_2260",
                                    "num_bnf_2270", "num_bnf_2280", "num_bnf_2285", "num_bnf_2290")
bnf_groupings$pis_stoma = c( "num_bnf_2305", "num_bnf_2310", "num_bnf_2315", "num_bnf_2320",
                             "num_bnf_2325", "num_bnf_2330", "num_bnf_2335", "num_bnf_2340",
                             "num_bnf_2345", "num_bnf_2346", "num_bnf_2350", "num_bnf_2355",
                             "num_bnf_2360", "num_bnf_2365", "num_bnf_2370", "num_bnf_2375",
                             "num_bnf_2380", "num_bnf_2385", "num_bnf_2390", "num_bnf_2392",
                             "num_bnf_2393", "num_bnf_2394", "num_bnf_2396", "num_bnf_2398")
bnf_groupings$pis_nutrition = c(
  "num_bnf_0904", "num_bnf_0908","num_bnf_0909","num_bnf_0911","num_bnf_0912"
)
bnf_groupings$pis_skin = c(
  "num_bnf_1301" ,"num_bnf_1302" ,"num_bnf_1303" ,"num_bnf_1304" ,
  "num_bnf_1305" ,"num_bnf_1308" ,"num_bnf_1310" ,"num_bnf_1311" ,
  "num_bnf_1313" ,"num_bnf_1314"
)
bnf_groupings$pis_supplements = c("num_bnf_0911","num_bnf_0912")
bnf_groupings$pis_vitamins=c("num_bnf_0905","num_bnf_0906")
bnf_groupings$pis_bandages=c("num_bnf_2002","num_bnf_2003")

bnf_groupings$pis_gut_motility=c("num_bnf_0102")
bnf_groupings$pis_antisecretory=c("num_bnf_0103")
bnf_groupings$pis_intestinal=c("num_bnf_0109")
bnf_groupings$pis_anticoagulant=c("num_bnf_0208")
bnf_groupings$pis_antifibrinolytic=c("num_bnf_0211")
bnf_groupings$pis_respiratory_corticosteroids=c("num_bnf_0302")
bnf_groupings$pis_bronco=c("num_bnf_0301")
bnf_groupings$pis_cromo=c("num_bnf_0303")
bnf_groupings$pis_mucolytics=c("num_bnf_0307")
bnf_groupings$pis_antiepileptic_Drugs=c("num_bnf_0408")
bnf_groupings$pis_parkinsonism=c("num_bnf_0409")
bnf_groupings$pis_sub_depend=c("num_bnf_0410")
bnf_groupings$pis_dementia=c("num_bnf_0411")
bnf_groupings$pis_antibacterial=c("num_bnf_0501")
bnf_groupings$pis_diabetes=c("num_bnf_0601")
bnf_groupings$pis_corticosteroids=c("num_bnf_0603")
bnf_groupings$pis_fluids=c("num_bnf_0902")
bnf_groupings$pis_Minerals=c("num_bnf_0905")
bnf_groupings$pis_rheumatic=c("num_bnf_1001")
bnf_groupings$pis_neuromuscular=c("num_bnf_1002")
bnf_groupings$pis_mydriatics=c("num_bnf_1105")
bnf_groupings$pis_catheters=c("num_bnf_2102")
bnf_groupings$pis_inotropic = c("num_bnf_0201")
bnf_groupings$pis_diuretics = c("num_bnf_0202")
bnf_groupings$pis_antiarrhythmics = c("num_bnf_0203")
bnf_groupings$pis_betablockers = c("num_bnf_0204")
bnf_groupings$pis_hypertensive_heart_failure = c("num_bnf_0205")
bnf_groupings$pis_antianginal = c("num_bnf_0206")
bnf_groupings$pis_antiplatelets = c("num_bnf_0209")
bnf_groupings$pis_lipid = c("num_bnf_0212")
bnf_groupings$pis_genitourinary = c("num_bnf_0704")
bnf_groupings$pis_cytotoxics = c("num_bnf_0801")
bnf_groupings$pis_immune = c("num_bnf_0802")
bnf_groupings$pis_sex_hormone_antagonists = c("num_bnf_0803")
bnf_groupings$pis_antianaemics = c("num_bnf_0901")
bnf_groupings$pis_topical_pain_relief = c("num_bnf_1003")
bnf_groupings$pis_antibacterial_eyes = c("num_bnf_1103")
bnf_groupings$pis_antiinflammatory_corticosteroids = c("num_bnf_1104")
bnf_groupings$pis_glaucoma = c("num_bnf_1106")
bnf_groupings$pis_local_anaesthetics = c("num_bnf_1107")
bnf_groupings$pis_ophthalmic = c("num_bnf_1108")
bnf_groupings$pis_ear = c("num_bnf_1201")
bnf_groupings$pis_nose = c("num_bnf_1202")
bnf_groupings$pis_oropharynx = c("num_bnf_1203")
bnf_groupings$pis_hosiery = c("num_bnf_2107")
bnf_groupings$pis_metabolic = c("num_bnf_0908")
bnf_groupings$pis_food = c("num_bnf_0909")


# BNF sections for the total count
# has various fields: individaul = all individual codes in the table on p31 of specs
# groups = codes that should be considered as a group
# excluded = the codes that should be excluded for ltc and fe
bnf_section_ref=list()
bnf_section_ref$individual=c(
	"num_bnf_0101","num_bnf_0102","num_bnf_0103","num_bnf_0105",
	"num_bnf_0106","num_bnf_0107","num_bnf_0109",
	
"num_bnf_0201" ,"num_bnf_0202" ,"num_bnf_0203" ,"num_bnf_0204" ,
	"num_bnf_0205" ,"num_bnf_0206" ,"num_bnf_0208" ,"num_bnf_0209" ,
	"num_bnf_0211" ,"num_bnf_0212",

"num_bnf_0301" ,"num_bnf_0302" ,"num_bnf_0303" ,"num_bnf_0304" ,
	"num_bnf_0306" ,"num_bnf_0307" ,"num_bnf_0309" ,"num_bnf_0310",

"num_bnf_0401" ,"num_bnf_0402" ,"num_bnf_0403" ,"num_bnf_0404" ,
	"num_bnf_0405" ,"num_bnf_0406" ,"num_bnf_0407" ,"num_bnf_0408" ,
	"num_bnf_0409" ,"num_bnf_0410" ,"num_bnf_0411",

"num_bnf_0501" ,"num_bnf_0502" ,"num_bnf_0503" ,"num_bnf_0504" ,
	"num_bnf_0505",

"num_bnf_0601" ,"num_bnf_0602" ,"num_bnf_0603" ,"num_bnf_0604" ,
	"num_bnf_0605" ,"num_bnf_0605" ,"num_bnf_0607",

"num_bnf_0704",

"num_bnf_0801" ,"num_bnf_0802" ,"num_bnf_0803",

"num_bnf_0901" ,"num_bnf_0902" ,"num_bnf_0904" ,"num_bnf_0905" ,
	"num_bnf_0906",

"num_bnf_1001" ,"num_bnf_1002" ,"num_bnf_1003",

"num_bnf_1103" ,"num_bnf_1104" ,"num_bnf_1105" ,"num_bnf_1106" ,
	"num_bnf_1107" ,"num_bnf_1108",

"num_bnf_1201" ,"num_bnf_1202" ,"num_bnf_1203",

"num_bnf_2002" ,"num_bnf_2003",

"num_bnf_2102" ,"num_bnf_2107"
)


# BNF sections for the total count
# has various fields: individaul = all individual codes in the table on p31 of specs
# groups = codes that should be considered as a group
bnf_section_ref=list()
bnf_section_ref$individual=c(
  "num_bnf_0101","num_bnf_0102","num_bnf_0103","num_bnf_0105",
  "num_bnf_0106","num_bnf_0107","num_bnf_0109",
  
  "num_bnf_0201" ,"num_bnf_0202" ,"num_bnf_0203" ,"num_bnf_0204" ,
  "num_bnf_0205" ,"num_bnf_0206" ,"num_bnf_0208" ,"num_bnf_0209" ,
  "num_bnf_0211" ,"num_bnf_0212",
  
  "num_bnf_0301" ,"num_bnf_0302" ,"num_bnf_0303" ,"num_bnf_0304" ,
  "num_bnf_0306" ,"num_bnf_0307" ,"num_bnf_0309" ,"num_bnf_0310",
  
  "num_bnf_0401" ,"num_bnf_0402" ,"num_bnf_0403" ,"num_bnf_0404" ,
  "num_bnf_0405" ,"num_bnf_0406" ,"num_bnf_0407" ,"num_bnf_0408" ,
  "num_bnf_0409" ,"num_bnf_0410" ,"num_bnf_0411",
  
  "num_bnf_0501" ,"num_bnf_0502" ,"num_bnf_0503" ,"num_bnf_0504" ,
  "num_bnf_0505",
  
  "num_bnf_0601" ,"num_bnf_0602" ,"num_bnf_0603" ,"num_bnf_0604" ,
  "num_bnf_0605" ,"num_bnf_0605" ,"num_bnf_0607",
  
  "num_bnf_0704",
  
  "num_bnf_0801" ,"num_bnf_0802" ,"num_bnf_0803",
  
  "num_bnf_0901" ,"num_bnf_0902" ,"num_bnf_0904" ,"num_bnf_0905" ,
  "num_bnf_0906",
  
  "num_bnf_1001" ,"num_bnf_1002" ,"num_bnf_1003",
  
  "num_bnf_1103" ,"num_bnf_1104" ,"num_bnf_1105" ,"num_bnf_1106" ,
  "num_bnf_1107" ,"num_bnf_1108",
  
  "num_bnf_1201" ,"num_bnf_1202" ,"num_bnf_1203",
  
  "num_bnf_2002" ,"num_bnf_2003",
  
  "num_bnf_2102" ,"num_bnf_2107"
)


bnf_section_ref$groups = list()
bnf_section_ref$groups$nutrition = c(
  "num_bnf_0908","num_bnf_0909","num_bnf_0911","num_bnf_0912"
)
bnf_section_ref$groups$skin = c(
  "num_bnf_1301" ,"num_bnf_1302" ,"num_bnf_1303" ,"num_bnf_1304" ,
  "num_bnf_1305" ,"num_bnf_1308" ,"num_bnf_1310" ,"num_bnf_1311" ,
  "num_bnf_1313" ,"num_bnf_1314"
)
bnf_section_ref$groups$combined_X = c(
  "num_bnf_2201" ,"num_bnf_2202" ,"num_bnf_2205" ,"num_bnf_2210" ,
  "num_bnf_2215" ,"num_bnf_2220" ,"num_bnf_2230" ,"num_bnf_2240" ,
  "num_bnf_2250" ,"num_bnf_2260" ,"num_bnf_2270" ,"num_bnf_2280" ,
  "num_bnf_2285" ,"num_bnf_2290"
)
bnf_section_ref$groups$combined_Y = c(
  "num_bnf_2305" ,"num_bnf_2310" ,"num_bnf_2315" ,"num_bnf_2320" ,
  "num_bnf_2325" ,"num_bnf_2330" ,"num_bnf_2335" ,"num_bnf_2340" ,
  "num_bnf_2345" ,"num_bnf_2346" ,"num_bnf_2350" ,"num_bnf_2355" ,
  "num_bnf_2360" ,"num_bnf_2365" ,"num_bnf_2370" ,"num_bnf_2375" ,
  "num_bnf_2380" ,"num_bnf_2385" ,"num_bnf_2390" ,"num_bnf_2392" ,
  "num_bnf_2393" ,"num_bnf_2394" ,"num_bnf_2396" ,"num_bnf_2398"
)




#**********************************************************#  
# ICD10 groupings --------------------------------
#**********************************************************#  

get_icd10_grouping_drug_alcohol_selfharm=function() {
  
  list(
    alcohol = c(
      "E244", "E512", "F100", "F101", "F102", "F103", "F104", 
      "F105", "F106", "F107", "F108", "F109", "G312", "G621", "G721", 
      "I426", "K292", "K700", "K701", "K702", "K703", "K704", "K705", 
      "K706", "K707", "K708", "K709", "K860", "O354", "P043", "Q860", 
      "R780", "T510", "T511", "T519", "X450", "X451", "X452", "X453", 
      "X454", "X455", "X456", "X457", "X458", "X459", "X650", "X651", 
      "X652", "X653", "X654", "X655", "X656", "X657", "X658", "X659", 
      "Y150", "Y151", "Y152", "Y153", "Y154", "Y155", "Y156", "Y157", 
      "Y158", "Y159", "Y573", "Y900", "Y902", "Y903", "Y904", "Y905", 
      "Y906", "Y907", "Y908", "Y909", "Y910", "Y911", "Y912", "Y913", 
      "Y914", "Y915", "Y916", "Y917", "Y918", "Y919", "Z502", "Z714", 
      "Z721"), 
    drug = c(
      "F110", "F111", "F112", "F113", "F114", "F115", 
      "F116", "F117", "F118", "F119", "F120", "F121", "F122", "F123", 
      "F124", "F125", "F126", "F127", "F128", "F129", "F130", "F131", 
      "F132", "F133", "F134", "F135", "F136", "F137", "F138", "F139", 
      "F140", "F141", "F142", "F143", "F144", "F145", "F146", "F147", 
      "F148", "F149", "F150", "F151", "F152", "F153", "F154", "F155", 
      "F156", "F157", "F158", "F159", "F160", "F161", "F162", "F163", 
      "F164", "F165", "F166", "F167", "F168", "F169", "F180", "F181", 
      "F182", "F183", "F184", "F185", "F186", "F187", "F188", "F189", 
      "F190", "F191", "F192", "F193", "F194", "F195", "F196", "F197", 
      "F198", "F199")
  )
}









#**********************************************************#  
# Patient information --------------------------------
#**********************************************************#  

#' transformer_patient_info
#' Returns basic information about each patient as a feature
#' That includes age, gender and SIMD_DECILE_2016_SCT
#' The latter is a categorised deprivation status
#' 
transformer_patient_info <- function(
  # Input data
  patients, episodes, list_of_data_tables,
  time_cutoff # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  
  # Further optional parameters (these have defaults)
  
){
  
  
  # Create the out_matrix
  out_matrix = patients %>% 
    transmute(
      id = id,
      gender = gender,
      age = as.integer((as.integer(as_date(as_datetime(time_cutoff))) - as.integer(as_date(date_of_birth)))/365),
      SIMD_DECILE_2016_SCT = SIMD_DECILE_2016_SCT
    ) #%>%
  #mutate(age = ifelse(age<0, 0, age))
  
  out_missing_default = as.list(rep(NA, length=length(names(out_matrix)[-1])))
  names(out_missing_default) <- names(out_matrix)[-1]
  
  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
  
}
