## Working directory
#setwd("/mnt/batch/tasks/shared/LS_root/mounts/clusters/sparratest/code/Users/james.liley/")

##' Generate mockups of EHR tables used in SPARRA score. Intended for testing of procedures.
##' 
##' @name ehr_mockup
##' @param ehr_directory file location to save tables
##' @param n size of each table
##' @param nid number of unique IDs
##' @param seed random seed
##' @returns NULL
ehr_mockup=function(ehr_directory="./",n=1000,nid=5000,seed=220) {

set.seed(seed)  
  
ehr_directory=paste0(ehr_directory,"/")  

## General rows
id_dob=as.Date.numeric(sample(80*365,nid,rep=T),origin="1930-01-01")
id_s=sample(c(1,2,0,9),nid,rep=T,prob=c(20,20,1,1))
id_simd=factor(sample(c("1=most deprived",2:9,"10=least deprived"),nid,rep=T),
                  levels=c("1=most deprived",2:9,"10=least deprived"))


## Lookup

UNIQUE_STUDY_ID=1:nid
GENDER=id_s
DOB=id_dob
SIMD_QUINTILE_2016_SCT=ceiling(as.numeric(id_simd)/2)
SIMD_DECILE_2016_SCT=as.numeric(id_simd)

dx=data.frame(UNIQUE_STUDY_ID,GENDER,DOB,SIMD_QUINTILE_2016_SCT,SIMD_DECILE_2016_SCT,
              stringsAsFactors=FALSE)


for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=format(dx[,i],"%d-%m-%Y")
haven::write_sav(dx,path=paste0(ehr_directory,"LOOKUP_SIMD_sex_dob_for_ATI.zsav"),compress=TRUE)




## AE2
TWELVE_MONTH_TIME_PERIOD=paste0("y",sample(4,n,rep=T))
ADMISSION_DATE=as.Date.numeric(sample(7*365,n,rep=T),origin="2013-01-01")
ADMISSION_TIME=paste0(sample(0:23,n,rep=T),":",sample(0:59,n,rep=T))
TRANSFER_DISCHARGE_DATE=as.Date(ADMISSION_DATE+sample(10,n,rep=T))
TRANSFER_DISCHARGE_TIME=paste0(sample(0:23,n,rep=T),":",sample(0:59,n,rep=T))
DIAGNOSIS_1=sample(99,n,rep=T)
DIAGNOSIS_2=sample(99,n,rep=T)
DIAGNOSIS_3=rep("",n)
DISCHARGE_DESTINATION=paste0("0",sample(9,n,rep=T),toupper(sample(letters,n,rep=T)))
LOCATION_CODE=paste0(toupper(sample(letters,n,rep=T)),sample(100:999,n,rep=T),toupper(sample(letters,n,rep=T)))
REFERRAL_SOURCE=paste0("0",sample(9,n,rep=T),toupper(sample(letters,n,rep=T)))
UNIQUE_STUDY_ID=sample(nid,n,rep=F)


dx=data.frame(
  TWELVE_MONTH_TIME_PERIOD,
  ADMISSION_DATE,
  ADMISSION_TIME,
  TRANSFER_DISCHARGE_DATE,
  TRANSFER_DISCHARGE_TIME,
  DIAGNOSIS_1,
  DIAGNOSIS_2,
  DIAGNOSIS_3,
  DISCHARGE_DESTINATION,
  LOCATION_CODE,
  REFERRAL_SOURCE,
  UNIQUE_STUDY_ID,
  stringsAsFactors=F
)

for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=format(dx[,i],"%d-%m-%Y")
haven::write_sav(dx,path=paste0(ehr_directory,"Final_AE2_extract_incl_UniqueStudyID.zsav"),compress=TRUE)


## Deaths
nd=sample(nid,100) # Number of recorded deaths
DATE_OF_DEATH=as.Date.numeric(sample(7*365,length(nd),rep=T),origin="2011-01-01")
PRIMARY_CAUSE_OF_DEATH=paste0(sample(letters,length(nd),rep=T),sample(800,length(nd),rep=T))
UNIQUE_STUDY_ID=nd

dx=data.frame(DATE_OF_DEATH,PRIMARY_CAUSE_OF_DEATH,UNIQUE_STUDY_ID,
              stringsAsFactors=FALSE)


for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=format(dx[,i],"%d-%b-%y")
haven::write_sav(dx,path=paste0(ehr_directory,"Final_deaths_extract_incl_UniqueStudyID_v2.zsav"),compress=TRUE)


## SMR00
id=sample(nid,n)

ATTENDANCE_FOLLOW_UP=sample(c("",1:5,8),n,replace=T)
CLINIC_ATTENDANCE=sample(c(1,5,8),n,replace=T)
CLINIC_DATE=as.Date.numeric(sample(7*365,n,rep=T),origin="2011-01-01")
DOB=id_dob[id]
DATE_OF_MAIN_OPERATION=as.Date(rep(NA,n)); s0=sample(n,round(0.3*n))
DATE_OF_MAIN_OPERATION[s0]=as.Date.numeric(sample(7*365,length(s0),rep=T),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_1=as.Date(rep(NA,n)); s1=sample(s0,round(length(s0)*0.3))
DATE_OF_OTHER_OPERATION_1[s1]=as.Date.numeric(sample(7*365,length(s1)),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_2=as.Date(rep(NA,n)); s2=sample(s1,round(length(s1)*0.3),rep=T)
DATE_OF_OTHER_OPERATION_2[s2]=as.Date.numeric(sample(7*365,length(s2),rep=T),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_3=as.Date(rep(NA,n)); s3=sample(s2,round(length(s2)*0.3))
DATE_OF_OTHER_OPERATION_3[s3]=as.Date.numeric(sample(7*365,length(s3),rep=T),origin="2011-01-01")
GP_PRACTICE_CODE=sample(10000,n,rep=T)
HBRES_CURRENT_DATE=paste0("S080000",sample(100,n,rep=T))
HBTREAT_CURRENT_DATE=paste0("S080000",sample(100,n,rep=T))
LOCATION=toupper(paste0(sample(letters,n,rep=T),sample(100,n,rep=T),sample(letters,n,rep=T)))
MAIN_CONDITION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(MAIN_CONDITION=="A"); MAIN_CONDITION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
OTHER_CONDITION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); OTHER_CONDITION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
OTHER_CONDITION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); OTHER_CONDITION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
OTHER_CONDITION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); OTHER_CONDITION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
OTHER_CONDITION_4=rep("",n); s4=sample(s3,round(length(s3)*0.5)); OTHER_CONDITION_1[s4]=toupper(paste0(sample(letters,length(s4),rep=T),sample(999,length(s4),rep=T)))
OTHER_CONDITION_5=rep("",n); s5=sample(s4,round(length(s4)*0.5)); OTHER_CONDITION_1[s5]=toupper(paste0(sample(letters,length(s5),rep=T),sample(999,length(s5),rep=T)))
MAIN_OPERATION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(MAIN_OPERATION=="A"); MAIN_OPERATION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
MODE_OF_CONTACT=sample(c("",1:3),n,rep=T)
OTHER_OPERATION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); OTHER_OPERATION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
OTHER_OPERATION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); OTHER_OPERATION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
OTHER_OPERATION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); OTHER_OPERATION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
PATIENT_CATEGORY=sample(c(2,3,4,5),n,rep=T)
REFERRAL_TYPE=sample(c(1:3),n,rep=T)
SEX=id_s[id]
SIGNIFICANT_FACILITY=sample(c("",11,31:39),n,rep=T)
SPECIALTY=toupper(paste0(sample(letters,n,rep=T),sample(99,n,rep=T)))
TWELVE_MONTH_TIME_PERIOD=paste0("y",sample(4,n,rep=T))
simd2016_sc_decile=id_simd[id]
simd2016_sc_quintile=ceiling(as.numeric(id_simd[id])/2)
simd2016_HB2014_decile=simd2016_sc_decile
simd2016_HB2014_quintile=simd2016_sc_quintile
simd2016tp15=as.numeric(simd2016_sc_decile)>8
simd2016bt15=as.numeric(simd2016_sc_decile)<2
UNIQUE_STUDY_ID=id

dx=data.frame(ATTENDANCE_FOLLOW_UP, 
                CLINIC_ATTENDANCE, 
                CLINIC_DATE, 
                DOB, 
                DATE_OF_MAIN_OPERATION, 
                DATE_OF_OTHER_OPERATION_1, 
                DATE_OF_OTHER_OPERATION_2, 
                DATE_OF_OTHER_OPERATION_3, 
                GP_PRACTICE_CODE, 
                HBRES_CURRENT_DATE, 
                HBTREAT_CURRENT_DATE, 
                LOCATION, 
                MAIN_CONDITION, 
                OTHER_CONDITION_1, 
                OTHER_CONDITION_2, 
                OTHER_CONDITION_3, 
                OTHER_CONDITION_4, 
                OTHER_CONDITION_5, 
                MAIN_OPERATION, 
                MODE_OF_CONTACT,
                OTHER_OPERATION_1, 
                OTHER_OPERATION_2, 
                OTHER_OPERATION_3, 
                PATIENT_CATEGORY, 
                REFERRAL_TYPE, 
                SEX, 
                SIGNIFICANT_FACILITY, 
                SPECIALTY, 
                TWELVE_MONTH_TIME_PERIOD, 
                simd2016_sc_decile, 
                simd2016_sc_quintile, 
                simd2016_HB2014_decile, 
                simd2016_HB2014_quintile, 
                simd2016tp15, 
                simd2016bt15, 
                UNIQUE_STUDY_ID,
              stringsAsFactors=FALSE)

  
for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=as.POSIXct(dx[,i]) #format(dx[,i],"%d-%m-%Y")
haven::write_sav(dx,path=paste0(ehr_directory,"Final_SMR00_extract_incl_UniqueStudyID.zsav"),compress=TRUE)



## SMR01
id=sample(nid,n)

DOB=id_dob[id]
ADMISSION_DATE=as.Date.numeric(sample(7*365,n,rep=T),origin="2011-01-01")
LENGTH_OF_STAY=sample(10,n,rep=T)
DISCHARGE_DATE=as.Date(ADMISSION_DATE+LENGTH_OF_STAY)
CIS_MARKER=sample(1500,n,rep=T)
HB_OF_RESIDENCE_NUMBER=sample(c("",1:15),n,rep=T)
MANAGEMENT_OF_PATIENT=sample(c("",1:8),n,rep=T)
CLINICAL_FACILITY_START=sample(c("","1G"),n,rep=T)
CLINICAL_FACILITY_END=rep("",n)
ADMISSION_TYPE=sample(c(36,11,33,30,18,31,10,35,38,32,12,39,20,22,34,19,21),n,rep=T)
ADMISSION_REASON=sample(c("11","","1C","13","1E","48","1M","10","15","1K","16","40","41","14","1B","1H","49","4A","1D","12","1F","1A","18","1J","42","19","43","47","44","48","17"),n,rep=T)
ADMISSION_TRANSFER_FROM_LOC=toupper(paste0(sample(letters,n,rep=T),sample(999,n,rep=T),sample(letters,n,rep=T)))
READY_FOR_DISCHARGE_DATE=as.Date(DISCHARGE_DATE-sample(2,n,rep=T))
DISCHARGE_TYPE=sample(c("",10:13,18:22,28:29,40:43),n,rep=T)
DISCHARGE_TRANSFER_TO=toupper(sample(c("","0",10:69,"DI",paste0("4",letters[1:8]),paste0("5",letters[1:8])),n,rep=T))
DISCHARGE_TRANSFER_TO_LOCATION=toupper(paste0(sample(letters,n,rep=T),sample(999,n,rep=T),sample(letters,n,rep=T)))
ADMISSION_TRANSFER_FROM=toupper(paste0(sample(letters,n,rep=T),sample(999,n,rep=T),sample(letters,n,rep=T)))
CHP=paste0("S030000",sample(100,n,rep=T))
INPATIENT_DAYCASE_IDENTIFIER=sample(c("I","D",""),n,rep=T)
HRG42=toupper(paste0(sample(letters,n,rep=T),sample(letters,n,rep=T),sample(99,n,rep=T),sample(letters,n,rep=T)))
DATE_OF_MAIN_OPERATION=as.Date(rep(NA,n)); s0=sample(n,round(0.3*n))
DATE_OF_MAIN_OPERATION[s0]=as.Date.numeric(sample(7*365,length(s0),rep=T),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_1=as.Date(rep(NA,n)); s1=sample(s0,round(length(s0)*0.3))
DATE_OF_OTHER_OPERATION_1[s1]=as.Date.numeric(sample(7*365,length(s1)),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_2=as.Date(rep(NA,n)); s2=sample(s1,round(length(s1)*0.3),rep=T)
DATE_OF_OTHER_OPERATION_2[s2]=as.Date.numeric(sample(7*365,length(s2),rep=T),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_3=as.Date(rep(NA,n)); s3=sample(s2,round(length(s2)*0.3))
DATE_OF_OTHER_OPERATION_3[s3]=as.Date.numeric(sample(7*365,length(s3),rep=T),origin="2011-01-01")
GP_PRACTICE_CODE=sample(10000,n,replace=T)
HBRES_CURRENT_DATE=paste0("S080000",sample(100,n,rep=T))
HBTREAT_CURRENT_DATE=paste0("S080000",sample(100,n,rep=T))
LOCATION=toupper(paste0(sample(letters,n,rep=T),sample(100,n,rep=T),sample(letters,n,rep=T)))
MAIN_CONDITION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(MAIN_CONDITION=="A"); MAIN_CONDITION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
OTHER_CONDITION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); OTHER_CONDITION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
OTHER_CONDITION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); OTHER_CONDITION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
OTHER_CONDITION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); OTHER_CONDITION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
OTHER_CONDITION_4=rep("",n); s4=sample(s3,round(length(s3)*0.5)); OTHER_CONDITION_1[s4]=toupper(paste0(sample(letters,length(s4),rep=T),sample(999,length(s4),rep=T)))
OTHER_CONDITION_5=rep("",n); s5=sample(s4,round(length(s4)*0.5)); OTHER_CONDITION_1[s5]=toupper(paste0(sample(letters,length(s5),rep=T),sample(999,length(s5),rep=T)))
MAIN_OPERATION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(MAIN_OPERATION=="A"); MAIN_OPERATION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
OTHER_OPERATION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); OTHER_OPERATION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
OTHER_OPERATION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); OTHER_OPERATION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
OTHER_OPERATION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); OTHER_OPERATION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
SEX=id_s[id]
SIGNIFICANT_FACILITY=sample(c("",11,31:39),n,rep=T)
SPECIALTY=toupper(paste0(sample(letters,n,rep=T),sample(99,n,rep=T)))
simd2016_sc_decile=id_simd[id]
simd2016_sc_quintile=ceiling(as.numeric(id_simd[id])/2)
simd2016_HB2014_decile=simd2016_sc_decile
simd2016_HB2014_quintile=simd2016_sc_quintile
simd2016tp15=as.numeric(simd2016_sc_decile)>8
simd2016bt15=as.numeric(simd2016_sc_decile)<2
UNIQUE_STUDY_ID=id

dx=data.frame(DOB, 
              ADMISSION_DATE, 
              LENGTH_OF_STAY, 
              DISCHARGE_DATE, 
              CIS_MARKER, 
              HB_OF_RESIDENCE_NUMBER, 
              MANAGEMENT_OF_PATIENT, 
              CLINICAL_FACILITY_START, 
              CLINICAL_FACILITY_END, 
              ADMISSION_TYPE, 
              ADMISSION_REASON, 
              ADMISSION_TRANSFER_FROM_LOC, 
              READY_FOR_DISCHARGE_DATE, 
              DISCHARGE_TYPE, 
              DISCHARGE_TRANSFER_TO, 
              DISCHARGE_TRANSFER_TO_LOCATION, 
              ADMISSION_TRANSFER_FROM, 
              CHP, 
              INPATIENT_DAYCASE_IDENTIFIER, 
              HRG42, 
              DATE_OF_MAIN_OPERATION, 
              DATE_OF_OTHER_OPERATION_1, 
              DATE_OF_OTHER_OPERATION_2, 
              DATE_OF_OTHER_OPERATION_3, 
              GP_PRACTICE_CODE, 
              HBRES_CURRENT_DATE, 
              HBTREAT_CURRENT_DATE, 
              LOCATION, 
              MAIN_CONDITION, 
              OTHER_CONDITION_1, 
              OTHER_CONDITION_2, 
              OTHER_CONDITION_3, 
              OTHER_CONDITION_4, 
              OTHER_CONDITION_5, 
              MAIN_OPERATION, 
              OTHER_OPERATION_1, 
              OTHER_OPERATION_2, 
              OTHER_OPERATION_3, 
              SEX, 
              SIGNIFICANT_FACILITY, 
              SPECIALTY, 
              simd2016_sc_decile, 
              simd2016_sc_quintile, 
              simd2016_HB2014_decile, 
              simd2016_HB2014_quintile, 
              simd2016tp15, 
              simd2016bt15, 
              UNIQUE_STUDY_ID,
              stringsAsFactors=FALSE)


for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=as.POSIXct(dx[,i]) #format(dx[,i],"%d-%m-%Y")
haven::write_sav(dx,path=paste0(ehr_directory,"Final_SMR01_extract_incl_UniqueStudyID.zsav"),compress=TRUE)


## SMR01E

id=sample(nid,n)
CARE_PACKAGE_IDENTIFIER=sample(c("",1:3000),n,rep=T)
DISCHARGE=sample(0:2,n,rep=T)
EPISODE_RECORD_KEY=sample(10^10,n,rep=F)
GLS_CIS_MARKER=sample(200,n,rep=T)
PATIENT_CATEGORY=sample(c(2,3,4,5),n,rep=T)
URBAN_RURAL_CODE=sample(c("",1:8),n,rep=T)
TWELVE_MONTH_TIME_PERIOD=paste0("y",sample(4,n,rep=T))
DOB=id_dob[id]
ADMISSION_DATE=as.Date.numeric(sample(7*365,n,rep=T),origin="2011-01-01")
LENGTH_OF_STAY=sample(10,n,rep=T)
DISCHARGE_DATE=as.Date(ADMISSION_DATE+LENGTH_OF_STAY)
CIS_MARKER=sample(1500,n,rep=T)
HB_OF_RESIDENCE_NUMBER=sample(c("",1:15),n,rep=T)
MANAGEMENT_OF_PATIENT=sample(c("",1:8),n,rep=T)
CLINICAL_FACILITY_START=sample(c("","1G"),n,rep=T)
CLINICAL_FACILITY_END=rep("",n)
DAYS_WAITING=sample(c(NA,1:200),n,rep=T)
ADMISSION_TYPE=sample(c(36,11,33,30,18,31,10,35,38,32,12,39,20,22,34,19,21),n,rep=T)
ADMISSION_REASON=sample(c("11","","1C","13","1E","48","1M","10","15","1K","16","40","41","14","1B","1H","49","4A","1D","12","1F","1A","18","1J","42","19","43","47","44","48","17"),n,rep=T)
ADMISSION_TRANSFER_FROM_LOC=toupper(paste0(sample(letters,n,rep=T),sample(999,n,rep=T),sample(letters,n,rep=T)))
READY_FOR_DISCHARGE_DATE=as.Date(DISCHARGE_DATE-sample(2,n,rep=T))
DISCHARGE_TYPE=sample(c("",10:13,18:22,28:29,40:43),n,rep=T)
DISCHARGE_TRANSFER_TO=toupper(sample(c("","0",10:69,"DI",paste0("4",letters[1:8]),paste0("5",letters[1:8])),n,rep=T))
DISCHARGE_TRANSFER_TO_LOCATION=toupper(paste0(sample(letters,n,rep=T),sample(999,n,rep=T),sample(letters,n,rep=T)))
ADMISSION_TRANSFER_FROM=toupper(paste0(sample(letters,n,rep=T),sample(999,n,rep=T),sample(letters,n,rep=T)))
CHP=paste0("S030000",sample(100,n,rep=T))
INPATIENT_DAYCASE_IDENTIFIER=sample(c("I","D",""),n,rep=T)
HRG42=toupper(paste0(sample(letters,n,rep=T),sample(letters,n,rep=T),sample(99,n,rep=T),sample(letters,n,rep=T)))
DATE_OF_MAIN_OPERATION=as.Date(rep(NA,n)); s0=sample(n,round(0.3*n))
DATE_OF_MAIN_OPERATION[s0]=as.Date.numeric(sample(7*365,length(s0),rep=T),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_1=as.Date(rep(NA,n)); s1=sample(s0,round(length(s0)*0.3))
DATE_OF_OTHER_OPERATION_1[s1]=as.Date.numeric(sample(7*365,length(s1)),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_2=as.Date(rep(NA,n)); s2=sample(s1,round(length(s1)*0.3),rep=T)
DATE_OF_OTHER_OPERATION_2[s2]=as.Date.numeric(sample(7*365,length(s2),rep=T),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_3=as.Date(rep(NA,n)); s3=sample(s2,round(length(s2)*0.3))
DATE_OF_OTHER_OPERATION_3[s3]=as.Date.numeric(sample(7*365,length(s3),rep=T),origin="2011-01-01")
GP_PRACTICE_CODE=sample(10000,n,rep=T)
HBRES_CURRENT_DATE=paste0("S080000",sample(100,n,rep=T))
HBTREAT_CURRENT_DATE=paste0("S080000",sample(100,n,rep=T))
LOCATION=toupper(paste0(sample(letters,n,rep=T),sample(100,n,rep=T),sample(letters,n,rep=T)))
MAIN_CONDITION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(MAIN_CONDITION=="A"); MAIN_CONDITION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
OTHER_CONDITION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); OTHER_CONDITION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
OTHER_CONDITION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); OTHER_CONDITION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
OTHER_CONDITION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); OTHER_CONDITION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
OTHER_CONDITION_4=rep("",n); s4=sample(s3,round(length(s3)*0.5)); OTHER_CONDITION_1[s4]=toupper(paste0(sample(letters,length(s4),rep=T),sample(999,length(s4),rep=T)))
OTHER_CONDITION_5=rep("",n); s5=sample(s4,round(length(s4)*0.5)); OTHER_CONDITION_1[s5]=toupper(paste0(sample(letters,length(s5),rep=T),sample(999,length(s5),rep=T)))
MAIN_OPERATION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(MAIN_OPERATION=="A"); MAIN_OPERATION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
OTHER_OPERATION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); OTHER_OPERATION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
OTHER_OPERATION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); OTHER_OPERATION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
OTHER_OPERATION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); OTHER_OPERATION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
SEX=id_s[id]
SIGNIFICANT_FACILITY=sample(c("",11,31:39),n,rep=T)
SPECIALTY=toupper(paste0(sample(letters,n,rep=T),sample(99,n,rep=T)))
simd2016_sc_decile=id_simd[id]
simd2016_sc_quintile=ceiling(as.numeric(id_simd[id])/2)
simd2016_HB2014_decile=simd2016_sc_decile
simd2016_HB2014_quintile=simd2016_sc_quintile
simd2016tp15=as.numeric(simd2016_sc_decile)>8
simd2016bt15=as.numeric(simd2016_sc_decile)<2
UNIQUE_STUDY_ID=id

dx=data.frame(CARE_PACKAGE_IDENTIFIER, 
              DISCHARGE, 
              EPISODE_RECORD_KEY, 
              GLS_CIS_MARKER, 
              PATIENT_CATEGORY, 
              URBAN_RURAL_CODE, 
              TWELVE_MONTH_TIME_PERIOD, 
              DOB, 
              ADMISSION_DATE, 
              LENGTH_OF_STAY, 
              DISCHARGE_DATE, 
              CIS_MARKER, 
              HB_OF_RESIDENCE_NUMBER, 
              MANAGEMENT_OF_PATIENT, 
              CLINICAL_FACILITY_START, 
              CLINICAL_FACILITY_END, 
              DAYS_WAITING, 
              ADMISSION_TYPE, 
              ADMISSION_REASON, 
              ADMISSION_TRANSFER_FROM_LOC, 
              READY_FOR_DISCHARGE_DATE, 
              DISCHARGE_TYPE, 
              DISCHARGE_TRANSFER_TO, 
              DISCHARGE_TRANSFER_TO_LOCATION, 
              ADMISSION_TRANSFER_FROM, 
              CHP, 
              INPATIENT_DAYCASE_IDENTIFIER, 
              HRG42, 
              DATE_OF_MAIN_OPERATION, 
              DATE_OF_OTHER_OPERATION_1, 
              DATE_OF_OTHER_OPERATION_2, 
              DATE_OF_OTHER_OPERATION_3, 
              GP_PRACTICE_CODE, 
              HBRES_CURRENT_DATE, 
              HBTREAT_CURRENT_DATE, 
              LOCATION, 
              MAIN_CONDITION, 
              OTHER_CONDITION_1, 
              OTHER_CONDITION_2, 
              OTHER_CONDITION_3, 
              OTHER_CONDITION_4, 
              OTHER_CONDITION_5, 
              MAIN_OPERATION, 
              OTHER_OPERATION_1, 
              OTHER_OPERATION_2, 
              OTHER_OPERATION_3, 
              SEX, 
              SIGNIFICANT_FACILITY, 
              SPECIALTY, 
              simd2016_sc_decile, 
              simd2016_sc_quintile, 
              simd2016_HB2014_decile, 
              simd2016_HB2014_quintile, 
              simd2016tp15, 
              simd2016bt15, 
              UNIQUE_STUDY_ID,
              stringsAsFactors=FALSE)

for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=as.POSIXct(dx[,i]) #format(dx[,i],"%d-%m-%Y")
haven::write_sav(dx,path=paste0(ehr_directory,"Final_SMR01E_extract_incl_UniqueStudyID.zsav"),compress=TRUE)




## SMR04

id=sample(nid,n)
ADMISSION_MAIN_CONDITION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(ADMISSION_MAIN_CONDITION=="A"); ADMISSION_MAIN_CONDITION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
ADMISSION_OTHER_CONDITION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); ADMISSION_OTHER_CONDITION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
ADMISSION_OTHER_CONDITION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); ADMISSION_OTHER_CONDITION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
ADMISSION_OTHER_CONDITION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); ADMISSION_OTHER_CONDITION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
ADMISSION_OTHER_CONDITION_4=rep("",n); s4=sample(s3,round(length(s3)*0.5)); ADMISSION_OTHER_CONDITION_1[s4]=toupper(paste0(sample(letters,length(s4),rep=T),sample(999,length(s4),rep=T)))
ADMISSION_REFERRAL_FROM=sample(c(1:9,toupper(letters[1:7])),n,rep=T)
ARRANGEMENTS_FOR_AFTERCARE_1=sample(c("",1:9,toupper(letters[1:3])),n,rep=T)
ARRANGEMENTS_FOR_AFTERCARE_2=sample(c("",1:9,toupper(letters[1:3])),prob=c(10,rep(1,12)),n,rep=T)
ARRANGEMENTS_FOR_AFTERCARE_3=sample(c("",1:9,toupper(letters[1:3])),prob=c(10,rep(1,12)),n,rep=T)
ARRANGEMENTS_FOR_AFTERCARE_4=sample(c("",1:9,toupper(letters[1:3])),prob=c(10,rep(1,12)),n,rep=T)
CARE_PLAN_ARRANGEMENTS=sample(c("",1:2),n,rep=T)
DAYS_WAITING=sample(c(NA,1:200),n,rep=T)
EPISODE_RECORD_KEY=sample(10^10,n,rep=F)
TWELVE_MONTH_TIME_PERIOD=paste0("y",sample(4,n,rep=T))
DOB=id_dob[id]
ADMISSION_DATE=as.Date.numeric(sample(7*365,n,rep=T),origin="2011-01-01")
ECT_1ST_TREATMENT_DATE=as.Date(ADMISSION_DATE + sample(20,n,rep=T)) 
ECT_TREATMENTS_NUM_THIS_EPISODE=sample(0:9,n,rep=T)
PREVIOUS_PSYCHIATRIC_CARE=sample(c(1:3,9),n,rep=T)
LENGTH_OF_STAY=sample(10,n,rep=T)
DISCHARGE_DATE=as.Date(ADMISSION_DATE+LENGTH_OF_STAY)
CIS_MARKER=sample(1500,n,rep=T)
HB_OF_RESIDENCE_NUMBER=sample(c("",1:15),n,rep=T)
MANAGEMENT_OF_PATIENT=sample(c("",1:8),n,rep=T)
CLINICAL_FACILITY_START=sample(c("","1G"),n,rep=T)
CLINICAL_FACILITY_END=rep("",n)
DAYS_WAITING=sample(c(NA,1:200),n,rep=T)
ADMISSION_TYPE=sample(c(36,11,33,30,18,31,10,35,38,32,12,39,20,22,34,19,21),n,rep=T)
ADMISSION_REASON=sample(c("11","","1C","13","1E","48","1M","10","15","1K","16","40","41","14","1B","1H","49","4A","1D","12","1F","1A","18","1J","42","19","43","47","44","48","17"),n,rep=T)
READY_FOR_DISCHARGE_DATE=as.Date(DISCHARGE_DATE-sample(2,n,rep=T))
DISCHARGE_TYPE=sample(c("",10:13,18:22,28:29,40:43),n,rep=T)
DISCHARGE_TRANSFER_TO=toupper(sample(c("","0",10:69,"DI",paste0("4",letters[1:8]),paste0("5",letters[1:8])),n,rep=T))
ADMISSION_TRANSFER_FROM=toupper(paste0(sample(letters,n,rep=T),sample(999,n,rep=T),sample(letters,n,rep=T)))
CHP=paste0("S030000",sample(100,n,rep=T))
HRG42=toupper(paste0(sample(letters,n,rep=T),sample(letters,n,rep=T),sample(99,n,rep=T),sample(letters,n,rep=T)))
DATE_OF_MAIN_OPERATION=as.Date(rep(NA,n)); s0=sample(n,round(0.3*n))
DATE_OF_MAIN_OPERATION[s0]=as.Date.numeric(sample(7*365,length(s0),rep=T),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_1=as.Date(rep(NA,n)); s1=sample(s0,round(length(s0)*0.3))
DATE_OF_OTHER_OPERATION_1[s1]=as.Date.numeric(sample(7*365,length(s1)),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_2=as.Date(rep(NA,n)); s2=sample(s1,round(length(s1)*0.3),rep=T)
DATE_OF_OTHER_OPERATION_2[s2]=as.Date.numeric(sample(7*365,length(s2),rep=T),origin="2011-01-01")
DATE_OF_OTHER_OPERATION_3=as.Date(rep(NA,n)); s3=sample(s2,round(length(s2)*0.3))
DATE_OF_OTHER_OPERATION_3[s3]=as.Date.numeric(sample(7*365,length(s3),rep=T),origin="2011-01-01")
GP_PRACTICE_CODE=sample(10000,n,rep=T)
HBRES_CURRENT_DATE=paste0("S080000",sample(100,n,rep=T))
HBTREAT_CURRENT_DATE=paste0("S080000",sample(100,n,rep=T))
LOCATION=toupper(paste0(sample(letters,n,rep=T),sample(100,n,rep=T),sample(letters,n,rep=T)))
MAIN_CONDITION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(MAIN_CONDITION=="A"); MAIN_CONDITION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
OTHER_CONDITION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); OTHER_CONDITION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
OTHER_CONDITION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); OTHER_CONDITION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
OTHER_CONDITION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); OTHER_CONDITION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
OTHER_CONDITION_4=rep("",n); s4=sample(s3,round(length(s3)*0.5)); OTHER_CONDITION_1[s4]=toupper(paste0(sample(letters,length(s4),rep=T),sample(999,length(s4),rep=T)))
OTHER_CONDITION_5=rep("",n); s5=sample(s4,round(length(s4)*0.5)); OTHER_CONDITION_1[s5]=toupper(paste0(sample(letters,length(s5),rep=T),sample(999,length(s5),rep=T)))
MAIN_OPERATION=sample(c("","A"),n,rep=T,prob=c(1,20)); s=which(MAIN_OPERATION=="A"); MAIN_OPERATION[s]=toupper(paste0(sample(letters,length(s),rep=T),sample(999,length(s),rep=T)))
OTHER_OPERATION_1=rep("",n); s1=sample(s,round(length(s)*0.5)); OTHER_OPERATION_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
OTHER_OPERATION_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); OTHER_OPERATION_1[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(999,length(s2),rep=T)))
OTHER_OPERATION_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); OTHER_OPERATION_1[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(999,length(s3),rep=T)))
SEX=id_s[id]
SIGNIFICANT_FACILITY=sample(c("",11,31:39),n,rep=T)
SPECIALTY=toupper(paste0(sample(letters,n,rep=T),sample(99,n,rep=T)))
simd2016_sc_decile=id_simd[id]
simd2016_sc_quintile=ceiling(as.numeric(id_simd[id])/2)
simd2016_HB2014_decile=simd2016_sc_decile
simd2016_HB2014_quintile=simd2016_sc_quintile
simd2016tp15=as.numeric(simd2016_sc_decile)>8
simd2016bt15=as.numeric(simd2016_sc_decile)<2
UNIQUE_STUDY_ID=id

dx=data.frame(ADMISSION_MAIN_CONDITION, 
              ADMISSION_OTHER_CONDITION_1, 
              ADMISSION_OTHER_CONDITION_2, 
              ADMISSION_OTHER_CONDITION_3, 
              ADMISSION_OTHER_CONDITION_4, 
              ADMISSION_REFERRAL_FROM, 
              ARRANGEMENTS_FOR_AFTERCARE_1, 
              ARRANGEMENTS_FOR_AFTERCARE_2, 
              ARRANGEMENTS_FOR_AFTERCARE_3, 
              ARRANGEMENTS_FOR_AFTERCARE_4, 
              CARE_PLAN_ARRANGEMENTS, 
              DAYS_WAITING, 
              EPISODE_RECORD_KEY, 
              TWELVE_MONTH_TIME_PERIOD, 
              DOB, 
              ADMISSION_DATE, 
              ECT_1ST_TREATMENT_DATE, 
              ECT_TREATMENTS_NUM_THIS_EPISODE, 
              PREVIOUS_PSYCHIATRIC_CARE, 
              LENGTH_OF_STAY, 
              DISCHARGE_DATE, 
              CIS_MARKER, 
              HB_OF_RESIDENCE_NUMBER, 
              MANAGEMENT_OF_PATIENT, 
              CLINICAL_FACILITY_START, 
              CLINICAL_FACILITY_END, 
              DAYS_WAITING, 
              ADMISSION_TYPE, 
              ADMISSION_REASON, 
              READY_FOR_DISCHARGE_DATE, 
              DISCHARGE_TYPE, 
              DISCHARGE_TRANSFER_TO, 
              ADMISSION_TRANSFER_FROM, 
              CHP, 
              HRG42, 
              DATE_OF_MAIN_OPERATION, 
              DATE_OF_OTHER_OPERATION_1, 
              DATE_OF_OTHER_OPERATION_2, 
              DATE_OF_OTHER_OPERATION_3, 
              GP_PRACTICE_CODE, 
              HBRES_CURRENT_DATE, 
              HBTREAT_CURRENT_DATE, 
              LOCATION, 
              MAIN_CONDITION, 
              OTHER_CONDITION_1, 
              OTHER_CONDITION_2, 
              OTHER_CONDITION_3, 
              OTHER_CONDITION_4, 
              OTHER_CONDITION_5, 
              MAIN_OPERATION, 
              OTHER_OPERATION_1, 
              OTHER_OPERATION_2, 
              OTHER_OPERATION_3, 
              SEX, 
              SIGNIFICANT_FACILITY, 
              SPECIALTY, 
              simd2016_sc_decile, 
              simd2016_sc_quintile, 
              simd2016_HB2014_decile, 
              simd2016_HB2014_quintile, 
              simd2016tp15, 
              simd2016bt15, 
              UNIQUE_STUDY_ID,
              stringsAsFactors=FALSE)

for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=as.POSIXct(dx[,i]) #format(dx[,i],"%d-%m-%Y")
haven::write_sav(dx,path=paste0(ehr_directory,"Final_SMR04_extract_v2_for_ATI_incl_prev_psych_care.sav"),compress=FALSE)




## PIS
nx=c(paste0("0",c(101:109,201:209,210:213,301:311,401:411,501:505,601:607,701:704,801:803,901:912)),
     c(1001:1003,1103:1109,1201:1203,1301:1315,1403:1405,1501:1502,1801,1803,1901:1914,2001:2018,2020,
       2101:2149,2201:2202,2205,2210,2215,2220,2230,2240,2250,2260,2270,2280,2285,2290,2305,2310,2315,
       2320,2325,2330,2335,2340,2345,2346,2350,2355,2360,2365,2370,2375,2380,2385,2390,2392,2393,2394,
       2396,2398))
cnames=paste0("NUM_BNF_",nx)


fnames=c("PIS_201305_to_201404_for_ATI_monthly_ID.zsav", 
         "PIS_201405_to_201504_for_ATI_monthly_ID.zsav", 
         "PIS_201505_to_201604_for_ATI_monthly_ID.zsav", 
         "PIS_201605_to_201704_for_ATI_monthly_ID.zsav", 
         "PIS_201705_to_201804_for_ATI_monthly_ID.zsav")

for (xyear in 1:5) {
  
pmat=c()
for (i in 1:length(cnames)) {
  presc=sample(c(0,1:20),n,rep=T,prob=c(10,dpois(1:20,lambda=2)))
  assign(cnames[i],presc)
  pmat=cbind(pmat,presc)
}
colnames(pmat)=cnames
NUMBER_OF_PAID_ITEMS=rowSums(pmat)
PAID_GIC_INCL_BB=rowSums(pmat>0.5)
UNIQUE_STUDY_ID=id
YEARMONTH=paste0("Y",xyear,"M",sample(12,n,rep=T))


dx=data.frame(UNIQUE_STUDY_ID,pmat,PAID_GIC_INCL_BB,NUMBER_OF_PAID_ITEMS,YEARMONTH,stringsAsFactors=F)
haven::write_sav(dx,path=paste0(ehr_directory,fnames[xyear]),compress=TRUE)
}



## Systemwatch

id=sample(nid,n)
TWELVE_MONTH_TIME_PERIOD=sample(paste0("y",1:5),n,rep=T)
LOCATION_CODE=paste0(toupper(sample(letters,n,rep=T)),sample(100:999,n,rep=T),toupper(sample(letters,n,rep=T)))
DATE_OF_ADMISSION=as.Date.numeric(sample(7*365,n,rep=T),origin="2011-01-01")
DATE_OF_DISCHARGE=as.Date(DATE_OF_ADMISSION+sample(0:20,n,rep=T))
SPECIALTY=sample(c(
  paste0("A",c(1:3,6:8,11,21,81:82)),
  paste0("A",c("A","B","C","D","F","FA","G","H","J","M","P","Q","R","V","W")),
  paste0("C",c(1:9,11:13,31,41:42,91)),
  paste0("C",c("A","B")),
  paste0("D",c(1,3:6,8)),
  paste0("E",11:12),
  paste0("F",c(1:3,31:32)),
  paste0("G",c(1:5,21:22)),
  paste0("H",1:2),
  paste0("J",4:5),
  paste0("R",c(1,2,4,5,9)),
  paste0("R",c("B","F","H","K","K1","K3","S")),
  paste0("T",c(2,21))),n,rep=T)
ADMISSION_TYPE=sample(c(36,11,33,30,18,31,10,35,38,32,12,39,20,22,34,19,21),n,rep=T)
ADMISSION_TRANSFER_FROM=toupper(paste0(sample(letters,n,rep=T),sample(999,n,rep=T),sample(letters,n,rep=T)))
DATE_OF_BIRTH=id_dob[id]
MANAGEMENT_OF_PATIENT=sample(c("",1:8),n,rep=T)
SEX=id_s[id]
SIGNIFICANT_FACILITY=sample(c("",11,31:39),n,rep=T)
DISCHARGE_TYPE=sample(c("",10:13,18:22,28:29,40:43),n,rep=T)
DISCHARGE_TRANSFER_TO=toupper(sample(c("","0",10:69,"DI",paste0("4",letters[1:8]),paste0("5",letters[1:8])),n,rep=T))
DIAGNOSIS_1=rep("",n); s1=sample(n,round(n*0.5)); DIAGNOSIS_1[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(999,length(s1),rep=T)))
DIAGNOSIS_2=rep("",n); s2=sample(s1,round(length(s1)*0.5)); DIAGNOSIS_2[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(100:999,length(s2),rep=T)))
DIAGNOSIS_3=rep("",n); s3=sample(s2,round(length(s2)*0.5)); DIAGNOSIS_3[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(100:999,length(s3),rep=T)))
DIAGNOSIS_4=rep("",n); s4=sample(s3,round(length(s3)*0.5)); DIAGNOSIS_4[s4]=toupper(paste0(sample(letters,length(s4),rep=T),sample(100:999,length(s4),rep=T)))
DIAGNOSIS_5=rep("",n); s5=sample(s4,round(length(s4)*0.5)); DIAGNOSIS_5[s5]=toupper(paste0(sample(letters,length(s5),rep=T),sample(100:999,length(s5),rep=T)))
DIAGNOSIS_6=rep("",n); s6=sample(s5,round(length(s5)*0.5)); DIAGNOSIS_6[s6]=toupper(paste0(sample(letters,length(s6),rep=T),sample(100:999,length(s6),rep=T)))
OPERATION_CODE_1A=rep("",n); s1=sample(n,round(n*0.3)); OPERATION_CODE_1A[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(100:999,length(s1),rep=T)))
OPERATION_CODE_1B=rep("",n); OPERATION_CODE_1B[s1]=toupper(paste0(sample(letters,length(s1),rep=T),sample(100:999,length(s1),rep=T)))
OPERATION_CODE_2A=rep("",n); s2=sample(s1,round(length(s1)*0.3)); OPERATION_CODE_2A[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(100:999,length(s2),rep=T)))
OPERATION_CODE_2B=rep("",n); OPERATION_CODE_2B[s2]=toupper(paste0(sample(letters,length(s2),rep=T),sample(100:999,length(s2),rep=T)))
OPERATION_CODE_3A=rep("",n); s3=sample(s2,round(length(s2)*0.3)); OPERATION_CODE_3A[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(100:999,length(s3),rep=T)))
OPERATION_CODE_3B=rep("",n); OPERATION_CODE_3B[s3]=toupper(paste0(sample(letters,length(s3),rep=T),sample(100:999,length(s3),rep=T)))
OPERATION_CODE_4A=rep("",n); s4=sample(s3,round(length(s3)*0.3)); OPERATION_CODE_4A[s4]=toupper(paste0(sample(letters,length(s4),rep=T),sample(100:999,length(s4),rep=T)))
OPERATION_CODE_4B=rep("",n); OPERATION_CODE_4B[s4]=toupper(paste0(sample(letters,length(s4),rep=T),sample(100:999,length(s4),rep=T)))
UNIQUE_STUDY_ID=id

dx=data.frame(TWELVE_MONTH_TIME_PERIOD,
              LOCATION_CODE,
              DATE_OF_ADMISSION,
              DATE_OF_DISCHARGE,
              SPECIALTY,
              ADMISSION_TYPE,
              ADMISSION_TRANSFER_FROM,
              DATE_OF_BIRTH,
              MANAGEMENT_OF_PATIENT,
              SEX,
              SIGNIFICANT_FACILITY,
              DISCHARGE_TYPE,
              DISCHARGE_TRANSFER_TO,
              DIAGNOSIS_1,
              DIAGNOSIS_2,
              DIAGNOSIS_3,
              DIAGNOSIS_4,
              DIAGNOSIS_5,
              DIAGNOSIS_6,
              OPERATION_CODE_1A,
              OPERATION_CODE_1B,
              OPERATION_CODE_2A,
              OPERATION_CODE_2B,
              OPERATION_CODE_3A,
              OPERATION_CODE_3B,
              OPERATION_CODE_4A,
              OPERATION_CODE_4B,
              UNIQUE_STUDY_ID,stringsAsFactors=FALSE)

for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=format(dx[,i],"%d-%m-%Y")
haven::write_sav(dx,path=paste0(ehr_directory,"Final_SystemWatch_extract_incl_UniqueStudyID_v2.zsav"),compress=TRUE)





## LTCs

id=sample(nid,n)
UNIQUE_STUDY_ID=id
ltc=c("ARTHRITIS","ASTHMA","ATRIAL_FIBRILLATION","COPD","CANCER","CEREBROVASCULAR_DISEASE",
      "CHRONIC_LIVER_DISEASE","DEMENTIA","DIABETES","EPILEPSY","HEART_DISEASE","HEART_FAILURE",
      "MULTIPLE_SCLEROSIS","PARKINSON_DISEASE","RENAL_FAILURE")
NUMBER_OF_LTCS=pmin(pmax(round(rchisq(n,df=2)),1),length(ltc))
xmat=matrix(NA,n,length(ltc)); for (i in 1:n) xmat[i,sample(length(ltc),NUMBER_OF_LTCS[i])]=1
dt0=as.Date("2018-05-01")

for (i in 1:length(ltc)) {
  dt=as.Date(-sample(13634,n,rep=T)*xmat[,i],origin=dt0)
  assign(paste0("FIRST_",ltc[i],"_EPISODE"),dt)
}


dx=data.frame(UNIQUE_STUDY_ID,
              NUMBER_OF_LTCS,stringsAsFactors=FALSE)
for (i in 1:length(ltc)) {
  d0=get(paste0("FIRST_",ltc[i],"_EPISODE"))
  dx=cbind(dx,d0)
  colnames(dx)[dim(dx)[2]]=paste0("FIRST_",ltc[i],"_EPISODE")
}

#for (i in 1:dim(dx)[2]) if (class(dx[1,i])=="Date") dx[,i]=format(dx[,i],"%d-%m-%Y")
haven::write_sav(dx,path=paste0(ehr_directory,"final_LTCs_using_new_ICD10codes_incl_uniquestudyid.zsav"),compress=TRUE)



## SPARRA v3 scores

UNIQUE_STUDY_ID=nid
times=outer(c("01","02","03","04","05","06","07","08","09","10","11","12"),2013:2018,
             function(x,y) paste0("SPARRA_RISK_SCORE_01_",x,"_",y))
times1=times[7:40]
times2=times[40:65]

dx1=data.frame(UNIQUE_STUDY_ID)
for (i in 1:length(times1)) {
  tx=round(100*(runif(length(nid))^3)); tx[which(runif(length(tx))<0.1)]=NA
  dx1=cbind(dx1,tx)
}
colnames(dx1)=c("UNIQUE_STUDY_ID",times1)

dx2=data.frame(UNIQUE_STUDY_ID,dx1[,dim(dx1)[2]])
for (i in 2:length(times2)) {
  tx=round(100*(runif(length(nid))^3)); tx[which(runif(length(tx))<0.1)]=NA
  dx2=cbind(dx2,tx)
}

colnames(dx2)=c("UNIQUE_STUDY_ID",times2)

for (i in 1:dim(dx1)[2]) if (class(dx1[1,i])=="Date") dx1[,i]=format(dx1[,i],"%d-%m-%Y")
for (i in 1:dim(dx2)[2]) if (class(dx2[1,i])=="Date") dx2[,i]=format(dx2[,i],"%d-%m-%Y")

haven::write_sav(dx1,path=paste0(ehr_directory,"SPARRAscores_1July13_1Apr16_ATI.sav"),compress=TRUE)
haven::write_sav(dx2,path=paste0(ehr_directory,"Final_SPARRA_extract_incl_UniqueStudyID.zsav"),compress=TRUE)

}





# All installed packages and versions
# To generate:
#  ip_windows=installed.packages()[,c(1,3)] # on windows
#   dput(ip_windows) (paste somewhere)
#  ip_linux=installed.packages()[,c(1,3)] # on linux
#   dput(ip_linux) (paste somewhere)

ip_linux=structure(c("assertthat", "autokeras", "backports", "base64enc", 
"BH", "bindr", "bindrcpp", "BiocManager", "bit", "bit64", "bitops", 
"blob", "brew", "broom", "callr", "caTools", "cellranger", "cli", 
"clipr", "clisymbols", "colorspace", "commonmark", "config", 
"covr", "crayon", "curl", "data.table", "DBI", "dbplyr", "desc", 
"devtools", "digest", "docopt", "dplyr", "dtplyr", "evaluate", 
"fansi", "feather", "foghorn", "forcats", "forge", "formatR", 
"fs", "fst", "generics", "ggplot2", "gh", "git2r", "glue", "gmailr", 
"gtable", "h2o", "haven", "highlight", "highr", "hms", "htmltools", 
"htmlwidgets", "httpuv", "httr", "hunspell", "igraph", "ini", 
"IRdisplay", "IRkernel", "jsonlite", "keras", "knitr", "labeling", 
"Lahman", "later", "lazyeval", "lintr", "littler", "lubridate", 
"magrittr", "markdown", "memoise", "microbenchmark", "mime", 
"mockery", "modelr", "munsell", "nycflights13", "openssl", "packrat", 
"parsedate", "pbdZMQ", "pillar", "pingr", "pkgbuild", "pkgconfig", 
"pkgdown", "pkgload", "plogr", "plyr", "praise", "prettyunits", 
"processx", "progress", "promises", "ps", "purrr", "r2d3", "R6", 
"rappdirs", "rcmdcheck", "RColorBrewer", "Rcpp", "RCurl", "readr", 
"readxl", "rematch", "rematch2", "remotes", "repr", "reprex", 
"reshape2", "reticulate", "rex", "rhub", "rlang", "rmarkdown", 
"RMySQL", "roxygen2", "RPostgreSQL", "rprojroot", "rsconnect", 
"RSQLite", "rstudioapi", "rversions", "rvest", "scales", "selectr", 
"sessioninfo", "shiny", "sourcetools", "sparklyr", "spelling", 
"stringdist", "stringi", "stringr", "tensorflow", "testit", "testthat", 
"tfruns", "tibble", "tidyr", "tidyselect", "tidyverse", "tinytex", 
"usethis", "utf8", "uuid", "viridisLite", "whisker", "whoami", 
"withr", "xfun", "XML", "xml2", "xopen", "xtable", "yaml", "zeallot", 
"base", "boot", "class", "cluster", "codetools", "compiler", 
"datasets", "foreign", "graphics", "grDevices", "grid", "KernSmooth", 
"lattice", "MASS", "Matrix", "methods", "mgcv", "nlme", "nnet", 
"parallel", "rpart", "spatial", "splines", "stats", "stats4", 
"survival", "tcltk", "tools", "utils", "0.2.0", "0.1.0", "1.1.3", 
"0.1-3", "1.66.0-1", "0.1.1", "0.2.2", "1.30.4", "1.1-14", "0.9-7", 
"1.0-6", "1.1.1", "1.0-6", "0.5.1", "3.1.0", "1.17.1.1", "1.1.0", 
"1.0.1", "0.4.1", "1.2.0", "1.3-2", "1.7", "0.3", "3.2.1", "1.3.4", 
"3.2", "1.11.8", "1.0.0", "1.2.2", "1.2.0", "2.0.1", "0.6.18", 
"0.6.1", "0.7.8", "0.0.2", "0.12", "0.4.0", "0.3.1", "1.0.2", 
"0.3.0", "0.1.0", "1.5", "1.2.6", "0.8.10", "0.0.2", "3.1.0", 
"1.0.1", "0.23.0", "1.3.0", "0.7.1", "0.2.0", "3.24.0.2", "2.0.0", 
"0.4.7.2", "0.7", "0.4.2", "0.3.6", "1.3", "1.4.5.1", "1.4.0", 
"3.0", "1.2.2", "0.3.1", "0.7.0", "0.8.14", "1.6", "2.2.4", "1.21", 
"0.3", "6.0-0", "0.7.5", "0.2.1", "1.0.3", "0.3.5", "1.7.4", 
"1.5", "0.9", "1.1.0", "1.4-6", "0.6", "0.4.1.1", "0.1.2", "0.5.0", 
"1.0.0", "1.1", "0.5.0", "1.1.3", "0.3-3", "1.3.1", "1.1.2", 
"1.0.2", "2.0.2", "1.3.0", "1.0.2", "0.2.0", "1.8.4", "1.0.0", 
"1.0.2", "3.2.1", "1.2.0", "1.0.1", "1.2.1", "0.2.5", "0.2.3", 
"2.3.0", "0.3.1", "1.3.2", "1.1-2", "1.0.0", "1.95-4.11", "1.3.0", 
"1.2.0", "1.0.1", "2.0.1", "2.0.2", "0.18", "0.2.1", "1.4.3", 
"1.10", "1.1.2", "1.0.2", "0.3.0.1", "1.11", "0.10.15", "6.1.1", 
"0.6-2", "1.3-2", "0.8.12", "2.1.1", "0.8", "1.0.3", "0.3.2", 
"1.0.0", "0.4-1", "1.1.1", "1.2.0", "0.1.7", "0.9.3", "2.0", 
"0.9.5.1", "1.2.4", "1.3.1", "1.10", "0.9", "2.0.1", "1.4", "1.4.2", 
"0.8.2", "0.2.5", "1.2.1", "0.9", "1.4.0", "1.1.4", "0.1-2", 
"0.3.0", "0.3-2", "1.2.0", "2.1.2", "0.4", "3.98-1.16", "1.2.0", 
"1.0.0", "1.8-3", "2.2.0", "0.1.0", "3.5.1", "1.3-20", "7.3-14", 
"2.0.7-1", "0.2-15", "3.5.1", "3.5.1", "0.8-70", "3.5.1", "3.5.1", 
"3.5.1", "2.23-15", "0.20-35", "7.3-50", "1.2-14", "3.5.1", "1.8-24", 
"3.1-137", "7.3-12", "3.5.1", "4.1-13", "7.3-11", "3.5.1", "3.5.1", 
"3.5.1", "2.42-3", "3.5.1", "3.5.1", "3.5.1"), .Dim = c(194L, 
2L), .Dimnames = list(c("assertthat", "autokeras", "backports", 
"base64enc", "BH", "bindr", "bindrcpp", "BiocManager", "bit", 
"bit64", "bitops", "blob", "brew", "broom", "callr", "caTools", 
"cellranger", "cli", "clipr", "clisymbols", "colorspace", "commonmark", 
"config", "covr", "crayon", "curl", "data.table", "DBI", "dbplyr", 
"desc", "devtools", "digest", "docopt", "dplyr", "dtplyr", "evaluate", 
"fansi", "feather", "foghorn", "forcats", "forge", "formatR", 
"fs", "fst", "generics", "ggplot2", "gh", "git2r", "glue", "gmailr", 
"gtable", "h2o", "haven", "highlight", "highr", "hms", "htmltools", 
"htmlwidgets", "httpuv", "httr", "hunspell", "igraph", "ini", 
"IRdisplay", "IRkernel", "jsonlite", "keras", "knitr", "labeling", 
"Lahman", "later", "lazyeval", "lintr", "littler", "lubridate", 
"magrittr", "markdown", "memoise", "microbenchmark", "mime", 
"mockery", "modelr", "munsell", "nycflights13", "openssl", "packrat", 
"parsedate", "pbdZMQ", "pillar", "pingr", "pkgbuild", "pkgconfig", 
"pkgdown", "pkgload", "plogr", "plyr", "praise", "prettyunits", 
"processx", "progress", "promises", "ps", "purrr", "r2d3", "R6", 
"rappdirs", "rcmdcheck", "RColorBrewer", "Rcpp", "RCurl", "readr", 
"readxl", "rematch", "rematch2", "remotes", "repr", "reprex", 
"reshape2", "reticulate", "rex", "rhub", "rlang", "rmarkdown", 
"RMySQL", "roxygen2", "RPostgreSQL", "rprojroot", "rsconnect", 
"RSQLite", "rstudioapi", "rversions", "rvest", "scales", "selectr", 
"sessioninfo", "shiny", "sourcetools", "sparklyr", "spelling", 
"stringdist", "stringi", "stringr", "tensorflow", "testit", "testthat", 
"tfruns", "tibble", "tidyr", "tidyselect", "tidyverse", "tinytex", 
"usethis", "utf8", "uuid", "viridisLite", "whisker", "whoami", 
"withr", "xfun", "XML", "xml2", "xopen", "xtable", "yaml", "zeallot", 
"base", "boot", "class", "cluster", "codetools", "compiler", 
"datasets", "foreign", "graphics", "grDevices", "grid", "KernSmooth", 
"lattice", "MASS", "Matrix", "methods", "mgcv", "nlme", "nnet", 
"parallel", "rpart", "spatial", "splines", "stats", "stats4", 
"survival", "tcltk", "tools", "utils"), c("Package", "Version"
)))
  
  
ip_windows=structure(c("abind", "ald", "askpass", "assertthat", "AUC", "backports", 
  "base64enc", "BBmisc", "benchmarkme", "benchmarkmeData", "BH", 
  "bitops", "brew", "broom", "CalibratR", "callr", "car", "carData", 
  "caret", "caTools", "cellranger", "checkmate", "classInt", "cli", 
  "clipr", "clisymbols", "coda", "colorspace", "commonmark", "corrplot", 
  "covr", "crayon", "crosstalk", "curl", "cvAUC", "cyclocomp", 
  "data.table", "DBI", "dbplyr", "deldir", "desc", "devtools", 
  "digest", "doParallel", "dplyr", "DT", "dtplyr", "e1071", "ellipsis", 
  "evaluate", "expm", "fansi", "fastmap", "fastmatch", "fBasics", 
  "fitdistrplus", "forcats", "foreach", "Formula", "fs", "fst", 
  "furrr", "future", "future.apply", "fuzzyjoin", "gamlss.dist", 
  "gdata", "generics", "geosphere", "ggalluvial", "ggplot2", "gh", 
  "git2r", "glmnet", "globals", "glue", "gmodels", "gower", "gplots", 
  "gridExtra", "groupdata2", "gss", "gtable", "gtools", "h2o", 
  "haven", "hexbin", "highr", "hms", "htmltools", "htmlwidgets", 
  "httpuv", "httr", "iml", "import", "ini", "ipred", "iterators", 
  "jsonlite", "kernlab", "knitr", "labeling", "labelled", "later", 
  "latex2exp", "lava", "lazyeval", "LearnBayes", "lifecycle", "lintr", 
  "listenv", "lme4", "lsei", "lsr", "lubridate", "magrittr", "maptools", 
  "maptpx", "markdown", "MatrixModels", "matrixStats", "memoise", 
  "Metrics", "mime", "minqa", "mlr", "ModelMetrics", "modelr", 
  "modeltools", "munsell", "nloptr", "NLP", "npsurv", "numbers", 
  "numDeriv", "nycflights13", "openssl", "openxlsx", "parallelMap", 
  "ParamHelpers", "pbkrtest", "pillar", "pkgbuild", "pkgconfig", 
  "pkgload", "plogr", "plotly", "plotrix", "plyr", "praise", "prediction", 
  "PReMiuM", "prettyunits", "pROC", "processx", "prodlim", "profvis", 
  "progress", "promises", "PRROC", "pryr", "ps", "purrr", "quantreg", 
  "R6", "randomForest", "ranger", "rcmdcheck", "RColorBrewer", 
  "Rcpp", "RcppEigen", "RCurl", "readr", "readxl", "recipes", "rematch", 
  "rematch2", "remotes", "reprex", "reshape2", "rex", "rgdal", 
  "rio", "rJava", "rlang", "rmarkdown", "ROCR", "roxygen2", "rprojroot", 
  "rstudioapi", "rversions", "rvest", "scales", "selectr", "sessioninfo", 
  "sf", "shiny", "slam", "sourcetools", "sp", "SparseM", "spData", 
  "spdep", "speedglm", "SQUAREM", "stabledist", "stringdist", "stringi", 
  "stringr", "styler", "sys", "testthat", "tibble", "tictoc", "tidyr", 
  "tidyselect", "tidyverse", "timeDate", "timeSeries", "tinytex", 
  "tm", "topicmodels", "units", "usethis", "utf8", "varhandle", 
  "vctrs", "viridisLite", "whisker", "withr", "xfun", "xgboost", 
  "xlsx", "xlsxjars", "XML", "xml2", "xmlparsedata", "xopen", "xtable", 
  "yaml", "zeallot", "zip", "base", "boot", "class", "cluster", 
  "codetools", "compiler", "datasets", "foreign", "graphics", "grDevices", 
  "grid", "KernSmooth", "lattice", "MASS", "Matrix", "methods", 
  "mgcv", "nlme", "nnet", "parallel", "rpart", "spatial", "splines", 
  "stats", "stats4", "survival", "tcltk", "tools", "translations", 
  "utils", "1.4-5", "1.2", "1.1", "0.2.1", "0.3.0", "1.1.5", "0.1-3", 
  "1.11", "1.0.2", "1.0.3", "1.69.0-1", "1.0-6", "1.0-6", "0.5.2", 
  "0.1.2", "3.3.2", "3.0-4", "3.0-2", "6.0-84", "1.17.1.2", "1.1.0", 
  "1.9.4", "0.4-2", "1.1.0", "0.7.0", "1.2.0", "0.19-3", "1.4-1", 
  "1.7", "0.84", "3.3.2", "1.3.4", "1.0.0", "4.2", "1.1.0", "1.1.0", 
  "1.12.6", "1.0.0", "1.4.2", "0.1-23", "1.2.0", "2.2.1", "0.6.22", 
  "1.0.15", "0.8.3", "0.9", "1.0.0", "1.7-2", "0.3.0", "0.14", 
  "0.999-4", "0.4.0", "1.0.1", "1.1-0", "3042.89", "1.0-14", "0.4.0", 
  "1.4.7", "1.2-3", "1.3.1", "0.9.0", "0.1.0", "1.16.0", "1.4.0", 
  "0.1.5", "5.1-5", "2.18.0", "0.0.2", "1.5-10", "0.10.0", "3.2.1", 
  "1.0.1", "0.26.1", "2.0-18", "0.12.5", "1.3.1", "2.18.1", "0.2.1", 
  "3.0.1.1", "2.3", "1.1.2", "2.1-10", "0.3.0", "3.8.1", "3.26.0.2", 
  "2.1.1", "1.27.3", "0.7", "0.5.2", "0.4.0", "1.5.1", "1.5.2", 
  "1.4.1", "0.10.0", "1.1.0", "0.3.1", "0.9-9", "1.0.12", "1.6", 
  "0.9-27", "1.25", "0.3", "2.2.1", "1.0.0", "0.4.0", "1.6.6", 
  "0.2.2", "2.15.1", "0.1.0", "2.0.0", "0.8.0", "1.1-21", "1.2-0", 
  "0.5", "1.7.4", "1.5", "0.9-8", "1.9-2", "0.8", "0.4-1", "0.55.0", 
  "1.1.0", "0.1.4", "0.5", "1.2.4", "2.15.0", "1.2.2", "0.1.5", 
  "0.2-22", "0.5.0", "1.2.1", "0.2-0", "0.4-0", "0.7-1", "2016.8-1.1", 
  "1.0.1", "1.4.1", "4.1.2", "1.4", "1.12", "0.4-7", "1.4.2", "1.0.6", 
  "2.0.3", "1.0.2", "0.2.0", "4.9.0", "3.7-6", "1.8.4", "1.0.0", 
  "0.3.14", "3.2.3", "1.0.2", "1.15.3", "3.4.1", "2019.10.13", 
  "0.3.6", "1.2.2", "1.1.0", "1.3.1", "0.1.4", "1.3.0", "0.3.3", 
  "5.51", "2.4.0", "4.6-14", "0.11.2", "1.3.3", "1.1-2", "1.0.2", 
  "0.3.3.5.0", "1.95-4.12", "1.3.1", "1.3.1", "0.1.7", "1.0.1", 
  "2.1.0", "2.1.0", "0.3.0", "1.4.3", "1.1.2", "1.4-7", "0.5.16", 
  "0.9-11", "0.4.1", "1.16", "1.0-7", "6.1.1", "1.3-2", "0.10", 
  "2.0.0", "0.3.4", "1.0.0", "0.4-1", "1.1.1", "0.8-0", "1.4.0", 
  "0.1-46", "0.1.7", "1.3-1", "1.77", "0.3.2", "1.1-3", "0.3-2", 
  "2017.10-1", "0.7-1", "0.9.5.5", "1.4.3", "1.3.1", "1.2.0", "3.3", 
  "2.2.1", "2.1.3", "1.0", "1.0.0", "0.2.5", "1.2.1", "3043.102", 
  "3042.102", "0.17", "0.7-6", "0.2-8", "0.6-5", "1.5.1", "1.1.4", 
  "2.0.4", "0.2.0", "0.3.0", "0.4", "2.1.2", "0.10", "0.90.0.2", 
  "0.6.1", "0.6.1", "3.98-1.20", "1.2.2", "1.0.3", "1.0.0", "1.8-4", 
  "2.1.19", "0.1.0", "2.0.4", "3.6.1", "1.3-22", "7.3-15", "2.1.0", 
  "0.2-16", "3.6.1", "3.6.1", "0.8-71", "3.6.1", "3.6.1", "3.6.1", 
  "2.23-15", "0.20-38", "7.3-51.4", "1.2-17", "3.6.1", "1.8-28", 
  "3.1-140", "7.3-12", "3.6.1", "4.1-15", "7.3-11", "3.6.1", "3.6.1", 
  "3.6.1", "2.44-1.1", "3.6.1", "3.6.1", "3.6.1", "3.6.1"), .Dim = c(272L, 
    2L), .Dimnames = list(c("abind", "ald", "askpass", "assertthat", 
      "AUC", "backports", "base64enc", "BBmisc", "benchmarkme", "benchmarkmeData", 
      "BH", "bitops", "brew", "broom", "CalibratR", "callr", "car", 
      "carData", "caret", "caTools", "cellranger", "checkmate", "classInt", 
      "cli", "clipr", "clisymbols", "coda", "colorspace", "commonmark", 
      "corrplot", "covr", "crayon", "crosstalk", "curl", "cvAUC", "cyclocomp", 
      "data.table", "DBI", "dbplyr", "deldir", "desc", "devtools", 
      "digest", "doParallel", "dplyr", "DT", "dtplyr", "e1071", "ellipsis", 
      "evaluate", "expm", "fansi", "fastmap", "fastmatch", "fBasics", 
      "fitdistrplus", "forcats", "foreach", "Formula", "fs", "fst", 
      "furrr", "future", "future.apply", "fuzzyjoin", "gamlss.dist", 
      "gdata", "generics", "geosphere", "ggalluvial", "ggplot2", "gh", 
      "git2r", "glmnet", "globals", "glue", "gmodels", "gower", "gplots", 
      "gridExtra", "groupdata2", "gss", "gtable", "gtools", "h2o", 
      "haven", "hexbin", "highr", "hms", "htmltools", "htmlwidgets", 
      "httpuv", "httr", "iml", "import", "ini", "ipred", "iterators", 
      "jsonlite", "kernlab", "knitr", "labeling", "labelled", "later", 
      "latex2exp", "lava", "lazyeval", "LearnBayes", "lifecycle", "lintr", 
      "listenv", "lme4", "lsei", "lsr", "lubridate", "magrittr", "maptools", 
      "maptpx", "markdown", "MatrixModels", "matrixStats", "memoise", 
      "Metrics", "mime", "minqa", "mlr", "ModelMetrics", "modelr", 
      "modeltools", "munsell", "nloptr", "NLP", "npsurv", "numbers", 
      "numDeriv", "nycflights13", "openssl", "openxlsx", "parallelMap", 
      "ParamHelpers", "pbkrtest", "pillar", "pkgbuild", "pkgconfig", 
      "pkgload", "plogr", "plotly", "plotrix", "plyr", "praise", "prediction", 
      "PReMiuM", "prettyunits", "pROC", "processx", "prodlim", "profvis", 
      "progress", "promises", "PRROC", "pryr", "ps", "purrr", "quantreg", 
      "R6", "randomForest", "ranger", "rcmdcheck", "RColorBrewer", 
      "Rcpp", "RcppEigen", "RCurl", "readr", "readxl", "recipes", "rematch", 
      "rematch2", "remotes", "reprex", "reshape2", "rex", "rgdal", 
      "rio", "rJava", "rlang", "rmarkdown", "ROCR", "roxygen2", "rprojroot", 
      "rstudioapi", "rversions", "rvest", "scales", "selectr", "sessioninfo", 
      "sf", "shiny", "slam", "sourcetools", "sp", "SparseM", "spData", 
      "spdep", "speedglm", "SQUAREM", "stabledist", "stringdist", "stringi", 
      "stringr", "styler", "sys", "testthat", "tibble", "tictoc", "tidyr", 
      "tidyselect", "tidyverse", "timeDate", "timeSeries", "tinytex", 
      "tm", "topicmodels", "units", "usethis", "utf8", "varhandle", 
      "vctrs", "viridisLite", "whisker", "withr", "xfun", "xgboost", 
      "xlsx", "xlsxjars", "XML", "xml2", "xmlparsedata", "xopen", "xtable", 
      "yaml", "zeallot", "zip", "base", "boot", "class", "cluster", 
      "codetools", "compiler", "datasets", "foreign", "graphics", "grDevices", 
      "grid", "KernSmooth", "lattice", "MASS", "Matrix", "methods", 
      "mgcv", "nlme", "nnet", "parallel", "rpart", "spatial", "splines", 
      "stats", "stats4", "survival", "tcltk", "tools", "translations", 
      "utils"), c("Package", "Version")))
  
  
  


##' Check package versions match that used on NSH
##' 
##' @name check_package_versions
##' @param linux packages are different on different machines, so need to check.
##' @returns list: correct=TRUE/FALSE, failed=list of packages with wrong version, version=current version, shouldbe=correct version
check_package_versions=function(linux) {
  ip=installed.packages() 
  if (linux) ip_standard=ip_linux else ip_standard=ip_windows
  ipx=intersect(ip[,1],ip_standard[,1])
  version_correct=ip_standard[match(ipx,ip_standard[,1]),2]
  version_is=ip[match(ipx,ip[,1]),3]
  correct=all(version_correct==version_is)
  w=which(version_correct!=version_is)
  if (length(w)==0) {
    failed=NULL; version=NULL; shouldbe=NULL
  } else {
    failed=ipx[w]; version=version_is[w]; shouldbe=version_correct[w]
  }
  return(list(correct=correct,failed=failed,version=version,shouldbe=shouldbe))
}