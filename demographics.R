#Demographics
library(tableone)
filenames <- list.files("C:/Users/wkayla/OneDrive for Business/Test/Macro", pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)

names(ldf)<-gsub("C:/Users/wkayla/OneDrive for Business/Test/Macro/","",filenames)

ldf$

Subject_ID=c("1027556",
"1027557",
"1027559",
"1027562",
"1027563",
"1027565",
"1027566",
"1027567",
"1027569",
"1027570",
"1027571",
"1027576",
"1027578",
"1027579",
"1027580",
"1027582")


demo=ldf$`_vapdemograph.csv`[ldf$`_vapdemograph.csv`$SubjectID%in%Subject_ID,]
demo=merge(demo,ldf$`_vaphospitali.csv`,"SubjectID")

demo$X_RACE=factor(demo$X_RACE,levels=c("1",   "2",   "3",   "4",   "4,5", "5",   "92"),
                      labels=c("American Indian or Alaska Native",
                               "Asian",
                               "Black or African American",
                               "Native Hawaiian or Other Pacific Islander",
                               "Native Hawaiian or Other Pacific Islander and White",
                               "White",
                               "Unknown"))

demo$BIRTHDATE=as.Date.character(demo$X_BIRTHDATE)
demo$HOSPADMITDATE=as.Date.character(demo$X_HOSPADMITDATE)


demo$`Age at Admit In Years`=as.numeric(as.character(round((demo$HOSPADMITDATE-demo$BIRTHDATE)/365,3)))

CreateTableOne(colnames(demo)[c(6:8,56)],data=demo)
