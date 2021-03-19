library(data.table)
library(haven)
library(dplyr)
library(survival)
## read the two datasets
dem<-fread("I:/2020-09-30/Demographics/dem_ctos_inv.dat")
outcome<-fread("I:/2020-09-30/Outcomes/outc_ct_os_inv.dat")
agreement<-read_sas('U:/crctable_vde.sas7bdat')
agreement_clean<-subset(agreement,agreement$crc_na !='X')
date <- read.csv("U:/capstone_perturbed_crcdate.csv")

## match with Roberta's dataset
dem_match<-dem[dem$ID %in% agreement_clean$id]
out_match<-outcome[outcome$ID %in% agreement_clean$id]



# 1) hormone trial
dem_ht <- subset(dem, HRTFLAG == 1)

# 2) perturb: WHI start date = CMS start date
colnames(date)[1] <- "ID"
date$whi_st <- as.Date(date$WHISTARTDT, format = "%m/%d/%Y")
date$medicare_st <- as.Date(date$CMSSTARTDT, format = "%m/%d/%Y")
date$medicare_end <- as.Date(date$CMSSTOPDT, format = "%m/%d/%Y")
date$medicare_diag <- as.Date(date$P_CMSCRCDT, format = "%m/%d/%Y")
date_same <- subset(date, date$whi_st == date$medicare_st)

# 1) + 2) + 3) found in roberta
colnames(agreement)[1] <- "ID"
agreement_ht <-subset(agreement,agreement$crc_na !='X' & 
                        agreement$ID %in% dem_ht$ID & 
                        agreement$ID %in% date_same$ID)
final_dt <-Reduce(function(x,y) merge(x=x,y=y,by='ID'), 
                  list(agreement_ht, date_same, dem)) #8,148 patients

# merge with COLORECTALDY
col_time <- out_match[,c("COLORECTAL", "COLORECTALDY","ID")]
all_dat <- Reduce(function(x,y) merge(x=x,y=y,by='ID'), 
                  list(final_dt, col_time))
head(all_dat, 10)

# still need to get start end dates for censored patients
# start is whi_st, already converted to as.date format
all_dat$whi_end <- ifelse(all_dat$HRTARM == "1" | all_dat$HRTARM == "2", "2004-03-31", "2002-07-31")
all_dat$time_to_end <- as.Date(all_dat$whi_end) - as.Date(all_dat$whi_st)
# if colorectal cancer, will have recorded event date, if not use study end date
all_dat$time_to_event <- ifelse(all_dat$COLORECTAL == 1, all_dat$COLORECTALDY, all_dat$time_to_end)
# if colorectal diagnosed after E+P study ends, censor them
all_dat$time_to_colo <- ifelse(all_dat$time_to_event >= all_dat$time_to_end, all_dat$time_to_end, all_dat$time_to_event)
head(test, 10)
all_dat$adj_colo <- ifelse(all_dat$time_to_event >= all_dat$time_to_end, 0, all_dat$COLORECTAL)
# pool the hormone therapy arms
all_dat$poolarm <- ifelse(all_dat$HRTARM == "1" | all_dat$HRTARM == "3", 1, 0)

# run the coxph
# whi
cox1 <- coxph(Surv(time_to_colo, adj_colo)~poolarm, data = all_dat)
mean(all_dat$time_to_colo)
summary(cox1)
# medicare
cox2 <- coxph(Surv(time, medicare_diag.x)~trt, data = joint_dt)


# run Yichen's code to get the Medicare dataset
# obtain the needed columns
medicare_dt <- hr_dt[,c("ID","time", "medicare_diag","trt")]
joint_dt <-Reduce(function(x,y) merge(x=x,y=y,by='ID'), 
                  list(medicare_dt,all_dat))
head(joint_dt)
library(boot)

# get the hr for whi
get.hr <- function(data, indices) {
  whi <- data[,c("time_to_colo", "adj_colo","poolarm")]
  d_whi <- whi[indices,]
  fit_whi <- coxph(Surv(time_to_colo, adj_colo)~poolarm, data=d_whi)
  exp(fit_whi$coefficient)
}
hr.whi <- boot(data=joint_dt, statistic=get.hr, R=5000)
boot.ci(hr.whi,type="norm")
# get hr for medicare
get.hr.cms <- function(data, indices) {
  medicare <- data[,c("time", "medicare_diag.x","trt")]
  d_medi <- medicare[indices,]
  fit_medi <- coxph(Surv(time, medicare_diag.x)~trt, data = d_medi)
  exp(fit_medi$coefficient)
}
hr.cms <- boot(data=joint_dt, statistic=get.hr.cms, R=5000)
boot.ci(hr.cms,type="norm")
#

compare.hr <- function(data, indices1, indices2) {
  medicare <- data[,c("time", "medicare_diag.x","trt")]
  whi <- data[,c("time_to_colo", "adj_colo","poolarm")]
  d_medi <- medicare[indices1,]
  d_whi <- whi[indices2,]
  fit_medi <- coxph(Surv(time, medicare_diag.x)~trt, data = d_medi)
  fit_whi <- coxph(Surv(time_to_colo, adj_colo)~poolarm, data=d_whi)
  as.numeric(exp(fit_medi$coefficient)/exp(fit_whi$coefficient))
}
results <- boot(data=joint_dt, statistic=compare.hr, R=5000)
boot.ci(results,type="norm")
mean(results$t < 1.140)

#####
# get the risk for the WHI diagnosis
# hormone
r.whi.h <- joint_dt %>% filter(poolarm == 1) %>% summarise(prop = mean(adj_colo))
# placebo
r.whi.p <- joint_dt %>% filter(poolarm == 0) %>% summarise(prop = mean(adj_colo))
# rd for WHI
r.whi.h - r.whi.p
# risk for the medicare diagnosis
r.med.h <- joint_dt %>% filter(trt == 1) %>% summarise(prop = mean(medicare_diag.x))
r.med.p <- joint_dt %>% filter(trt == 0) %>% summarise(prop = mean(medicare_diag.x))
r.med.h - r.med.p
# bootstrap
compare.rd <- function(data, indices1, indices2) {
  medicare <- data[,c("time", "medicare_diag.x","trt")]
  whi <- data[,c("time_to_colo", "adj_colo","poolarm")]
  d_medi <- medicare[indices1,]
  d_whi <- whi[indices2,]
  # get the risk for the WHI diagnosis
  # hormone
  r.whi.h <- d_whi %>% filter(poolarm == 1) %>% summarise(prop = mean(adj_colo)) %>% .$prop
  # placebo
  r.whi.p <- d_whi %>% filter(poolarm == 0) %>% summarise(prop = mean(adj_colo)) %>% .$prop
  # rd for WHI
  rd.whi <- r.whi.h - r.whi.p
  r.med.h <- d_medi %>% filter(trt == 1) %>% summarise(prop = mean(medicare_diag.x)) %>% .$prop
  r.med.p <- d_medi %>% filter(trt == 0) %>% summarise(prop = mean(medicare_diag.x)) %>% .$prop
  rd.med <- r.med.h - r.med.p
  rd.med - rd.whi
}
results.rd <- boot(data=joint_dt, statistic=compare.rd, R=5000)
boot.ci(results.rd,type="norm")
mean(results.rd$t < 0.00102)



# rd for whi
get.rd <- function(data, indices) {
  whi <- data[,c("time_to_colo", "adj_colo","poolarm")]
  d_whi <- whi[indices,]
  # hormone
  r.whi.h <- d_whi %>% filter(poolarm == 1) %>% summarise(prop = mean(adj_colo)) %>% .$prop
  # placebo
  r.whi.p <- d_whi %>% filter(poolarm == 0) %>% summarise(prop = mean(adj_colo)) %>% .$prop
  # rd for WHI
  r.whi.h - r.whi.p
}
rd.whi <- boot(data=joint_dt, statistic=get.rd, R=5000)
boot.ci(rd.whi,type="norm")


#rd for medicare
get.rd.cms <- function(data, indices) {
  medicare <- data[,c("time", "medicare_diag.x","trt")]
  d_medi <- medicare[indices,]
  # get the risk for the WHI diagnosis
  r.med.h <- d_medi %>% filter(trt == 1) %>% summarise(prop = mean(medicare_diag.x)) %>% .$prop
  r.med.p <- d_medi %>% filter(trt == 0) %>% summarise(prop = mean(medicare_diag.x)) %>% .$prop
  rd.med <- r.med.h - r.med.p
  rd.med
}
rd.cms <- boot(data=joint_dt, statistic=get.rd.cms, R=5000)
boot.ci(rd.cms,type="norm")

######
joint_dt$yy <- ifelse(joint_dt$crc_yy == "X", 1, 0)
joint_dt$yn <- ifelse(joint_dt$crc_yn == "X", 1, 0)
joint_dt$ny <- ifelse(joint_dt$crc_ny == "X", 1, 0)
joint_dt$nn <- ifelse(joint_dt$crc_nn == "X", 1, 0)

hormone <- joint_dt %>% filter(poolarm == 1) 
nrow(hormone)
placebo <- joint_dt %>% filter(poolarm == 0) 
nrow(placebo)

# whi yes and medicare yes
sum(placebo$yy)
sum(hormone$yy)
# whi no medicare yes
sum(placebo$ny)
sum(hormone$ny)
# whi yes medicare no
sum(placebo$yn)
sum(hormone$yn)
# whi no medicare no
sum(placebo$nn)
sum(hormone$nn)

#calculate HR for all the participants in hormone trials
enrollment <- read.table("U:/ppt_rand_enroll_dates.dat", header = TRUE, fill = TRUE)
enrollment$enroll.char <- as.character(enrollment$ENRLDATE)
# paste 0 before ENRLDATE
enrollment$enroll0 <- sapply(enrollment$enroll.char, function(x) ifelse(nchar(x) == 7, paste0("0",x), x))
enrollment$start <- as.Date(enrollment$enroll0, format = "%m%d%Y")
head(enrollment)
# dem_ht hormone treatment indicator
# outcome outcome dataset
outcome.col <- outcome[,c("COLORECTAL", "COLORECTALDY","ID")]
dem_ht.ID <- dem_ht[,c("HRTARM","ID")]
ht_dat <- Reduce(function(x,y) merge(x=x,y=y,by='ID'), 
                 list(outcome.col, dem_ht.ID, enrollment))
head(ht_dat)
nrow(ht_dat)
# end date for the two trials
ht_dat$end <- ifelse(ht_dat$HRTARM == "1" | ht_dat$HRTARM == "2", "2004-03-31", "2002-07-31")
ht_dat$followup <- as.Date(ht_dat$end) - as.Date(ht_dat$start)
ht_dat$time_to_event <- ifelse(ht_dat$COLORECTAL == 1, ht_dat$COLORECTALDY, ht_dat$followup)
# if colorectal diagnosed after E+P study ends, censor them
ht_dat$time_to_cancer <- ifelse(ht_dat$time_to_event >= ht_dat$followup, ht_dat$followup, ht_dat$time_to_event)
ht_dat$adj_colo <- ifelse(ht_dat$time_to_event > ht_dat$followup, 0, ht_dat$COLORECTAL)
# pool the hormone therapy arms
ht_dat$poolarm <- ifelse(ht_dat$HRTARM == "1" | ht_dat$HRTARM == "3", 1, 0)
# hr for all
hr.full <- coxph(Surv(time_to_cancer, adj_colo) ~ poolarm, data = ht_dat)
summary(hr.full)
library(dplyr)
EP <- ht_dat %>% filter(HRTARM == "3" | HRTARM == "4")
EP$arm <- ifelse(EP$HRTARM == "3", 1,0)
hr.ep <- coxph(Surv(time_to_cancer, adj_colo) ~ arm, data = EP)
summary(hr.ep)
