# COVID-19 HGI – clinical value project analysis plan 

This contains the analysis description for four sets of analyses for COVID-19 clinical value project.
If you have any question, please contact `tomoko.nakanishi`at`mail.mcgill.ca`.
Although this code is based on R language, we do not restrict each cohort to run the analysis on python or any other languages.
Please do the following analysis both using **ALL COVID+ individuals and ONLY hospitalized individuals.**

1. Survival analysis for chr3 variant
2. Associations between COVID-19 severity/complications and chr3 variant.



## 1. Survival analysis for chr3 variant
### 1-1.	Cox proportional hazards model for all-cause mortality from the date of diagnosis.

Survival analysis adjusting for age, sex, age*sex, study, and genetic PCs1:10.

```{r}
library(coxph)
library(tidyr)
library(dplyr)
data_EUR <- data_EUR %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))
cox <- coxph(Surv(time, death) ~ snp + sex + age_at_diagnosis*sex + age2 + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data_EUR)

#beta for snp
summary(cox)$coefficients[1,1]
#se for snp
summary(cox)$coefficients[1,3]
#pvalue for snp
summary(cox)$coefficients[1,5]

```
Here, data_EUR is a data frame which contains the following variables.

`time` Days from the date of COVID-19 diagnosis. (if missing, use the date of hospitalization)

`death` all–cause death

`snp` carrier status of chr3:45823240:T:C_C allele (rs10490770)

`age_at_diagnosis` age

`age2` age squared (age^2)

`study` If only one study/cohort included, please ignore.

`PC` genetic PCs

### 1-2. Competitive risk model for covid-19 related mortality.

Competitive risk model adjusting for age, sex, age*sex, study, and genetic PCs1:10

covid-19 related death was defined as doctor diagnosed or death with cause of death as ICD10 codes of `U71` or `U72` .

```{r}
library(cmprsk)
#assign covid-19 related death as 2
final <- final %>% mutate(death = ifelse(cause_of_death == 1 & death == 1, 2, death))
covs1 <- model.matrix(~ snp + sex*age_at_diagnosis + sex + age2 + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data_EUR)[, -1]
shr_fit <- 
  crr(
    ftime = data_EUR$time,
    fstatus = data_EUR$death,
    cov1 = covs1,
    cencode = 0, failcode = 2
  )
res <- summary.crr(shr_fit, conf.int = 0.95)
#beta for snp
res$coef[1,1]
#se for snp
res$coef[1,3]
#pvalue for snp
res$coef[1,5]
```

* Please consider removing some of the covariates espectially when the number of event (death) is small.
* Please run the analysis if total carrier count is greater than 10 in your subgroup. Otherwise, please provide us only the counts of  death_carrier, surviver_carrier, surviver_noncarrier, death_noncarrier as described below.
* First, stratify by genetically determined ancestry (continental-wise, EUR, AMR, SAS, EAS, AFR) and perform the above analysis seperately.
* Consider removing relatives, since we could not account for relatedness in this analysis.


### 1-3. File format for submission

Please submit the statistics file in `tsv` format.
File name should be `{Cohort}_{Date}_survival.tsv`.

example
`UKB_20200110_survival.tsv`


| study |  beta |  se  | pvalue | pop | type | death_carrier | surviver_carrier | surviver_noncarrier | death_noncarrier | covariates_used | agegroup |
----|----|----|----|----|----|----|----|----|----|----|----|
| UKB |  0.2090753 | 0.1190136 | 0.07896314 | EUR  |  all_cause | 91 | 523 | 2439 | 388 | age,age2,sex,age*sex,PC1:10,study | ALL |
| UKB |  0.9862004 | 0.6539996 | 0.13156667 | SAS  |  all_cause | 6 | 51 | 60 | 4 | age,sex | age≤60 |
| UKB |  0.4004562 | 0.5081928 | 0.43069645 | AMR  |  all_cause | 5 | 16 | 126 | 18 | age,sex | age>60 | 
| UKB |  NA | NA | NA | AFR |  all_cause | 0 | 2 | 162 | 44 | NA | ALL |
| UKB |  NA | NA | NA | EAS |  all_cause | 0 | 1 | 52 | 2 | NA | age≤60 |
| UKB |  0.3258524 | 0.1297687 | 0.01200000 | EUR |  covid_related | 78 | 370 | 1932 | 309 | age,age2,sex,age*sex,PC1:10,study | age>60 |
| UKB |  1.3078685 | 0.6820407 | 0.05500000 | SAS |  covid_related | 6 | 41 | 57 | 3 | age,sex | ALL |
| UKB |  NA | NA | NA | AMR |  covid_related | 1 | 5 | 36 | 0 | age,sex | age≤60 | 
| UKB |  NA | NA | NA | AFR |  covid_related | 0 | 1 | 146 | 33 | age,sex | age>60 | 
| UKB |  NA | NA | NA | EAS |  covid_related | 0 | 1 | 45 | 0 | NA | ALL |


## 2. Associations between COVID-19 severity/complications and chr3 variant.

### 2-1. Analysis details.

We would test the association between chr3:45823240:T:C_C allele dosage (rs10490770) and the following eight binary outcomes.  

All of these event are only counted if these occured within 30 days from the date of diagnosis (if missing use date of hospitalization instead).

* severity 
1. `icu_admit` ICU admission. **For ALL COVID+ analysis, code those with non-hospitalized as 0, those with hospitalization but without ICU admission as missing.
2. `high_who_score` highest WHO score >= 6 according to the score below.

```
-1 unknown
0 Uninfected, no viral RNA detected
1 Asymptomatic, viral RNA detected
2 Symtomatic, independent
3 Symptomatic, assistance needed
4 Hospitalized; no oxygen therapy
5 Hospitalised; oxygen by mask or nasal prongs
6 Hospitalized; oxygen by NIV or high flow
7 Intubation and mechanical ventilation, P/F ≥ 150 or SpO2/FiO2 ≥ 200
8 Mechanical ventilation PF < 150 (SpO2/FiO2 < 200) or vasopressors
9 Mechanical ventilation PF < 150 and vasopressors, dialysis or ECMO
10 Dead
```

* complications
3. `resp_severe` Need for mechanical ventilation (including oxygen by NIV or high flow) or ICD-10 codes of following (`J80`,`J9600`,`J9609`,`Z991`) or OPCS4 code of following (`E851`,`E852`). **For ALL COVID+ analysis, code those with no-oxygen as 0, those with only oxygen supplement as missing.**
4. `aki` Renal complication: doctor-diagnosed acute renal injury (AKI), highest creatinine > 1.5xULN, or ICD-10 codes of AKI (`N17*`). 
5. `hepatic` Hepatic complication: doctor-diagnosed hepatic complications, highest ALT > 3x upper limit of normal (ULN), or ICD-10 codes of acute hepatic failure (`K720`)
6. `cardiovascular` Cardiovascular complication: doctor-diagnosed acute myocardial infarction (AMI) or stroke, highest troponin T or troponin I > ULN, or ICD-10 codes of AMI (`I21*`) or stroke (`I61`,`I62`, `I63`, `I64`,`I65`,`I66*`)
7. `vte` doctor-diagnosed venous thromboembolism (VTE: pulmonary embolism or deep venous thromboembolism), or ICD-10 codes of VTE (`I81`, `I82*`, `I26*`)

As a comparison, it would be great to have association analyses with `AgeGroup`, `Sex`, `smoking status`, `BMI` if possible.
Please use following groupings.

```{r}
final <- final %>% mutate(Agegroup = case_when(age_at_diagnosis < 40 ~ "Age<40",
                                               age_at_diagnosis >= 40 & age_at_diagnosis < 50 ~ "Age40-50",
                                               age_at_diagnosis >= 50 & age_at_diagnosis < 60 ~ "0Age50-60",
                                               age_at_diagnosis >= 60 & age_at_diagnosis < 70 ~ "Age60-70",
                                               age_at_diagnosis >= 70 & age_at_diagnosis < 80 ~ "Age70-80",
                                               age_at_diagnosis >= 80 ~ "Age>80"),
                          Male = ifelse(sex == 0, 1, 0))

final <- final %>% mutate(smoke = case_when(smoking == 0 ~ "current",
                                            smoking == 1 ~ "ex",
                                            smoking == 2 ~ "0none",
                                            TRUE ~ "NA"),
                          BMI = case_when(weight/((height/100)^2) >= 40 ~ "class4",
                                          weight/((height/100)^2) < 40 & weight/((height/100)^2) >=35 ~ "class3",
                                          weight/((height/100)^2) < 35 & weight/((height/100)^2) >=30 ~ "class2",
                                          weight/((height/100)^2) < 30 ~ "0none",
                                          TRUE ~ "NA"))
```



example code 

```{r}
#association analysis for i th outcome and j th snp (here chr3:45823240:T:C_C allele dosage)
LM <- glm(outcome[i]," ~ `",snps[j],"` + age_at_diagnosis + age2 + sex + age_at_diagnosis*sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=data_EUR, family ="binomial")

#beta for snp
summary(LM)$coefficient[2,1]
#se for snp
summary(LM)$coefficient[2,2]
#pvalue for snp
summary(LM)$coefficient[2,4]



```

### 2-2. File format for submission

Please submit the statistics file in `tsv` format.
File name should be `{Cohort}_{Date}_complications.tsv`.

example
`UKB_20200110_complications.tsv`

| study |  outcome | beta |  se  | pvalue | pop | N_case | N_control | risk_factor | covariates_used | hospitalized_only
|----|----|----|----|----|----|----|----|----|----|----|
| UKB |  resp_severe | 0.2090753 | 0.1190136 | 0.07896314 | EUR | 270 | 607 | chr3:45823240:T:C_C | age,age2,sex,age*sex,PC1:10,study | yes |
| UKB |  resp_severe | 0.49685214 | 0.2533789 | 0.0498899920 | EUR | 578 | 593  | chr3:45823240:T:C_C | age,age2,sex,age*sex,PC1:10,study | no |
| UKB |  resp_severe | 0.49685214 | 0.2533789 | 0.0498899920 | EUR | 578 | 593  | past_smoker | age,age2,sex,age*sex,PC1:10,study | yes |

* Please describe the snp as following coding `chr3:45823240:T:C_C` (build 38) and make sure you test the association with the dosage of `C` allele.



