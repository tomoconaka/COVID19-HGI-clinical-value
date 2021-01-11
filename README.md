# COVID-19 HGI – clinical value project analysis plan 

## 1-1.	Cox proportional hazards model for all-cause mortality within 30 days from the date of diagnosis.

```{r}
library(coxph)
cox <- coxph(Surv(time, death) ~ snp + sex + age_at_diagnosis*sex + age2 + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = data_EUR)
```

`time` Days from the date of COVID-19 diagnosis. (if missing, use the date of hospitalization)

`death` all–cause death within 30 days from the starting date (COVID-19 diagnosis or the date of hospitalization)

`snp` carrier status of chr3:45823240:T:C_C allele (rs10490770)

```{r}
final <- final %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))
```

`age2` age squared (age^2)

`study` If only one study/cohort included, please ignore.

`PC` genetic PCs

* Please consider removing some of the covariates espectially when the number of event (death) is small.
* First, stratify by genetically determined ancestry (continental-wise, EUR, AMR, SAS, EAS, AFR) and perform the above analysis seperately.

## 1-2. Competitive risk model for covid-19 related mortality within 30 days from the date of diagnosis.



## 1.3 submitting file format


| study |  beta |  se  | pvalue | pop | type | death_carrier | surviver_carrier | surviver_noncarrier | death_noncarrier
----|----|----|----|----|----|----|----|----|----|
| UKB |  0.2090753 | 0.1190136 | 0.07896314 | EUR  |  all_cause | 91 | 523 | 2439 | 388 |
| UKB |  0.9862004 | 0.6539996 | 0.13156667 | SAS  |  all_cause | 6 | 51 | 60 | 4 |
| UKB |  0.4004562 | 0.5081928 | 0.43069645 | AMR  |  all_cause | 5 | 16 | 126 | 18 |
| UKB |  NA | NA | NA | AFR |  all_cause | 0 | 2 | 162 | 44 |
| UKB |  NA | NA | NA | EAS |  all_cause | 0 | 1 | 52 | 2 |




