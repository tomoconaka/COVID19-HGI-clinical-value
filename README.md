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

### submitting file format

