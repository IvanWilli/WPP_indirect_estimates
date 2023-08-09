library(tidyverse)
library(DemoTools)
library(lubridate)
library(data.table)
options(tibble.print_min=50)
source("R/intercensal_survival.R")

# replicate some publications ---------------------------------------------------------------

# match: PANAMA 1960-1970 from UNX (Table 174)
panama_match <- intercensal_survival(
  c1 = c(90071,76598,63635,54431,45634,37818,32179,28724,23974,20618,15068,11999,10283,6737,5242,6756),
  c2 = c(114017,106944,85253,73381,63010,50924,40885,36115,29409,25360,21775,17632,13004,10061,6690,9873),
  date1 = "1960/12/11",
  date2 = "1970/05/10",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "match")

# in UNX removes outliers from ages 10, 20, 35, 55 and 65 (implying the use of the initial OAG) 
# and get the mean 16.1 which is between 57.5 and 60 of e0.
# by default the function only removes age 10 because nSx and give a result a higher one.
# But is possible to reconstruct the result using the output where is found the best level for each age:
panama_match$rank_match %>% mutate(color = ifelse(age %in% c(10, 20, 35, 55, 65), "outlier", "used")) %>% 
  ggplot(aes(age, e0_interp, color=color)) + geom_point(size = 2)
panama_match$rank_match %>% filter(!age %in% c(10, 20, 35, 55, 65)) %>% summarise(mean(e0_interp))
# to replicate this during the process, parameters ages_fit or e0_range_accept should be used
intercensal_survival(
  c1 = c(90071,76598,63635,54431,45634,37818,32179,28724,23974,20618,15068,11999,10283,6737,5242,6756),
  c2 = c(114017,106944,85253,73381,63010,50924,40885,36115,29409,25360,21775,17632,13004,10061,6690,9873),
  date1 = "1960/12/11",
  date2 = "1970/05/10",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "match",
  ages_fit = seq(0,75,5)[!seq(0,75,5) %in% c(10, 20, 35, 55, 65)])$selected_level

# projection. Indonesia Males 1980-1990, SA pattern (Preston, 2001; Box 11.4) - OK
india_fproj <- intercensal_survival(
  c1 = c(10815974, 10832383, 9131871, 7512541, 5978576, 5612684, 4022625, 4190944, 3644053, 3012756, 2717883, 1720501, 1559230, 811113, 689074, 688422),
  c2 = c(10760859,11928095,11044127,9520440,7583305,7457150,6584325,5788441,4010254,3723922,3289190,2321621,2219069,1329162,945876,867636),
  date1 = "1980/07/01",
  date2 = "1990/07/01",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "m",
  mlt_family = "UN_South_Asian",
  method = "fproj",
  e0_accept = c(51,57),
  ages_fit = seq(0,40,5))

# feeney: japan females 1960-1970. OK - level at age 5, 68.7 (page 27)
# level 72.5 has ~ e80=6.06. It´s not exact but very accurate
MLTlookup %>% filter(age == 80, type == "CD_West", sex == 2)
japan_feeney <- intercensal_survival(
  c1 <- c(3831870,4502304,5397061,4630775,4193184,4114704,3770907,3274822,
          2744786,2559755,2160716,1839025,1494043,1133409,870238,577972,313781, 131547),
  c2 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
          3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191, 53116),
  date1 = "1960/08/18",
  date2 = "1970/08/18",
  age1 = seq(0,85,5),
  age2 = seq(0,85,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "feeney",
  ages_fit = 5,
  mlt_e0_logit_feeney = 72.5)
japan_feeney$rank_match

# Preston (1983) logit - r. Indian Females 1961-1971 - OK e5 = 51
load("data/all_censuses.rda")
all_censuses %>% filter(LocName == "India", TimeLabel %in% c(1961,1971)) %>% distinct(TimeMid)
c1 <- all_censuses %>% 
  filter(LocName == "India", TimeLabel == 1961, SexName == "Female") %>% 
  select(AgeStart, AgeSpan, DataValue)
c2 <- all_censuses %>% 
  filter(LocName == "India", TimeLabel == 1971, SexName == "Female") %>% 
  select(AgeStart, AgeSpan, DataValue) %>% 
  mutate(AgeStart = pmin(AgeStart, 70)) %>% 
  group_by(AgeStart) %>% summarise(DataValue = sum(DataValue))
intercensal_surv_Preston_logit_var_r(c1 = c1$DataValue,
                                     c2 = c2$DataValue,
                                     date1 = 1961.164,
                                     date2 = 1971.249,
                                     age1 = seq(0,70,5),
                                     age2 = seq(0,70,5),
                                     sex = "f",
                                     mlt_family = "CD_South",
                                     mlt_e0 = 50,
                                     ages_fit = seq(5,65,5),
                                     q0_5 = 1 - .776)

# variable-r (Preston-Bennet): PANAMA 1960-1970 from UNX (Table 186). - OK
# It´s picked ages_fit 10:30 as Manual choose, because growing too much below that
# promedio e(10) 17.1 ~ 60.3, en comparación con 60.1
panama_var_r <- intercensal_survival(
                    c1 = c(90071,76598,63635,54431,45634,37818,32179,28724,23974,20618,15068,11999,10283,6737,5242,6756),
                    c2 = c(114017,106944,85253,73381,63010,50924,40885,36115,29409,25360,21775,17632,13004,10061,6690,9873),
                    date1 = "1960/12/11",
                    date2 = "1970/05/10",
                    age1 = seq(0,75,5),
                    age2 = seq(0,75,5),
                    sex = "f",
                    mlt_family = "CD_West",
                    method = "var-r",
                    ages_fit = 10:30)
panama_var_r$selected_level$mean_e0

# others ------------------------------------------------------------------

# feeney: zimbawe 1982-1992 females. level 65 has ~ e70=9.7. BUT He uses Brass MLT in Anex II.3
zimbawe_feeney <- intercensal_survival(
  c1 <- c(666513,620383,519647,413331,364837,281551,207121,170467,139774,110583,
          91039,60906,65374,38928,30553,46842),
  c2 <- c(798430,835296,734331,634658,524836,377773,327407,260436,190152,143928,
          147839,87023,84499,51075,62691,68635),
  date1 = "1982/08/18",
  date2 = "1992/08/18",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "feeney",
  ages_fit = 5,
  mlt_e0_logit_feeney = 65)

#projection - level 17.4 (e0=61) - OK
panama_match <- intercensal_survival(
  c1 = c(90071,76598,63635,54431,45634,37818,32179,28724,23974,20618,15068,11999,10283,6737,5242,6756),
  c2 = c(114017,106944,85253,73381,63010,50924,40885,36115,29409,25360,21775,17632,13004,10061,6690,9873),
  date1 = "1960/12/11",
  date2 = "1970/05/10",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "fproj",
  ages_fit = c(5, 10, 15, 20))
