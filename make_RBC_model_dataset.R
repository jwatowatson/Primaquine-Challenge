library(tidyverse)

## load trial data
load('Data/PQdata.rds')
PQdat$Day = as.numeric(PQdat$Day)

PQdat = PQdat %>% filter(label != 'ADPQ 18', # subject 18 not included in analysis
                         Day >= 0, Day <= 30) %>%
  group_by(label, Day) %>%
  mutate(Study_Day = unique(Day),
         Haemocue_hb = mean(Haemocue_hb,na.rm=T),
         CBC_hb = mean(CBC_hb, na.rm = T),
         CBC_retic = mean(CBC_retic, na.rm=T),
         Manual_retic = mean(Manual_retic, na.rm = T),
         dose = sum(PQdose, na.rm = T)) %>%
  ungroup()%>%
  distinct(label, Day, .keep_all = T) %>%
  arrange(label, Day) %>% group_by(label, Day) %>%
  mutate(dosemgkg = dose/weight,
         Mean_retic = mean(c(Manual_retic, CBC_retic), na.rm = T),
         Mean_hb = mean(c(Haemocue_hb, CBC_hb), na.rm = T))
PQdat$ID = ifelse(!is.na(PQdat$Ascending_ID), PQdat$Ascending_ID, PQdat$label)
PQdat$ID2 = PQdat$label

cols = c('ID', 'ID2', 'study', 'Study_Day','G6PD_variant','Haemocue_hb', 
         'CBC_hb', 'CBC_retic', 'Manual_retic','Mean_hb','Mean_retic',
         'dose', 'dosemgkg', 'weight')

PQdat = PQdat[, cols]
PQdat$CBC_hb = ifelse(is.nan(PQdat$CBC_hb), NA, PQdat$CBC_hb)
PQdat$CBC_retic = ifelse(is.nan(PQdat$CBC_retic), NA, PQdat$CBC_retic)
PQdat$Mean_retic = ifelse(is.nan(PQdat$Mean_retic), NA, PQdat$Mean_retic)


writeLines('Unique individuals:')
unique(sort(PQdat$ID))

writeLines('Unique IDs with repeats:')
unique(sort(PQdat$ID2))


save(PQdat, file = 'Data/RBC_model_data.RData')


## Setup models


