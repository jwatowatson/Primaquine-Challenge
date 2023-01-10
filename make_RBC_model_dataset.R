## load trial data
Ascending_dat = read.csv('PQAscending.csv')
Ascending_dat = Ascending_dat[!Ascending_dat$label=='PQ-mrdt01-018', ] # exclude subject 18
Ascending_dat$ID_combined = Ascending_dat$label

PQ45_dat = read.csv('PQ45.csv')
PQ45_dat$ID_combined = PQ45_dat$Ascending_ID
PQ45_dat$ID_combined[is.na(PQ45_dat$ID_combined)]=PQ45_dat$label[is.na(PQ45_dat$ID_combined)]

cols = c('label', 'study', 'visit', 'Haemocue_hb','G6PD_variant',
         'CBC_hb','Mean_Retic','sddose','weight', 'ID_combined')

PQ_dat_all = rbind(Ascending_dat[, cols], PQ45_dat[, cols])
PQ_dat_all$ID = as.numeric(as.factor(PQ_dat_all$label))
PQ_dat_all$ID2 = as.numeric(as.factor(PQ_dat_all$ID_combined))
PQ_dat_all$Study_Day = PQ_dat_all$visit

ID_map = PQ_dat_all[!duplicated(PQ_dat_all$label), c('label','ID','ID2', 'ID_combined')]

all_na = is.na(PQ_dat_all$Haemocue_hb) & is.na(PQ_dat_all$CBC_hb) & is.na(PQ_dat_all$Mean_Retic)
PQ_dat_all = PQ_dat_all[!all_na, ]

writeLines('Unique IDs:')
unique(sort(PQ_dat_all$ID))

PQ_dat_all = PQ_dat_all %>% filter(Study_Day<30)
IDs = unique(PQ_dat_all$ID)

PQ_dat_all = dplyr::arrange(PQ_dat_all, ID, Study_Day)
PQ_dat_all$dosemgkg = PQ_dat_all$sddose/PQ_dat_all$weight


save(PQ_dat_all, file = 'RBC_model_data.RData')