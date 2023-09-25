####################### - Fit RBC model to data - #############################
## load trial data
library(tidyverse)
load('Data/RBC_model_data.RData')
source('master_functions.R')

### make all datasets
dat_stan_list = list()
dat_stan_list[[1]] = make_stan_dataset(my_data = PQdat) # individuals recruited into both have same fundamental parameters
dat_stan_list[[2]] = make_stan_dataset(my_data = PQdat,ID_subject = 'ID2') # individuals recruited into both have different fundamental parameters

### Only ascending
PQ_ascending = PQdat %>% filter(study=='Part1')
PQ_pred = PQdat%>% filter(ID=='PQ45 13')
dat_stan_list[[3]] = make_stan_dataset(my_data = PQ_ascending,data_pred = PQ_pred) # individuals recruited into both have same fundamental parameters


### Only single dose
PQ_single = PQdat %>% filter(study=='Part2')
PQ_pred = PQdat%>% filter(study=='Part1', ID=='ADPQ 22')
dat_stan_list[[4]] = make_stan_dataset(my_data = PQ_single, data_pred = PQ_pred) # individuals recruited into both have same fundamental parameters

### Small dataset for testing
PQ_test = PQdat %>% filter(study=='Part1',
                           ID2 %in% c("ADPQ 1",
                                      "ADPQ 10",
                                      "ADPQ 11", 
                                      "ADPQ 12",
                                      "ADPQ 13",
                                      "ADPQ 14",
                                      "ADPQ 15", 
                                      "ADPQ 16" ,
                                      "ADPQ 17"))
PQ_test_pred = PQdat %>% filter(study=='Part1', ID2=='ADPQ 20')

dat_stan_list[[5]] = make_stan_dataset(my_data = PQ_test, data_pred = PQ_test) # individuals recruited into both have same fundamental parameters

dat_stan_list[[5]]$drug_regimen_pred = rep(0.5, 50)
dat_stan_list[[5]]$N_pred = 50


### Single individual for testing
PQ_test = PQdat %>% filter(study=='Part1',
                           ID2 %in% c("ADPQ 13"))
PQ_test_pred = PQdat %>% filter(study=='Part1', ID2=='ADPQ 20')

dat_stan_list[[6]] = make_stan_dataset(my_data = PQ_test, data_pred = PQ_test) # individuals recruited into both have same fundamental parameters

dat_stan_list[[6]]$drug_regimen_pred = rep(0.5, 50)
dat_stan_list[[6]]$N_pred = 50


K_add=length(dat_stan_list)
# leave one out
for(kk in 1:length(unique(PQdat$ID2))){
  id2 = unique(PQdat$ID2)[kk]
  PQ_fit = PQdat %>% filter(ID2 != id2)
  PQ_pred = PQdat %>% filter(ID2 == id2)
  
  dat_stan_list[[kk+K_add]] = 
    make_stan_dataset(my_data = PQ_fit,
                      ID_subject = 'ID2',
                      data_pred = PQ_pred
                      ) 
}
save(dat_stan_list, file = 'Rout/stan_data_list.RData')