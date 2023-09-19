####################### - Fit RBC model to data - #############################
## load trial data
load('Data/RBC_model_data.RData')
source('master_functions.R')

### make all datasets
dat_stan_list = list()
dat_stan_list[[1]] = make_stan_dataset(my_data = PQdat) # individuals recruited into both have same fundamental parameters
dat_stan_list[[2]] = make_stan_dataset(my_data = PQdat,ID_subject = 'ID2') # individuals recruited into both have different fundamental parameters
# leave one out
for(kk in 1:length(unique(PQdat$ID2))){
  id2 = unique(PQdat$ID2)[kk]
  PQ_fit = PQdat %>% filter(ID2 != id2)
  PQ_pred = PQdat %>% filter(ID2 == id2)
  
  dat_stan_list[[kk+2]] = 
    make_stan_dataset(my_data = PQ_fit,
                      ID_subject = 'ID2',
                      data_pred = PQ_pred
                      ) 
}
save(dat_stan_list, file = 'Rout/stan_data_list.RData')