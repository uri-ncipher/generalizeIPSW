# generalizeIPSW

#Function that computes the IPSW estimator, its variance and 95% CI
#data: combined trial (already subsetted so outcome data complete)
#and cohort data with indicator for trial participation (s), covariates, outcome and predictors
#nt: size of target population
#selvar: vector of variables for sampling score model in combined trial and cohort data
#trt: treatment or exposure variable in the trial
#outcome: outcome variable in the trial
