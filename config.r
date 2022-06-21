#############
####config###
#############
#input file
wd_data = "./data/"  
wd_functions = "./fucntions/"  
input_data_file = "sample_data.csv"

#algorithm and loss function
set_algorithm  = "exhaustive" #"exhaustive" #"greedy"

size_var = "paid" #name of the financial amount variable for reserving

#For exhaustive
#Only categorical are allowed.
#Only one variable is allowed now, if combination of several variables is required, combinations must be prepared in database.
#Indicators of all possible partitionings for each variable one hot encoded 
if(set_algorithm == "exhaustive"){
  categ_var_exhaustive = "LoBcnt" 
}

#For greedy algorithm
#assign categorical and metric variables
if(set_algorithm == "greedy"){
  metric_var = c("age")#list of metric vars, no encoding is done - only used for greedy
  categ_var_enc = c("inj_part", "cc") #list of categorical vars to be "development factor encoded"- only used for greedy
  categ_var = c("LoBcnt") #list of categorical variables to be one hot encoded
}
#subset of variables to be used for optimization #e.g. c("LoBcnt", "cc", "inj_part", "age" ) 
use_set_var =  c("LoBcnt") #Variables from all categories may be used for greedy alg. #c("LoBcnt","inj_part","age")


trunc_penalty_mean = 0.00 #Truncation for penalty term mean. 
penalty_fun = "median" #or use median to estimate mean penalty if numerically unstable (trunc_penalty_mean is then redundant) If not defined, truncated mean is used.
variance_power = 1 #power parameter for variance function (Var = sigma * mu^variance_power)
ignore_dev = 2 #how many of the corner values should be omitted in triangles

prob_grid = seq(from = 0.2, to = 0.8, by = 0.1) #grid for numerical variables for greedy algorithm
min_volume = 0.0000000001 #stopping criterion do not split if node volume is smaller
calc_penalty = T #Include also penalty term?
calc_future = T #Calculate also error for future predictions. Lower right triangle must be included in the data.