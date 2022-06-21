library(stringr)
library(ChainLadder)
library(partitions)
library(fastDummies)
library(data.tree)
library(plyr) #for mapvalues used in fact_enc (encoding of numerical variables)

rm(list = ls()) #clear
setwd("<workdir>")
source("config.r")

#######################
###initial processes###
#######################
#load functions
source(paste0(wd_functions, "all_rules_fcn_v01.R"))
source(paste0(wd_functions, "optim_tri_partition_functions_v11.R"))
source(paste0(wd_functions, "triangle_fcn_v01.R"))
source(paste0(wd_functions, "evaluation_fcn_v04.R"))
source(paste0(wd_functions, "tree_build_fnc_v07.R"))#for greedy specific

#load data
dat = read.csv(paste0(wd_data, input_data_file))
dat$dev_yr = as.factor(dat$dev_yr) #dev_yr are factors 

#encoding categorical variables for greedy -> creates indicators for all categories as well as its combinations
if(set_algorithm == "greedy"){
  split_variables_all =c()
  for(i_catvar in use_set_var){ #for all variables to be used for optimization
    if(i_catvar%in%categ_var){ #if variable is categorical
      dsgn = comb_catvar(catvar = i_catvar, dat = dat)
      dat = cbind(dat,dsgn)
      rm(dsgn)
      split_variables_all =c(split_variables_all, colnames(dat)[str_detect(string = colnames(dat), pattern = paste0("^",i_catvar))])
    }#if categorical
    #drop original categorical variables
    split_variables_all = split_variables_all[!(split_variables_all %in% categ_var)]
  }#for cat vars
  for(i_catvar in use_set_var){ #for all variables to be used for optimization
    if(i_catvar%in%categ_var_enc){#if categ var is to be encoded as a metric variable
      new_var_name = paste0(i_catvar,"_enc")
      #encoding for metric variables (used only for greedy)
      num_var_enc = fact_enc_plain(enc_var = i_catvar, dat.node = dat)
      dat[new_var_name] = num_var_enc
      split_variables_all = c(split_variables_all, new_var_name)
    }
  }
  for(i_catvar in use_set_var){ #for all variables to be used for optimization
    if(i_catvar%in%metric_var){#add all numeric (no encoding is performed)
      split_variables_all = c(split_variables_all, i_catvar)
    }
  }
  

} #if greedy


#encoding categorical variable for exhaustive search -> creates indicators for all categories (one hot encoding)
if(set_algorithm == "exhaustive"){
  dat <- dummy_cols(dat, select_columns = categ_var_exhaustive)
  split_variables_all =colnames(dat)[str_detect(string = colnames(dat), pattern = paste0(categ_var_exhaustive,"_[0-9]"))]
}


#for exhaustive
if(set_algorithm  == "exhaustive"){

  #setup rules list
  rules_lst = create_all_rules(base_set = split_variables_all) #creates rules for all possible partitions


  n_rules_lst = length(rules_lst)
  rule_elements = extract_rule_elements(rules_lst)
  
  
  #set up result containers
  n_sim = 1 #internal
  i_sim = 1 #internal
  err = matrix(NA, nrow = n_sim, ncol = n_rules_lst)
  loss = err
  likelihood = err
  penalty = err
  
  #initiate the algorithm
  sig_fixed = getDeviance(rules = "T", tri_data = dat, ignore_dev=ignore_dev, sig_fixed =  NULL, calc_penalty=F,  precalc = F, rule_elements = NA)$sig_fixed 
  pre_triangles = precalc_triangles(in.d=dat, rule_elements)
  
  #evaluate errors for all rules

  for(i_rules in 1:n_rules_lst){
    rules = rules_lst[[i_rules]] #extract set of rules
    node_list = seq(length(rules)) #fake labels
    if(calc_future == T){  #calculate future error
      err[1, i_rules] =  getTriErr(rules,node_list, dat, size_var)[1]
    }
    #res = getDeviance(rules = rules, tri_data = dat.node_sampled, ignore_dev=ignore_dev, sig_fixed = sig_fixed, calc_penalty=T,  precalc = F, rule_elements = NA)
    res = getDeviance(rules = rules, tri_data = pre_triangles, ignore_dev=ignore_dev, sig_fixed = sig_fixed, calc_penalty=T,  precalc = T, rule_elements = rule_elements)
    loss[1, i_rules] = res$crit
    likelihood[1, i_rules] = res$likelihood
    penalty[1, i_rules] = res$penalty_term
  }

  #Evaluate results stored in err, loss, likelihood, penalty
  best_report = create_best_report(err = err, loss= loss, likelihood=likelihood, penalty = penalty, use_penalty = T)
  print('Best partitioning using exhaustive search:')
  print(best_report)
} #if exhaustive


#for greedy
if(set_algorithm  == "greedy"){
  #initiate 
  skip_to_next <- FALSE
  
  
  tbl_tree = treeBuilder(dat.node.orig = dat, min_volume, ignore_dev, c(metric_var, paste0(categ_var_enc, "_enc")), split_variables_all,prob_grid)
  print(tbl_tree)
  plot_tree(tbl_tree)
  #####################
  ###evaluate errors###
  #####################
  rules = tbl_tree[tbl_tree$action == "leaf","rule"]
  node_list = tbl_tree[tbl_tree$action == "leaf","node"]
  
  #extract penalty
  loss_output = data.frame()
  loss_output[1,"rule"] = paste(rules, collapse ="," )
  min_loss = min(tbl_tree[tbl_tree$action == "leaf","node_loss"] )
  loss_output[1,c("loss", "likelihood", "penalty")] = tbl_tree[tbl_tree$action == "internal" &tbl_tree$split_loss == min_loss,c("split_loss", "likelihood", "penalty")] 
  
  
  
  
  if(calc_future == T){  #calculate future error
    loss_output["err"] = (getTriErr(rules,node_list, dat, size_var)[1])^2 #future IBNR error
    loss_output = loss_output[,c(1,5,2,3,4)]
  }
  print(loss_output)
}
  
