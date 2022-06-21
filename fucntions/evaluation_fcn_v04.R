create_best_report = function(err, loss, likelihood, penalty, use_penalty = T){
#Evaluate results stored in err, loss, likelihood, penalty
#Finds rule with lowest loss in each simulation (row)
#Creates data frame with rule, err, loss, likelihood, penalty for best rule on each simulation
  n_sim = nrow(err) #number of simulations - merged batches
  ######pick the best in each simulation##########
  if(use_penalty == T){
    best_rule_index = apply(loss, 1, FUN = which.min) #which rule was the best in each row
  }else{ #if penalty is not used, just likelihood is used to pick the rule
    best_rule_index = apply(likelihood, 1, FUN = which.min) #which rule was the best in each row
  }
  tmp_init = matrix(NA, nrow = n_sim)
  
  best_report = data.frame("rule" = "", "err" = tmp_init, "loss" = tmp_init, "likelihood"=tmp_init, "penalty"=tmp_init, "leaf_cnt" = tmp_init, stringsAsFactors = F)
  
  
  for(i_sim in 1:n_sim){
    best_report[i_sim,"rule"] =paste(rules_lst[[best_rule_index[i_sim]]], collapse = " , " ) 
    best_report[i_sim,"err"] = err[i_sim,best_rule_index[i_sim]]^2
    best_report[i_sim,"loss"] = loss[i_sim,best_rule_index[i_sim]]
    best_report[i_sim,"likelihood"] = likelihood[i_sim,best_rule_index[i_sim]]
    best_report[i_sim,"penalty"] = penalty[i_sim,best_rule_index[i_sim]]
    best_report[i_sim,"leaf_cnt"] = length(rules_lst[[best_rule_index[i_sim]]])
  }
  
  return(best_report)
}


create_mean_report = function(err, loss, likelihood, penalty){
#creates report for average errors and losses of all rules
#ordered by mean error of the rules ascending (top are best performing rules)
#Also trim mean is created if limited amount of simulations is performed
  not.na.flt = !is.na(err[,1]) #just for sure?
  tmp_init = matrix(NA,nrow = n_rules_lst)
  res_report = data.frame("rule" = "", "mean_err" = tmp_init, "mean_loss" = tmp_init, "mean_likelihood"=tmp_init, "mean_penalty"=tmp_init, "mean_err_trim" = tmp_init, "leaf_cnt" = tmp_init , "orig_index" = tmp_init, stringsAsFactors =F)
  res_report$mean_err = apply((err[not.na.flt,])^2,2,mean)^sqrt_err
  res_report$mean_loss = apply(loss[not.na.flt,],2,mean)
  res_report$mean_likelihood = apply(likelihood[not.na.flt,],2,mean)
  res_report$mean_penalty = apply(penalty[not.na.flt,],2,mean) #res_report$mean_loss - res_report$mean_likelihood
  res_report$mean_err_trim = apply((err[not.na.flt,])^2,2,mean, trim = 0.05)^sqrt_err
  res_report$orig_index = seq(n_rules_lst) #save original indexing
  res_report$leaf_cnt = lengths(rules_lst) #number of leafs
  #fill in rules
  for(i_rules in 1:n_rules_lst){
    res_report[i_rules, "rule"] =paste(rules_lst[[i_rules]], collapse = " & " )
  }
  
  #add relative error to no partitioning error
  
  
  #order res_report by real ERROR
  res_report <- res_report[order(res_report$mean_err),]
  return(res_report)
}

create_top_res = function(mean_report, best_means, top_incl = 5){
#Merges mean results for algorithm selecting rule with best (min) loss and top n results of all rules
#First row is always best rule algorithm, other n rows are ordered (just a copy of res_report)

  top_res = mean_report[1:top_incl,] #filter top n results
  
  top_res = rbind(top_res[1,], top_res) #duplicate first row and replace with best partitioning
  top_res[1,"rule"] = "best rule"
  top_res[1,"orig_index"] = NA
  top_res[1,c("mean_err", "mean_loss", "mean_likelihood", "mean_penalty", "mean_err_trim",  "leaf_cnt")] = best_means[c("err","loss","likelihood","penalty","err_trim",  "leaf_cnt")]
  top_res[,"order"] = rank(top_res$mean_err)
  return(top_res)
}


create_selected_res = function(mean_report, best_means, selected_ind, selected_labels, err){
  #Merges mean results for algorithm selecting rule with best (min) loss and top n results of all rules
  #First row is always best rule algorithm, other n rows are ordered (just a copy of res_report)
  
  top_res = mean_report[mean_report$orig_index%in%selected_ind,] #filter top n results
  
  top_res = rbind(top_res[1,], top_res) #duplicate first row and replace with best partitioning
  top_res[1,"rule"] = "best rule"
  top_res[1,"orig_index"] = NA
  top_res[1,c("mean_err", "mean_loss", "mean_likelihood", "mean_penalty", "leaf_cnt")] = best_means[c("err","loss","likelihood","penalty", "leaf_cnt")]
  top_res[,"order"] = rank(top_res$mean_err)
  
  tmp = cbind(best_report$err, m_err[,selected_ind]^2)
  colnames(tmp) = selected_labels
  boxplot(tmp)
  summary(tmp)
  
  return(top_res)
}
