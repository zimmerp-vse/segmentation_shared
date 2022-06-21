
treeBuilder = function(dat.node.orig, min_volume, ignore_dev, metric_var, split_variables_all,prob_grid){

nr_claims_orig = length(unique(dat.node.orig$id))
#dat.node = dat.node.orig 
tri_size = length(levels(dat.node.orig[,"dev_yr"])) #dev_yr are factors

tbl_tree = data.frame("node" = 1,
                      "rule" = "T", #T for root
                      "action" = "split",
                      "split_variable" = NA,
                      "node_loss" = NA,
                      "likelihood" = NA,
                      "penalty" = NA,
                      "parental_node" = "",
                      "thr" = NA,
                      "thr_quantile" = NA,
                      "level" = 1,
                      "split_loss"= NA,
                      stringsAsFactors = F
)



#results for the root
i_node = 1
res.total = getDeviance(rules = "T", tri_data = dat.node.orig, ignore_dev=ignore_dev, sig_fixed = NULL, calc_penalty) 
tbl_tree[i_node, "node_loss"] = res.total[["crit"]]
tbl_tree[i_node, "likelihood"] = res.total[["likelihood"]]
tbl_tree[i_node, "penalty"] = res.total[["penalty_term"]]
tbl_tree[i_node, "volume"] = 1 #100 % for root

#Set variables to split for the root node 
split_var_grid_lst = list()
split_var_grid_lst[[1]] = update_split_var(rule =tbl_tree[1,"rule"] , dat.node.orig, split_variables_all, metric_var = metric_var, tri_size, prob_grid) #parent rule is identiacal for each split variable  


total_volume = res.total[["tri_volume"]] #size of all to calculate percentage size of nodes 

sig_fixed = res.total$sig_fixed #set NULL to recalculate sigma each tim

perform_splits = T #initial


while(perform_splits){
  

  nodes_to_split = tbl_tree[tbl_tree[,"action"]=="split","node"]
  res_nodes = list() #save temproary results for all actual nodes to split
  crit_nodes = matrix(NA, nrow = length(nodes_to_split), ncol = 2) #save crit for each actual node
  crit_nodes[,1] = nodes_to_split
  colnames(crit_nodes) = c("node", "crit")
  
  for(i_node in nodes_to_split){ #for each node with action == "split", i_node is not coiunter but number of the node
    #split_var_grid_lst = update_split_var(rule =tbl_tree[i_node,"rule"] , dat.node.orig, split_variables_all, size_var = size_var, tri_size, prob_grid) #parent rule is identiacal for each split variable
    
    node_split_variables = split_var_grid_lst[[i_node]][["split_variables"]]
    #thr.grid = split_var_grid_lst[[i_node]][["thr.grid"]]
    #set up node optimization table
    #print(i_node)
    #print(node_split_variables)
    tbl_node = data.frame("node" = i_node,
                          "split_variable" = node_split_variables,
                          "thr" = NA,
                          "thr_quantile" = NA,
                          "loss" = NA,
                          "likelihood" = NA,
                          "penalty" = NA,
                          "parent_rule" = tbl_tree[i_node,"rule"],
                          "child_rule1" =NA,
                          "child_rule2" =NA,
                          "parent_volume" = tbl_tree[i_node,"volume"],
                          "child_volume1" =NA,
                          "child_volume2" =NA,
                          stringsAsFactors = F
    )

    for(split_var in node_split_variables){
      #print(c(i_node, split_var))
      if(split_var %in% metric_var){ #is split by size
        #set up grid to be optimized
        #thr.grid = quantile(dat.node[dat.node$dev_yr==tri_size,"paid"], c(seq(from = 0.35, to = 0.85, by = 0.05), seq(from = 0.87, to = 0.95, by = 0.01)))
        
        #thr.grid = quantile(dat.node[dat.node$dev_yr==tri_size,"paid"], c(seq(from = 0.6, to = 0.9, by = 0.05)))  
        #thr.grid = thr.grid[thr.grid>0]
        thr.grid = split_var_grid_lst[[i_node]][["thr_grid_lst"]][[split_var]]
        n.thr = length(thr.grid)
        #results matrix for grid
        LLRes = matrix(NA, nrow = n.thr)
        likelihood = matrix(NA, nrow = n.thr)
        penalty = matrix(NA, nrow = n.thr)
        volume = matrix(NA, nrow = n.thr, ncol = 2)
        for(i.grid in 1:n.thr){#for grid
          thr =thr.grid[i.grid] #i.grid=1
          #create rules for split by thr
          parent_rule = tbl_tree[i_node, "rule"] 
          new_rules = c(paste(parent_rule,paste(split_var, "<=", thr,sep=""),sep = "&"), paste(parent_rule,paste(split_var,">", thr,sep=""),sep = "&"))
          #rules_tmp = tbl_tree[nodes_to_split[nodes_to_split!=i_node] ,"rule"] #all other rules than the splited
          rules_tmp = tbl_tree[(tbl_tree$action == "leaf"|tbl_tree$action == "split")&tbl_tree$node!=i_node,"rule"] #all other rules than the splited that are to be splitted or leaf
          rules_tmp = c(rules_tmp, new_rules)
          #calculates results for the split
          res = getDeviance(rules = rules_tmp, tri_data = dat.node.orig, ignore_dev=ignore_dev, sig_fixed = sig_fixed, calc_penalty) #res.total[["sig_dev"]]
          LLRes[i.grid] =res[["crit"]]
          likelihood[i.grid] = res[["likelihood"]]
          penalty[i.grid] = res[["penalty_term"]]
          #save volumes for new_rules
          volume[i.grid,c(1,2)] =  res[["tri_volume"]][rules_tmp%in%new_rules]
          
        }#for i.grid
        
        #smooth the results
        smoothed.LL = tryCatch(predict(loess(LLRes~thr.grid)),error=function(e) e, warning=function(w) w)  #smoothed line 
        if(is(smoothed.LL,"warning") |is(smoothed.LL,"error") ){#in case not enough to smooth
          smoothed.LL = as.numeric(LLRes)
          #print(LLRes)
        } 
        #tt <- tryCatch(ddd(5),error=function(e) e, warning=function(w) w)
        
        #save results for optimal thr
        opt_i_grid = which.min(smoothed.LL) 
        tbl_node[tbl_node$split_variable == split_var,"thr"] = thr.grid[opt_i_grid]
        tbl_node[tbl_node$split_variable == split_var,"thr_quantile"] = names(thr.grid)[opt_i_grid]
        tbl_node[tbl_node$split_variable == split_var,"loss"] = smoothed.LL[opt_i_grid]
        tbl_node[tbl_node$split_variable == split_var,"likelihood"] = likelihood[opt_i_grid]   
        tbl_node[tbl_node$split_variable == split_var,"penalty"] = penalty[opt_i_grid] 
        tbl_node[tbl_node$split_variable == split_var,"child_volume1"] = volume[opt_i_grid,1]/total_volume
        tbl_node[tbl_node$split_variable == split_var,"child_volume2"] = volume[opt_i_grid,2]/total_volume
        #create new rules for opt tree
        new_rules = c(paste(parent_rule,paste(split_var, "<=", thr.grid[opt_i_grid],sep=""),sep = "&"), paste(parent_rule,paste(split_var,">", thr.grid[opt_i_grid],sep=""),sep = "&"))
        
        #block boundary solutions
        if(tbl_node[tbl_node$split_variable == split_var,"thr"]== min(thr.grid)|tbl_node[tbl_node$split_variable == split_var,"thr"]== max(thr.grid)){
          tbl_node[tbl_node$split_variable == split_var,"loss"] = tbl_node[tbl_node$split_variable == split_var,"loss"]*10000
        }
        
      }else{#if categorical
        parent_rule = tbl_tree[i_node, "rule"] 
        new_rules = c(paste(parent_rule,paste(split_var,"==0", sep = ""),sep = "&"), paste(parent_rule,paste(split_var,"==1", sep = ""),sep = "&"))
        #rules_tmp = tbl_tree[nodes_to_split[nodes_to_split!=i_node] ,"rule"] #all other rules than the splited
        rules_tmp = tbl_tree[(tbl_tree$action == "leaf"|tbl_tree$action == "split")&tbl_tree$node!=i_node,"rule"] #all other rules than the splited that are to be splitted or leaf
        rules_tmp = c(rules_tmp, new_rules)
        
        res = getDeviance(rules = rules_tmp,tri_data = dat.node.orig, ignore_dev=ignore_dev, sig_fixed = sig_fixed, calc_penalty)  #res.total[["sig_dev"]]
        
        tbl_node[tbl_node$split_variable == split_var,"loss"] = res[["crit"]]
        tbl_node[tbl_node$split_variable == split_var,"likelihood"] = res[["likelihood"]]
        tbl_node[tbl_node$split_variable == split_var,"penalty"] = res[["penalty_term"]]
        volume = res[["tri_volume"]][rules_tmp%in%new_rules]
        tbl_node[tbl_node$split_variable == split_var,"child_volume1"] = volume[1]/total_volume
        tbl_node[tbl_node$split_variable == split_var,"child_volume2"] = volume[2]/total_volume
        
      }
      #add child rules 
      tbl_node[tbl_node$split_variable == split_var,"child_rule1"] = new_rules[1]
      tbl_node[tbl_node$split_variable == split_var,"child_rule2"] = new_rules[2]
      
    }#for each split_var
    #print(tbl_node) 
    #Stop crit based on volume of the children - >if after split too small, set crit above total crit
    #if(sum(tbl_node$child_rule1<min_volume | tbl_node$child_rule2<min_volume))
    tbl_node[tbl_node$child_volume1<min_volume | tbl_node$child_volume2<min_volume,"loss"] = tbl_tree[1, "node_loss"] + 100000 #set crit above total loss + large intiger
    #boundary solution blocking within the i_var loop
    
    res_nodes[[i_node]]=tbl_node #save the results of spliting for each actual node
    #find the best split variable and save for optimization between nodes
    best_split_var = which.min(tbl_node[,"loss"])
    crit_nodes[,"crit"] = tbl_node[best_split_var, "loss"]
  }#for i_node
  #find best node to split
  best_split_node = crit_nodes[which.min(crit_nodes[,"crit"]),"node"] #id of the best node
  tbl_node = res_nodes[[best_split_node]] #delete all other node tables (keep only optimal)
  
  ###Update table of  nodes###
  best_split_var = which.min(tbl_node[,"loss"])
  
  ###Stop criterion that stops the growth for all nodes
  #If best split of all nodes worse than best previous, stop
  stop_crit_all = tbl_node[best_split_var,"loss"]>min(tbl_tree[, "node_loss"])
  #stop_crit_all = F #!!!!
  ###Stop criterion - do not split specific node
  #stop_crit = T
  
  
  
  if(stop_crit_all == T){##STOP ALL!!!
    tbl_tree[nodes_to_split, "action"] = "leaf" #all split nodes set to leaf
    
    #}else if(stop_crit == T){#do not perfor split
    #  tbl_tree[i_node, "action"] = "leaf"
  }else{#perform split
    tbl_tree[best_split_node, "action"] = "internal"
    tbl_tree[best_split_node, "split_variable"] = tbl_node[best_split_var, "split_variable"]
    tbl_tree[best_split_node, "split_loss"] = tbl_node[best_split_var, "loss"]
    tbl_tree[best_split_node, "likelihood"] = tbl_node[best_split_var, "likelihood"]
    tbl_tree[best_split_node, "penalty"] = tbl_node[best_split_var, "penalty"]
    tbl_tree[best_split_node, "thr"] = tbl_node[best_split_var, "thr"]
    tbl_tree[best_split_node, "thr_quantile"] = tbl_node[best_split_var, "thr_quantile"]
    
    ###add child nodes to current node
    #assign numbers to the children
    max_existing_node = max(tbl_tree$node)
    tbl_tree[max_existing_node+1, "node"] = max_existing_node+1
    tbl_tree[max_existing_node+2, "node"] = max_existing_node+2
    #add child rules
    tbl_tree[max_existing_node+1, "rule"] = tbl_node[best_split_var, "child_rule1"]
    tbl_tree[max_existing_node+2, "rule"] = tbl_node[best_split_var, "child_rule2"]
    #add node loss
    tbl_tree[max_existing_node+1, "node_loss"]=tbl_node[best_split_var, "loss"]
    tbl_tree[max_existing_node+2, "node_loss"]=tbl_node[best_split_var, "loss"]
    #tbl_tree[max_existing_node+1, "likelihood"]=tbl_node[best_split_var, "likelihood"]
    #tbl_tree[max_existing_node+2, "likelihood"]=tbl_node[best_split_var, "likelihood"]
    #tbl_tree[max_existing_node+1, "penalty"]=tbl_node[best_split_var, "penalty"]
    #tbl_tree[max_existing_node+2, "penalty"]=tbl_node[best_split_var, "penalty"]
    
    #add level
    tbl_tree[max_existing_node+1, "level"]=tbl_tree[best_split_node, "level"] + 1
    tbl_tree[max_existing_node+2, "level"]=tbl_tree[best_split_node, "level"] + 1
    #add node volume
    tbl_tree[max_existing_node+1, "volume"]=tbl_node[best_split_var, "child_volume1"]#/total_volume
    tbl_tree[max_existing_node+2, "volume"]=tbl_node[best_split_var, "child_volume2"]#/total_volume
    #add split action
    tbl_tree[max_existing_node+1, "action"]="split"
    tbl_tree[max_existing_node+2, "action"]="split"
    
    #update variables to split for child nodes
    split_var_grid_lst[[max_existing_node+1]] = update_split_var(rule =tbl_tree[max_existing_node+1,"rule"] , dat.node.orig, split_variables_all, metric_var = metric_var, tri_size, prob_grid) #parent rule is identiacal for each split variable  
    split_var_grid_lst[[max_existing_node+2]] = update_split_var(rule =tbl_tree[max_existing_node+2,"rule"] , dat.node.orig, split_variables_all, metric_var = metric_var, tri_size, prob_grid) #parent rule is identiacal for each split variable  
    
    #check, if no variable to split the node
    if(length(split_var_grid_lst[[max_existing_node+1]][["split_variables"]])==0){
      tbl_tree[max_existing_node+1, "action"] = "leaf"
    }
    if(length(split_var_grid_lst[[max_existing_node+2]][["split_variables"]])==0){
      tbl_tree[max_existing_node+2, "action"] = "leaf"
    }
    
  }#if stop perform split    
  #print(tbl_tree) 
  if(sum(tbl_tree$action =="split" )==0){perform_splits=F} #if nothing to split, stop  
  
  
}#whle perform splits

return(tbl_tree)

} #function




update_split_var = function(rule, dat.node.orig, split_variables_all, metric_var, tri_size, prob_grid){
  dat.node = getNodeData(flt.char = rule, in.data = dat.node.orig )
  flt = logical(length = length(split_variables_all))
  thr_grid_lst = list()
  for(i_var in split_variables_all){
    if(i_var %in% metric_var){ #size_var is always present
      
      thr.grid = quantile(dat.node[dat.node$dev_yr==tri_size,i_var], prob_grid )  
      thr.grid = thr.grid[thr.grid>0] #only nonzero threshold is taken
      thr.grid = thr.grid[!duplicated(thr.grid, fromLast = T)] #only unique values are taken
      thr.grid = thr.grid[!(thr.grid == max(dat.node[,i_var])|thr.grid ==min(dat.node[,i_var]))] #drop boundary values
      if(length(thr.grid)>=2){flt[ split_variables_all == i_var] =T}else{flt[ split_variables_all == i_var] =F}
      thr_grid_lst[[i_var]] = thr.grid
    }else{ #for other categorial variables
      flt[ split_variables_all == i_var] = sum(table(dat.node[,i_var])!=0)>1
    }
  }
  if(!exists(x = "thr.grid")){thr.grid = NULL} #case size_var is not used 
  return(list("split_variables" = split_variables_all[flt], "thr_grid_lst" = thr_grid_lst))
}



#####plot the tree#####
plot_tree = function(tbl_tree){
  path_str = tbl_tree$rule
  path_str[1] = "Total"
  #path_str = gsub(" ", "", path_str)  #remove white spaces
  
  
  my_frame = data.frame(
    "levelName" = tbl_tree$node,
    "pathString" = path_str,
    "p" = 100*round(tbl_tree$volume,2),
    "cost" = tbl_tree$node_loss
  )
  my_tree = as.Node(my_frame, pathDelimiter = "&")
  print(my_tree)
  
  GetNodeLabel = function(node){
    pthStr = strsplit(node$pathString, "/")
    pthStr = pthStr[[1]][length(pthStr[[1]])]
    return(paste(pthStr, paste0(", size ",node$p, " %")))
  }
  
  SetNodeStyle(my_tree, label = GetNodeLabel)
  plot(my_tree)
}




###############extract results table for all simulations ######
extract_res_tab = function(err, res_lst){
  
  res_tab = data.frame("i_sim" = NA,
                       "loss_root" = NA,
                       "loss" = NA,
                       "likelihood" = NA,
                       "penalty" = NA,
                       "err_splitted" = NA,
                       "err_total" = NA,
                       "err_lob" = NA,
                       "nr_leafs" = NA,
                       stringsAsFactors = F
  )
  
  
  
  for(i_sim in 1:nrow(err)){
    tbl_tree = res_lst[[i_sim]][["tbl_tree"]]
    res_tab[i_sim, "i_sim"] = i_sim
    res_tab[i_sim, "loss_root"] = tbl_tree[1,"node_loss"]
    
    #last split
    min_loss = min(tbl_tree[tbl_tree$action == "leaf","node_loss"]) #leaf just for sure
    if(nrow(tbl_tree)==1){ #in case no split
      last_split = 1
    }else{
      last_split= which(tbl_tree[,"split_loss"] == min_loss)  #this is the row of the last split containing final splitting  likelihood and penalty
    }
    res_tab[i_sim, "loss"] = min_loss
    res_tab[i_sim, "likelihood"] = tbl_tree[last_split,"likelihood"]
    res_tab[i_sim, "penalty"] = tbl_tree[last_split,"penalty"]
    res_tab[i_sim, "nr_leafs"] = sum(tbl_tree[,"action"]=="leaf")
    res_tab[i_sim, "err_splitted"] = err[i_sim,"err_splitted"] 
    res_tab[i_sim, "err_total"] = err[i_sim,"err_total"]
    res_tab[i_sim, "err_lob"] = err[i_sim,"err_lob"]
  }
  return(res_tab)
}




