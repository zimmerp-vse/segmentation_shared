
getDeviance = function(rules, tri_data, ignore_dev, sig_fixed = NULL, calc_penalty = T, precalc = F, rule_elements = NA){
 
  #tri_data either data directly (in.d) or precalc_triangles
  
  #get triangles for a split
  if(!precalc){ #directly from data
    nodes_tri = get_tri_list(rules, tri_data = tri_data, sig_fixed = sig_fixed, precalc = F)
  }else{#from precalculated
    nodes_tri = get_tri_list(rules, tri_data=tri_data, sig_fixed = sig_fixed, precalc = T, rule_elements=rule_elements)
  }
  tri_size =  dim(nodes_tri[["D"]][[1]])[1]  #length(levels(in.d[,"dev_yr"])) #dev_yr are factors
  n_triangles = length(rules)
  #aggregate results
  nodes_tri[["T"]] = sum_tri(nodes_tri[["D"]])
  nodes_tri[["T_pred"]] = sum_tri(nodes_tri[["D_pred"]])
  nodes_tri[["TF"]] = as.triangle(ata(Triangle = nodes_tri[["T"]], NArow.rm = TRUE, colname.sep = "-")) #ge-to-age factors
  nodes_tri[["TF_pred"]] = as.triangle(ata(Triangle = nodes_tri[["T_pred"]], NArow.rm = TRUE, colname.sep = "-")) #ge-to-age factors
  
  nodes_tri[["vars_tri_T"]] = getVars_tri(tri=nodes_tri[["T"]],df=nodes_tri[["TF"]],tri_preds =nodes_tri[["T_pred"]], DF_pred = nodes_tri[["TF_pred"]], sig_fixed = sig_fixed, ignore_dev = ignore_dev, n_triangles= n_triangles)#  fix_sig_dev = NULL => calculate sigmas and fix them for now
  nodes_tri[["vars_tri_T"]] = addDeviance(vars_tri_T = nodes_tri[["vars_tri_T"]], sig_fixed = sig_fixed)
  
  flt.deviance = is.finite(nodes_tri[["vars_tri_T"]]$deviance)&!is.na(nodes_tri[["vars_tri_T"]]$deviance)
  n.obs_total = sum(flt.deviance) 
  dev.in = sum(nodes_tri[["vars_tri_T"]][flt.deviance, "deviance"])/n.obs_total
  #print(nodes_tri[["vars_tri_T"]][flt.deviance, "deviance"])
  
  ###get covariances between two triangles
  #instantiate 2D list of combinations of all triangles
  rhos <- as.list(matrix(NA,n_triangles^2))
  dim(rhos) <- c(n_triangles,n_triangles)
  cov_long = rhos #initiate covariances 
  
  #place correlations of two triangles for all devs
  #a_tri = 1
  #b_tri = 2
  for(a_tri in 1:n_triangles){
    for(b_tri in 1:a_tri){
      if(a_tri>b_tri){
        rhos[[a_tri,b_tri]] = getRhoTri(vars_tri_list=nodes_tri[["vars_tri_list"]],a_tri=a_tri, b_tri=a_tri, ignore_dev) #vector of cov for each dev
        cov_long[[a_tri,b_tri]] = addCov_long(vars_tri_list=nodes_tri[["vars_tri_list"]], rhos_ab = rhos[[a_tri,b_tri]], a_tri=a_tri,b_tri=a_tri)  #covariance of couples of triangles for each dev aqnd occ
      } 
    }#a_tri
  }#b_tri
  
  #calculate the penalty
  if(calc_penalty == T){
    penalty = getPenalty(nodes_tri = nodes_tri, cov_long = cov_long, sig_fixed = sig_fixed, tri_size =tri_size)
  }else{
    penalty = 0 
  }
  #calculate volume of diagonals after split
  tri_volume = getTri_volume(nodes_tri[["D"]])
  return(list("crit" =dev.in+2*penalty, "likelihood" = dev.in, "penalty_term"=2*penalty, "sig_fixed"=nodes_tri[["vars_tri_T"]]["sig"], "tri_volume" = tri_volume, "pseudoresids" = nodes_tri[["vars_tri_T"]][flt.deviance, "deviance"]))
}


getPenalty = function(nodes_tri, cov_long, sig_fixed, tri_size){
  vars_tri = nodes_tri[["vars_tri_T"]]
  n_triangles = length(nodes_tri[["D"]])
  
  sum_w_cov_AB = matrix(0,nrow = (tri_size-1)^2) #instantiate with fake 0
  for(a_tri in 1:n_triangles){
    for(b_tri in 1:n_triangles){
      if(a_tri == b_tri){
        sum_w_cov_AB = cbind(sum_w_cov_AB,as.matrix(nodes_tri[["vars_tri_list"]][[a_tri]]["w_occ_yr"]*nodes_tri[["vars_tri_list"]][[a_tri]]["Var_D"]))
      }else{
        if(a_tri>b_tri){c_AB = cov_long[[a_tri,b_tri]]}else{c_AB = cov_long[[b_tri,a_tri]]} #mirroring of cov(A_ij, B_ij)
        sum_w_cov_AB = cbind(sum_w_cov_AB, as.matrix(nodes_tri[["vars_tri_list"]][[a_tri]]["w_occ_yr"]*c_AB))
      }#if
    }#b_tri
  }#a_tri
  #get rid of the first column of zeros created only to define the size
  sum_w_cov_AB = sum_w_cov_AB[,2:ncol(sum_w_cov_AB)]
  
  if(n_triangles>1){ #for total (n_triangles == 1) nothing to aggregate
    vars_tri[,"covD_D_hat"] = rowSums(sum_w_cov_AB)
  }else{
    vars_tri[,"covD_D_hat"] =sum_w_cov_AB
  }
  
  vars_tri[,"covf_f_hat"] = vars_tri[,"covD_D_hat"] / (vars_tri[,"D"]*vars_tri[,"D_pred"])

  #fixed sigma???
  if(is.null(sig_fixed)){
    vars_tri[,"covf_p"] = vars_tri[,"covf_f_hat"]/vars_tri[,"Var_f"] #replace varf_ij by varf_fixed if fixed is required
  }else{
    vars_tri[,"covf_p"] = vars_tri[,"covf_f_hat"]/vars_tri[,"Var_f_fixed"] #replace varf_ij by varf_fixed if fixed is required
  }
  flt.cov = is.finite(vars_tri[,"covf_p"])&!is.na(vars_tri[,"covf_p"])
  if(penalty_fun == "median"){
      return(median(vars_tri[flt.cov,"covf_p"]))
    }else{
      return(mean(vars_tri[flt.cov,"covf_p"], trim = trunc_penalty_mean))
    }
  #return(median(vars_tri[flt.cov,"covf_p"]))
  #return(median(vars_tri[flt.cov,"covD_D_hat"]))
}

addDeviance = function(vars_tri_T, sig_fixed = NULL){#sig_fixed only to indicate which to use
  if(is.null(sig_fixed)){ #if fixed not required
    denom = vars_tri_T$Var_f
  }else{#add fixed sig if needed
    #vars_tri_T["var_f_fixed"] = merge(x=vars_tri_T,y=fix_sig_dev, by= 'dev',all.x=TRUE)[,"sig.y"]
    denom = vars_tri_T$Var_f_fixed
  }
  
  #add deviance to total
  vars_tri_T["deviance"] = (vars_tri_T$df-vars_tri_T$df_pred)^2/denom
  return(vars_tri_T)
  #vars_tri_T = getVars_tri(D=D,DF=DF,D_pred = D_pred, DF_pred=DF_pred, fix_sig_dev = NULL, i_tri = i_tri, ignore_dev = 3, n_triangles= n_triangles)
}
  
addOccYrSums = function(vars_tri){
  #abc = c("A","B", "D")
  #work_col = abc[i_tri]
  #if(i_tri==3){work_col = "A"}elseif{work_col = "B"}
  flt.acc = !is.na(vars_tri[,"D_pred"]) & !is.nan(vars_tri[,"D_pred"])
  sum_ij = aggregate(vars_tri[flt.acc,"D"],
                     by = list(vars_tri[flt.acc,"dev"]),
                     FUN = sum)
  names(sum_ij) = c("dev", "sum_ij" )
  sum_ij = merge(x = vars_tri, y= sum_ij, by = "dev")
  sum_ij["w_ij"] = sum_ij[,"D_pred"]/sum_ij[,"sum_ij"]
  return(sum_ij[,"w_ij"])
}



  
get_tri_list=function(rules, tri_data=NA, sig_fixed, precalc=F, rule_elements=NA ){
#calculate list of triangles from given rules either directly form data or from precalculated triangles
  if(!precalc){ #directly from data
    D = build_tri(rules = rules, in.d = tri_data, node_ids = NULL)
  }else{ #from precalculated
    D = build_tri_precalc(rules = rules, pre_triangles = tri_data, rule_elements, node_ids = NULL)
  }
  
  tri_size = dim(D[[1]])[1] #length(levels(in.d[,"dev_yr"])) #dev_yr are factors
  n_triangles = length(rules)
  
  D_pred = predict_D(D) #calculate separate ch-l predicitions before aggregation
  #Triangles of development factors
  DF= list()
  for(i_tri in 1:n_triangles){
    tri_tmp = as.triangle(ata(Triangle = D[[i_tri]], NArow.rm = FALSE, colname.sep = "-")) #age-to-age factors
    #get rid of last row full of NA
    tri_tmp = tri_tmp[-nrow(tri_tmp),]

    DF[[i_tri]] = as.triangle(tri_tmp)
  }
  
  DF_pred = list()
  for(i_tri in 1:n_triangles){
    tri_tmp = as.triangle(ata(Triangle = D_pred[[i_tri]], NArow.rm = FALSE, colname.sep = "-")) #ge-to-age factors
    #get rid of last row full of NA
    tri_tmp = tri_tmp[-nrow(tri_tmp),]
    DF_pred[[i_tri]]= as.triangle(tri_tmp)
  }
  
  
  #calculate variances
  vars_tri_list = list()
  for(i_tri in 1:n_triangles){
    #vars_tri_list[[i_tri]] = getVars_tri(D=D,DF=DF,D_pred = D_pred, DF_pred=DF_pred, fix_sig_dev = NULL, i_tri = i_tri, ignore_dev = 3, n_triangles= n_triangles)#  fix_sig_dev = NULL => calculate sigmas and fix them for now
    vars_tri_list[[i_tri]] = getVars_tri(tri=D[[i_tri]],df=DF[[i_tri]],tri_preds = D_pred[[i_tri]], DF_pred=DF_pred[[i_tri]], ignore_dev = 3, n_triangles= n_triangles,sig_fixed = sig_fixed)#  fix_sig_dev = NULL => calculate sigmas and fix them for now
  }                                      #tri,          df,               tri_preds,                preds,                  DF_pred, ignore_dev, n_triangles, sig_fixed

return(list("D" = D,"D_pred" = D_pred, "DF" = DF, "DF_pred" =DF_pred,"vars_tri_list" = vars_tri_list))
}


###caluclate covariances
#estimate correlation coefficient for each occ and dec for a couple of triangles a_tri, b_tri
getRhoTri = function(vars_tri_list,a_tri, b_tri, ignore_dev){
  tri_size = ncol(vars_tri_list["D"][[1]]) #expect all tri with the the samelength#max(vars_tri_list[[3]][,"dev"]+1)  
  vars_tri = vars_tri_list[[1]][,c("origin", "dev")] #local template 
  vars_tri[,"rho_varA_varB"] = (vars_tri_list[[a_tri]][,"df"] -vars_tri_list[[a_tri]][,"df_pred"] )*sqrt(vars_tri_list[[a_tri]][,"D"]^(-2+variance_power))*(vars_tri_list[[b_tri]][,"df"] -vars_tri_list[[b_tri]][,"df_pred"] )*sqrt(vars_tri_list[[b_tri]][,"D"]^(-2+variance_power))
  vars_tri[,"rho_ij"] = vars_tri[,"rho_varA_varB"] / (sqrt(vars_tri_list[[a_tri]][,"Var_D"])*sqrt(vars_tri_list[[b_tri]][,"Var_D"])  )
  
  flt.acc =  !is.na(vars_tri_list[[1]][,"df_pred"]) #filter upper tri
  rho_j = aggregate(vars_tri[flt.acc,"rho_ij"],
                    by = list(vars_tri[flt.acc,"dev"]),
                    FUN = sum)/
    (aggregate(vars_tri[flt.acc,"rho_ij"],
               by = list(vars_tri[flt.acc,"dev"]),
               FUN = length)-1)
  rho_j[,1] = seq(1:nrow(rho_j))
  colnames(rho_j) = c("dev","rho")
  
  rho_j[(nrow(rho_j)-ignore_dev+1):nrow(rho_j), "rho"]=0
  #plot(rho_j)
  return(rho_j)
}


addCov_long = function(vars_tri_list, rhos_ab, a_tri,b_tri){
  #leftjoin of correlations for each occ and dec for a couple of triangles a_tri, b_tri
  #creates covariances
  vars_tri = vars_tri_list[[a_tri]] #any of the vars_tri
  vars_tri = merge(x=vars_tri,y=rhos_ab, by="dev",all.x=TRUE)
  vars_tri[["covAB"]] = vars_tri[["rho"]]*sqrt(vars_tri_list[[a_tri]][,"Var_D"])*sqrt(vars_tri_list[[b_tri]][,"Var_D"])
  return(vars_tri[["covAB"]])
}



####Calculate variances#####
#D, DF, D_pred, DF_pred, fix_sig_dev, i_tri, ignore_dev, n_triangles
#getVars_tri = function(tri, df, tri_preds, preds, DF_pred, fix_sig_dev, i_tri, ignore_dev, n_triangles){
  
getVars_tri = function(tri, df, tri_preds, DF_pred, ignore_dev, n_triangles, sig_fixed){
    #tri = D[[i_tri]]
    #df = DF[[i_tri]]
    #preds = DF_pred[[i_tri]]
    #tri_preds = D_pred[[i_tri]]
    
    sig_dev = getSigTri(tri, df, DF_pred,ignore_dev) #for last 3 dev yrs, sig are set to 0 and then deleed 
    
    #fix_sig_dev
  
  #plot(sig_dev[1,])
  vars_tri = suppressWarnings(as.data.frame(DF_pred))
  names(vars_tri) =c("origin", "dev", "value")
  vars_tri = merge(x=vars_tri,y=sig_dev, by= 'dev',all.x=TRUE)
  cum_val = as.data.frame(as.triangle(tri[-ncol(tri),-ncol(tri)]))
  vars_tri[,"D"] = cum_val[,"value"]
  cum_pred_val = as.data.frame(as.triangle(tri_preds[-ncol(tri_preds),-ncol(tri_preds)]))
  vars_tri[,"D_pred"] = cum_pred_val[,"value"]
  vars_tri = suppressWarnings(cbind(vars_tri[,c("origin", "dev", "value", "D", "D_pred", "sig")], as.data.frame(df)[,"value"]))
  names(vars_tri) = c("origin", "dev","df_pred", "D", "D_pred", "sig","df")
  #add variance
  vars_tri[,"Var_D"] = vars_tri[,"sig"]*(vars_tri[,"D_pred"]^variance_power)
  vars_tri[,"Var_f"] = vars_tri[,"sig"]/(vars_tri[,"D_pred"]^(2-variance_power))
  vars_tri[,"w_occ_yr"] = addOccYrSums(vars_tri)
  if(!is.null(sig_fixed)){vars_tri[,"Var_f_fixed"] = sig_fixed/(vars_tri[,"D_pred"]^(2-variance_power))}
  return(vars_tri)
  #add deviance
  #vars_tri$deviance = ((vars_tri$df-vars_tri$pred)/vars_tri$std)^2
  #return(list("vars_tri" = vars_tri,"sig_dev" = sig_dev))
}


getSigTri = function(tri, df, DF_pred, ignore_dev){
  tri_size = ncol(tri)
  DIF =as.triangle(matrix(NA, nrow=tri_size, ncol = tri_size) )
  sig_dev = matrix(NA, nrow=1, ncol = tri_size)
  for (i_dev in 2:tri_size){
    tri_rows = tri_size-i_dev+1
    DIF[1:tri_rows,i_dev] = (df[1:tri_rows,i_dev-1] - DF_pred[1:tri_rows,i_dev-1])^2*tri[1:tri_rows,i_dev-1]*(tri[1:tri_rows,i_dev-1]^(-variance_power+1) )
    sig_dev[1,i_dev-1] = sum(DIF[!is.na(DIF[,i_dev]),i_dev]) / (sum(!is.na(DIF[,i_dev])) -1) #n-1 in denominator
    #sig_dev[1,i_dev-1] = mean(DIF[!is.na(DIF[,i_dev]),i_dev])
  }
  #set the last sigma (with only 1 observation) to the value of last but one
  #sig_dev[1,tri_size] = sig_dev[1,tri_size-1]
  
  #set last ignore.dev values to 0 and consequently ignore
  sig_dev[1,(-tri_size+ignore_dev):0] = 0
  
  #plot(sig_dev[1,])
  #transform
  sig_dev =  data.frame("dev" = seq(tri_size),"sig" = t(sig_dev))
  return(sig_dev)
}





###evaluation functions for simulation

getTriErr = function(rules,node_list, dat.node.orig, size_var){
  #to calculate the true error for a given set of rules
  D= build_tri(rules = rules, in.d = dat.node.orig, node_ids = node_list)
  tri_size = length(levels(dat.node.orig[,"dev_yr"])) #dev_yr are factors
  n_triangles = length(rules)
  
  ibnr_model = c()
  latest = c()
  cellEst = c()
  rowEst = c()
  row_latest = c()
  for(i_tri in 1:n_triangles){
    #print(i_tri)
    tri_model = suppressWarnings(MackChainLadder(Triangle = D[[i_tri]], est.sigma = "Mack"))
    latest = c(latest,summary(tri_model)$Totals[1,1]) #sum of diagonal
    ibnr_model = c(ibnr_model,summary(tri_model)$Totals[4,1]) # 
    
    #for testing MSE
    #last_col = dim(aaa)[1]
    cellEst = c(cellEst, tri_model$FullTriangle[tri_size,3]) #one particular cell estimate
    rowEst = c(rowEst, summary(tri_model)$ByOrigin$IBNR[tri_size-1]) # one but last #=tri_model$FullTriangle[9,tri_size]
    row_latest =c(row_latest, summary(tri_model)$ByOrigin$Latest[tri_size-1]) #cell on diagonal
  }
  
  ibnr_real = sum(dat.node.orig[dat.node.orig$dev_yr==tri_size,size_var]) -sum(latest)
  err = ibnr_real-sum(ibnr_model)
  
  #for testing MSE
  cell_real = sum(dat.node.orig[dat.node.orig$dev_yr==3 & dat.node.orig$occ_yr== 2004,size_var]) 
  err_cell =  cell_real - sum(cellEst) #error of a selected cell
  
  row_real = sum(dat.node.orig[dat.node.orig$dev_yr==tri_size & dat.node.orig$occ_yr== 2002,size_var]) -sum(row_latest)
  err_row = row_real - sum(rowEst) #error of selected row
  
  #tri_model$FullTriangle
  return(c(err, sum(ibnr_model), ibnr_real, err_cell, err_row))
}
