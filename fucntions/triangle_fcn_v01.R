
#calculate total tri
#to aggregate all triangles in list D
sum_tri = function(D){
  n_triangles = length(D)
  T = D[[1]]
  if(n_triangles>1){
    for(i_tri in 2:n_triangles){
      T = T + D[[i_tri]]
    }
  }#if
  return(T)
}



#to create a list of triangles for a set of rules
build_tri = function(rules, in.d=tri_data, node_ids){
  D = list()
  #number of triangles
  n_triangles = length(rules)
  #length of the triangle
  tri_size = length(levels(in.d[,"dev_yr"])) #dev_yr are factors 
  
  #in.d.art = activateAllCells(in.d, n_triangles)
  #in.d = rbind(in.d,in.d.art)
  min_occ = min(in.d$occ_yr)
  max_occ = max(in.d$occ_yr)
  
  for(i_tri in 1:n_triangles){
    
    D[[i_tri]] = get_tri(in.d, rule=rules[i_tri],  min_occ, max_occ,tri_size)
    #dat_tri = getNodeData(in.data=in.d, flt.char=rules[i_tri])
    #dat_tri_art = activateAllCells(in.data = dat_tri,min_occ = min_occ ,max_occ = max_occ, n_triangles)
    #dat_tri = rbind(dat_tri,dat_tri_art)
    #cumTri <- xtabs(paid~occ_yr+dev_yr,data=dat_tri)
    #cumTri[,tri_size:1][lower.tri(cumTri)]<-NA
    #cumTri <- as.data.frame(t(cumTri))
    #cumTri <- as.triangle(cumTri,
    #                      origin="occ_yr",
    #                      dev="dev_yr",
    #                      value="Freq")
    #cumTri[cumTri==0] = 0.000001 # to prevent dividing 0 and having different nr. of occurrence years
    #D[[i_tri]] <- cumTri
    
  }
  if(!is.null(node_ids)){
    names(D) = node_ids #name triangles by node ids
  }
  return(D)
}

#to create a list of triangles for a set of rules
build_tri_precalc = function(rules, pre_triangles, rule_elements, node_ids){
  D = list()
  #number of triangles
  n_triangles = length(rules)
  #length of the triangle
  #tri_size = length(levels(in.d[,"dev_yr"])) #dev_yr are factors 
  #min_occ = min(in.d$occ_yr)
  #max_occ = max(in.d$occ_yr)
  
  for(i_rule in rules){
    #rule=rules[i_tri]
    which_rule  = which(rule_elements == i_rule) #search for given rule amongst precalculated triangles
    D[[i_rule]] = pre_triangles[[which_rule]] 
      #get_tri(in.d = dat_tri, rule=rules[i_tri],  min_occ, max_occ,tri_size)
    
  }
  if(!is.null(node_ids)){
    names(D) = node_ids #name triangles by node ids
  }
  return(D)
}



get_tri = function(in.d, rule,  min_occ, max_occ,tri_size){ 
  #to create a single triangle based on a rule
  #min_occ = min(in.d$occ_yr)
  #max_occ = max(in.d$occ_yr)
  #tri_size = length(levels(in.d[,"dev_yr"])) #dev_yr are factors
  
  dat_tri = getNodeData(in.data=in.d, flt.char=rule)
  dat_tri_art = activateAllCells(in.data = dat_tri,min_occ = min_occ ,max_occ = max_occ, n_triangles)
  dat_tri = rbind(dat_tri,dat_tri_art)
  cumTri <- xtabs(paid~occ_yr+dev_yr,data=dat_tri)
  cumTri[,tri_size:1][lower.tri(cumTri)]<-NA
  cumTri <- as.data.frame(t(cumTri))
  cumTri <- as.triangle(cumTri,
                        origin="occ_yr",
                        dev="dev_yr",
                        value="Freq")
  cumTri[cumTri==0] = 0.000001 # to prevent dividing 0 and having different nr. of occurrence years
  return(cumTri)
}



####Triangle predictions
#filtr node data
create_rules_binom = function(split_variable, thr = NULL){
  if(split_variable == "size"){#if split done by size
    #rules = c(paste(split_variable,"<=",thr, sep = ""), paste(split_variable,">",thr, sep = ""))
    rules = c(paste("paid<=",thr, sep = ""), paste("paid>",thr, sep = ""))
  }else{
    rules = c(paste(split_variable,"==0", sep = ""), paste(split_variable,"==1", sep = ""))
  }
  return(rules)
}

getNodeData = function(in.data, flt.char){
  newdata <- subset(in.data, eval(parse(text =flt.char )))
  return(newdata)
}


predict_D = function(D){
  preds = list()
  for(i_tri in 1:length(D)){
    preds[[i_tri]] = getPreds_cum(D[[i_tri]])
  }
  return(preds)
}


#calculate predictions from triangles
getPreds_cum = function(tri){
  factors.train1 = getFactors(tri) #estimate factors of small tri
  preds = as.triangle(getFitted(cumTri = tri,devFactors = factors.train1))
  return(preds)
}

getFactors <- function(trojuhelnik)
{
  trojuhelnik[is.na(trojuhelnik)]<-0L
  r <- dim(trojuhelnik)[1]
  s <- dim(trojuhelnik)[2]
  rozdil <- r - s
  pomocTroj=trojuhelnik
  if(rozdil>0) diag(pomocTroj[-c(0:rozdil),s:1])=0L
  else diag(pomocTroj[,s:1])=0L
  zobrazeni=rep(1,dim(trojuhelnik)[1])
  faktory=(t(trojuhelnik[,-1]) %*% zobrazeni)/(t(pomocTroj[,-s]) %*% zobrazeni)
  return(faktory)
}

#calculates fitted values for a triangle
getFitted=function(cumTri,devFactors)
{
  delka=dim(cumTri)[1]
  pomocMat=matrix(0,delka,delka)
  pomocMat[1:delka,]=log(c(devFactors,1))
  pomocMat=t(pomocMat)
  pomocMat[,delka:1][lower.tri(pomocMat)]=0L
  diag(pomocMat[,delka:1])=0L
  pomocMat2=matrix(1,delka,delka)
  pomocMat2[upper.tri(pomocMat2)]=0L
  pomocMat=exp(pomocMat%*%pomocMat2)
  pomocMat2[1:delka,]=diag(cumTri[,delka:1])
  trojuhelnik_fitted=pomocMat2/pomocMat
  trojuhelnik_fitted[,delka:1][lower.tri(trojuhelnik_fitted)]=0L
  trojuhelnik_fitted[trojuhelnik_fitted==0] <- NA
  return(trojuhelnik_fitted)
}



activateAllCells = function(in.data, min_occ, max_occ, n_triangles){
  #creates dummy file . All occurrence years appear. 0 in first dev year.
  oc.yrs.all = seq(from = min_occ, to = max_occ)
  in.d.art =  in.data[1,] #copy first row
  #in.d.art[1:((n.thr.p+1)*length(oc.yrs.all)),] = NA
  in.d.art[1:length(oc.yrs.all),] = NA
  in.d.art[,"occ_yr"] = oc.yrs.all
  in.d.art[,"paid"] = 0
  in.d.art[,"dev_yr"] = 1
  return(in.d.art)
}


getTri_volume = function(D){
  tri_volume = c()
  for(i_tri in 1:length(D)){
    tri_diag = getLatestCumulative(D[[i_tri]])
    tri_volume = c(tri_volume, sum(tri_diag))
  }
  return(tri_volume)
}



###development encoding for data generator
fact_enc_plain = function(enc_var, dat.node){ #for encoding in working example
  #enc_var = "inj_part"
  if(is.factor(dat.node[,enc_var])){values_encoding = levels(dat.node[,enc_var])}
  if(is.integer(dat.node[,enc_var])){values_encoding = unique(dat.node[,enc_var])}
  rules = paste(enc_var,"==",values_encoding, sep = "")
  D= build_tri(rules = rules, in.d = dat.node, node_ids = values_encoding)
  factors_tri = data.frame(values_encoding)
  factors_tri[,"ult_f"] = NA
  
  for(i_tri in 1:length(D)){
    factors_tri[i_tri,"ult_f"] = prod(getFactors(D[[i_tri]])) 
  }
  return(mapvalues(x = dat.node[,enc_var], from = factors_tri[,"values_encoding"], to= factors_tri[,"ult_f"]))
}  

fact_enc = function(enc_var, dat.node){ #speific for simulation
  #enc_var = "inj_part"
  if(is.factor(dat.node[,enc_var])){values_encoding = levels(dat.node[,enc_var])}
  if(is.integer(dat.node[,enc_var])){values_encoding = unique(dat.node[,enc_var])}
  rules = paste(enc_var,"==",values_encoding, sep = "")
  D= build_tri(rules = rules, in.d = dat.node, node_ids = values_encoding)
  factors_tri = data.frame(values_encoding)
  factors_tri[,"ult_f"] = NA
  
  for(i_tri in 1:length(D)){
    factors_tri[i_tri,"ult_f"] = prod(getFactors(D[[i_tri]])) 
  }
  dat.node[,paste(enc_var, "_enc", sep = "")] = mapvalues(x = dat.node[,enc_var], from = factors_tri[,"values_encoding"], to= factors_tri[,"ult_f"])
  return(as.numeric(levels( dat.node[,paste(enc_var, "_enc", sep = "")]))[ dat.node[,paste(enc_var, "_enc", sep = "")]])  
}
