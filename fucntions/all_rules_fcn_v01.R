#############create all possible partitons rules##################

extract_rule = function(part_tmp){
  #helper function for create_all_rules function
  #
  rule_tmp =c()
  for(i_part in 1:length(part_tmp)){
    rule_tmp =c(rule_tmp, paste(paste(part_tmp[[i_part]], "==1", sep = ""), collapse = "|"))
  }
  return(rule_tmp)
}


create_all_rules = function(base_set){
  #to create rules for all possible partitions
  #base_set = c(split_variables_all)
  parts <- listParts(length(base_set))
  #print(parts)
  out <- rapply(parts, function(ii) base_set[ii], how="replace")
  #print(out)
  
  rules_lst = list()
  
  for(i_rule in 1:length(out)){#i_rule = 1
    rules_lst[[i_rule]] = extract_rule(part_tmp = out[[i_rule]])
  }#for i_rule
  
  return(rules_lst)
}

extract_rule_elements = function(rules_lst){
  #To find all elements (triangles) of all rules
  #Normally set of all subsets of base set, but is more general
  rule_elements = c()
  for(i_rule in 1:length(rules_lst)){
    rule_elements = c(rule_elements,rules_lst[[i_rule]])
    rule_elements = unique(rule_elements)
  }
  return(rule_elements)
}

precalc_triangles = function(in.d, rule_elements){
#this is to precalculate triangles for all elements in rule_elements  
  
  by_elements = list() #list of elements separated into vectors of elementary elements
  n_elements = length(rule_elements)
  for(i_element in 1:n_elements){
    by_elements[i_element] = strsplit(rule_elements[i_element], fixed = T, split = "|")
  }
  element_legths = lengths(by_elements) #number of elementary elements in element
  
  #initialize list for triangles for all elements
  init = 0 
  precalc_triangles = rep(list(init),n_elements) 
  
  tri_size = length(levels(in.d[,"dev_yr"])) #dev_yr are factors 
  min_occ = min(in.d$occ_yr)
  max_occ = max(in.d$occ_yr)
  
  #get triangles for all elementary elements first
  for(i_element in 1:n_elements){ #i_element = 3
    
    if(element_legths[i_element] == 1){ #if elementary element, get triangle from data
      precalc_triangles[[i_element]] = get_tri(in.d, rule =by_elements[[i_element]], min_occ, max_occ,tri_size)
    }
  }
  
  #for non-elementary elements, aggregate elementary elements (elementary elements must be calculated first!)
  for(i_element in 1:n_elements){ #i_element = 1
    
    if(element_legths[i_element] != 1){ #if non-elementary element
      rule =by_elements[[i_element]]
      #tri_tmp = 0 #must be deleted after each element!!! 
      for(i_rule_el in rule){ #i_rule_el = rule[1]
        which_elementary  = which(rule_elements == i_rule_el)
        if(!exists("tri_tmp")){ #if first tri was not yet created, get the first element
          tri_tmp = precalc_triangles[[which_elementary]]
        }else{ #if partial aggregate already exists
          tri_tmp =tri_tmp+ precalc_triangles[[which_elementary]]
        } #if first tri_tmp
        
      }#i_rule_el
      precalc_triangles[[i_element]] = tri_tmp #save aggregated tri_tmp
      rm(tri_tmp) #reset tri_tmp
    }#if non elementary
    
  }#i_element
  return(precalc_triangles)
  
}



comb_catvar=function(catvar, dat){
  n_all_lobs =length(unique(dat[,catvar]))#length(levels(data.source[,"LoB"])) #to how many lobs split
  #create indicators for LoBs
  #create all combinations
  if((n_all_lobs %% 2) == 0){#even
    half_n_lob = n_all_lobs/2
  }else{#odd
    half_n_lob = ceiling(n_all_lobs/2)
  }
  combs = list()
  for(i_n_lob in 1:half_n_lob){
    combs[[i_n_lob]] = combn(x = seq(n_all_lobs), simplify = T, m=(i_n_lob))
  }
  #print(combs)
  
  
  #create LoB design matrix
  prev_combs = c()
  comb_dsgn = data.frame(as.integer(dat[,catvar]))
  names(comb_dsgn) =  catvar 
  for(i_n_lob in 1:half_n_lob){
    comb_tmp = combs[[i_n_lob]]
    #lob_tmp = as.integer(dat.node1$LoB)
    for(i_com in 1:ncol(comb_tmp)){
      nr_lob_name = paste0(comb_tmp[,i_com], collapse = '') #combination name (numeric part)
      diff_name = paste0(setdiff(seq(n_all_lobs),comb_tmp[,i_com]), collapse = '') # difference to all Lobs set
      if(!(diff_name%in%prev_combs)){ #if difference not previously added, add this combination
        comb_name = paste0(catvar,"_", nr_lob_name)
        comb_dsgn[,comb_name] = comb_dsgn[,catvar]%in%comb_tmp[,i_com]
        prev_combs = c(prev_combs,nr_lob_name)
      }
    }
    
  }
  return(comb_dsgn[-1])
}


