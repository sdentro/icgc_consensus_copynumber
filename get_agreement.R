
#' Issues
#'  - grep TODO in this script



# get_rounded_agreement = function(samplename, method_segmentsfile, method_purityfile) {
#   
#   breakpoints = read.table(paste0(samplename, "_consensus_breakpoints.txt"), header=T, stringsAsFactors=F)
#   segments = breakpoints2segments(breakpoints)
#   
#   res = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=T)
#   map_dkfz = res$map_dkfz
#   map_mustonen = res$map_mustonen
#   map_peifer = res$map_peifer
#   map_vanloowedge = res$map_vanloowedge
#   map_broad = res$map_broad
#   all_maps = list(map_broad=map_broad, map_dkfz=map_dkfz, map_mustonen=map_mustonen, map_vanloowedge=map_vanloowedge, map_peifer=map_peifer)
#   
#   res = get_frac_genome_agree(all_maps, segments)
#   frac_agree = res$frac_agree
#   seg_agree = res$agree
#   cn_states = res$cn_states
# }

#' @param min_methods_agree The minimum number of methods that is required to agree
#' @param min_methods_with_call_on_segment The minimum number of methods with a call for a segment to be considered for agreement
#' @param method_overruled A data frame with a single row and a column per method. Each cell contains TRUE if the method is to be overruled
get_frac_genome_agree = function(samplename, all_data, segments, min_methods_agree=0, min_methods_agree_x=0, min_methods_agree_y=0, min_methods_with_call_on_segment=2, min_methods_with_call_on_segment_x=2, min_methods_with_call_on_segment_y=2, method_overruled=NA, allowed_methods_x_female=c("dkfz", "mustonen", "vanloowedge", "jabba", "broad"), allowed_methods_x_male=c("dkfz", "mustonen", "jabba", "broad"), allowed_methods_y=c("dkfz", "jabba")) {
  # breakpoints = read.table(paste0(samplename, "_consensus_breakpoints.txt"), header=T, stringsAsFactors=F)
  # segments = breakpoints2segments(breakpoints)

  # res = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=mustonen_has_header)
  map_dkfz = all_data$map_dkfz
  map_mustonen = all_data$map_mustonen
  map_peifer = all_data$map_peifer
  map_vanloowedge = all_data$map_vanloowedge
  map_broad = all_data$map_broad
  map_jabba = all_data$map_jabba
  all_maps = list(map_broad=map_broad, map_dkfz=map_dkfz, map_mustonen=map_mustonen, map_vanloowedge=map_vanloowedge, map_peifer=map_peifer, map_jabba=map_jabba)
  # combined_status = data.frame(segments, dkfz=map_dkfz$status, mustonen=map_mustonen$status, peifer=map_peifer$status, vanloowedge=map_vanloowedge$status, broad=map_broad$status)
  combined_status = get_combined_status(segments, map_vanloowedge, map_dkfz, map_mustonen, map_peifer, map_broad, map_jabba)
  
  # Order both combined_status and method_overruled to have the same order - if method_overruled is supplied
  if (!is.na(method_overruled)) {
    method_overruled = method_overruled[colnames(combined_status)[4:9]]
  }
  
  
  agree = rep(F, nrow(segments))
  cn_states = list()
  num_methods = rep(0, nrow(segments))
  for (i in 1:nrow(segments)) {

    ###########################################
    # Do different things when addressing X and Y because there are fewer methods reporting
    ###########################################
    # Order of methods is: dkfz, mustonen, peifer, vanloowedge, broad, jabba
    if (segments$chromosome[i]=="X") {
      if (sex=="female") { allowed_methods_x = allowed_methods_x_female 
      } else { allowed_methods_x = allowed_methods_x_male }
      selection = !is.na(combined_status[i,4:9]) & colnames(combined_status)[4:9] %in% allowed_methods_x
      # These can be overwritten because X and Y are always the last chromosomes to be considered
      min_methods_with_call_on_segment = min_methods_with_call_on_segment_x
      min_methods_agree = min_methods_agree_x
      
    } else if (segments$chromosome[i]=="Y") {
      selection = !is.na(combined_status[i,4:9]) & colnames(combined_status)[4:9] %in% allowed_methods_y
      # These can be overwritten because X and Y are always the last chromosomes to be considered
      min_methods_with_call_on_segment = min_methods_with_call_on_segment_y
      min_methods_agree = min_methods_agree_y
      
    } else {
      selection = !is.na(combined_status[i,4:9])
    }
    
    # Check if method overruling is requested
    if (!is.na(method_overruled)) {
      selection = selection & !method_overruled
    }
    # Finally select the indices that correspond the methods with results
    methods_with_result = (4:9)[selection]
    
    ###########################################
    # If no methods report a result, skip
    ###########################################
    if (length(methods_with_result)==0) {
      next
    }
    
    ###########################################
    # All methods agree
    ###########################################
    all_methods_agree = length(methods_with_result) >= min_methods_with_call_on_segment &
                        sum(combined_status[i,methods_with_result]=="clonal", na.rm=T) >= min_methods_agree &
                        all(combined_status[i,methods_with_result]=="clonal")
    
    if (!is.na(method_overruled)) {
      # All overruled methods agree - but it must be more than 50% of methods and at least the minimum number
      all_methods_agree_no_overrule = length(methods_with_result) >= min_methods_with_call_on_segment && 
        (sum(combined_status[i,methods_with_result]=="clonal")/length(methods_with_result) > 0.5) && 
        all(combined_status[i,methods_with_result]=="clonal") 
    } else {
      # Not required to exclude overruled methods, so set this to true
      all_methods_agree_no_overrule = T
    }
    
    ###########################################
    # Accept a segment if all methods agree
    ###########################################
    if (all_methods_agree & all_methods_agree_no_overrule) {
      
      inventory = data.frame()
      for (j in 1:length(all_maps)) {
        map = all_maps[[j]]
        
        if (!is.na(map) && length(map$cn_states) >= i && !is.null(map$cn_states[[i]])) {
          method_name = gsub("map_", "", names(all_maps)[j])
          if (method_name %in% names(combined_status)[methods_with_result]) {
            seg = map$cn_states[[i]][[1]]
            inventory = rbind(inventory, data.frame(method=gsub("map_", "", names(all_maps)[j]), major_cn=seg$major_cn, minor_cn=seg$minor_cn))
          }
        }
      }
      
      cn_states[[i]] = inventory
      inventory = na.omit(inventory)
      agree[i] = sum(inventory$major_cn==inventory$major_cn[1] & inventory$minor_cn==inventory$minor_cn[1]) >= min_methods_agree & 
        length(methods_with_result) >= min_methods_with_call_on_segment &
        nrow(inventory)>0 & inventory$major_cn[1] > -1 & inventory$minor_cn[1] > -2
      num_methods[i] = nrow(inventory)
    } else {
      cn_states[[i]] = NA
    }
  }
  
  segments$size = segments$end - segments$start
  frac_genome_agree = round(sum(segments$size[agree]/1000) / sum(segments$size/1000), 2)
  
  return(list(frac_agree=data.frame(samplename=samplename, frac_genome_agree=frac_genome_agree), segments=segments, agree=agree, cn_states=cn_states, num_methods_agree=num_methods))
}


get_all_cn_fits = function(all_data, segment_index, allowed_methods) {
  vanloowedge = NULL
  if (!is.na(all_data$map_vanloowedge)) {
    if (!is.na(all_data$map_vanloowedge$cn_states[[segment_index]])  && !is.na(all_data$map_vanloowedge$cn_states[[segment_index]]) & "vanloowedge" %in% allowed_methods & length(all_data$map_vanloowedge$cn_states) >= segment_index) {
      vanloowedge = all_data$map_vanloowedge$cn_states[[segment_index]]
    }
  }
  
  mustonen = NULL
  if (!is.na(all_data$map_mustonen)) {
    if (!is.na(all_data$map_mustonen$cn_states[[segment_index]])  && !is.na(all_data$map_mustonen$cn_states[[segment_index]]) && "mustonen" %in% allowed_methods & length(all_data$map_mustonen$cn_states) >= segment_index) {
      mustonen = all_data$map_mustonen$cn_states[[segment_index]]
    }
  }
  
  peifer = NULL
  if (!is.na(all_data$map_peifer)) {
    if (!is.na(all_data$map_peifer$cn_states[[segment_index]]) && !is.na(all_data$map_peifer$cn_states[[segment_index]]) && "peifer" %in% allowed_methods & length(all_data$map_peifer$cn_states) >= segment_index) {
      peifer = all_data$map_peifer$cn_states[[segment_index]]
    }
  }
  
  broad = NULL
  if (!is.na(all_data$map_broad)) {
    if (!is.na(all_data$map_broad$cn_states[[segment_index]]) && !is.na(all_data$map_broad$cn_states[[segment_index]]) && "broad" %in% allowed_methods & length(all_data$map_broad$cn_states) >= segment_index) {
      broad = all_data$map_broad$cn_states[[segment_index]]
    }
  } 

  dkfz = NULL  
  if (!is.na(all_data$map_dkfz)) {
    if (!is.na(all_data$map_dkfz$cn_states[[segment_index]]) && !is.na(all_data$map_dkfz$cn_states[[segment_index]]) && "dkfz" %in% allowed_methods & length(all_data$map_dkfz$cn_states) >= segment_index) {
      dkfz = all_data$map_dkfz$cn_states[[segment_index]]
    }
  } 
  
  jabba = NULL  
  if (!is.na(all_data$map_jabba)) {
    if (!is.na(all_data$map_jabba$cn_states[[segment_index]]) && !is.na(all_data$map_jabba$cn_states[[segment_index]]) && "jabba" %in% allowed_methods & length(all_data$map_jabba$cn_states) >= segment_index) {
      jabba = all_data$map_jabba$cn_states[[segment_index]]
    }
  } 
  
  check_missing = function(cn_state) if (!is.null(cn_state)) { return(cn_state[[1]][, c("chromosome", "start", "end", "copy_number", "major_cn", "minor_cn")]) } else { return(NA) } # || !is.na(cn_state)
  
  return(list(vanloowedge=check_missing(vanloowedge),
               mustonen=check_missing(mustonen),
               peifer=check_missing(peifer),
               broad=check_missing(broad),
               dkfz=check_missing(dkfz),
               jabba=check_missing(jabba)))
}



get_frac_genome_agree_maj_vote = function(samplename, all_data, segments, allowed_methods_x_female, allowed_methods_x_male, allowed_methods_y, min_methods_agree=0, min_methods_agree_x=0, min_methods_agree_y=0, min_methods_with_call_on_segment=2, min_methods_with_call_on_segment_x=2, min_methods_with_call_on_segment_y=2, method_overruled=NA) {
  

  # Iterate over the segments and fetch a majority of the methods agrees on the CN state


  # res = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=mustonen_has_header)
  map_dkfz = all_data$map_dkfz
  map_mustonen = all_data$map_mustonen
  map_peifer = all_data$map_peifer
  map_vanloowedge = all_data$map_vanloowedge
  map_broad = all_data$map_broad
  map_jabba = all_data$map_jabba
  all_maps = list(map_broad=map_broad, map_dkfz=map_dkfz, map_mustonen=map_mustonen, map_vanloowedge=map_vanloowedge, map_peifer=map_peifer, map_jabba=map_jabba)
  # combined_status = data.frame(segments, dkfz=map_dkfz$status, mustonen=map_mustonen$status, peifer=map_peifer$status, vanloowedge=map_vanloowedge$status, broad=map_broad$status)
  combined_status = get_combined_status(segments, map_vanloowedge, map_dkfz, map_mustonen, map_peifer, map_broad, map_jabba)

  
  agree = rep(F, nrow(segments))
  num_methods = rep(0, nrow(segments))
  cn_states = list()
  for (i in 1:nrow(segments)) {

    # establish all CN states for this segment - depending on the chromosome
    if (segments$chromosome[i]=="X" & sex=="female") {
      all_states = do.call(rbind, get_all_cn_fits(all_data, i, allowed_methods=allowed_methods_x_female))  
    } else if (segments$chromosome[i]=="X" & sex=="male") {
      all_states = do.call(rbind, get_all_cn_fits(all_data, i, allowed_methods=allowed_methods_x_male))  
    } else if (segments$chromosome[i]=="Y") {
      all_states = do.call(rbind, get_all_cn_fits(all_data, i, allowed_methods=allowed_methods_y))  
    } else {
      all_states = do.call(rbind, get_all_cn_fits(all_data, i, allowed_methods=colnames(combined_status)))  
    }

    # If a method is overruled, remove it from the list here
    if (!is.na(method_overruled)) {
      method_overruled_name = colnames(method_overruled)[unlist(method_overruled[1,])]
      all_states = all_states[!(row.names(all_states) %in% method_overruled_name),]
    }
    
    # Work out what is the majority vote and by how many samples - for only the clonal fits
    methods_with_clonal_solution = colnames(combined_status)[!is.na(combined_status[i,]) & combined_status[i,]=="clonal"]
    # Check that we've got enough methods with a solution - we'll check for enough methods agreeing below
    if (length(methods_with_clonal_solution) >= min_methods_with_call_on_segment) {
      clonal_solutions = all_states[row.names(all_states) %in% methods_with_clonal_solution, ]
      
      # Choose, if a consensus is reached
      combined_states = table(paste(clonal_solutions$major_cn, clonal_solutions$minor_cn, sep="_"))
      # more than 50% of the methods must agree
      min_methods_agreement = length(methods_with_clonal_solution) * 0.5
      majority_vote = which(combined_states > min_methods_agreement & combined_states >= min_methods_with_call_on_segment)
      
      if (length(majority_vote) == 1) {
        majorit_vote_states = as.numeric(unlist(strsplit(names(combined_states)[majority_vote], "_")))
        majority_vote_major = majorit_vote_states[1]
        majority_vote_minor = majorit_vote_states[2]
        
        new_state = clonal_solutions[1,c("major_cn", "minor_cn"), drop=F]
        new_state$major_cn = majority_vote_major
        new_state$minor_cn = majority_vote_minor
        new_state$num_methods = length(methods_with_clonal_solution)
        new_state$num_methods_agree = combined_states[[majority_vote]]
        
        agree[i] = T
        num_methods[i] = combined_states[[majority_vote]]
        cn_states[[i]] = new_state
        
      } else {
        # No suitable calls to create a majority vote
      }
    }
  }

  segments$size = segments$end - segments$start
  frac_genome_agree = round(sum(segments$size[agree]/1000) / sum(segments$size/1000), 2)

  return(list(frac_agree=data.frame(samplename=samplename, frac_genome_agree=frac_genome_agree), segments=segments, agree=agree, num_methods_agree=num_methods, cn_states=cn_states))
}


#' @param min_methods_agree The minimum number of methods that is required to agree
#' @param min_methods_with_call_on_segment The minimum number of methods with a call for a segment to be considered for agreement
#' @param method_overruled A data frame with a single row and a column per method. Each cell contains TRUE if the method is to be overruled
get_frac_genome_agree_rounded = function(samplename, all_data_clonal, all_data_1, all_data_2, all_data_3, segments, do_majority_vote=FALSE, min_methods_agree=0, min_methods_agree_x=0, min_methods_agree_y=0, min_methods_with_call_on_segment=2, min_methods_with_call_on_segment_x=2, min_methods_with_call_on_segment_y=2, method_overruled=NA, allowed_methods_x_female=c("dkfz", "mustonen", "vanloowedge", "jabba", "broad"), allowed_methods_x_male=c("dkfz", "mustonen", "jabba", "broad"), allowed_methods_y=c("dkfz", "jabba")) {
  
  get_all_states = function(all_maps_1, methods_with_result) {
    inventory_1 = data.frame()
    for (j in 1:length(all_maps_1)) {
      map = all_maps_1[[j]]
      
      if (!is.na(map) && length(map$cn_states) >= i && !is.null(map$cn_states[[i]])) {
        method_name = gsub("map_", "", names(all_maps_1)[j])
        if (method_name %in% names(combined_status)[methods_with_result]) {
          seg = map$cn_states[[i]][[1]]
          inventory_1 = rbind(inventory_1, data.frame(method=gsub("map_", "", names(all_maps_1)[j]), major_cn=seg$major_cn, minor_cn=seg$minor_cn))
        }
      }
    }
    return(inventory_1)
  }
  
  
  greedy_get_allele = function(inventory_1, inventory_2, min_methods_agree, min_methods_with_call_on_segment, do_majority_vote) {
    # Greedy algorithm to find the single most prevalent state and make that into a copy number call
    allele_inventory = table(c(inventory_1$major_cn, inventory_1$minor_cn, inventory_2$major_cn[! inventory_2$method %in% clonal_votes], inventory_2$minor_cn[! inventory_2$method %in% clonal_votes]))
    
    num_methods = length(inventory_1$method)
    hits_inventory = data.frame(allele=as.numeric(names(allele_inventory)))
    scoring = as.data.frame(matrix(0, nrow=length(allele_inventory), ncol=num_methods))
    colnames(scoring) = inventory_1$method
    hits_inventory = data.frame(hits_inventory, scoring)
    for (j in 1:length(allele_inventory)) {
      allele = hits_inventory$allele[j]
      for (method in inventory_1$method) {
        
        inv_1_method = inventory_1[inventory_1$method==method,]
        inv_2_method = inventory_2[inventory_2$method==method,]
        cn_state_1 = inv_1_method$major_cn==allele | inv_1_method$minor_cn==allele
        cn_state_2 = inv_2_method$major_cn[!method %in% clonal_votes]==allele | inv_2_method$minor_cn[!method %in% clonal_votes]==allele
        
        if (method %in% clonal_votes) {
          hits_inventory[j, method] = as.numeric(cn_state_1)
        } else {
          hits_inventory[j, method] = as.numeric(cn_state_1 | cn_state_2)
        }
      }
    }
    
    # Get the qualified hits (i.e. supported by more than 50% of the methods) and take the best one
    num_hits = rowSums(hits_inventory[,-1])
    # There are two options: By passing the minimum number of methods or by majority vote with a minimum number of methods reporting
    if (!do_majority_vote) {
      qualified_hits = num_hits > min_methods_agree
    } else {
      qualified_hits = num_hits > ceiling(num_methods*0.5) & num_methods > min_methods_with_call_on_segment
    }
    
    if (!any(qualified_hits)) {
      best_hit = NA
    } else {
      best_hit = max(num_hits[qualified_hits])
    }
    
    if (!is.na(best_hit) && sum(num_hits==best_hit)==1) {
      # One result, save that as allele 1 and now look for 2
      allele_1 = hits_inventory[which(num_hits==best_hit), 1]
      
      # We have one state shared by all methods, fix this and find the other
      inv_1_hit = inventory_1[inventory_1$major_cn==allele_1 | inventory_1$minor_cn==allele_1,]
      inv_2_hit = inventory_2[(inventory_2$major_cn==allele_1 | inventory_2$minor_cn==allele_1) & !inventory_2$method %in% clonal_votes,]
      
      other_allele_options = c()
      if (nrow(inv_1_hit) > 0) {
        for (j in 1:nrow(inv_1_hit)) {
          if (inv_1_hit$major_cn[j]==allele_1) {
            other_allele_options = c(other_allele_options, inv_1_hit$minor_cn[j])
          } else {
            other_allele_options = c(other_allele_options, inv_1_hit$major_cn[j])
          }
        }
      }
      
      if (nrow(inv_2_hit) > 0) {
        for (j in 1:nrow(inv_2_hit)) {
          if (inv_2_hit$major_cn[j]==allele_1) {
            other_allele_options = c(other_allele_options, inv_2_hit$minor_cn[j])
          } else {
            other_allele_options = c(other_allele_options, inv_2_hit$major_cn[j])
          }
        }
      }
      
      other_allele_options = table(other_allele_options)
      best_hit_other = which.max(other_allele_options)
      if (sum(other_allele_options==other_allele_options[best_hit_other])==1) {
        allele_2 = other_allele_options[best_hit_other]
      } else {
        allele_2 = NA
      }
      
    } else {
      # There is no single best allele 1
      allele_1 = NA
      allele_2 = NA
    }
    return(list(allele_1=allele_1, allele_2=allele_2))
  }
  
  
  
  
  # breakpoints = read.table(paste0(samplename, "_consensus_breakpoints.txt"), header=T, stringsAsFactors=F)
  # segments = breakpoints2segments(breakpoints)
  
  # res = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=mustonen_has_header)
  map_dkfz_1 = all_data_1$map_dkfz
  map_mustonen_1 = all_data_1$map_mustonen
  map_peifer_1 = all_data_1$map_peifer
  map_vanloowedge_1 = all_data_1$map_vanloowedge
  map_broad_1 = all_data_1$map_broad
  map_jabba_1 = all_data_1$map_jabba
  all_maps_1 = list(map_broad=map_broad_1, map_dkfz=map_dkfz_1, map_mustonen=map_mustonen_1, map_vanloowedge=map_vanloowedge_1, map_peifer=map_peifer_1, map_jabba=map_jabba_1)
  
  map_dkfz_2 = all_data_2$map_dkfz
  map_mustonen_2 = all_data_2$map_mustonen
  map_peifer_2 = all_data_2$map_peifer
  map_vanloowedge_2 = all_data_2$map_vanloowedge
  map_broad_2 = all_data_2$map_broad
  map_jabba_2 = all_data_2$map_jabba
  all_maps_2 = list(map_broad=map_broad_2, map_dkfz=map_dkfz_2, map_mustonen=map_mustonen_2, map_vanloowedge=map_vanloowedge_2, map_peifer=map_peifer_2, map_jabba=map_jabba_2)

  # Use the non-rounded data to build the status
  combined_status = get_combined_status(segments, all_data_clonal$map_vanloowedge, all_data_clonal$map_dkfz, all_data_clonal$map_mustonen, all_data_clonal$map_peifer, all_data_clonal$map_broad, all_data_clonal$map_jabba)  
  # Order both combined_status and method_overruled to have the same order - if method_overruled is supplied
  if (!is.na(method_overruled)) {
    method_overruled = method_overruled[colnames(combined_status)[4:9]]
  }
  
  agree = rep(F, nrow(segments))
  cn_states = list()
  num_methods = rep(0, nrow(segments))
  one_allele_saved = rep(NA, nrow(segments))
  for (i in 1:nrow(segments)) {
    clonal_votes = colnames(combined_status)[which(combined_status[i,] == "clonal")]
    ###########################################
    # Do different things when addressing X and Y because there are fewer methods reporting
    ###########################################
    # Order of methods is: dkfz, mustonen, peifer, vanloowedge, broad, jabba
    if (segments$chromosome[i]=="X") {
      if (sex=="female") { allowed_methods_x = allowed_methods_x_female
      } else { allowed_methods_x = allowed_methods_x_male }
      selection = !is.na(combined_status[i,4:9]) & colnames(combined_status)[4:9] %in% allowed_methods_x
      # These can be overwritten because X and Y are always the last chromosomes to be considered
      min_methods_with_call_on_segment = min_methods_with_call_on_segment_x
      min_methods_agree = min_methods_agree_x
      
    } else if (segments$chromosome[i]=="Y") {
      selection = !is.na(combined_status[i,4:9]) & colnames(combined_status)[4:9] %in% allowed_methods_y
      # These can be overwritten because X and Y are always the last chromosomes to be considered
      min_methods_with_call_on_segment = min_methods_with_call_on_segment_y
      min_methods_agree = min_methods_agree_y
      
    } else {
      selection = !is.na(combined_status[i,4:9])
    }
    
    # Check if method overruling is requested
    if (!is.na(method_overruled)) {
      selection = selection & !method_overruled
    }
    # Finally select the indices that correspond the methods with results
    methods_with_result = (4:9)[selection]
    
    #' If all methods think it's clonal there is no rounded majority vote
    if (length(methods_with_result)==0 || all(combined_status[i,methods_with_result]=="clonal", na.rm=T)) {
      cn_states[[i]] = NA
      next
    }
    
    ###########################################
    # Fetch all the CN states for this segment
    ###########################################    
    inventory_1 = get_all_states(all_maps_1, methods_with_result)
    inventory_2 = get_all_states(all_maps_2, methods_with_result)
    
    ###########################################
    # Perform the voting - option 1
    ###########################################    
    major_configs = paste(inventory_1$major_cn[! inventory_1$method %in% clonal_votes], inventory_1$minor_cn[! inventory_1$method %in% clonal_votes], sep="_")
    minor_configs = paste(inventory_2$major_cn[! inventory_2$method %in% clonal_votes], inventory_2$minor_cn[! inventory_2$method %in% clonal_votes], sep="_")
    clonal_configs = paste(inventory_1$major_cn[inventory_1$method %in% clonal_votes], inventory_1$minor_cn[inventory_1$method %in% clonal_votes], sep="_")
    votes = table(c(major_configs, minor_configs, clonal_configs))
    best_vote = max(votes)
    if (!do_majority_vote & best_vote < min_methods_agree) {
      # Did not meet threshold for enough methods agree
      best_vote = NA
    } else if (do_majority_vote & nrow(inventory_1) < min_methods_with_call_on_segment) {
      # Not enough methods reporting for a reliable majority vote
      best_vote = NA
    }
    
    if (sum(votes==best_vote, na.rm=T) == 1) {    
      # We've found a majority vote, save it and done
      maj = as.numeric(unlist(stringr::str_split(names(votes)[votes==best_vote], "_"))[1])
      min = as.numeric(unlist(stringr::str_split(names(votes)[votes==best_vote], "_"))[2])
      agree[i] = T
      combination_method = "vote_combined_states"
      
    } else {
      ###########################################
      # A tie - Perform the voting - option 2
      ###########################################    
      res = greedy_get_allele(inventory_1, inventory_2, min_methods_agree, min_methods_with_call_on_segment, do_majority_vote)
      
      # If both alleles are now defined we have consensus
      if (!is.na(res$allele_1) & !is.na(res$allele_2)) {
        agree[i] = T
        maj = max(c(res$allele_1, res$allele_2))
        min = min(c(res$allele_1, res$allele_2))
        combination_method = "vote_single_state"
      } else if (!is.na(res$allele_1)) {
        one_allele_saved[i] = res$allele_1
      }
    }
    
    
    if (agree[i]) {
      new_state = inventory_1[1,c("method", "major_cn", "minor_cn"), drop=F]
      new_state$method = combination_method
      new_state$major_cn = maj
      new_state$minor_cn = min
      new_state$num_methods_agree = sum((inventory_1$major_cn==maj & inventory_1$minor_cn==min) | (inventory_2$major_cn==maj & inventory_2$minor_cn==min))
      
      
      num_methods[i] = new_state$num_methods_agree
      cn_states[[i]] = new_state
    } else {
      cn_states[[i]] = NA
    }
  }  
  
  segments$size = segments$end - segments$start
  frac_genome_agree = round(sum(segments$size[agree]/1000) / sum(segments$size/1000), 2)
  
  return(list(frac_agree=data.frame(samplename=samplename, frac_genome_agree=frac_genome_agree), segments=segments, agree=agree, cn_states=cn_states, num_methods_agree=num_methods, one_allele_saved=one_allele_saved))
}

test_purities = function(purities, consensus_profile) {
  
  # rBacktransform = function(rho, nA, nB, psi) {
  #   return(gamma_param*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2))
  # }
  
  bBacktransform = function(rho, nA, nB) {
    return((1-rho+rho*nB)/(2-2*rho+rho*(nA+nB)))
  }
  
  # Calculates a confidence in the scale of 0-100 with 100 being best. Deviation from the expected baf is literally rescaled to this range. 
  # Would normally expect confidence values in the high 90s for clonal aberrations
  bConf = function(btsm, b) {
    return(ifelse(btsm!=0.5, pmin(100, pmax(0, ifelse(b==0.5, 100, 100*(1-abs(btsm-b)/abs(b-0.5))))), NA))
  }
  
  
  # Test all star 3 segments
  star3_segments = which(consensus_profile$star==3 & !is.na(consensus_profile$major_cn) & !is.na(consensus_profile$minor_cn))
  if (length(star3_segments)==0) {
    return(NA)
  } 

  # Calc expected baf given the fit and purity and calculate the confidence
  bhat = do.call(cbind, lapply(purities[1,], function(purity) { 1-bBacktransform(purity, consensus_profile$major_cn[star3_segments], consensus_profile$minor_cn[star3_segments]) }))
  bConf_vlw = do.call(cbind, lapply(1:ncol(bhat), function(i) { round(bConf(bhat[,i], consensus_profile$vanloowedge_baf[star3_segments]), 4) }))
  bConf_broad = do.call(cbind, lapply(1:ncol(bhat), function(i) { round(bConf(bhat[,i], consensus_profile$broad_baf[star3_segments]), 4) }))
  
  output = matrix(NA, 1, ncol(bConf_vlw)*2)
  output[1,] = c(apply(bConf_vlw, 2, median, na.rm=T), apply(bConf_broad, 2, median, na.rm=T))
  output = as.data.frame(output)
  method_names = unlist(lapply(colnames(purities), function(x) unlist(strsplit(x, "_"))[2]))
  colnames(output) = c(paste(method_names, "bafConfVanloowedge", sep="_"),
                       paste(method_names, "bafConfBroad", sep="_"))
  output$numsegments_bafConf = length(star3_segments)
  return(output)
}

get_ploidy = function(segments, map, broad=F) {
  if (!is.na(map)) {
    cn_bb = collapse2bb(segments=segments, cn_states=map$cn_states, broad=broad)
    return(list(ploidy=round(calc_ploidy(cn_bb), 4), status=get_ploidy_status(cn_bb)))
  } else {
    return(list(ploidy=NA, status=NA))
  }
}

#####################################################################
# Original agreement
#####################################################################
# setwd("~/Documents/Projects/icgc/consensus_subclonal_copynumber/final_run_testing/")

library(readr)
source("~/repo/icgc_consensus_copynumber_final/util.R")
max.plot.cn=4
num_threads=1


args = commandArgs(T)
samplename = args[1]
outdir = args[2]
sex = args[3]

# setwd("/Users/sd11/Documents/Projects/icgc/consensus_subclonal_copynumber/6aa00162-6294-4ce7-b6b7-0c3452e24cd6")
# outdir = "./"
# samplename = "6aa00162-6294-4ce7-b6b7-0c3452e24cd6"

# setwd("/Users/sd11/Documents/Projects/icgc/consensus_subclonal_copynumber/final_run_testing")
# samplename = "005e85a3-3571-462d-8dc9-2babfc7ace21"
# sex = "female"
# outdir = "output"


# samplename = "04b9837e-9ab5-4eb9-9a9c-ef49e3a62662"
# sex = "male"


breakpoints_file = file.path("consensus_bp", paste0(samplename, ".txt"))
# expected_ploidy_file = "consensus.20161103.purity.ploidy.txt.gz" # Removed after ploidy has been reinferred after fixes
expected_ploidy_file = "icgc_pcawg_reference_ploidy_final_alpha.txt"
marked_unknown_file = "icgc_marked_uknown.lst"
purity_based_overrulings_file = "icgc_purity_and_ploidy_overrulings.txt"
# the reference ploidy is multiplied by this factor to determine how much of a deviation is tolerated
max_expected_ploidy_diff_factor = 0.25
allowed_methods_x_female = c("dkfz", "mustonen", "vanloowedge", "jabba", "broad")
allowed_methods_x_male = c("dkfz", "mustonen", "broad")
allowed_methods_y = c("dkfz", "jabba", "broad")

# Table with overrulings 
# overrulings_pivot = readr::read_tsv("~/Documents/Projects/icgc/consensus_subclonal_copynumber/manual_review_overrulings_pivot_table.txt")
overrulings_pivot = readr::read_tsv("manual_review_overrulings_pivot_table.txt")
overrulings_pivot = overrulings_pivot[overrulings_pivot$samplename==samplename,]

if (nrow(overrulings_pivot)==1 & sum(!is.na(overrulings_pivot)) > 1) {
  method_overruled = as.data.frame(!is.na(overrulings_pivot[,2:6]))
} else {
  method_overruled = NA
}

if (file.exists(breakpoints_file)) {
  print("Reading in and mapping data...")
  
  breakpoints = read.table(file.path("consensus_bp", paste0(samplename, ".txt")), header=T, stringsAsFactors=F)
  segments = breakpoints2segments(breakpoints)
  
  dkfz_segmentsfile = paste0("dkfz/", samplename, "_segments.txt")
  dkfz_purityfile = "purity_ploidy_dkfz.txt"
  vanloowedge_segmentsfile = paste0("vanloowedge/", samplename, "_segments.txt")
  vanloowedge_purityfile = "purity_ploidy_vanloowedge.txt"
  peifer_segmentsfile = paste0("peifer/", samplename, "_segments.txt")
  peifer_purityfile = "purity_ploidy_peifer.txt"
  mustonen_segmentsfile = paste0("mustonen/", samplename, "_segments.txt")
  mustonen_purityfile = "purity_ploidy_mustonen.txt"
  broad_segmentsfile = paste0("broad/", samplename, "_segments.txt")
  broad_purityfile = "purity_ploidy_broad.txt"
  jabba_segmentsfile = paste0("jabba/", samplename, "_segments.txt")
  jabba_purityfile = "purity_ploidy_jabba.txt"
  vanloowedge_baflogrfile = paste0("vanloowedge_baflogr/", samplename, "_baflogr.txt")
  broad_baflogrfile = paste0("broad_baflogr/", samplename, "_baflogr.txt")
  
  method_segmentsfile = list(dkfz=dkfz_segmentsfile,
                             vanloowedge=vanloowedge_segmentsfile,
                             peifer=peifer_segmentsfile,
                             mustonen=mustonen_segmentsfile,
                             broad=broad_segmentsfile,
                             jabba=jabba_segmentsfile)
  
  method_purityfile = list(dkfz=dkfz_purityfile,
                           vanloowedge=vanloowedge_purityfile,
                           peifer=peifer_purityfile,
                           mustonen=mustonen_purityfile,
                           broad=broad_purityfile,
                           jabba=jabba_purityfile)
  
  method_baflogr = list(vanloowedge=vanloowedge_baflogrfile,
                        broad=broad_baflogrfile)
  
  all_data_clonal = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, method_baflogr, sex=sex, mustonen_has_header=F, num_threads=num_threads)
  
  #####################################################################
  # Overrule methods
  #####################################################################
  # Calc ploidy of all profiles and overrule those that are not concordant
  ploidy_vanloowedge = get_ploidy(segments, all_data_clonal$map_vanloowedge)
  ploidy_broad = get_ploidy(segments, all_data_clonal$map_broad, broad=T)
  ploidy_peifer = get_ploidy(segments, all_data_clonal$map_peifer)
  ploidy_dkfz = get_ploidy(segments, all_data_clonal$map_dkfz)
  ploidy_mustonen = get_ploidy(segments, all_data_clonal$map_mustonen)
  ploidy_jabba = get_ploidy(segments, all_data_clonal$map_jabba)
  
  expected_ploidy = read.table(expected_ploidy_file, header=T, stringsAsFactors=F)
  expected_ploidy = expected_ploidy[expected_ploidy$samplename==samplename, "ploidy"]
  max_expected_ploidy_diff = expected_ploidy*max_expected_ploidy_diff_factor
  overrulings = list(broad=F, mustonen=F, dkfz=F, peifer=F, vanloowedge=F, jabba=F)
  # Compare to expected
  if (length(expected_ploidy) > 0 && !is.na(expected_ploidy)) {
    if (!is.na(ploidy_vanloowedge$ploidy) && abs(ploidy_vanloowedge$ploidy-expected_ploidy) > max_expected_ploidy_diff) {
      print("Overruling Battenberg ploidy")
      all_data_clonal$map_vanloowedge = NA
      all_data_clonal$dat_vanloowedge = NA
      overrulings$vanloowedge = T
    }
    
    if (!is.na(ploidy_broad$ploidy) && abs(ploidy_broad$ploidy-expected_ploidy) > max_expected_ploidy_diff) {
      print("Overruling ABSOLUTE ploidy")
      all_data_clonal$map_broad = NA
      all_data_clonal$dat_broad = NA
      overrulings$broad = T
    }
    
    if (!is.na(ploidy_dkfz$ploidy) && abs(ploidy_dkfz$ploidy-expected_ploidy) > max_expected_ploidy_diff) {
      print("Overruling ACEseq ploidy")
      all_data_clonal$map_dkfz = NA
      all_data_clonal$dat_dkfz = NA
      overrulings$dkfz = T
    }
    
    if (!is.na(ploidy_peifer$ploidy) && abs(ploidy_peifer$ploidy-expected_ploidy) > max_expected_ploidy_diff) {
      print("Overruling Sclust ploidy")
      all_data_clonal$map_peifer = NA
      all_data_clonal$dat_peifer = NA
      overrulings$peifer = T
    }
    
    if (!is.na(ploidy_mustonen$ploidy) && abs(ploidy_mustonen$ploidy-expected_ploidy) > max_expected_ploidy_diff) {
      print("Overruling CloneHD ploidy")
      all_data_clonal$map_mustonen = NA
      all_data_clonal$dat_mustonen = NA
      overrulings$mustonen = T
    }
    
    if (!is.na(ploidy_jabba$ploidy) && abs(ploidy_jabba$ploidy-expected_ploidy) > max_expected_ploidy_diff) {
      print("Overruling JaBbA ploidy")
      all_data_clonal$map_jabba = NA
      all_data_clonal$dat_jabba = NA
      overrulings$jabba = T
    }
  }
  method_overruled = data.frame(t(data.frame(unlist(overrulings))), stringsAsFactors=F)
  row.names(method_overruled) = NULL
  
  # Add the purity based overrulings
  purity_based_overrulings = readr::read_tsv(purity_based_overrulings_file)
  purity_based_overrulings = purity_based_overrulings[purity_based_overrulings$samplename==samplename, names(method_overruled)]
  method_overruled = data.frame(method_overruled | purity_based_overrulings)
  
  if (all(as.logical(method_overruled[1,]))) {
    print("All methods overruled, quitting now")
    print(data.frame(ploidy_vanloowedge=ploidy_vanloowedge$ploidy, ploidy_broad=ploidy_broad$ploidy, ploidy_dkfz=ploidy_dkfz$ploidy, ploidy_peifer=ploidy_peifer$ploidy, ploidy_mustonen=ploidy_mustonen$ploidy, ploidy_jabba=ploidy_jabba$ploidy, reference=expected_ploidy))
    q(save="no")
  }
  
  # Remove those that are overruled
  if (any(as.logical(method_overruled))) {
    for (method in names(method_overruled)[as.logical(method_overruled)]) {
      all_data_clonal[grepl(method, names(all_data_clonal))] = NA
    }
  }
  
  #####################################################################
  # Clonal complete agreement
  #####################################################################
  print("Getting clonal agreement...")
  agreement_clonal = get_frac_genome_agree(samplename, 
                                           all_data_clonal, 
                                           segments, 
                                           min_methods_agree=6, 
                                           min_methods_agree_x=ifelse(sex=="male", length(allowed_methods_x_male),  length(allowed_methods_x_female)),
                                           min_methods_agree_y=length(allowed_methods_y), 
                                           allowed_methods_x_female=allowed_methods_x_female, 
                                           allowed_methods_x_male=allowed_methods_x_male, 
                                           allowed_methods_y=allowed_methods_y)
  
  #####################################################################
  # Agreement exclude 1
  #####################################################################
  print("Getting exclude 1 agreement...")
  agreement_clonal_exclude_1 = get_frac_genome_agree(samplename, 
                                                     all_data_clonal, 
                                                     segments, 
                                                     min_methods_agree=5, 
                                                     min_methods_agree_x=ifelse(sex=="male", length(allowed_methods_x_male)-1,  length(allowed_methods_x_female)-1),
                                                     min_methods_agree_y=length(allowed_methods_y), 
                                                     allowed_methods_x_female=allowed_methods_x_female, 
                                                     allowed_methods_x_male=allowed_methods_x_male, 
                                                     allowed_methods_y=allowed_methods_y)
  
  #####################################################################
  # Agreement after excluding overruled profiles
  #####################################################################
  print("Getting exclude overruled agreement...")
  # Check that there is an entry and that there is at least one method overruled
  # In this case we exclude the overruled methods, and we accept if all others agree if there are a required
  # minimum number of methods with a call
  # if (nrow(overrulings_pivot)==1 & sum(!is.na(overrulings_pivot)) > 1) {
  if (any(unlist(overrulings))) {
    agreement_clonal_overrule = get_frac_genome_agree(samplename, 
                                                      all_data_clonal, 
                                                      segments, 
                                                      method_overruled=method_overruled, 
                                                      min_methods_with_call_on_segment=3, 
                                                      min_methods_agree=sum(!is.na(method_overruled)),
                                                      min_methods_agree_x=2, # Keep X and Y steady
                                                      min_methods_agree_y=2,
                                                      allowed_methods_x_female=allowed_methods_x_female, 
                                                      allowed_methods_x_male=allowed_methods_x_male, 
                                                      allowed_methods_y=allowed_methods_y)
  } else {
    # No methods have been overruled for this sample - so clonal agreement it is
    agreement_clonal_overrule = agreement_clonal
  }
  
  #####################################################################
  # Agreement after rounding
  #####################################################################
  print("Loading rounded...")
  dkfz_segmentsfile = file.path(outdir, "dkfz_rounded_clonal", paste0(samplename, "_segments.txt"))
  vanloowedge_segmentsfile = file.path(outdir, "vanloowedge_rounded_clonal", paste0(samplename, "_segments.txt"))
  peifer_segmentsfile = file.path(outdir, "peifer_rounded_clonal", paste0(samplename, "_segments.txt"))
  mustonen_segmentsfile = file.path(outdir, "mustonen_rounded_clonal", paste0(samplename, "_segments.txt"))
  broad_segmentsfile = file.path(outdir, "broad_rounded_clonal", paste0(samplename, "_segments.txt"))
  jabba_segmentsfile = file.path(outdir, "jabba_rounded_clonal", paste0(samplename, "_segments.txt"))
  
  method_segmentsfile = list(dkfz=ifelse(!method_overruled$dkfz, dkfz_segmentsfile, NA),
                             vanloowedge=ifelse(!method_overruled$vanloowedge, vanloowedge_segmentsfile, NA),
                             peifer=ifelse(!method_overruled$peifer, peifer_segmentsfile, NA),
                             mustonen=ifelse(!method_overruled$mustonen, mustonen_segmentsfile, NA),
                             broad=ifelse(!method_overruled$broad, broad_segmentsfile, NA),
                             jabba=ifelse(!method_overruled$jabba, jabba_segmentsfile, NA))
  
  all_data_rounded = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, method_baflogr=NULL, sex=sex, mustonen_has_header=T, num_threads=num_threads)

  print("Loading alt rounded...")
  dkfz_segmentsfile = file.path(outdir, "dkfz_rounded_alt_clonal/", paste0(samplename, "_segments.txt"))
  vanloowedge_segmentsfile = file.path(outdir, "vanloowedge_rounded_alt_clonal/", paste0(samplename, "_segments.txt"))
  peifer_segmentsfile = file.path(outdir, "peifer_rounded_alt_clonal/", paste0(samplename, "_segments.txt"))
  mustonen_segmentsfile = file.path(outdir, "mustonen_rounded_alt_clonal/", paste0(samplename, "_segments.txt"))
  broad_segmentsfile = file.path(outdir, "broad_rounded_alt_clonal/", paste0(samplename, "_segments.txt"))
  jabba_segmentsfile = file.path(outdir, "jabba_rounded_alt_clonal", paste0(samplename, "_segments.txt"))
  
  method_segmentsfile = list(dkfz=ifelse(!method_overruled$dkfz, dkfz_segmentsfile, NA),
                             vanloowedge=ifelse(!method_overruled$vanloowedge, vanloowedge_segmentsfile, NA),
                             peifer=ifelse(!method_overruled$peifer, peifer_segmentsfile, NA),
                             mustonen=ifelse(!method_overruled$mustonen, mustonen_segmentsfile, NA),
                             broad=ifelse(!method_overruled$broad, broad_segmentsfile, NA),
                             jabba=ifelse(!method_overruled$jabba, jabba_segmentsfile, NA))
  
  all_data_rounded_alt = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, method_baflogr=NULL, sex=sex, mustonen_has_header=T, num_threads=num_threads)
  
  print("Getting rounded agreement...")
  agreement_rounded = get_frac_genome_agree_rounded(samplename, 
                                                    all_data_clonal, 
                                                    all_data_rounded, 
                                                    all_data_rounded_alt, 
                                                    NULL, 
                                                    segments, 
                                                    min_methods_agree=4, 
                                                    min_methods_agree_x=2, 
                                                    min_methods_agree_y=2, 
                                                    method_overruled=method_overruled, 
                                                    allowed_methods_x_female=allowed_methods_x_female, 
                                                    allowed_methods_x_male=allowed_methods_x_male, 
                                                    allowed_methods_y=allowed_methods_y)
  
  #####################################################################
  # Agreement with accepting majority vote
  #####################################################################
  print("Getting majority vote agreement...")
  agreement_clonal_majority_vote = get_frac_genome_agree_maj_vote(samplename, 
                                                                   all_data_clonal, 
                                                                   segments, 
                                                                   method_overruled=method_overruled,
                                                                   min_methods_agree=4,
                                                                   min_methods_agree_x=2, 
                                                                   min_methods_agree_y=2, 
                                                                   allowed_methods_x_female=allowed_methods_x_female, 
                                                                   allowed_methods_x_male=allowed_methods_x_male, 
                                                                   allowed_methods_y=allowed_methods_y)
  
  agreement_rounded_majority_vote = get_frac_genome_agree_rounded(samplename, 
                                                                all_data_clonal, 
                                                                all_data_rounded, 
                                                                all_data_rounded_alt, 
                                                                NULL, 
                                                                segments, 
                                                                method_overruled=method_overruled, 
                                                                allowed_methods_x_female=allowed_methods_x_female, 
                                                                allowed_methods_x_male=allowed_methods_x_male, 
                                                                allowed_methods_y=allowed_methods_y,
                                                                do_majority_vote=TRUE)
  
  #####################################################################
  # Minority agreement - better than a single method
  #####################################################################
  # agreement_clonal_minority = get_frac_genome_agree(samplename, 
  #                                          all_data_clonal, 
  #                                          segments, 
  #                                          min_methods_agree=2, 
  #                                          min_methods_agree_x=ifelse(sex=="male", length(allowed_methods_x_male),  length(allowed_methods_x_female)),
  #                                          min_methods_agree_y=length(allowed_methods_y), 
  #                                          allowed_methods_x_female=allowed_methods_x_female, 
  #                                          allowed_methods_x_male=allowed_methods_x_male, 
  #                                          allowed_methods_y=allowed_methods_y)
  # 
  # agreement_rounded_minority = get_frac_genome_agree_rounded(samplename, 
  #                                                   all_data_clonal, 
  #                                                   all_data_rounded, 
  #                                                   all_data_rounded_alt, 
  #                                                   NULL, 
  #                                                   segments, 
  #                                                   min_methods_agree=2, 
  #                                                   min_methods_agree_x=2, 
  #                                                   min_methods_agree_y=2, 
  #                                                   method_overruled=method_overruled, 
  #                                                   allowed_methods_x_female=allowed_methods_x_female, 
  #                                                   allowed_methods_x_male=allowed_methods_x_male, 
  #                                                   allowed_methods_y=allowed_methods_y)
  
  #####################################################################
  # Piece together a complete agreement profile
  #####################################################################
  create_consensus_profile = function(segments, agreement_clonal, agreement_clonal_exclude_1, agreement_clonal_overrule, agreement_clonal_majority_vote, agreement_rounded, agreement_rounded_majority_vote, map_broad_baflogr, map_vanloowedge_baflogr) {
    consensus_profile = data.frame()
    r = agreement_rounded$cn_states
    # for (i in 1:length(r)) {
    for (i in 1:nrow(segments)) {
      if (agreement_clonal$agree[i]) {
        # if clonal agree, choose that and assign 3*
        new_entry = agreement_clonal$cn_states[[i]][1,2:3]
        new_entry$star = 3
        new_entry$level = "a"
        new_entry$methods_agree = agreement_clonal$num_methods_agree[i]
        
      } else if (agreement_clonal_exclude_1$agree[i]) {
        # if clonal agree except 1, choose that and assign 3*
        new_entry = agreement_clonal_exclude_1$cn_states[[i]][1,2:3]
        new_entry$star = 3
        new_entry$level = "b" 
        new_entry$methods_agree = agreement_clonal_exclude_1$num_methods_agree[i]
        
      } else if (agreement_clonal_overrule$agree[i]) {
        # pivot table
        new_entry = agreement_clonal_overrule$cn_states[[i]][1,2:3]
        if (agreement_clonal_overrule$num_methods_agree[i] %in% c(5,6)) {
          new_entry$star = 3
          new_entry$level = "c"
        } else {
          new_entry$star = 2
          new_entry$level = "c"
        }
        new_entry$methods_agree = agreement_clonal_overrule$num_methods_agree[i]
        
      } else if (agreement_clonal_majority_vote$agree[i]) {
        # else if rounded clonal agree, choose that and assign 2*
        new_entry = agreement_clonal_majority_vote$cn_states[[i]][1,1:2]
        new_entry$star = 2
        new_entry$level = "d"
        new_entry$methods_agree = agreement_clonal_majority_vote$num_methods_agree[i]  
        
      } else if (agreement_rounded$agree[i]) {
        # else if rounded clonal agree, choose that and assign 2*
        new_entry = agreement_rounded$cn_states[[i]][1,2:3]
        new_entry$star = 2
        new_entry$level = "e"
        new_entry$methods_agree = agreement_rounded$num_methods_agree[i]
        
      } else if (agreement_rounded_majority_vote$agree[i]) {
        # majority vote by > 50% of the methods 
        new_entry = agreement_rounded_majority_vote$cn_states[[i]][1,2:3]
        new_entry$star = 2
        new_entry$level = "f"
        row.names(new_entry) = NULL
        new_entry$methods_agree = agreement_rounded_majority_vote$num_methods_agree[i]
        
      } else {
        # else no solution for now, below will select the best method for this sample
        new_entry = data.frame(major_cn=NA, minor_cn=NA)
        new_entry$star = 1
        new_entry$level = "g"
        new_entry$methods_agree = NA
      }
      
      
      if (!is.na(map_broad_baflogr) && !is.null(map_broad_baflogr$cn_states[[i]]) && !is.na(map_broad_baflogr$cn_states[[i]])) {
        new_entry$broad_baf = map_broad_baflogr$cn_states[[i]][[1]][1,4]
        new_entry$broad_logr = map_broad_baflogr$cn_states[[i]][[1]][1,5]
      } else {
        new_entry$broad_baf = NA
        new_entry$broad_logr = NA
      }
      
      if (!is.na(map_vanloowedge_baflogr) && !is.null(map_vanloowedge_baflogr$cn_states[[i]]) && !is.na(map_vanloowedge_baflogr$cn_states[[i]])) {
        new_entry$vanloowedge_baf = map_vanloowedge_baflogr$cn_states[[i]][[1]][1,4]
        new_entry$vanloowedge_logr = map_vanloowedge_baflogr$cn_states[[i]][[1]][1,5]
      } else {
        new_entry$vanloowedge_baf = NA
        new_entry$vanloowedge_logr = NA
      }
      consensus_profile = rbind(consensus_profile, new_entry)
    }
    return(consensus_profile)
  }
  
  print("Building initial consensus profile...")
  consensus_profile = create_consensus_profile(segments, 
                                               agreement_clonal, 
                                               agreement_clonal_exclude_1, 
                                               agreement_clonal_overrule, 
                                               agreement_clonal_majority_vote,
                                               agreement_rounded, 
                                               agreement_rounded_majority_vote, 
                                               all_data_clonal$map_broad_baflogr, 
                                               all_data_clonal$map_vanloowedge_baflogr)
  
  print("Making figure input")
  profile_bb = collapseRoundedClonal2bb(data.frame(segments, consensus_profile))
  print("Making figure")
  p = plot_profile(profile_bb, "Consensus - after rounding", max.plot.cn=max.plot.cn)
  png(file.path(outdir, "figures", paste0(samplename, "_consensus_rounded.png")), height=300, width=1300)
  print(p)
  dev.off()
  print("Made figure")
  
  #####################################################################
  # Calc agreement with both the created consensus per method
  #####################################################################
  
  calc_method_agreement = function(all_maps, segments, consensus_profile, segment_status) {
    agreement = list(dkfz=0, vanloowedge=0, peifer=0, mustonen=0, broad=0, jabba=0)
    # collect_mustonen_2 = c()
    for (i in 1:nrow(consensus_profile)) {
      # Exclude X and Y as not all methods report it and it would skew the metric
      if (!(segments$chromosome[i] %in% c("X", "Y"))) {
        # Iterate over all the segment mappings
        for (j in which(grepl("map", names(all_maps)) & !grepl("baflogr", names(all_maps)))) {
          method = gsub("map_", "", names(all_maps)[j])
          
          if (!is.na(all_maps[[j]]) && !is.na(all_maps[[j]]$status[i]) && all_maps[[j]]$status[i]==segment_status) {
            
            # Check if agreement
            # TODO bugfix : missing value where TRUE/FALSE needed in if statement
            method_segment = all_maps[[j]]$cn_states[[i]][[1]]
            if (!is.null(method_segment) && !is.na(consensus_profile$major_cn[i]) && (!is.na(method_segment$major_cn) | !is.na(method_segment$minor_cn)) &&
                consensus_profile$major_cn[i]==method_segment$major_cn & consensus_profile$minor_cn[i]==method_segment$minor_cn) {
              # Add the extra agreement to the tally of this method
              # if (method=="mustonen") { collect_mustonen_2 = c(collect_mustonen_2, i) }
              agreement[[method]] = agreement[[method]] + (segments$end[i]/1000-segments$start[i]/1000)
            }
          }
        }
      }
    }
    
    genome_size = sum(segments$end/1000 - segments$start/1000)
    frac_agreement = lapply(agreement, function(x) x / genome_size)
    return(frac_agreement)
  }
  print("Calculating method agreements with profile so far...")
  frac_agreement_clonal = calc_method_agreement(all_data_clonal, segments, consensus_profile, "clonal")
  clonal_ranking = sort(unlist(frac_agreement_clonal), decreasing=T)
  if (all(frac_agreement_clonal==0)) {
    # If there is no agreement, i.e. no best method, make sure the NOT overruled methods are ranked highest
    clonal_ranking = c(colnames(method_overruled)[method_overruled==FALSE], colnames(method_overruled)[method_overruled==TRUE])
    names(clonal_ranking) = clonal_ranking
  }
  
  frac_agreement_rounded = calc_method_agreement(all_data_rounded, segments, consensus_profile, "clonal")
  rounded_ranking = sort(unlist(frac_agreement_rounded), decreasing=T)
  if (all(frac_agreement_clonal==0)) {
    # If there is no agreement, i.e. no best method, make sure the NOT overruled methods are ranked highest
    rounded_ranking = c(colnames(method_overruled)[method_overruled==FALSE], colnames(method_overruled)[method_overruled==TRUE])
    names(rounded_ranking) = rounded_ranking
  }
  
  #####################################################################
  # Create the consensus using the method that is most often agreeing with the consensus so far
  #####################################################################
  update_consensus_profile = function(consensus_profile, rounded_ranking, all_data_clonal, segments, identified_preferences, identified_preferences_maj_vote, method_overruled, allowed_methods_x_female, allowed_methods_x_male) {
    # identified_preferences contains alleles that have been identified during rounding, but the other allele count not be found
  
    #' Function that looks for whether a method is allowed to make a call. Some methods do not
    #' participate on X and Y, while single methods are restricted from adding large homozygous
    #' deletions
    check_if_allowed = function(segments, closest_method, method_name, i, allowed_methods_x_female, allowed_methods_x_male) {
      # Don't allow not to be used methods for X and Y segments
      if (segments$chromosome[i] %in% c("X", "Y") & sex=="female" & !method_name %in% allowed_methods_x_female) {
        print(paste0("Excluding ", method_name, " from adding call on, chrom: ", segments$chromosome[i]))
        return(FALSE)
      }
      
      if (segments$chromosome[i] %in% c("X", "Y") & sex=="male" & !method_name %in% allowed_methods_x_male) {
        print(paste0("Excluding ", method_name, " from adding call on, chrom: ", segments$chromosome[i]))
        return(FALSE)
      }
      
      if (!is.null(closest_method[[i]]) && !is.na(closest_method[[i]]) && nrow(closest_method[[i]][[1]])>0 && !is.na(closest_method[[i]][[1]]$minor_cn) && !is.na(closest_method[[i]][[1]]$major_cn)) {
        
        # Exclude a single method from adding a large hom del
        maj = closest_method[[i]][[1]]$major_cn[1]
        min = closest_method[[i]][[1]]$minor_cn[1]
        
        if (maj==0 & min==0 & (segments$end[i]/1000000-segments$start[i]/1000000) > 20) {
          print(paste0("Excluding ", method_name, " from adding in a large hom del, chrom: ", segments$chromosome[i], " size: ", (segments$end[i]/1000000-segments$start[i]/1000000), "Mb"))
          return(FALSE)
        }
      }
        
      # No other clause triggered, therefore allowed
      return(TRUE)
    }
    
    method_reports_identified_allele = function(identified_preferences, closest_method_profile, i) {
      identified_allele = identified_preferences[i]
      if (is.na(identified_allele) | is.na(closest_method_profile[[i]])) {
        return(FALSE)
      } else {
        contains_identified_allele = any(closest_method_profile[[i]][[1]]$major_cn==identified_allele | closest_method_profile[[i]][[1]]$minor_cn==identified_allele, na.rm=T)
        return(contains_identified_allele)
      }
    }
    
    get_maj_min_call = function(contains_identified_allele, contains_identified_allele_maj_vote, identified_preferences, identified_preferences_maj_vote, closest_method_profile, i, method_name) {
      entry = closest_method_profile[[i]][[1]]
      if (contains_identified_allele) {
        # Pick the best solution with bias from the rounded consensus
        index = which(entry$major_cn==identified_preferences[i] | entry$minor_cn==identified_preferences[i])
        index = which.max(entry$ccf[index])
      } else if (contains_identified_allele_maj_vote) {
        # Pick the best solution with bias from the rounded consensus - taken by majority vote
        index = which(entry$major_cn==identified_preferences_maj_vote[i] | entry$minor_cn==identified_preferences_maj_vote[i])
        index = which.max(entry$ccf[index])
      } else {
        # Pick the highest CCF call
        index = which.max(entry$ccf)
      }
      # Round the CN states - in case of DKFZ selected
      maj = round(entry$major_cn[index])
      min = round(entry$minor_cn[index])
      return(list(maj=maj, min=min))
    }
    
    closest_method = names(rounded_ranking)[1]
    # Take the best method and exclude baflogr entries here
    closest_method_index = which(grepl(paste0("map_", closest_method), names(all_data_clonal))  & !grepl("baflogr", names(all_data_clonal)))
    closest_method_profile = all_data_clonal[[closest_method_index]]$cn_states
    for (i in 1:nrow(consensus_profile)) {

      if (is.na(consensus_profile$major_cn[i]) && is.na(consensus_profile$minor_cn[i])) {
        
        # Check if the method is allowed to make a call singlehandidly on this segment
        is_allowed = check_if_allowed(segments, closest_method_profile, closest_method, i, allowed_methods_x_female, allowed_methods_x_male)
        
        # Check if a previously identified allele is reported by this method
        contains_identified_allele = method_reports_identified_allele(identified_preferences, closest_method_profile, i)
        contains_identified_allele_maj_vote = method_reports_identified_allele(identified_preferences_maj_vote, closest_method_profile, i)
        
        # Make flags for whether this segment has preference and whether there is a previously identified allele
        has_preference = !is.na(identified_preferences[i]) | !is.na(identified_preferences_maj_vote[i])
        not_has_preference_or_contains_identified = ((contains_identified_allele | contains_identified_allele_maj_vote) | !has_preference)
        
        if (has_preference & !contains_identified_allele) {
          print(paste0(i, " : ", closest_method, " not right state ", identified_preferences[i], " or ", identified_preferences_maj_vote[i]))
        }
        
        # Take the fit of the method that is most often agreeing with the consensus
        if (is_allowed & not_has_preference_or_contains_identified & !is.null(closest_method_profile[[i]]) && !is.na(closest_method_profile[[i]]) && nrow(closest_method_profile[[i]][[1]])>0 && !is.na(closest_method_profile[[i]][[1]]$minor_cn) && !is.na(closest_method_profile[[i]][[1]]$major_cn)) {
          res = get_maj_min_call(contains_identified_allele, contains_identified_allele_maj_vote, identified_preferences, identified_preferences_maj_vote, closest_method_profile, i, closest_method)
          consensus_profile$major_cn[i] = res$maj
          consensus_profile$minor_cn[i] = res$min
          
        } else {
          # The most often agreeing method does not make a call, iterate over the others until we find one
          for (j in which(grepl("map_", names(all_data_clonal)) & !grepl("baflogr", names(all_data_clonal)))) {
            if (!is.na(all_data_clonal[[j]])) {
              other_closest_method = all_data_clonal[[j]]$cn_states
              method_name = unlist(stringr::str_split(names(all_data_clonal)[j], "_"))[2]  
              if (method_overruled[[method_name]]) { next }
                
              # Check if the method is allowed to make a call singlehandidly on this segment
              is_allowed = check_if_allowed(segments, other_closest_method, method_name, i, allowed_methods_x_female, allowed_methods_x_male)
              # Check if a previously identified allele is reported by this method
              contains_identified_allele = method_reports_identified_allele(identified_preferences, other_closest_method, i)
              contains_identified_allele_maj_vote = method_reports_identified_allele(identified_preferences_maj_vote, other_closest_method, i)
              
              # Make flags for whether this segment has preference and whether there is a previously identified allele
              has_preference = !is.na(identified_preferences[i]) | !is.na(identified_preferences_maj_vote[i])
              # not_has_preference_or_contains_identified = ((contains_identified_allele | contains_identified_allele_maj_vote) & has_preference)
              
              if (!is_allowed | ((has_preference & !contains_identified_allele) & (has_preference & !contains_identified_allele_maj_vote))) { next }
              
              if (!is.null(other_closest_method[[i]]) && !is.na(other_closest_method[[i]]) && nrow(other_closest_method[[i]][[1]])>0 && !is.na(other_closest_method[[i]][[1]]$minor_cn) && !is.na(other_closest_method[[i]][[1]]$major_cn)) {
                print(paste0("Segment ", i))
                print(paste0("Filling in with ", method_name))
                res = get_maj_min_call(contains_identified_allele, contains_identified_allele_maj_vote, identified_preferences, identified_preferences_maj_vote, other_closest_method, i, closest_method)
                print(res)
                consensus_profile$major_cn[i] = res$maj
                consensus_profile$minor_cn[i] = res$min
                consensus_profile$star[i] = 1
                if (((contains_identified_allele | contains_identified_allele_maj_vote) & has_preference)) {
                  # At least other methods agree on one of the alleles
                  consensus_profile$level[i] = "g"
                } else {
                  # No agreement at all
                  consensus_profile$level[i] = "h"
                }
                consensus_profile$methods_agree[i] = 1
                break # stop the loop as we've found a match
              }
            }
          }
        }
      }
    }
    return(consensus_profile)
  }
  print("Filling in remaining segments with best method...")
  consensus_profile = update_consensus_profile(consensus_profile, 
                                               rounded_ranking, 
                                               all_data_clonal, 
                                               segments, 
                                               agreement_rounded$one_allele_saved, 
                                               agreement_rounded_majority_vote$one_allele_saved,
                                               method_overruled,
                                               allowed_methods_x_female, 
                                               allowed_methods_x_male)
  
  #####################################################################
  # Mark ploidy unknown samples
  #####################################################################
  # If this sample is marked as unknown, we overrule all the segment levels and stars and assign level g to make it stand out
  marked_unknown = readr::read_tsv(marked_unknown_file, col_names=F)
  if (samplename %in% marked_unknown[[1]]) {
    consensus_profile$star = 1
    consensus_profile$level = "i"
  }
  
  #####################################################################
  # Cleanup and write profile
  #####################################################################
  # Pad empty entries if there are no calls for the last segment(s). This can occur for the Y chromosome
  if (nrow(consensus_profile) < nrow(segments)) {
    print(paste0("Diff before: ", nrow(consensus_profile), " ", nrow(segments)))
    template_entry = consensus_profile[1,,drop=F]
    template_entry[1, 1:ncol(template_entry)] = NA
    for (i in nrow(consensus_profile):nrow(segments)) {
      print("Adding template entry")
      consensus_profile = rbind(consensus_profile, template_entry)
    }
  }
  print(paste0("Diff after: ", nrow(consensus_profile), " ", nrow(segments)))
  
  consensus_profile = data.frame(segments, consensus_profile)
  
  is_major_negative = consensus_profile$major_cn < 0
  is_minor_negative = consensus_profile$minor_cn < 0
  is_major_negative[is.na(is_major_negative)] = F
  is_minor_negative[is.na(is_minor_negative)] = F
  if (any(is_major_negative | is_minor_negative)) {
    print(paste0(samplename, " contains segments with CN state lower than 1"))
  }
  
  write.table(consensus_profile, file=file.path(outdir, "consensus_profile", paste0(samplename, "_consensus_profile.txt")), quote=F, sep="\t", row.names=F)
  
  # TODO bugfix : Is creating the full data.frame here needed? it's done above already?
  profile_bb = collapseRoundedClonal2bb(data.frame(segments, consensus_profile))
  p = plot_profile(profile_bb, "Consensus - after rounding and using best method", max.plot.cn=max.plot.cn)
  png(file.path(outdir, "figures", paste0(samplename, "_consensus_rounded_bestMethod.png")), height=300, width=1300)
  print(p)
  dev.off()
  
  #####################################################################
  # Write some summary stats
  #####################################################################
  
  save.image(file.path(outdir, "saves", paste0(samplename, "_agreement_save.RData")))
  
  print("Building summary stats file...")
  res = parse_all_purities(samplename, method_purityfile)
  purity_dkfz = res$dkfz
  purity_vanloowedge = res$vanloowedge
  purity_broad = res$broad
  purity_peifer = res$peifer
  purity_mustonen = res$mustonen
  purity_jabba = res$jabba
  purities = data.frame(purity_dkfz=purity_dkfz, purity_vanloowedge=purity_vanloowedge, purity_peifer=purity_peifer, purity_mustonen=purity_mustonen, purity_broad=purity_broad, purity_jabba=purity_jabba)
  
  purities_tested = test_purities(purities, consensus_profile)
  
  ploidy_vanloowedge = get_ploidy(segments, all_data_clonal$map_vanloowedge)
  ploidy_broad = get_ploidy(segments, all_data_clonal$map_broad, broad=T)
  ploidy_peifer = get_ploidy(segments, all_data_clonal$map_peifer)
  ploidy_dkfz = get_ploidy(segments, all_data_clonal$map_dkfz)
  ploidy_mustonen = get_ploidy(segments, all_data_clonal$map_mustonen)
  ploidy_jabba = get_ploidy(segments, all_data_clonal$map_jabba)
  ploidy_consensus = round(calc_ploidy(profile_bb), 4)
  
  ploidies = data.frame(ploidy_vanloowedge=ploidy_vanloowedge$ploidy, ploidy_broad=ploidy_broad$ploidy, ploidy_peifer=ploidy_peifer$ploidy, ploidy_dkfz=ploidy_dkfz$ploidy, ploidy_mustonen=ploidy_mustonen$ploidy, ploidy_jabba=ploidy_jabba$ploidy, ploidy_consensus=ploidy_consensus,
                        status_vanloowedge=ploidy_vanloowedge$status, status_broad=ploidy_broad$status, status_peifer=ploidy_peifer$status, status_dkfz=ploidy_dkfz$status, status_mustonen=ploidy_mustonen$status, status_jabba=ploidy_jabba$status)
  
  
  agreement_summary = as.data.frame(t(data.frame(
    c(
      round(unlist(frac_agreement_clonal),4), 
      unlist(names(clonal_ranking)),
      round(unlist(frac_agreement_rounded),4),
      unlist(names(rounded_ranking)))
    )))
  
  
  colnames(agreement_summary) = c(paste0("agree_clonal_", names(frac_agreement_clonal)), 
                                   paste0("rank_clonal_", 1:length(names(clonal_ranking))),
                                   paste0("agree_rounded_", names(frac_agreement_rounded)), 
                                   paste0("rank_rounded_", 1:length(names(rounded_ranking))))
  
  # See if there are clonal aberrations
  has_clonal_aberration = any(consensus_profile$major_cn!=1 | consensus_profile$minor_cn!=1)
  
  seg_lengths = round(consensus_profile$end/1000000 - consensus_profile$start/1000000)
  has_largish_clonal_aberration = any((consensus_profile$major_cn!=1 | consensus_profile$minor_cn!=1)[seg_lengths > 10])
  has_large_clonal_aberration = any((consensus_profile$major_cn!=1 | consensus_profile$minor_cn!=1)[seg_lengths > 50])
  
  if (is.na(method_overruled)) {
    overruling = data.frame(exclude_broad=FALSE,
                            exclude_mustonen=FALSE,
                            exclude_dkfz=FALSE,
                            exclude_peifer=FALSE,
                            exclude_vanloowedge=FALSE,
                            exclude_jabba=FALSE)
  } else {
    overruling = data.frame(exclude_broad=method_overruled$broad,
                            exclude_mustonen=method_overruled$mustonen,
                            exclude_dkfz=method_overruled$dkfz,
                            exclude_peifer=method_overruled$peifer,
                            exclude_vanloowedge=method_overruled$vanloowedge,
                            exclude_jabba=method_overruled$jabba)
  }
  
  summary_data = data.frame(samplename=samplename, 
                            overruling,
                            agreement_level_a=agreement_clonal$frac_agree[,2], 
                            agreement_level_b=agreement_clonal_overrule$frac_agree[,2], 
                            agreement_level_c=agreement_rounded$frac_agree[,2], 
                            agreement_level_d=agreement_rounded_majority_vote$frac_agree[,2], 
                            has_clonal_aberration=has_clonal_aberration, 
                            has_largish_clonal_aberration=has_largish_clonal_aberration, 
                            has_large_clonal_aberration=has_large_clonal_aberration,
                            agreement_summary,
                            purities,
                            ploidies,
                            purities_tested)
  write.table(summary_data, file=file.path(outdir, "summary_stats", paste0(samplename, "_summary_stats.txt")), quote=F, sep="\t", row.names=F)
} else {
  print("No breakpoints file found")
}

q(save="no")
