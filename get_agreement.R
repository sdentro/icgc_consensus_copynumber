

#' 1) Load in a profile and map the segments to the reference
#' 2) across all segments see which we agree on
#' 3) calculate the fraction of the genome on which we agree
#' 

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

get_frac_genome_agree = function(samplename, all_data, segments) {
  # breakpoints = read.table(paste0(samplename, "_consensus_breakpoints.txt"), header=T, stringsAsFactors=F)
  # segments = breakpoints2segments(breakpoints)

  # res = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=mustonen_has_header)
  map_dkfz = all_data$map_dkfz
  map_mustonen = all_data$map_mustonen
  map_peifer = all_data$map_peifer
  map_vanloowedge = all_data$map_vanloowedge
  map_broad = all_data$map_broad
  all_maps = list(map_broad=map_broad, map_dkfz=map_dkfz, map_mustonen=map_mustonen, map_vanloowedge=map_vanloowedge, map_peifer=map_peifer)
  # combined_status = data.frame(segments, dkfz=map_dkfz$status, mustonen=map_mustonen$status, peifer=map_peifer$status, vanloowedge=map_vanloowedge$status, broad=map_broad$status)
  combined_status = get_combined_status(segments, map_vanloowedge, map_dkfz, map_mustonen, map_peifer, map_broad)
  
  
  agree = rep(F, nrow(segments))
  cn_states = list()
  num_methods = rep(0, nrow(segments))
  for (i in 1:nrow(segments)) {
    
    methods_with_result = (5:8)[!is.na(combined_status[i,5:8])]
    if (length(methods_with_result) > 0 && all(combined_status[i,methods_with_result]=="clonal")) {
      
      inventory = data.frame()
      for (j in 1:length(all_maps)) {
        map = all_maps[[j]]
        
        if (!is.na(map) && !is.null(map$cn_states[[i]])) {
          seg = map$cn_states[[i]][[1]]
          inventory = rbind(inventory, data.frame(method=names(all_maps)[j], major_cn=seg$major_cn, minor_cn=seg$minor_cn))
        }
      }
      cn_states[[i]] = inventory
      inventory = na.omit(inventory)
      agree[i] = all(inventory$major_cn==inventory$major_cn[1]) & all(inventory$minor_cn==inventory$minor_cn[1]) & nrow(inventory)>0
      num_methods[i] = nrow(inventory)
    }
  }
  
  segments$size = segments$end - segments$start
  frac_genome_agree = round(sum(segments$size[agree]/1000) / sum(segments$size/1000), 2)
  
  return(list(frac_agree=data.frame(samplename=samplename, frac_genome_agree=frac_genome_agree), segments=segments, agree=agree, cn_states=cn_states, num_methods=num_methods))
}

#####################################################################
# Original agreement
#####################################################################
# setwd("~/Documents/Projects/icgc/consensus_subclonal_copynumber/")

source("~/repo/icgc_consensus_copynumber/util.R")
max.plot.cn=4


args = commandArgs(T)
samplename = args[1]
outdir = args[2]

# outdir = "output"
# samplename = "6aa00162-6294-4ce7-b6b7-0c3452e24cd6"

breakpoints = read.table(file.path("consensus_bp", paste0(samplename, ".txt")), header=T, stringsAsFactors=F)
segments = breakpoints2segments(breakpoints)

dkfz_segmentsfile = paste0("dkfz/", samplename, "_segments.txt")
dkfz_purityfile = "purity_ploidy_dkfz.txt"
vanloowedge_segmentsfile = paste0("vanloowedge/", samplename, "_segments.txt")
vanloowedge_purityfile = "purity_ploidy_vanloowedge.txt"
peifer_segmentsfile = paste0("peifer/", samplename, "_segments.txt")
peifer_purityfile = "purity_ploidy_peifer.txt"
mustonen_segmentsfile = paste0("mustonen/", samplename, ".penalty0.95_segments.txt")
mustonen_purityfile = "purity_ploidy_mustonen.txt"
broad_segmentsfile = paste0("broad/", samplename, "_segments.txt")
broad_purityfile = NA

method_segmentsfile = list(dkfz=dkfz_segmentsfile,
                           vanloowedge=vanloowedge_segmentsfile,
                           peifer=peifer_segmentsfile,
                           mustonen=mustonen_segmentsfile,
                           broad=broad_segmentsfile)

method_purityfile = list(dkfz=dkfz_purityfile,
                         vanloowedge=vanloowedge_purityfile,
                         peifer=peifer_purityfile,
                         mustonen=mustonen_purityfile,
                         broad=broad_purityfile)

all_data_clonal = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=F)
agreement_clonal = get_frac_genome_agree(samplename, all_data_clonal, segments)

#####################################################################
# Agreement after rounding
#####################################################################
dkfz_segmentsfile = file.path(outdir, "dkfz_rounded_clonal", paste0(samplename, "_segments.txt"))
vanloowedge_segmentsfile = file.path(outdir, "vanloowedge_rounded_clonal", paste0(samplename, "_segments.txt"))
peifer_segmentsfile = file.path(outdir, "peifer_rounded_clonal", paste0(samplename, "_segments.txt"))
mustonen_segmentsfile = file.path(outdir, "mustonen_rounded_clonal", paste0(samplename, "_segments.txt"))
broad_segmentsfile = file.path(outdir, "broad_rounded_clonal", paste0(samplename, "_segments.txt"))

method_segmentsfile = list(dkfz=dkfz_segmentsfile,
                           vanloowedge=vanloowedge_segmentsfile,
                           peifer=peifer_segmentsfile,
                           mustonen=mustonen_segmentsfile,
                           broad=broad_segmentsfile)

all_data_rounded = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=T)

agreement_rounded = get_frac_genome_agree(samplename, all_data_rounded, segments)

#####################################################################
# Agreement on median cn states
#####################################################################
# No longer used as this has unwanted side affects when methods are disagreeing and can create non-integer solutions
# 
# r = agreement_rounded$cn_states
# cn_states = list()
# num_methods = rep(0, length(r))
# agree = rep(F, length(r))
# for (i in 1:length(r)) {
#   
#   if (!agreement_rounded$agree[i]) {
#     major = median(r[[i]]$major_cn, na.rm=T)
#     minor = median(r[[i]]$minor_cn, na.rm=T)
#     num_methods[i] = sum(!is.na(r[[i]]$major_cn))
#     agree[i] = T
#   } else {
#     # There are no na's in this data.frame and there is agreement on the cn state, so its safe to take the first entry here
#     major = r[[i]]$major_cn[1]
#     minor = r[[i]]$minor_cn[1]
#     num_methods[i] = 1
#     agree[i] = T
#   }
#   cn_states[[i]] = data.frame(method="median", major_cn=major, minor_cn=minor)
# }
# frac_agree = NA
# agreement_median = list(frac_agree=frac_agree, segments=r$segments, agree=agree, cn_states=cn_states)


#####################################################################
# Piece together a complete agreement profile
#####################################################################
create_consensus_profile = function(agreement_clonal, agreement_rounded) {
  consensus_profile = data.frame()
  r = agreement_rounded$cn_states
  for (i in 1:length(r)) {
    
    if (agreement_clonal$agree[i]) {
      # if clonal agree, choose that and assign 1*
      new_entry = agreement_clonal$cn_states[[i]][1,2:3]
      new_entry$star = 1
      
      
    } else if (agreement_rounded$agree[i]) {
      # else if rounded clonal agree, choose that and assign 2*
      new_entry = agreement_rounded$cn_states[[i]][1,2:3]
      new_entry$star = 2
        
    } else {
      # else take the median solution and assign 3*
      # new_entry = agreement_median$cn_states[[i]][1,2:3]
      new_entry = data.frame(major_cn=NA, minor_cn=NA)
      new_entry$star = 3
    }
    
    consensus_profile = rbind(consensus_profile, new_entry)
  }
  return(consensus_profile)
}

consensus_profile = create_consensus_profile(agreement_clonal, agreement_rounded)

profile_bb = collapseRoundedClonal2bb(data.frame(segments, consensus_profile))
p = plot_profile(profile_bb, "Consensus - after rounding", max.plot.cn=max.plot.cn)
png(file.path(outdir, "figures", paste0(samplename, "_consensus_rounded.png")), height=300, width=1300)
print(p)
dev.off()

#####################################################################
# Calc agreement with both the created consensus per method
#####################################################################

calc_method_agreement = function(all_maps, segments, consensus_profile, segment_status) {
  agreement = list(dkfz=0, vanloowedge=0, peifer=0, mustonen=0, broad=0)
  
  for (i in 1:nrow(consensus_profile)) {
    # Iterate over all the segment mappings
    for (j in which(grepl("map", names(all_maps)))) {
      method = gsub("map_", "", names(all_maps)[j])
      
      if (!is.na(all_maps[[j]]) && !is.na(all_maps[[j]]$status[i]) && all_maps[[j]]$status[i]==segment_status) {
        
        # Check if agreement
        method_segment = all_maps[[j]]$cn_states[[i]][[1]]
        if (!is.null(method_segment) && !is.na(consensus_profile$major_cn[i]) && (!is.na(method_segment$major_cn) | !is.na(method_segment$minor_cn)) &&
            consensus_profile$major_cn[i]==method_segment$major_cn & consensus_profile$minor_cn[i]==method_segment$minor_cn) {
          # Add the extra agreement to the tally of this method
          agreement[[method]] = agreement[[method]] + (segments$end[i]/1000-segments$start[i]/1000)
        }
      }
    }
  }
  
  genome_size = sum(segments$end/1000 - segments$start/1000)
  frac_agreement = lapply(agreement, function(x) x / genome_size)
  return(frac_agreement)
}

frac_agreement_clonal = calc_method_agreement(all_data_clonal, segments, consensus_profile, "clonal")
clonal_ranking = sort(unlist(frac_agreement_clonal), decreasing=T)

frac_agreement_rounded = calc_method_agreement(all_data_rounded, segments, consensus_profile, "clonal")
rounded_ranking = sort(unlist(frac_agreement_rounded), decreasing=T)


#####################################################################
# Create the consensus using the method that is most often agreeing with the consensus so far
#####################################################################
update_consensus_profile = function(consensus_profile, rounded_ranking, all_data_rounded) {

  closest_method = names(rounded_ranking)[1]
  closest_method_index = which(grepl(paste0("map_", closest_method), names(all_data_rounded)))
  closest_method_profile = all_data_rounded[[closest_method_index]]$cn_states
  for (i in 1:nrow(consensus_profile)) {
    if (is.na(consensus_profile$major_cn[i]) && is.na(consensus_profile$minor_cn[i])) {
      
      # Take the fit of the method that is most often agreeing with the consensus
      if (!is.null(closest_method_profile[[i]])) {
        consensus_profile$major_cn[i] = closest_method_profile[[i]][[1]]$major_cn[1]
        consensus_profile$minor_cn[i] = closest_method_profile[[i]][[1]]$minor_cn[1]
      } else {
        
      # The most often agreeing method does not make a call, iterate over the others until we find one
        for (j in which(grepl("map_", names(all_data_rounded)))) {
          if (!is.na(all_data_rounded[[j]])) {
            other_closest_method = all_data_rounded[[j]]$cn_states
            if (!is.null(other_closest_method[[i]])) {
              consensus_profile$major_cn[i] = other_closest_method[[i]][[1]]$major_cn[1]
              consensus_profile$minor_cn[i] = other_closest_method[[i]][[1]]$minor_cn[1]
              consensus_profile$star[i] = 4
            }
          }
        }
      }
    }
  }
  return(consensus_profile)
}

consensus_profile = update_consensus_profile(consensus_profile, rounded_ranking, all_data_rounded)
consensus_profile = data.frame(segments, consensus_profile)
write.table(consensus_profile, file=file.path(outdir, "consensus_profile", paste0(samplename, "_consensus_profile.txt")), quote=F, sep="\t", row.names=F)

profile_bb = collapseRoundedClonal2bb(data.frame(segments, consensus_profile))
p = plot_profile(profile_bb, "Consensus - after rounding and using best method", max.plot.cn=max.plot.cn)
png(file.path(outdir, "figures", paste0(samplename, "_consensus_rounded_bestMethod.png")), height=300, width=1300)
print(p)
dev.off()

#####################################################################
# Write some summary stats
#####################################################################

get_ploidy = function(segments, map) {
  if (!is.na(map)) {
    cn_bb = collapse2bb(segments=segments, cn_states=map$cn_states)
    return(list(ploidy=calc_ploidy(cn_bb), status=get_ploidy_status(cn_bb)))
  } else {
    return(list(ploidy=NA, status=NA))
  }
}

ploidy_vanloowedge = get_ploidy(segments, all_data_clonal$map_vanloowedge)
ploidy_broad = get_ploidy(segments, all_data_clonal$map_broad)
ploidy_peifer = get_ploidy(segments, all_data_clonal$map_peifer)
ploidy_dkfz = get_ploidy(segments, all_data_clonal$map_dkfz)
ploidy_mustonen = get_ploidy(segments, all_data_clonal$map_mustonen)
ploidies = data.frame(ploidy_vanloowedge=ploidy_vanloowedge$ploidy, ploidy_broad=ploidy_broad$ploidy, ploidy_peifer=ploidy_peifer$ploidy, ploidy_dkfz=ploidy_dkfz$ploidy, ploidy_mustonen=ploidy_mustonen$ploidy,
                      status_vanloowedge=ploidy_vanloowedge$status, status_broad=ploidy_broad$status, status_peifer=ploidy_peifer$status, status_dkfz=ploidy_dkfz$status, status_mustonen=ploidy_mustonen$status)


agreement_summary = as.data.frame(t(data.frame(
  c(
    unlist(frac_agreement_clonal), 
    unlist(names(clonal_ranking)),
    unlist(frac_agreement_rounded),
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

summary_data = data.frame(samplename=samplename, 
                           agreement_clonal=agreement_clonal$frac_agree[,2], 
                           agreement_rounded=agreement_rounded$frac_agree[,2], 
                           has_clonal_aberration=has_clonal_aberration, 
                           has_largish_clonal_aberration=has_largish_clonal_aberration, 
                           has_large_clonal_aberration=has_large_clonal_aberration,
                           agreement_summary,
                           ploidies)
write.table(summary_data, file=file.path(outdir, "summary_stats", paste0(samplename, "_summary_stats.txt")), quote=F, sep="\t", row.names=F)
