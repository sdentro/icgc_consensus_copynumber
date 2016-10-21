

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

get_frac_genome_agree = function(samplename, method_segmentsfile, method_purityfile, mustonen_has_header) {
  breakpoints = read.table(paste0(samplename, "_consensus_breakpoints.txt"), header=T, stringsAsFactors=F)
  segments = breakpoints2segments(breakpoints)

  res = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=mustonen_has_header)
  map_dkfz = res$map_dkfz
  map_mustonen = res$map_mustonen
  map_peifer = res$map_peifer
  map_vanloowedge = res$map_vanloowedge
  map_broad = res$map_broad
  all_maps = list(map_broad=map_broad, map_dkfz=map_dkfz, map_mustonen=map_mustonen, map_vanloowedge=map_vanloowedge, map_peifer=map_peifer)
  combined_status = data.frame(segments, dkfz=map_dkfz$status, mustonen=map_mustonen$status, peifer=map_peifer$status, vanloowedge=map_vanloowedge$status, broad=map_broad$status)
  
  agree = rep(F, nrow(segments))
  cn_states = list()
  num_methods = rep(0, nrow(segments))
  for (i in 1:nrow(segments)) {
    
    methods_with_result = (5:8)[!is.na(combined_status[i,5:8])]
    if (length(methods_with_result) > 0 && all(combined_status[i,methods_with_result]=="clonal")) {
      
      inventory = data.frame()
      for (j in 1:length(all_maps)) {
        map = all_maps[[j]]
        
        if (!is.null(map$cn_states[[i]])) {
          seg = map$cn_states[[i]][[1]]
          inventory = rbind(inventory, data.frame(method=names(all_maps)[j], major_cn=seg$major_cn, minor_cn=seg$minor_cn))
        }
      }
      cn_states[[i]] = inventory
      inventory = na.omit(inventory)
      agree[i] = all(inventory$major_cn==inventory$major_cn[1]) & all(inventory$minor_cn==inventory$minor_cn[1])
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
dkfz_segmentsfile = paste0(samplename, "_dkfz/", samplename, "_segments.txt")
dkfz_purityfile = paste0(samplename, "_dkfz/puritiesPloidies_20160906.txt")
vanloowedge_segmentsfile = paste0(samplename, "_vanloowedge/", samplename, "_segments.txt.gz")
vanloowedge_purityfile = paste0(samplename, "_vanloowedge/", samplename, "_purity_ploidy.txt")
peifer_segmentsfile = paste0(samplename, "_peifer/", samplename, "_segments.txt")
peifer_purityfile = paste0(samplename, "_peifer/", samplename, "_purity_ploidy.txt")
mustonen_segmentsfile = paste0(samplename, "_mustonen/", samplename, ".penalty0.95_segments.txt")
mustonen_purityfile = paste0(samplename, "_mustonen/", samplename, ".penalty0.95.purity.ploidy.txt")
broad_segmentsfile = paste0(samplename, "_broad/", samplename, "_segments.txt")
broad_purityfile = paste0(samplename, "_broad/", samplename, "_purity_ploidy.txt")

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

outdir = "output"
samplename = "6aa00162-6294-4ce7-b6b7-0c3452e24cd6"
agreement_clonal = get_frac_genome_agree(samplename, method_segmentsfile, method_purityfile, mustonen_has_header=F)

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

agreement_rounded = get_frac_genome_agree(samplename, method_segmentsfile, method_purityfile, mustonen_has_header=T)

#####################################################################
# Agreement on median cn states
#####################################################################

r = agreement_rounded$cn_states
cn_states = list()
num_methods = rep(0, length(r))
agree = rep(F, length(r))
for (i in 1:length(r)) {
  
  if (!agreement_rounded$agree[i]) {
    major = median(r[[i]]$major_cn, na.rm=T)
    minor = median(r[[i]]$minor_cn, na.rm=T)
    num_methods[i] = sum(!is.na(r[[i]]$major_cn))
    agree[i] = T
  } else {
    # There are no na's in this data.frame and there is agreement on the cn state, so its safe to take the first entry here
    major = r[[i]]$major_cn[1]
    minor = r[[i]]$minor_cn[1]
    num_methods[i] = 1
    agree[i] = T
  }
  cn_states[[i]] = data.frame(method="median", major_cn=major, minor_cn=minor)
}
frac_agree = NA
agreement_median = list(frac_agree=frac_agree, segments=r$segments, agree=agree, cn_states=cn_states)

#####################################################################
# Piece together a complete agreement profile
#####################################################################
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
    new_entry = agreement_median$cn_states[[i]][1,2:3]
    new_entry$star = 3
  }
  
  consensus_profile = rbind(consensus_profile, new_entry)
}
write.table(consensus_profile, file=paste0(samplename, "_"), quote=F, sep="\t", row.names=F)
