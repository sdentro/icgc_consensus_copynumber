source("~/repo/icgc_consensus_copynumber/util.R")


#####################################################################
# Rounded CN main method
#####################################################################
create_rounded_copynumber = function(samplename, segments, outdir, method_segmentsfile, method_purityfile, max.plot.cn=4) {
  
  res = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile)
  map_dkfz = res$map_dkfz
  map_mustonen = res$map_mustonen
  map_peifer = res$map_peifer
  map_vanloowedge = res$map_vanloowedge
  map_broad = res$map_broad
  
  combined_status = get_combined_status(segments, map_vanloowedge, map_dkfz, map_mustonen, map_peifer, map_broad)
  
  #####################################################################
  # Make inventory
  #####################################################################
  
  all_clonal = apply(combined_status[,c("dkfz", "peifer", "vanloowedge")], 1, function(x) { all(x=="clonal") })
  all_clonal[is.na(all_clonal)] = F
  all_subclonal = apply(combined_status[,c("dkfz", "peifer", "vanloowedge")], 1, function(x) { all(x=="subclonal") })
  all_subclonal[is.na(all_subclonal)] = F
  any_na = apply(combined_status[,c("dkfz", "peifer", "vanloowedge")], 1, function(x) { any(is.na(x)) })
  one_subclonal = apply(combined_status[,c("dkfz", "peifer", "vanloowedge")], 1, function(x) { sum(x=="subclonal")==1 })
  one_subclonal[is.na(one_subclonal)] = F
  two_subclonal = apply(combined_status[,c("dkfz", "peifer", "vanloowedge")], 1, function(x) { sum(x=="subclonal")==2 })
  two_subclonal[is.na(two_subclonal)] = F
  three_subclonal = apply(combined_status[,c("dkfz", "peifer", "vanloowedge")], 1, function(x) { sum(x=="subclonal")==3 })
  three_subclonal[is.na(three_subclonal)] = F
  
  # sanity check
  if (any((all_clonal+all_subclonal+any_na+one_subclonal+two_subclonal+three_subclonal) > 1)) {
    print("Found segment with multiple status?")
  }
  
  status_inventory = rep(NA, nrow(combined_status))
  status_inventory[all_clonal] = "all_clonal"
  status_inventory[all_subclonal] = "all_subclonal"
  status_inventory[any_na] = "any_na"
  status_inventory[one_subclonal] = "one_subclonal"
  status_inventory[two_subclonal] = "two_subclonal"
  status_inventory[three_subclonal] = "three_subclonal"
  
  combined_status$inventory = status_inventory
  write.table(combined_status, file=file.path(outdir, "status_inventory", paste0(samplename, "_status_inventory.txt")), row.names=F, sep="\t", quote=F)
  
  #####################################################################
  # Make plots
  #####################################################################
  print("Making figures")
  print("Battenberg")
  if (!is.na(map_vanloowedge)) {
    cn_bb_vanloowedge = collapse2bb(segments=segments, cn_states=map_vanloowedge$cn_states)
  } else {
    cn_bb_vanloowedge = parse_dummy_cn_profile()
  }
  plot_vanloowedge = plot_profile(cn_bb_vanloowedge, "Battenberg", max.plot.cn=max.plot.cn)
  
  print("ABSOLUTE")
  if (!is.na(map_broad)) {
    cn_bb_broad = collapse2bb(segments=segments, cn_states=map_broad$cn_states, broad=T)
  } else {
    cn_bb_broad = parse_dummy_cn_profile()
  }
  plot_broad = plot_profile(cn_bb_broad, "ABSOLUTE", max.plot.cn=max.plot.cn)
  
  print("ACEseq")
  if (!is.na(map_dkfz)) {
    cn_bb_dkfz = collapse2bb(segments=segments, cn_states=map_dkfz$cn_states)
  } else {
    cn_bb_dkfz = parse_dummy_cn_profile()
  }
  plot_dkfz = plot_profile(cn_bb_dkfz, "ACEseq", max.plot.cn=max.plot.cn)
  
  print("CloneHD")
  if (!is.na(map_mustonen)) {
    cn_bb_mustonen = collapse2bb(segments=segments, cn_states=map_mustonen$cn_states)
  } else {
    cn_bb_mustonen = parse_dummy_cn_profile()
  }
  plot_mustonen = plot_profile(cn_bb_mustonen, "CloneHD", max.plot.cn=max.plot.cn)
  
  print("Sclust")
  if (!is.na(map_peifer)) {
    cn_bb_peifer = collapse2bb(segments=segments, cn_states=map_peifer$cn_states)
  } else {
    cn_bb_peifer = parse_dummy_cn_profile()
  }
  plot_peifer = plot_profile(cn_bb_peifer, "Sclust", max.plot.cn=max.plot.cn)
  print("Done plot profile")
  
  png(file.path(outdir, "figures", paste0(samplename, "_copynumber_complete.png")), height=1500, width=1300)
  grid.arrange(arrangeGrob(plot_broad, plot_dkfz, plot_vanloowedge, plot_mustonen, plot_peifer, ncol=1), 
               top=textGrob(samplename, gp = gpar(fontsize=25, face=2, col="black")))
  dev.off()
  
  #####################################################################
  # Make rounded clonal copy number
  #####################################################################
  print("Looking at segments")
  rounded_clonal = list(broad=data.frame(), dkfz=data.frame(), vanloowedge=data.frame(), peifer=data.frame(), mustonen=data.frame())
  for (i in 1:nrow(segments)) {
    # TODO: at the moment this does not reset the CP/CCFs

    if (!is.na(map_broad) && length(map_broad$cn_states) >= i) {
      rounded_clonal$broad = rbind(rounded_clonal$broad, round_broad(map_broad, i))
    } else {
      temp_entry = get_dummy_cn_entry(segments[i,,drop=F])
      temp_entry$historically_clonal = 0
      rounded_clonal$broad = rbind(rounded_clonal$broad, temp_entry)
    }
    
    if (!is.na(map_vanloowedge) && length(map_vanloowedge$cn_states) >= i) {
      rounded_clonal$vanloowedge = rbind(rounded_clonal$vanloowedge, round_vanloo_wedge(map_vanloowedge, i))
    } else {
      temp_entry = get_dummy_cn_entry(segments[i,,drop=F])
      colnames(temp_entry)[7] = "cellular_prevalence"
      temp_entry$ccf = NA
      rounded_clonal$vanloowedge = rbind(rounded_clonal$vanloowedge, temp_entry)
    }

    if (!is.na(map_peifer) && length(map_peifer$cn_states) >= i) {
      rounded_clonal$peifer = rbind(rounded_clonal$peifer, round_peifer(map_peifer, i))
    } else {
      temp_entry = get_dummy_cn_entry(segments[i,,drop=F])
      colnames(temp_entry)[7] = "cellular_prevalence"
      temp_entry$ccf = NA
      rounded_clonal$peifer = rbind(rounded_clonal$peifer, temp_entry)
    }

    if (!is.na(map_mustonen) && length(map_mustonen$cn_states) >= i) {
      rounded_clonal$mustonen = rbind(rounded_clonal$mustonen, round_mustonen(map_mustonen, i))
    } else {
      temp_entry = get_dummy_cn_entry(segments[i,,drop=F])
      colnames(temp_entry)[7] = "cellular_prevalence"
      temp_entry$ccf = NA
      rounded_clonal$mustonen = rbind(rounded_clonal$mustonen, temp_entry)
    }

    if (!is.na(map_dkfz) && length(map_dkfz$cn_states) >= i) {
      rounded_clonal$dkfz = rbind(rounded_clonal$dkfz, round_dkfz(map_dkfz, i))
    } else {
      temp_entry = get_dummy_cn_entry(segments[i,,drop=F])
      colnames(temp_entry)[7] = "cellular_prevalence"
      temp_entry$dh = 0
      temp_entry$covRatio = 0
      temp_entry$ccf = NA
      rounded_clonal$dkfz = rbind(rounded_clonal$dkfz, temp_entry)
    }
  }
  print("Calling unique")
  for (i in 1:length(rounded_clonal)) {
    rounded_clonal[[i]] = unique(rounded_clonal[[i]])
  }
  
  print("Writing data to file")
  # write out the complete profiles
  for (i in 1:length(rounded_clonal)) {
    cn_method = names(rounded_clonal)[i]
    if (all(rounded_clonal[[i]]$copy_number==-2)) {
      next
    } else {
      write.table(rounded_clonal[[i]], file=file.path(outdir, paste0(cn_method, "_rounded_clonal"), paste0(samplename, "_segments.txt")), quote=F, sep="\t", row.names=F)
    }
  }
  
  #####################################################################
  # Plot rounded clonal
  #####################################################################
  print("Making figures - rounded")
  print("Battenberg")
  cn_bb_vanloowedge = collapseRoundedClonal2bb(cn_states=rounded_clonal$vanloowedge)
  plot_vanloowedge = plot_profile(cn_bb_vanloowedge, "Battenberg", max.plot.cn=max.plot.cn)
  
  print("ABSOLUTE")
  cn_bb_broad = collapseRoundedClonal2bb(cn_states=rounded_clonal$broad)
  plot_broad = plot_profile(cn_bb_broad, "ABSOLUTE", max.plot.cn=max.plot.cn)
  
  print("ACEseq")
  cn_bb_dkfz = collapseRoundedClonal2bb(cn_states=rounded_clonal$dkfz)
  plot_dkfz = plot_profile(cn_bb_dkfz, "ACEseq", max.plot.cn=max.plot.cn)
  
  print("CloneHD")
  cn_bb_mustonen = collapseRoundedClonal2bb(cn_states=rounded_clonal$mustonen)
  plot_mustonen = plot_profile(cn_bb_mustonen, "CloneHD", max.plot.cn=max.plot.cn)
  
  print("Sclust")
  cn_bb_peifer = collapseRoundedClonal2bb(cn_states=rounded_clonal$peifer)
  plot_peifer = plot_profile(cn_bb_peifer, "Sclust", max.plot.cn=max.plot.cn)
  
  png(file.path(outdir, "figures", paste0(samplename, "_copynumber_rounded_clonal.png")), height=1500, width=1300)
  grid.arrange(arrangeGrob(plot_broad, plot_dkfz, plot_vanloowedge, plot_mustonen, plot_peifer, ncol=1), 
               top=textGrob(samplename, gp = gpar(fontsize=25, face=2, col="black")))
  dev.off()
}


#####################################################################
# Script
#####################################################################
#' TODO
#'  - CP/CCF is not reset
#'  - No purity values for Broad included, maybe not needed
#'  - File paths at the top of the create_rounded_copynumber could be improved, maybe move all the input into an input directory

args = commandArgs(T)
samplename = args[1]
outdir = args[2]

# setwd("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/consensus_copynumber/2016_09_consensus_breakpoints_fullruns/20161024_rounding/")
# samplename = "6aa00162-6294-4ce7-b6b7-0c3452e24cd6"
# outdir = "output"


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


# outdir = "output"
# samplename = "6aa00162-6294-4ce7-b6b7-0c3452e24cd6"
breakpoints_file = file.path("consensus_bp", paste0(samplename, ".txt"))
if (file.exists(breakpoints_file)) {
  breakpoints = read.table(breakpoints_file, header=T, stringsAsFactors=F)
  segments = breakpoints2segments(breakpoints)
  create_rounded_copynumber(samplename, segments, outdir, method_segmentsfile, method_purityfile, max.plot.cn=4)
}

q(save="no")







