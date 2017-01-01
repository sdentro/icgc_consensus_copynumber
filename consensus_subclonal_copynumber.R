source("~/repo/icgc_consensus_copynumber/util.R")


#####################################################################
# Rounded CN main method
#####################################################################
create_rounded_copynumber = function(samplename, segments, outdir, method_segmentsfile, method_purityfile, method_baflogr, sex, max.plot.cn=4, rounding_up=T, make_figures=T) {
  
  #####################################################################
  # Read in and map all the raw data
  #####################################################################
  res = parse_all_profiles(samplename, segments, method_segmentsfile, method_purityfile, method_baflogr, sex)
  map_dkfz = res$map_dkfz
  map_mustonen = res$map_mustonen
  map_peifer = res$map_peifer
  map_vanloowedge = res$map_vanloowedge
  map_broad = res$map_broad
  
  # Calc ploidy of all profiles
  ploidy_vanloowedge = get_ploidy(segments, res$map_vanloowedge)
  ploidy_broad = get_ploidy(segments, res$map_broad, broad=T)
  ploidy_peifer = get_ploidy(segments, res$map_peifer)
  ploidy_dkfz = get_ploidy(segments, res$map_dkfz)
  ploidy_mustonen = get_ploidy(segments, res$map_mustonen)
  
  ploidies = data.frame(samplename=samplename,
                        ploidy_vanloowedge=ploidy_vanloowedge$ploidy, ploidy_broad=ploidy_broad$ploidy, ploidy_peifer=ploidy_peifer$ploidy, ploidy_dkfz=ploidy_dkfz$ploidy, ploidy_mustonen=ploidy_mustonen$ploidy,
                        status_vanloowedge=ploidy_vanloowedge$status, status_broad=ploidy_broad$status, status_peifer=ploidy_peifer$status, status_dkfz=ploidy_dkfz$status, status_mustonen=ploidy_mustonen$status)
  write.table(ploidies, file=file.path(outdir, "raw_ploidy", paste0(samplename, "_raw_ploidy.txt")), sep="\t", quote=F, row.names=F)
  
  res = parse_all_purities(samplename, method_purityfile)
  purity_dkfz = res$dkfz
  purity_vanloowedge = res$vanloowedge
  purity_broad = res$broad
  purity_peifer = res$peifer
  purity_mustonen = res$mustonen
  
  combined_status = get_combined_status(segments, map_vanloowedge, map_dkfz, map_mustonen, map_peifer, map_broad)
  
  #####################################################################
  # Make inventory
  #####################################################################
  
  all_clonal = apply(combined_status[,c("broad", "dkfz", "peifer", "vanloowedge")], 1, function(x) { all(x=="clonal") })
  all_clonal[is.na(all_clonal)] = F
  all_subclonal = apply(combined_status[,c("broad", "dkfz", "peifer", "vanloowedge")], 1, function(x) { all(x=="subclonal") })
  all_subclonal[is.na(all_subclonal)] = F
  any_na = apply(combined_status[,c("broad", "dkfz", "peifer", "vanloowedge")], 1, function(x) { any(is.na(x)) })
  one_subclonal = apply(combined_status[,c("broad", "dkfz", "peifer", "vanloowedge")], 1, function(x) { sum(x=="subclonal")==1 })
  one_subclonal[is.na(one_subclonal)] = F
  two_subclonal = apply(combined_status[,c("broad", "dkfz", "peifer", "vanloowedge")], 1, function(x) { sum(x=="subclonal")==2 })
  two_subclonal[is.na(two_subclonal)] = F
  three_subclonal = apply(combined_status[,c("broad", "dkfz", "peifer", "vanloowedge")], 1, function(x) { sum(x=="subclonal")==3 })
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
  if (make_figures) {
    print("Making figures")
    
    make_plot = function(mapping, segments, plot_title, max.plot.cn) {
      print(plot_title)
      if (!is.na(mapping)) {
        cn_bb = collapse2bb(segments=segments, cn_states=mapping$cn_states)
      } else {
        cn_bb = parse_dummy_cn_profile()
      }
      plot = plot_profile(cn_bb, plot_title, max.plot.cn=max.plot.cn)
      return(plot)
    }
    
    plot_vanloowedge = make_plot(map_vanloowedge, segments, "Battenberg", max.plot.cn)
    plot_broad = make_plot(map_broad, segments, "ABSOLUTE", max.plot.cn)
    plot_dkfz = make_plot(map_dkfz, segments, "ACEseq", max.plot.cn)
    plot_mustonen = make_plot(map_mustonen, segments, "CloneHD", max.plot.cn)
    plot_peifer = make_plot(map_peifer, segments, "Sclust", max.plot.cn)
    print("Done plot profile")
    
    png(file.path(outdir, "figures", paste0(samplename, "_copynumber_complete.png")), height=1500, width=1300)
    grid.arrange(arrangeGrob(plot_broad, plot_dkfz, plot_vanloowedge, plot_mustonen, plot_peifer, ncol=1), 
                 top=textGrob(samplename, gp = gpar(fontsize=25, face=2, col="black")))
    dev.off()
  }
  
  #####################################################################
  # Make rounded clonal copy number
  #####################################################################
  print("Rounding segments")
  rounded_clonal = list(broad=data.frame(), dkfz=data.frame(), vanloowedge=data.frame(), peifer=data.frame(), mustonen=data.frame())
  for (i in 1:nrow(segments)) {
    
    if (!is.na(map_broad) && length(map_broad$cn_states) >= i) {
      rounded_clonal$broad = rbind(rounded_clonal$broad, round_broad(map_broad, i, rounding_up=rounding_up))
    } else {
      temp_entry = get_dummy_cn_entry(segments[i,,drop=F])
      temp_entry$historically_clonal = 0
      rounded_clonal$broad = rbind(rounded_clonal$broad, temp_entry)
    }
    
    if (!is.na(map_vanloowedge) && length(map_vanloowedge$cn_states) >= i) {
      rounded_clonal$vanloowedge = rbind(rounded_clonal$vanloowedge, round_vanloo_wedge(map_vanloowedge, i, purity_vanloowedge, rounding_up=rounding_up))
    } else {
      temp_entry = get_dummy_cn_entry(segments[i,,drop=F])
      colnames(temp_entry)[7] = "cellular_prevalence"
      temp_entry$ccf = NA
      rounded_clonal$vanloowedge = rbind(rounded_clonal$vanloowedge, temp_entry)
    }

    if (!is.na(map_peifer) && length(map_peifer$cn_states) >= i) {
      rounded_clonal$peifer = rbind(rounded_clonal$peifer, round_peifer(map_peifer, i, purity_peifer, rounding_up=rounding_up))
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
      rounded_clonal$dkfz = rbind(rounded_clonal$dkfz, round_dkfz(map_dkfz, i, purity_dkfz, rounding_up=rounding_up))
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

      if (rounding_up) {
        dir_postfix = paste0(cn_method, "_rounded_clonal")
      } else {
        dir_postfix = paste0(cn_method, "_rounded_alt_clonal")
      }
      write.table(rounded_clonal[[i]], file=file.path(outdir, dir_postfix, paste0(samplename, "_segments.txt")), quote=F, sep="\t", row.names=F)
  }
  
  #####################################################################
  # Plot rounded clonal
  #####################################################################
  if (make_figures) {
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
}


#####################################################################
# Script
#####################################################################
#' TODO
#'  - File paths at the top of the create_rounded_copynumber could be improved, maybe move all the input into an input directory

args = commandArgs(T)
samplename = args[1]
outdir = args[2]
rounding_up = as.logical(args[3])
make_figures = as.logical(args[4])
sex = args[5]

# setwd("/Users/sd11/Documents/Projects/icgc/consensus_subclonal_copynumber/6aa00162-6294-4ce7-b6b7-0c3452e24cd6")
# setwd("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/consensus_copynumber/2016_09_consensus_breakpoints_fullruns/20161024_rounding/")
# samplename = "6aa00162-6294-4ce7-b6b7-0c3452e24cd6"

# setwd("/Users/sd11/Documents/Projects/icgc/consensus_subclonal_copynumber/final_run_testing")
# samplename = "003819bc-c415-4e76-887c-931d60ed39e7"
# sex = "female"
# outdir = "output"


dkfz_segmentsfile = paste0("dkfz/", samplename, "_segments.txt")
dkfz_purityfile = "purity_ploidy_dkfz.txt"
vanloowedge_segmentsfile = paste0("vanloowedge/", samplename, "_segments.txt")
vanloowedge_purityfile = "purity_ploidy_vanloowedge.txt"
peifer_segmentsfile = paste0("peifer/", samplename, "_segments.txt")
peifer_purityfile = "purity_ploidy_peifer.txt"
#mustonen_segmentsfile = paste0("mustonen/", samplename, ".penalty0.95_segments.txt")
mustonen_segmentsfile = paste0("mustonen/", samplename, "_segments.txt")
mustonen_purityfile = "purity_ploidy_mustonen.txt"
broad_segmentsfile = paste0("broad/", samplename, "_segments.txt")
broad_purityfile = "purity_ploidy_broad.txt"
vanloowedge_baflogrfile = paste0("vanloowedge_baflogr/", samplename, "_baflogr.txt")
broad_baflogrfile = paste0("broad_baflogr/", samplename, "_baflogr.txt")

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

method_baflogr = list(vanloowedge=vanloowedge_baflogrfile,
                      broad=broad_baflogrfile)

breakpoints_file = file.path("consensus_bp", paste0(samplename, ".txt"))
if (file.exists(breakpoints_file)) {
  breakpoints = read.table(breakpoints_file, header=T, stringsAsFactors=F)
  segments = breakpoints2segments(breakpoints)
  create_rounded_copynumber(samplename, 
                            segments, 
                            outdir, 
                            method_segmentsfile, 
                            method_purityfile, 
                            method_baflogr, 
                            sex,
                            max.plot.cn=4, 
                            rounding_up=rounding_up, 
                            make_figures=make_figures)
}

q(save="no")







