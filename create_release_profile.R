#' Function that adds a line of NAs for missing (i.e. NULL) mapped entries
#' this convenience function is used to make the cn_states list have the same
#' number of entries as the segments it was mapped to. do.call(rbind,) does
#' not insert a line when the list entry is NULL
padd_empty_entries = function(map, methodname) {
  if (!is.null(map$cn_states[[1]])) {
    template = map$cn_states[[1]][[1]]
  } else if (!is.null(map$cn_states[[2]])) {
    template = map$cn_states[[2]][[1]]
  } else {
    stop(paste0("Could not find template for ", methodname))
  }
  template[1,1:ncol(template)] = NA
  
  cn_states = list()
  for (i in 1:length(map$cn_states)) {
    if (is.null(map$cn_states[[i]])) {
      cn_states[[i]] = template
    } else {
      cn_states[[i]] = map$cn_states[[i]][[1]]
    }
  }
  return(cn_states)
}

reset_overruled_annotations = function(anno, overrulings_pivot, methodid) {
  if (nrow(overrulings_pivot) > 0 & !is.na(overrulings_pivot[1, methodid])) {
    anno[1:nrow(anno), 1:ncol(anno)] = NA
  }
  return(anno)
}

combine_all_annotations = function(all_annotations, overrulings_pivot) {
  if (all(unlist(lapply(all_annotations$map_vanloowedge$cn_states, function(x) nrow(x[[1]]))) == 1)) {
    anno_vanloowedge = do.call(rbind, padd_empty_entries(all_annotations$map_vanloowedge, "Battenberg"))
    anno_vanloowedge = anno_vanloowedge[,c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A", "SDfrac_A", "SDfrac_A_BS", "frac1_A_0.025", "frac1_A_0.975")]
    colnames(anno_vanloowedge) = paste0("battenberg_", colnames(anno_vanloowedge))
  } else {
    print("Found too many annotations for some segments from Battenberg")
  }
  
  if (all(unlist(lapply(all_annotations$map_broad$cn_states, function(x) nrow(x[[1]]))) == 1)) {
    # anno_broad = do.call(rbind, lapply(all_annotations$map_broad$cn_states, function(x) x[[1]]))
    anno_broad = do.call(rbind, padd_empty_entries(all_annotations$map_broad, "ABSOLUTE"))
    anno_broad = anno_broad[,c("broad_major_cn", "broad_minor_cn", "broad_het_error", "broad_cov_error")]
    colnames(anno_broad) = paste0("absolute_", colnames(anno_broad))
  } else {
    print("Found too many annotations for some segments from ABSOLUTE")
  }
  
  if (all(unlist(lapply(all_annotations$map_dkfz$cn_states, function(x) nrow(x[[1]]))) == 1)) {
    # anno_dkfz = do.call(rbind, lapply(all_annotations$map_dkfz$cn_states, function(x) x[[1]]))
    anno_dkfz = do.call(rbind, padd_empty_entries(all_annotations$map_dkfz, "ACEseq"))
    anno_dkfz = anno_dkfz[,c("copy_number", "minor_cn", "major_cn", "ccf", "dh", "covRatio")]
    colnames(anno_dkfz) = paste0("aceseq_", colnames(anno_dkfz))
  } else {
    print("Found too many annotations for some segments from ACEseq")
  }
  
  if (all(unlist(lapply(all_annotations$map_mustonen$cn_states, function(x) nrow(x[[1]]))) == 1)) {
    # anno_mustonen = do.call(rbind, lapply(all_annotations$map_mustonen$cn_states, function(x) x[[1]]))
    anno_mustonen = do.call(rbind, padd_empty_entries(all_annotations$map_mustonen, "CloneHD"))
    anno_mustonen = anno_mustonen[,c("copy_number", "minor_cn", "major_cn")]
    anno_mustonen$ccf = 1
    colnames(anno_mustonen) = paste0("clonehd_", colnames(anno_mustonen))
  } else {
    print("Found too many annotations for some segments from CloneHD")
  }
  
  if (all(unlist(lapply(all_annotations$map_peifer$cn_states, function(x) nrow(x[[1]]))) == 1)) {
    # anno_peifer = do.call(rbind, lapply(all_annotations$map_peifer$cn_states, function(x) x[[1]]))
    anno_peifer = do.call(rbind, padd_empty_entries(all_annotations$map_peifer, "Sclust"))
    anno_peifer = anno_peifer[,c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]
    colnames(anno_peifer) = paste0("sclust_", colnames(anno_peifer))
  } else {
    print("Found too many annotations for some segments from Sclust")
  }
  
  # Check for whether a method has been overruled and reset annotations if needed
  anno_broad = reset_overruled_annotations(anno_broad, overrulings_pivot, "broad")
  anno_dkfz = reset_overruled_annotations(anno_dkfz, overrulings_pivot, "dkfz")
  anno_vanloowedge = reset_overruled_annotations(anno_vanloowedge, overrulings_pivot, "vanloo_wedge")
  anno_mustonen = reset_overruled_annotations(anno_mustonen, overrulings_pivot, "mustonen")
  anno_peifer = reset_overruled_annotations(anno_peifer, overrulings_pivot, "peifer")
  
  return(data.frame(anno_broad, anno_dkfz, anno_vanloowedge, anno_mustonen, anno_peifer))
}


# samplename = "6aa00162-6294-4ce7-b6b7-0c3452e24cd6"
source("~/repo/icgc_consensus_copynumber/util.R")

args = commandArgs(T)
samplename = args[1]
outdir = args[2]

cons_profile_file = file.path("output/consensus_profile", paste0(samplename, "_consensus_profile.txt"))
breakpoints_file = file.path("consensus_bp", paste0(samplename, ".txt"))
overrulings_pivot = readr::read_tsv("manual_review_overrulings_pivot_table.txt")
overrulings_pivot = overrulings_pivot[overrulings_pivot$samplename==samplename,]

if (file.exists(cons_profile_file) & file.exists(breakpoints_file)) {
  print("Reading in and mapping data...")
  dat = read.table(cons_profile_file, header=T, stringsAsFactors=F)  
  breakpoints = read.table(file.path("consensus_bp", paste0(samplename, ".txt")), header=T, stringsAsFactors=F)
  segments = breakpoints2segments(breakpoints)

  dkfz_segmentsfile = paste0("dkfz_annotations/", samplename, "_annotations.txt")
  dkfz_purityfile = "purity_ploidy_dkfz.txt"
  vanloowedge_segmentsfile = paste0("vanloowedge_annotations/", samplename, "_annotations.txt")
  vanloowedge_purityfile = "purity_ploidy_vanloowedge.txt"
  peifer_segmentsfile = paste0("peifer_annotations/", samplename, "_annotations.txt")
  peifer_purityfile = "purity_ploidy_peifer.txt"
  mustonen_segmentsfile = paste0("mustonen_annotations/", samplename, "_annotations.txt")
  mustonen_purityfile = "purity_ploidy_mustonen.txt"
  broad_segmentsfile = paste0("broad_annotations/", samplename, "_annotations.txt")
  broad_purityfile = "purity_ploidy_broad.txt"
  
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
  
  dat = na.omit(dat)
  dat$total_cn = dat$major_cn+dat$minor_cn
  dat = dat[, c("chromosome", "start", "end", "total_cn", "major_cn", "minor_cn", "star", "level")]
  
  # Map the annotations against the loaded consensus profile
  all_annotations = parse_all_profiles(samplename, dat, method_segmentsfile, method_purityfile, method_baflogr=NULL, mustonen_has_header=F)  
  combined_annotations = combine_all_annotations(all_annotations)
  
  # PCAWG11 profile with full annotations
  dat = data.frame(dat, combined_annotations)
  write.table(dat, file=file.path(outdir, "pcawg11_consensus_profile", paste0(samplename, ".consensus.20161103.somatic.cna.annotated.txt")), quote=F, row.names=F, sep="\t")
  
  # PCAWG-wide profile
  dat = dat[, c("chromosome", "start", "end", "total_cn", "major_cn", "minor_cn", "star")]
  write.table(dat, file=file.path(outdir, "consensus_profile_final", paste0(samplename, ".consensus.20161103.somatic.cna.txt")), quote=F, row.names=F, sep="\t")
}