suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(parallel))

# Static for plot
RECT_HEIGHT = 0.07

breakpoints2segments = function(breakpoints) {
  segments = data.frame()
  for (chrom in unique(breakpoints$chromosome)) {
    selection = breakpoints$chromosome==chrom
    breakpoints_chrom = breakpoints[selection,]
    
    segments = rbind(segments, 
                     data.frame(chromosome=chrom, start=breakpoints_chrom$position[1:(nrow(breakpoints_chrom)-1)], end=breakpoints_chrom$position[2:nrow(breakpoints_chrom)]-1))
  }
  return(segments)
}

parse_dkfz = function(segmentsfile, purityfile, samplename, dkfz_subclonality_cutoff=0.1, perform_rounding=T) {
  if (!is.na(segmentsfile) && file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
    
    purity = parse_dkfz_purity(purityfile, samplename)
    dat$ccf = dat$cellular_prevalence / purity
    
    # Check for X and Y in males as they don't have allele specific CN in the given files
    sel = (dat$chromosome==23 | dat$chromosome==24) & !is.na(dat$copy_number) & is.na(dat$major_cn) & is.na(dat$minor_cn)
    if (any(sel)) {
      # Set as major allele as we don't expect a minor allele in males
      dat$major_cn[sel] = dat$copy_number[sel]
      dat$minor_cn[sel] = 0
    }
    
    # If the CN does not deviate from an integer by a supplied cutoff the segment should be considered clonal. Here we round those values that are supposedly clonal
    cn_deviation = abs(dat$copy_number-round(dat$copy_number))
    cn_deviation = ifelse(is.na(cn_deviation), 0, cn_deviation)
    major_deviation = abs(dat$major_cn-round(dat$major_cn))
    major_deviation = ifelse(is.na(major_deviation), 0, major_deviation)
    minor_deviation = abs(dat$minor_cn-round(dat$minor_cn))
    minor_deviation = ifelse(is.na(minor_deviation), 0, minor_deviation)
    # Make an inventory
    below_subclonal_threshold = cn_deviation <= dkfz_subclonality_cutoff & major_deviation <= dkfz_subclonality_cutoff & minor_deviation <= dkfz_subclonality_cutoff
    # Round where appropriate
    if (perform_rounding & any(below_subclonal_threshold)) {
      dat$major_cn[below_subclonal_threshold] = round(dat$major_cn[below_subclonal_threshold])
      dat$minor_cn[below_subclonal_threshold] = round(dat$minor_cn[below_subclonal_threshold])
      dat$copy_number[below_subclonal_threshold] = dat$major_cn[below_subclonal_threshold] + dat$minor_cn[below_subclonal_threshold]
    }
    
    # Replace 23 and 24 with X and Y
    if (23 %in% dat$chromosome) {
      dat$chromosome[dat$chromosome==23] = "X"
    }
    
    if (24 %in% dat$chromosome) {
      dat$chromosome[dat$chromosome==24] = "Y"
    }
    
    maj_too_low = dat$major_cn < 0
    maj_too_low[is.na(maj_too_low)] = F
    min_too_low = dat$minor_cn < 0
    min_too_low[is.na(min_too_low)] = F
    if (any(maj_too_low | min_too_low)) {
      dat$major_cn[maj_too_low | min_too_low] = NA
      dat$minor_cn[maj_too_low | min_too_low] = NA
    }
    
    return(dat)
  } else {
    return(NA)
  }
}

parse_dkfz_purity = function(purityfile, samplename) {
  purity = read.table(purityfile, header=F, stringsAsFactors=F)
  purity = purity[purity$V1==samplename,3]
  if (length(purity)==0) {
    purity = NA
  }
  return(purity)
}

parse_vanloowedge = function(segmentsfile, purityfile, samplename, sex) {
  if (!is.na(segmentsfile) && file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
    purity = parse_vanloowedge_purity(purityfile, samplename)
    if ("clonal_frequency" %in% colnames(dat)) {
      dat$ccf = dat$clonal_frequency / purity
      colnames(dat)[7] = "cellular_prevalence"
    } else if ("cellular_prevalence" %in% colnames(dat)) {
      dat$ccf = dat$cellular_prevalence / purity
      colnames(dat)[7] = "cellular_prevalence"
    } else {
      # Annotations don't need any adjustments
    }
    
    # Remove calls for X as it is a combination of X and Y
    if (sex=="male") {
      dat = dat[dat$chromosome != "X",]
    }
    
    return(dat)
  } else {
    return(NA)
  }
}

parse_vanloowedge_purity = function(purityfile, samplename) {
  purity = read.table(purityfile, header=T, stringsAsFactors=F)
  purity = unique(purity[purity$sample==samplename,]$purity)
  if (length(purity) > 1) {
    print(paste0("parse_vanloowedge_purity - found multiple purities for sample ", samplename))
  }
  if (length(purity)==0) {
    purity = NA
  }
  return(purity)
}

parse_peifer = function(segmentsfile, purityfile, samplename) {
  if (!is.na(segmentsfile) && file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
    # Remove negative size segments
    dat = dat[(dat$end-dat$start) > 0, ]
    if ("cellular_prevalence" %in% colnames(dat)) {
      purity = parse_peifer_purity(purityfile, samplename)
      # What should be CP is encoded as CCF
      dat$ccf = dat$cellular_prevalence
      dat$cellular_prevalence = dat$ccf * purity
    } else {
      # Annotations don't need any adjustments
    }
    return(dat)
  } else {
    return(NA)
  }
}

parse_peifer_purity = function(purityfile, samplename) {
  purity = read.table(purityfile, header=T, stringsAsFactors=F)
  purity = unique(purity[purity$sample==samplename,]$purity)
  if (length(purity) > 1) {
    print(paste0("parse_peifer_purity - found multiple purities for sample ", samplename))
  }
  if (length(purity)==0) {
    purity = NA
  }
  return(purity)
}

parse_mustonen = function(segmentsfile, purityfile, samplename, has_header=F) {
  if (!is.na(segmentsfile) && file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=has_header, stringsAsFactors=F)
    if (!has_header) {
      colnames(dat) = c("chromosome", "start", "end", "copy_number", "major_cn", "minor_cn", "cellular_prevalence")
      # Remove 1bp segments
      dat = dat[dat$end-dat$start > 0,]
      # Add 1 so that segments do not overlap
      dat$end = dat$end - 1
    }
    purity = parse_mustonen_purity(purityfile, samplename)
    dat$ccf = dat$cellular_prevalence / purity
    
    # Apply a few filters to get rid of artifacts
    dat = dat[!(dat$start==dat$end),]
    dat = dat[!(dat$chromosome=="13" & dat$start==1),]
    dat = dat[!(dat$chromosome=="14" & dat$start==1),]
    dat = dat[!(dat$chromosome=="15" & dat$start==1),]
    dat = dat[!(dat$chromosome=="21" & dat$start==1),]
    
    # Replace 23 and 24 with X and Y
    if (23 %in% dat$chromosome) {
      dat$chromosome[dat$chromosome==23] = "X"
    }
    
    if (24 %in% dat$chromosome) {
      dat$chromosome[dat$chromosome==24] = "Y"
    }
    
    return(dat)
  } else {
    return(NA)
  }
}

parse_mustonen_purity = function(purityfile, samplename) {
  purity = read.table(purityfile, header=T, stringsAsFactors=F)#[1,2]
  purity = unique(purity[purity$sample==samplename,]$purity)
  # if (length(purity)==0) {
  #   purity = NA
  # }
  return(purity)
}

parse_broad = function(segmentsfile, purityfile, samplename) {
  if (!is.na(segmentsfile) && file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
    # Offset the start by 1 to make sure it does not overlap with the previous segment
    if ("ccf" %in% colnames(dat)) {
      dat$end = dat$end - 1
      dat = dat[!is.na(dat$copy_number) & !is.na(dat$major_cn) & !is.na(dat$minor_cn),]
    } else {
      # Annotations don't need any adustments
    }
    #' Data already in CCF supplied, no need to convert
    # colnames(dat) = c("chromosome", "start", "end", "copy_number", "major_cn", "minor_cn", "ccf")
    # purity = read.table(purityfile, header=F, stringsAsFactors=F)
    # dat$ccf = dat$cellular_prevalence / purity[1,2]
    # dat$cellular_prevalence = dat$ccf * purity
    return(dat)
  } else {
    return(NA)
  }
}

parse_broad_purity = function(purityfile, samplename) {
  purity = read.table(purityfile, header=T, stringsAsFactors=F)
  purity = unique(purity[purity$sample==samplename,]$purity)
  if (length(purity) > 1) {
    print(paste0("parse_broad_purity - found multiple purities for sample ", samplename))
  }
  if (length(purity)==0) {
    purity = NA
  }
  return(purity)
}

parse_jabba = function(segmentsfile) {
  if (!is.na(segmentsfile) && file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
    dat$ccf = 1 # only clonal CN
    
    if (any(dat$chromosome=="Y")) {
      dat$major_cn[dat$chromosome=="Y"] = dat$copy_number[dat$chromosome=="Y"]
      dat$minor_cn[dat$chromosome=="Y"] = 0
    }
    
    return(dat[,c("chromosome", "start", "end", "copy_number", "major_cn", "minor_cn", "cellular_prevalence", "ccf")])
  } else {
    return(NA)
  }
}

parse_jabba_purity = function(purityfile, samplename) {
  purity = read.table(purityfile, header=T, stringsAsFactors=F)#[1,2]
  if (samplename %in% purity$sample) {
    purity = unique(purity[purity$sample==samplename,]$purity)
    return(purity)
  } else {
    return(NA)
  }
}

parse_all_profiles = function(samplename, segments, method_segmentsfile, method_purityfile, method_baflogr, sex, mustonen_has_header=F, cn_round_dkfz=T, num_threads=1) {

  do_mapping = function(index, all_dat, method_names, segments) {
    dat = all_dat[[index]]
    method_name = method_names[index]
    if (!is.na(dat)) {
      dat_map = mapdata(segments, dat, is_dkfz=method_name=="dkfz", is_broad=method_name=="broad")
    } else {
      dat_map = NA
    }
    return(list(dat=dat, map=dat_map))
  }
  
  dat_dkfz = parse_dkfz(method_segmentsfile[["dkfz"]], method_purityfile[["dkfz"]], samplename, perform_rounding=cn_round_dkfz)
  dat_vanloowedge = parse_vanloowedge(method_segmentsfile[["vanloowedge"]], method_purityfile[["vanloowedge"]], samplename, sex)
  dat_peifer = parse_peifer(method_segmentsfile[["peifer"]], method_purityfile[["peifer"]], samplename)
  dat_mustonen = parse_mustonen(method_segmentsfile[["mustonen"]], method_purityfile[["mustonen"]], samplename, has_header=mustonen_has_header)
  dat_broad = parse_broad(method_segmentsfile[["broad"]], method_purityfile[["broad"]], samplename)
  dat_jabba = parse_jabba(method_segmentsfile[["jabba"]])
  
  res = mclapply(1:6,
               do_mapping,
               list(dat_dkfz, dat_vanloowedge, dat_peifer, dat_mustonen, dat_broad, dat_jabba),
               c("dkfz", "vanloowedge", "peifer", "mustonen", "broad", "jabba"),
               segments,
               mc.cores=num_threads)
  
  map_dkfz = res[[1]]$map
  map_vanloowedge = res[[2]]$map
  map_peifer = res[[3]]$map
  map_mustonen = res[[4]]$map
  map_broad = res[[5]]$map
  map_jabba = res[[6]]$map
  
  if (!is.null(method_baflogr)) {
    if (file.exists(method_baflogr$vanloowedge)) {
      baflogr_vanloowedge = read.table(method_baflogr$vanloowedge, header=T, stringsAsFactors=F)
      map_vanloowedge_baflogr = mapdata(segments, baflogr_vanloowedge)
      
      # Pad empty data in case last segments / chromosomes were not reported on
      if (length(map_vanloowedge_baflogr$cn_states) < nrow(segments)) {
        for (i in length(map_vanloowedge_baflogr$cn_states):nrow(segments)) {
          map_vanloowedge_baflogr$cn_states[[i]] = NA
          map_vanloowedge_baflogr$status[i] = NA
        }
      }
      
    } else { 
      map_vanloowedge_baflogr = NA
    }
    
    if (file.exists(method_baflogr$broad)) {
      baflogr_broad = read.table(method_baflogr$broad, header=T, stringsAsFactors=F)
      map_broad_baflogr = mapdata(segments, baflogr_broad)
      
      # Pad empty data in case last segments / chromosomes were not reported on
      if (length(map_broad_baflogr$cn_states) < nrow(segments)) {
        for (i in length(map_broad_baflogr$cn_states):nrow(segments)) {
          map_broad_baflogr$cn_states[[i]] = NA
          map_broad_baflogr$status[i] = NA
        }
      }
    } else { 
      map_broad_baflogr = NA
    }
    
  } else {
    map_vanloowedge_baflogr = NA
    map_broad_baflogr = NA
  }
  
  return(list(dat_dkfz=dat_dkfz, map_dkfz=map_dkfz,
              dat_vanloowedge=dat_vanloowedge, map_vanloowedge=map_vanloowedge,
              dat_peifer=dat_peifer, map_peifer=map_peifer,
              dat_mustonen=dat_mustonen, map_mustonen=map_mustonen,
              dat_broad=dat_broad, map_broad=map_broad,
              dat_jabba=dat_jabba, map_jabba=map_jabba,
              map_vanloowedge_baflogr=map_vanloowedge_baflogr, map_broad_baflogr=map_broad_baflogr))
}

parse_all_purities = function(samplename, method_purityfile) {
  purity_broad = parse_broad_purity(method_purityfile[["broad"]], samplename)
  purity_dkfz = parse_dkfz_purity(method_purityfile[["dkfz"]], samplename)
  purity_mustonen = parse_mustonen_purity(method_purityfile[["mustonen"]], samplename)
  purity_peifer = parse_peifer_purity(method_purityfile[["peifer"]], samplename)
  purity_vanloowedge = parse_vanloowedge_purity(method_purityfile[["vanloowedge"]], samplename)
  purity_jabba = parse_jabba_purity(method_purityfile[["jabba"]], samplename)
  return(list(broad=purity_broad, dkfz=purity_dkfz, mustonen=purity_mustonen, peifer=purity_peifer, vanloowedge=purity_vanloowedge, jabba=purity_jabba))
}

parse_dummy_cn_profile = function(nmaj=NA, nmin=NA) {
  temp = read.table("segmentation_chrom_arms_full.txt", header=T, stringsAsFactors=F)
  temp$nMaj1_A = nmaj
  temp$nMin1_A = nmin
  return(temp)
}

get_dummy_cn_entry = function(segment) {
  return(data.frame(chromosome=segment$chromosome, start=segment$start, end=segment$end, copy_number=NA, major_cn=NA, minor_cn=NA, ccf=NA))
}

#' DKFZ does not make separate calls, a deviation from integer should be used to detect subclonality
asses_dkfz_clonal_status = function(cn_segments, i, dkfz_subclonality_cutoff) {
  cn_deviation = abs(cn_segments$copy_number[i]-round(cn_segments$copy_number[i]))
  cn_deviation = ifelse(is.na(cn_deviation), 0, cn_deviation)
  
  major_deviation = abs(cn_segments$major_cn[i]-round(cn_segments$major_cn[i]))
  major_deviation = ifelse(is.na(major_deviation), 0, major_deviation)
  
  minor_deviation = abs(cn_segments$minor_cn[i]-round(cn_segments$minor_cn[i]))
  minor_deviation = ifelse(is.na(minor_deviation), 0, minor_deviation)
  
  if (any(c(cn_deviation, major_deviation, minor_deviation) > dkfz_subclonality_cutoff)) {
    # enough deviation, so subclonal
    status = "subclonal"
  } else {
    # not enough deviation, therefore clonal
    status = "clonal"
  }
  return(status)
}

#' Map reported cn segments to the given bp segments
#' @return Yields a list with two fields: status (with clonal/subclonal/NA classifications) and cn_states (with the assigned cn states)
mapdata = function(bp_segments, cn_segments, is_dkfz=F, dkfz_subclonality_cutoff=0.1, is_broad=F) {
  
  merge_broad_segments = function(cn_segments, overlap, bp_segment) {
    # Perform merging of clonal segments
    temp_segs = cn_segments[queryHits(overlap),,drop=F]
    
    # Case 1: Only a single segment left
    if (nrow(temp_segs)==1) { return(temp_segs) }
    
    # Case 2: Only a single clonal segment left
    if (sum(temp_segs$historically_clonal==1)==1 & sum(temp_segs$historically_clonal==0) > 0) { return(temp_segs) }
    
    # Case 3: All elements the same state, just merge
    if (isTRUE(all.equal(max(temp_segs$major_cn), min(temp_segs$major_cn))) & isTRUE(all.equal(max(temp_segs$minor_cn), min(temp_segs$minor_cn)))) {
      merged_entry = temp_segs[1,,drop=F]
      merged_entry$end[1] = temp_segs$end[nrow(temp_segs)]
      return(merged_entry)
    }
    
    # Case 4: Attempt to merge after removing small segments that fall within the consensus segment
    # Remove small segments that fall completely within the consensus segment - if the consensus segment is large enough
    is_to_small = (temp_segs$start > bp_segment$start & temp_segs$end < bp_segment$end & (temp_segs$end-temp_segs$start) < 1000000 & (bp_segment$end-bp_segment$start) > 3000000)
    temp_segs = temp_segs[!is_to_small,]
    
    merged = T
    while (merged & nrow(temp_segs) > 1) {
      merged = F
      merged_temp_segs = data.frame()
      prev = NULL
      for (j in 2:nrow(temp_segs)) {
        if (temp_segs$minor_cn[j-1]==temp_segs$minor_cn[j] & temp_segs$major_cn[j-1]==temp_segs$major_cn[j] & temp_segs$ccf[j]==1 & temp_segs$ccf[j-1]) {
          # Can merge as major/minor states the same
          merged_entry = temp_segs[j-1,]
          merged_entry$end = temp_segs$end[j]
          merged_temp_segs = rbind(merged_temp_segs, merged_entry)
          merged = T
        } else if (j==nrow(temp_segs)) {
          # Last step and cannot merge the final segment, so add the last two segments to complete
          merged_temp_segs = rbind(merged_temp_segs, temp_segs[j-1,])
          merged_temp_segs = rbind(merged_temp_segs, temp_segs[j,])
        } else {
          # No merging done, just add one segment
          merged_temp_segs = rbind(merged_temp_segs, temp_segs[j-1,])
        }
      }
      if (merged) {
        temp_segs = merged_temp_segs
      }
    }
    return(temp_segs)
  }
  
  map_segment = function(i, cns_gr, bps_gr, bp_segments, cn_segments, is_dkfz, dkfz_subclonality_cutoff, is_broad) {
    overlap = findOverlaps(cns_gr, bps_gr[i,])
    
    # No overlap, no call
    if (length(overlap)==0) {
      return(list(cn_states=NA, status=NA))
      
      # One segment overlaps, but could be subclonal in DKFZ output
    } else if (length(overlap)==1) {
      
      if (is_dkfz) {
        status = asses_dkfz_clonal_status(cn_segments, queryHits(overlap), dkfz_subclonality_cutoff)
      } else {
        # clonal
        status = "clonal"
      }
      cn_states = list(cn_segments[queryHits(overlap),])
      
      # More than one segment overlaps, this could be subclonal, but it could also be that another segment slightly  
      # overlaps with the breakpoint defined segment and there are not really two calls that overlap
    } else { 
      # get only segments that overlap at least 50%
      overlap = findOverlaps(cns_gr, bps_gr[i,], minoverlap=round((bp_segments$end[i]-bp_segments$start[i])*0.5))
      
      # No segments overlap 50%
      if (length(overlap)==0 & is_broad) {
        # Get segments that overlap at least 1% of the consensus segment to get rid of those that overlap just one bp
        overlap = findOverlaps(cns_gr, bps_gr[i,], minoverlap=round((bp_segments$end[i]-bp_segments$start[i])*0.01))
        
        # Sometimes the next segment just bleeds in, so here remove entries that overlap substantially with the next or previous segment
        # This is not ideal as it may be that this segment was merged with the next and we're removing legit signal here. However there must
        # be more than one segment aligning
        if (i < nrow(bp_segments) & length(overlap) > 1) {
          overlap_next = findOverlaps(cns_gr, bps_gr[i+1,], minoverlap=round((bp_segments$end[i+1]-bp_segments$start[i+1])*0.8))
          overlap = setdiff(overlap, overlap_next)
        }
        
        if (i > 1 & length(overlap) > 1) {
          # Remove segments that substantially overlap with the previous segment
          overlap_prev = findOverlaps(cns_gr, bps_gr[i-1,], minoverlap=round((bp_segments$end[i-1]-bp_segments$start[i-1])*0.8))
          overlap = setdiff(overlap, overlap_prev)
        }
        
        if (length(overlap)==0) {
          print(paste0("mapdata broad - found multiple clonal segments that overlap, but also with other segments ", i))
          return(list(cn_states=NA, status=NA))
        }
        
        # temp_segs = merge_broad_segments(cn_segments, overlap, bp_segments[i,,drop=F])
        temp_segs = cn_segments[queryHits(overlap),,drop=F]
        if (nrow(temp_segs) == 1) {
          status = "clonal"
          cn_states = list(temp_segs)
          return(list(cn_states=NA, status=NA)) # skip the rest of the loop as it attempts to store the unmerged segments
        } else if (sum(temp_segs$historically_clonal==1)==1 & sum(temp_segs$historically_clonal==0)>=1) {
          # Subclonality is encoded as one historical and at least one not state
          status = "subclonal"
          cn_states = list(temp_segs)
        } else {
          print(bp_segments[i,])
          print(temp_segs)
          print(paste0("mapdata broad - found multiple clonal segments that cannot be merged for single consensus segment ", i))
          return(list(cn_states=NA, status=NA))
        }
        
        
        # No segments overlap 50% - other pipelines
      } else if (length(overlap)==0) {
        # In this case the cn segment is smaller than the breakpoint given segment
        # find the breakpoint segment that overlaps with 95% of the given cn segment and make that mapping
        segment_overlaps_bp = round(cn_segments$end[queryHits(overlap)]-cn_segments$start[queryHits(overlap)])
        
        # Check for zero sized segments
        if (sum(segment_overlaps_bp > 0)==0) {
          # No more candidates left, so no overlap
          return(list(cn_states=NA, status=NA))
        }
        
        # Remove segments that are 1bp long
        smallest_segment_overlaps = min(segment_overlaps_bp[segment_overlaps_bp > 0])
        overlap = findOverlaps(cns_gr, bps_gr[i,], minoverlap=round(smallest_segment_overlaps*0.95))
        
        # There are no segments that overlap considerably with the given segment, no call
        if (length(overlap)==0) {
          return(list(cn_states=NA, status=NA))
          
          # One segment overlaps considerably, therefore clonal
        } else if (length(overlap)==1) {
          if (is_dkfz) {
            status = asses_dkfz_clonal_status(cn_segments, queryHits(overlap), dkfz_subclonality_cutoff)
          } else {
            # clonal
            status = "clonal"
          }
          
          # Two segments overlap considerably, therefore subclonal
        } else {
          status = "subclonal"
        }
        
        # clonal - because only one segment overlaps considerably
      } else if (length(overlap)==1) {
        if (is_dkfz) {
          status = asses_dkfz_clonal_status(cn_segments, queryHits(overlap), dkfz_subclonality_cutoff)
        } else {
          # clonal
          status = "clonal"
        }
        
        # subclonal - because multiple segments overlap considerably
      } else {
        status = "subclonal"
      }
      cn_states = list(cn_segments[queryHits(overlap),])
    }
    return(list(cn_states=cn_states, status=status))
  }
  
  
  bps_gr = makeGRangesFromDataFrame(bp_segments)
  cns_gr = makeGRangesFromDataFrame(cn_segments, keep.extra.columns=T)
  
  status = rep(NA, length(bps_gr))
  cn_states = list()
  
  # res = mclapply(1:nrow(bp_segments), map_segment, cns_gr, bps_gr, bp_segments, cn_segments, is_dkfz, dkfz_subclonality_cutoff, is_broad, mc.cores=2)
  res = lapply(1:nrow(bp_segments), map_segment, cns_gr, bps_gr, bp_segments, cn_segments, is_dkfz, dkfz_subclonality_cutoff, is_broad)
  
  for (i in 1:length(res)) {
    if (!is.null(res[[i]]$status) & !is.null(res[[i]]$cn_states)) {
      status[[i]] = res[[i]]$status
      cn_states[[i]] = res[[i]]$cn_states
    }
  }
  
  # # for (i in 1:length(overlap)) {
  # for (i in 1:nrow(bp_segments)) {
  #   res = map_segment(i, cns_gr, bps_gr, bp_segments, cn_segments, is_dkfz, dkfz_subclonality_cutoff, is_broad)
  #   if (!is.null(res$status) & !is.null(res$cn_states)) {
  #     status[[i]] = res$status
  #     cn_states[[i]] = res$cn_states
  #   }
  # }
  return(list(status=status, cn_states=cn_states)) 
}

#' Collapse the cn states inventory into a BB profile that can be plotted
collapse2bb = function(segments, cn_states, broad=F) {
  bb_template = parse_bb_template()
  cn_bb = data.frame()
  for (i in 1:length(cn_states)) {

    new_bb_seg = bb_template
    new_bb_seg$chr[1] = segments$chromosome[i]
    new_bb_seg$startpos[1] = segments$start[i]
    new_bb_seg$endpos[1] = segments$end[i]
    
    cn_states_i = cn_states[[i]][[1]]
    
    # Take the largest segment if there are multiple clonal segments from the broad data
    if (broad && !is.null(cn_states_i) && !is.na(cn_states_i) && nrow(cn_states_i)>0) {
      are_clonal = which(cn_states_i$ccf==1)
      if (length(are_clonal) > 1) {
        largest = are_clonal[which.max(cn_states_i$end[are_clonal] - cn_states_i$start[are_clonal])]
        cn_states_i = cn_states_i[c(largest, which(cn_states_i$ccf!=1)),]
      }
    }
    
    if (is.null(cn_states_i) || is.na(cn_states_i) || nrow(cn_states_i)==0) {
      # Do not add a segment for where there is no call
      next
    } else if (nrow(cn_states_i)==0) {
      # Do not add a segment for where there is no call
      next
    } else if (nrow(cn_states_i)==1) {
      # Clonal copy number
      new_bb_seg$nMaj1_A[1] = cn_states_i$major_cn
      new_bb_seg$nMin1_A[1] = cn_states_i$minor_cn
      new_bb_seg$frac1_A = 1
      cn_bb = rbind(cn_bb, new_bb_seg)
    } else if (nrow(cn_states_i)==2 & !broad) {
      # Subclonal copy number
      new_bb_seg$nMaj1_A[1] = cn_states_i$major_cn[1]
      new_bb_seg$nMin1_A[1] = cn_states_i$minor_cn[1]
      new_bb_seg$frac1_A[1] = cn_states_i$ccf[1]
      new_bb_seg$nMaj2_A[1] = cn_states_i$major_cn[2]
      new_bb_seg$nMin2_A[1] = cn_states_i$minor_cn[2]
      new_bb_seg$frac2_A[1] = cn_states_i$ccf[2]
      cn_bb = rbind(cn_bb, new_bb_seg)
    } else if (nrow(cn_states_i)==2 & broad) {
      # Subclonal copy number - 2 states only
      ancestral = which(cn_states_i$ccf==1)
      descendant = (1:2)[(1:2)!=ancestral]
      
      new_bb_seg$nMaj1_A[1] = cn_states_i$major_cn[ancestral]
      new_bb_seg$nMin1_A[1] = cn_states_i$minor_cn[ancestral]
      new_bb_seg$frac1_A[1] = cn_states_i$ccf[ancestral] - cn_states_i$ccf[descendant]
      
      new_bb_seg$nMaj2_A[1] = cn_states_i$major_cn[descendant]
      new_bb_seg$nMin2_A[1] = cn_states_i$minor_cn[descendant]
      new_bb_seg$frac2_A[1] = cn_states_i$ccf[descendant]
      
      cn_bb = rbind(cn_bb, new_bb_seg)
    } else if (nrow(cn_states_i)>=2 & broad) {
      # Subclonal copy number - more than 2 states
      # TODO: averaging the ancestral with each subclone may be better
      ancestral = which(cn_states_i$ccf==1)
      descendant = which(cn_states_i$ccf!=1)
      
      new_bb_seg$nMaj1_A[1] = cn_states_i$major_cn[ancestral]
      new_bb_seg$nMin1_A[1] = cn_states_i$minor_cn[ancestral]
      new_bb_seg$frac1_A[1] = cn_states_i$ccf[ancestral] - sum(cn_states_i$ccf[descendant])
      
      # Combine the subclones into a single
      new_bb_seg$nMaj2_A[1] = sum(sapply(descendant, function(x) cn_states_i$major_cn[x]*cn_states_i$ccf[x]))
      new_bb_seg$nMin2_A[1] = sum(sapply(descendant, function(x) cn_states_i$minor_cn[x]*cn_states_i$ccf[x]))
      new_bb_seg$frac2_A[1] = sum(cn_states_i$ccf[descendant])
      cn_bb = rbind(cn_bb, new_bb_seg)
      # print(paste0("Too many fits, cannot put into data format. segment: ", i))
    } else {
      print(paste0("Too many fits, cannot put into data format. segment: ", i))
    }
  }
  return(cn_bb)
}

#' Collapse data in rounded clonal format to BB native
collapseRoundedClonal2bb = function(cn_states) {
  bb_template = parse_bb_template()
  cn_bb = data.frame()
  
  for (i in 1:nrow(cn_states)) {
    
    new_bb_seg = bb_template
    new_bb_seg$chr = cn_states$chromosome[i]
    new_bb_seg$startpos = cn_states$start[i]
    new_bb_seg$endpos = cn_states$end[i]
    
    new_bb_seg$nMaj1_A[1] = cn_states$major_cn[i]
    new_bb_seg$nMin1_A[1] = cn_states$minor_cn[i]
    new_bb_seg$frac1_A = 1
    cn_bb = rbind(cn_bb, new_bb_seg)
    
  }
  return(cn_bb)
}

createPng = function(p, filename, height, width) { png(filename, height=height, width=width); print(p); dev.off() }

plot_profile = function(subclones, method_name, max.plot.cn=3) {
  subclones$chr = factor(subclones$chr, levels=mixedsort(unique(subclones$chr)))
  subclones$plot_maj = ifelse(is.na(subclones$frac2_A), subclones$nMaj1_A, subclones$nMaj1_A*subclones$frac1_A + subclones$nMaj2_A*subclones$frac2_A)
  subclones$plot_min = ifelse(is.na(subclones$frac2_A), subclones$nMin1_A, subclones$nMin1_A*subclones$frac1_A + subclones$nMin2_A*subclones$frac2_A)
  subclones$plot_cn_total = subclones$plot_maj + subclones$plot_min
  p1 = ggplot(subclones) + 
    geom_hline(data=data.frame(y=seq(0,max.plot.cn,0.5)), mapping=aes(slope=0, yintercept=y), colour="black", alpha=0.3) +
    # Not plotting the major allele
    # geom_rect(mapping=aes(xmin=startpos, xmax=endpos, ymin=plot_maj, ymax=(plot_maj+RECT_HEIGHT)), fill="purple") +
    # Removed shifting to prevent overplotting
    # geom_rect(mapping=aes(xmin=startpos, xmax=endpos, ymin=plot_min-RECT_HEIGHT, ymax=(plot_min)), fill="#2f4f4f") +
    # geom_rect(mapping=aes(xmin=startpos, xmax=endpos, ymin=plot_cn_total, ymax=(plot_cn_total+RECT_HEIGHT)), fill="#E69F00") +
    geom_rect(mapping=aes(xmin=startpos, xmax=endpos, ymin=(plot_min-RECT_HEIGHT), ymax=(plot_min+RECT_HEIGHT)), fill="#2f4f4f") +
    geom_rect(mapping=aes(xmin=startpos, xmax=endpos, ymin=(plot_cn_total-RECT_HEIGHT), ymax=(plot_cn_total+RECT_HEIGHT)), fill="#E69F00") +
    facet_grid(~chr, scales="free_x", space = "free_x") +
    scale_y_continuous(breaks=0:max.plot.cn, limits=c(-0.2,max.plot.cn), name="Copy Number") +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y = element_text(colour="black",size=22,face="plain"),  
                       axis.title.y = element_text(colour="black",size=24,face="plain"),
                       strip.text.x = element_text(colour="black",size=24,face="plain"),
                       plot.title = element_text(colour="black",size=24,face="plain")) + 
    ggtitle(method_name)
  return(p1)
}

parse_bb_template = function() {
  bb_template = read.table("segmentation_chrom_arms_full.txt", header=T, stringsAsFactors=F)[1,,drop=F]
  bb_template$chr = NA
  bb_template$startpos = NA
  bb_template$endpos = NA
  return(bb_template)
}

#####################################################################
# Round subclonal CNAs
#####################################################################
round_vanloo_wedge = function(map, i, purity, rounding_up=T) {
  if (!is.null(map$cn_states[[i]]) && !is.na(map$cn_states[[i]]) && nrow(map$cn_states[[i]][[1]]) > 1) {
    dat = map$cn_states[[i]][[1]]
    if (rounding_up) {
      index_major_clone = which.max(dat$ccf)
    } else {
      index_major_clone = which.min(dat$ccf)
    }
    dat = map$cn_states[[i]][[1]][index_major_clone,,drop=F]
    dat$cellular_prevalence = purity
    dat$ccf = 1
    return(dat)
  } else if (is.null(map$cn_states[[i]]) || is.na(map$cn_states[[i]])) {
    return(data.frame())
  } else {
    return(map$cn_states[[i]][[1]])
  }
}

round_peifer = function(map, i, purity, rounding_up=T) {
  if (!is.null(map$cn_states[[i]]) && !is.na(map$cn_states[[i]]) && nrow(map$cn_states[[i]][[1]]) > 1) {
    dat = map$cn_states[[i]][[1]]
    if (rounding_up) {
      index_major_clone = which.max(dat$ccf)
    } else {
      index_major_clone = which.min(dat$ccf)
    }
    dat = map$cn_states[[i]][[1]][index_major_clone,,drop=F]
    dat$cellular_prevalence = purity
    dat$ccf = 1
    return(dat)
  } else if (is.null(map$cn_states[[i]]) || is.na(map$cn_states[[i]])) {
    return(data.frame())
  } else {
    return(map$cn_states[[i]][[1]])
  }
}

round_mustonen = function(map, i) {
  if (!is.na(map$cn_states[[i]][[1]])) {
    return(map$cn_states[[i]][[1]][1,,drop=F])
  } else {
    return(data.frame())
  }
}

round_jabba = function(map, i) {
  if (!is.na(map$cn_states[[i]][[1]])) {
    return(map$cn_states[[i]][[1]][1,,drop=F])
  } else {
    return(data.frame())
  }
}

round_dkfz = function(map, i, purity, rounding_up=T) {
  if (!is.null(map$cn_states[[i]]) && !is.na(map$cn_states[[i]]) && nrow(map$cn_states[[i]][[1]]) == 1) {
    temp = map$cn_states[[i]][[1]]
    if (rounding_up) {
      temp$minor_cn = ceiling(temp$minor_cn)
      temp$major_cn = ceiling(temp$major_cn)
    } else {
      temp$minor_cn = floor(temp$minor_cn)
      temp$major_cn = floor(temp$major_cn)
    }
    temp$copy_number = temp$minor_cn + temp$major_cn
    temp$cellular_prevalence = purity
    temp$ccf = 1
    return(temp)
  } else if (is.null(map$cn_states[[i]]) || is.na(map$cn_states[[i]])) {
    return(data.frame())
  } else {
    return(data.frame())
    print(paste0("round_dkfz - fit contains multiple states for segment ", i))
  }
}

round_broad = function(map, i, rounding_up=T) {
  if (!is.null(map$cn_states[[i]]) && !is.na(map$cn_states[[i]]) && nrow(map$cn_states[[i]][[1]]) > 1) {
    dat = map$cn_states[[i]][[1]]
    if (rounding_up & sum(dat$historically_clonal==1)==1) {
      # Rounding up means taking the historically clonal state
      dat = dat[dat$historically_clonal==1,,drop=F]
    } else if (!rounding_up & sum(dat$historically_clonal==1)==1) {
      # Rounding down means taking the major subclonal state
      dat = dat[which.max(dat$ccf[dat$historically_clonal < 1]),,drop=F]
    } else {
      print(paste0("round_broad - multiple segments, pick largest for segment ", i))
      dat = dat[which.max(dat$end-dat$start),,drop=F]
    }
    dat$ccf = 1
    return(dat)
  } else if (is.null(map$cn_states[[i]]) || is.na(map$cn_states[[i]])) {
    return(data.frame())
  } else {
    return(map$cn_states[[i]][[1]])
  }
}

get_combined_status = function(segments, map_vanloowedge, map_dkfz, map_mustonen, map_peifer, map_broad, map_jabba) {
  
  get_status = function(map) {
    if (is.na(map)) {
      status = rep(NA, nrow(segments))
    } else {
      status = map$status
    }
    return(status)
  }
  
  vanloowedge = get_status(map_vanloowedge)
  dkfz = get_status(map_dkfz)  
  mustonen = get_status(map_mustonen)
  peifer = get_status(map_peifer)
  broad = get_status(map_broad)
  jabba = get_status(map_jabba)  
  
  combined_status = data.frame(segments, dkfz=dkfz, mustonen=mustonen, peifer=peifer, vanloowedge=vanloowedge, broad=broad, jabba=jabba)
  return(combined_status)
}

calc_ploidy = function(subclones) {
  subclones$length = round((subclones$endpos-subclones$startpos)/1000)
  cn_state_one = (subclones$nMaj1_A+subclones$nMin1_A)*subclones$frac1_A
  cn_state_two = ifelse(!is.na(subclones$frac2_A), (subclones$nMaj2_A+subclones$nMin2_A)*subclones$frac2_A, 0)
  ploidy = sum((cn_state_one+cn_state_two) * subclones$length, na.rm=T) / sum(subclones$length, na.rm=T)
  return(ploidy)
}



get_ploidy_status = function(subclones, min_frac_genome_state=0.2) {
  seg_length = (subclones$endpos/1000)-(subclones$startpos/1000)
  is_clonal = is.na(subclones$frac2_A) & !is.na(subclones$nMaj1_A)
  dipl = subclones$nMin1_A==1 & subclones$nMin1_A==1
  tetrpl = subclones$nMin1_A==2 & subclones$nMin1_A==2
  dipl_frac_genome = sum(seg_length[subclones$nMin1_A==1 & subclones$nMin1_A==1], na.rm=T) / sum(seg_length)
  tetrpl_frac_genome = sum(seg_length[subclones$nMin1_A==2 & subclones$nMin1_A==2], na.rm=T) / sum(seg_length)
  
  if (sum(seg_length[dipl & is_clonal], na.rm=T) >= sum(seg_length[tetrpl & is_clonal], na.rm=T) & dipl_frac_genome > min_frac_genome_state) {
    status = "diploid"
  } else if (tetrpl_frac_genome > min_frac_genome_state) {
    status = "tetraploid"
  } else {
    status = "other"
  }
  return(status)
}

get_ploidy = function(segments, map, broad=F) {
  if (!is.na(map)) {
    cn_bb = collapse2bb(segments=segments, cn_states=map$cn_states, broad=broad)
    return(list(ploidy=round(calc_ploidy(cn_bb), 4), status=get_ploidy_status(cn_bb)))
  } else {
    return(list(ploidy=NA, status=NA))
  }
}

