library(GenomicRanges)
library(gtools)
library(ggplot2)
library(grid)
library(gridExtra)

# Static for plot
RECT_HEIGHT = 0.15

breakpoints2segments = function(breakpoints) {
  segments = data.frame()
  for (chrom in unique(breakpoints$chromosome)) {
    selection = breakpoints$chromosome==chrom
    breakpoints_chrom = breakpoints[selection,]
    
    segments = rbind(segments, 
                     data.frame(chromosome=chrom, start=breakpoints_chrom$position[1:(nrow(breakpoints_chrom)-1)], end=breakpoints_chrom$position[2:nrow(breakpoints_chrom)]))
  }
  return(segments)
}

parse_dkfz = function(segmentsfile, purityfile, samplename, dkfz_subclonality_cutoff=0.1) {
  if (file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
    
    purity = read.table(purityfile, header=F, stringsAsFactors=F)
    dat$ccf = dat$cellular_prevalence / purity[purity$V1==samplename,4]
    
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
    if (any(below_subclonal_threshold)) {
      dat$major_cn[below_subclonal_threshold] = round(dat$major_cn[below_subclonal_threshold])
      dat$minor_cn[below_subclonal_threshold] = round(dat$minor_cn[below_subclonal_threshold])
      dat$copy_number[below_subclonal_threshold] = dat$major_cn[below_subclonal_threshold] + dat$minor_cn[below_subclonal_threshold]
    }
    return(dat)
  } else {
    return(NA)
  }
}

parse_vanloowedge = function(segmentsfile, purityfile, samplename) {
  if (file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
    purity = read.table(purityfile, header=T, stringsAsFactors=F)
    if ("clonal_frequency" %in% colnames(dat)) {
      dat$ccf = dat$clonal_frequency / purity[purity$sample==samplename,]$purity
    } else {
      dat$ccf = dat$cellular_prevalence / purity[purity$sample==samplename,]$purity
    }
    colnames(dat)[7] = "cellular_prevalence"
    return(dat)
  } else {
    return(NA)
  }
}

parse_peifer = function(segmentsfile, purityfile, samplename) {
  if (file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
    purity = read.table(purityfile, header=T, stringsAsFactors=F)
    dat$ccf = dat$cellular_prevalence
    dat$cellular_prevalence = dat$ccf * purity[purity$sample==samplename,]$purity
    return(dat)
  } else {
    return(NA)
  }
}

parse_mustonen = function(segmentsfile, purityfile, samplename, has_header=F) {
  if (file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=has_header, stringsAsFactors=F)
    if (!has_header) {
      colnames(dat) = c("chromosome", "start", "end", "copy_number", "major_cn", "minor_cn", "cellular_prevalence")
    }
    purity = read.table(purityfile, header=F, stringsAsFactors=F)
    dat$ccf = dat$cellular_prevalence / purity[1,2]
    return(dat)
  } else {
    return(NA)
  }
}

parse_broad = function(segmentsfile, purityfile, samplename) {
  if (file.exists(segmentsfile)) {
    dat = read.table(segmentsfile, header=T, stringsAsFactors=F)
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

parse_all_profiles = function(samplename, segments, method_segmentsfile, method_purityfile, mustonen_has_header=F) {
  
  dat_dkfz = parse_dkfz(method_segmentsfile[["dkfz"]], method_purityfile[["dkfz"]], samplename)
  if (!is.na(dat_dkfz)) {
    map_dkfz = mapdata(segments, dat_dkfz, is_dkfz=T)
  } else {
    map_dkfz = NA
  }
  
  dat_vanloowedge = parse_vanloowedge(method_segmentsfile[["vanloowedge"]], method_purityfile[["vanloowedge"]], samplename)
  if (!is.na(dat_vanloowedge)) {
    map_vanloowedge = mapdata(segments, dat_vanloowedge)
  } else {
    map_vanloowedge = NA
  }
  
  dat_peifer = parse_peifer(method_segmentsfile[["peifer"]], method_purityfile[["peifer"]], samplename)
  if (!is.na(dat_peifer)) {
    map_peifer = mapdata(segments, dat_peifer)
  } else {
    map_peifer = NA
  }
  
  dat_mustonen = parse_mustonen(method_segmentsfile[["mustonen"]], method_purityfile[["mustonen"]], samplename, has_header=mustonen_has_header)
  if (!is.na(dat_mustonen)) {
    map_mustonen = mapdata(segments, dat_mustonen)
  } else {
    map_mustonen = NA
  }
  
  dat_broad = parse_broad(method_segmentsfile[["broad"]], method_purityfile[["broad"]], samplename)
  if (!is.na(dat_broad)) {
    map_broad = mapdata(segments, dat_broad)
  } else {
    map_broad = NA
  }
  
  return(list(dat_dkfz=dat_dkfz, map_dkfz=map_dkfz,
              dat_vanloowedge=dat_vanloowedge, map_vanloowedge=map_vanloowedge,
              dat_peifer=dat_peifer, map_peifer=map_peifer,
              dat_mustonen=dat_mustonen, map_mustonen=map_mustonen,
              dat_broad=dat_broad, map_broad=map_broad))
}

parse_dummy_cn_profile = function(nmaj=-1, nmin=-1) {
  temp = read.table("segmentation_chrom_arms_full.txt", header=T, stringsAsFactors=F)
  temp$nMaj1_A = nmaj
  temp$nMin1_A = nmin
  return(temp)
}

get_dummy_cn_entry = function(segment) {
  return(data.frame(chromosome=segment$chromosome, start=segment$start, end=segment$end, copy_number=-2, major_cn=-1, minor_cn=-1, ccf=NA))
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
mapdata = function(bp_segments, cn_segments, is_dkfz=F, dkfz_subclonality_cutoff=0.1) {
  bps_gr = makeGRangesFromDataFrame(bp_segments)
  cns_gr = makeGRangesFromDataFrame(cn_segments, keep.extra.columns=T)
  overlap = findOverlaps(cns_gr, bps_gr)
  
  status = rep(NA, length(bps_gr))
  cn_states = list()
  
  # for (i in 1:length(overlap)) {
  for (i in 1:nrow(bp_segments)) {
    overlap = findOverlaps(cns_gr, bps_gr[i,])
    
    # No overlap, no call
    if (length(overlap)==0) {
      next
      
      # One segment overlaps, but could be subclonal in DKFZ output
    } else if (length(overlap)==1) {
      
      if (is_dkfz) {
        status[i] = asses_dkfz_clonal_status(cn_segments, queryHits(overlap), dkfz_subclonality_cutoff)
      } else {
        # clonal
        status[i] = "clonal"
      }
      cn_states[[i]] = list(cn_segments[queryHits(overlap),])
      
      # More than one segment overlaps, this could be subclonal, but it could also be that another segment slightly  
      # overlaps with the breakpoint defined segment and there are not really two calls that overlap
    } else { 
      # get only segments that overlap at least 50%
      overlap = findOverlaps(cns_gr, bps_gr[i,], minoverlap=round((bp_segments$end[i]-bp_segments$start[i])*0.5))
      
      # No segments overlap 50%
      if (length(overlap)==0) {
        overlap = findOverlaps(cns_gr, bps_gr[i,])
        
        # In this case the cn segment is smaller than the breakpoint given segment
        # find the breakpoint segment that overlaps with 95% of the given cn segment and make that mapping
        segment_overlaps_bp = round(cn_segments$end[queryHits(overlap)]-cn_segments$start[queryHits(overlap)])
        
        # Check for zero sized segments
        if (sum(segment_overlaps_bp > 0)==0) {
          # No more candidates left, so no overlap
          next
        }
        
        # Remove segments that are 1bp long
        smallest_segment_overlaps = min(segment_overlaps_bp[segment_overlaps_bp > 0])
        overlap = findOverlaps(cns_gr, bps_gr[i,], minoverlap=round(smallest_segment_overlaps*0.95))
        
        # There are no segments that overlap considerably with the given segment, no call
        if (length(overlap)==0) {
          next
          
          # One segment overlaps considerably, therefore clonal
        } else if (length(overlap)==1) {
          if (is_dkfz) {
            status[i] = asses_dkfz_clonal_status(cn_segments, queryHits(overlap), dkfz_subclonality_cutoff)
          } else {
            # clonal
            status[i] = "clonal"
          }
          
          # Two segments overlap considerably, therefore subclonal
        } else {
          status[i] = "subclonal"
        }
        
        # clonal - because only one segment overlaps considerably
      } else if (length(overlap)==1) {
        if (is_dkfz) {
          status[i] = asses_dkfz_clonal_status(cn_segments, queryHits(overlap), dkfz_subclonality_cutoff)
        } else {
          # clonal
          status[i] = "clonal"
        }
        
        # subclonal - because multiple segments overlap considerably
      } else {
        status[i] = "subclonal"
      }
      cn_states[[i]] = list(cn_segments[queryHits(overlap),])
    }
  }
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
    if (is.null(cn_states_i)) {
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
      # Subclonal copy number
      ancestral = which(cn_states_i$ccf==1)
      descendant = (1:2)[(1:2)!=ancestral]
      
      new_bb_seg$nMaj1_A[1] = cn_states_i$major_cn[ancestral]
      new_bb_seg$nMin1_A[1] = cn_states_i$minor_cn[ancestral]
      new_bb_seg$frac1_A[1] = cn_states_i$ccf[ancestral] - cn_states_i$ccf[descendant]
      new_bb_seg$nMaj2_A[1] = cn_states_i$major_cn[descendant]
      new_bb_seg$nMin2_A[1] = cn_states_i$minor_cn[descendant]
      new_bb_seg$frac2_A[1] = cn_states_i$ccf[descendant]
      cn_bb = rbind(cn_bb, new_bb_seg)
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
    geom_rect(mapping=aes(xmin=startpos, xmax=endpos, ymin=plot_maj, ymax=(plot_maj+RECT_HEIGHT)), fill="purple") +
    geom_rect(mapping=aes(xmin=startpos, xmax=endpos, ymin=plot_min-RECT_HEIGHT, ymax=(plot_min)), fill="#2f4f4f") +
    geom_rect(mapping=aes(xmin=startpos, xmax=endpos, ymin=plot_cn_total, ymax=(plot_cn_total+RECT_HEIGHT)), fill="#E69F00") +
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
round_vanloo_wedge = function(map, i) {
  if (!is.null(map$cn_states[[i]]) && nrow(map$cn_states[[i]][[1]]) > 1) {
    dat = map$cn_states[[i]][[1]]
    index_major_clone = which.max(dat$ccf)
    return(map$cn_states[[i]][[1]][index_major_clone,,drop=F])
  } else if (is.null(map$cn_states[[i]])) {
    return(data.frame())
  } else {
    return(map$cn_states[[i]][[1]])
  }
}

round_peifer = function(map, i) {
  if (!is.null(map$cn_states[[i]]) && nrow(map$cn_states[[i]][[1]]) > 1) {
    dat = map$cn_states[[i]][[1]]
    index_major_clone = which.max(dat$ccf)
    return(map$cn_states[[i]][[1]][index_major_clone,,drop=F])
  } else if (is.null(map$cn_states[[i]])) {
    return(data.frame())
  } else {
    return(map$cn_states[[i]][[1]])
  }
}

round_mustonen = function(map, i) {
  return(map$cn_states[[i]][[1]][1,,drop=F])
}

round_dkfz = function(map, i) {
  if (!is.null(map$cn_states[[i]]) && nrow(map$cn_states[[i]][[1]]) == 1) {
    temp = map$cn_states[[i]][[1]]
    temp$minor_cn = round(temp$minor_cn)
    temp$major_cn = round(temp$major_cn)
    temp$copy_number = temp$minor_cn + temp$major_cn
    return(temp)
  } else if (is.null(map$cn_states[[i]])) {
    return(data.frame())
  } else {
    stop(paste0("DKFZ fit contains multiple states for segment ", i))
  }
}

round_broad = function(map, i) {
  if (is.na(map)) {
    dat = parse_bb_template()[1,,drop=F]
    dat$n
  } else if (!is.null(map$cn_states[[i]]) && nrow(map$cn_states[[i]][[1]]) > 1) {
    dat = map$cn_states[[i]][[1]]
    return(dat[dat$historically_clonal==1,,drop=F])
  } else if (is.null(map$cn_states[[i]])) {
    return(data.frame())
  } else {
    return(map$cn_states[[i]][[1]])
  }
}

get_combined_status = function(segments, map_vanloowedge, map_dkfz, map_mustonen, map_peifer, map_broad) {
  if (is.na(map_vanloowedge)) {
    vanloowedge = rep(NA, nrow(segments))
  } else {
    vanloowedge = map_vanloowedge$status
  }
  
  if (is.na(map_dkfz)) {
    dkfz = rep(NA, nrow(segments))
  } else {
    dkfz = map_dkfz$status
  }
  
  if (is.na(map_mustonen)) {
    mustonen = rep(NA, nrow(segments))
  } else {
    mustonen = map_mustonen$status
  }
  
  if (is.na(map_peifer)) {
    peifer = rep(NA, nrow(segments))
  } else {
    peifer = map_peifer$status
  }
  
  if (is.na(map_broad)) {
    broad = rep(NA, nrow(segments))
  } else {
    broad = map_broad$status
  }
  
  combined_status = data.frame(segments, dkfz=dkfz, mustonen=mustonen, peifer=peifer, vanloowedge=vanloowedge, broad=broad)
  return(combined_status)
}

calc_ploidy = function(subclones) {
  subclones$length = round((subclones$endpos-subclones$startpos)/1000)
  cn_state_one = (subclones$nMaj1_A+subclones$nMin1_A)*subclones$frac1_A
  cn_state_two = ifelse(!is.na(subclones$frac2_A), (subclones$nMaj2_A+subclones$nMin2_A)*subclones$frac2_A, 0)
  ploidy = sum((cn_state_one+cn_state_two) * subclones$length, na.rm=T) / sum(subclones$length, na.rm=T)
  return(ploidy)
}

get_ploidy_status = function(subclones) {
  seg_length = (subclones$endpos/1000)-(subclones$startpos/1000)
  is_clonal = is.na(subclones$frac2_A)
  dipl = subclones$nMin1_A==1 & subclones$nMin1_A==1
  tetrpl = subclones$nMin1_A==2 & subclones$nMin1_A==2
  
  if (sum(seg_length[dipl & is_clonal], na.rm=T) >= sum(seg_length[tetrpl & is_clonal], na.rm=T)) {
    status = "diploid"
  } else {
    status = "tetraploid"
  }
  return(status)
}

