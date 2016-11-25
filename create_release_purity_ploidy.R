#' args = commandArgs(T)
#' 
#' This script was used to create the preliminary consensus. The purity_ploidy stats files did not include the ploidy of the 
#' final consensus profile and the mustonen ploidy was incorrect. So they are added in within this script. The script also
#' creates a series of figures that show how the agreement level increases as more levels are added and about the purities.
#' 
#' The final ploidy calculation is added to the get_agreement pipeline. The mustonen bug needs to still be fixed
#' 
#' 
#' 

library(ggplot2)
library(reshape2)
createPng = function(p, filename, height, width) { png(filename, width=width, height=height); print(p); dev.off(); }

make_plots = F

# setwd("/Users/sd11/Documents/Projects/icgc/consensus_subclonal_copynumber/complete_run/purity_ploidy")
###############################################################################
# Load the data
###############################################################################
dat = read.table("purity_ploidy_table/purity_ploidy.txt", header=T, stringsAsFactors=F)
ploidy = read.table("purity_ploidy_table/ploidy.txt", header=T, stringsAsFactors=F)
dat$ploidy_consensus = ploidy$consensus_ploidy
purity_mustonen = read.table("purity_ploidy_mustonen.txt", header=F, stringsAsFactors=F)
colnames(purity_mustonen) = c("samplename", "purity", "ploidy")

# Add in mustonen correct purities
for (i in 1:nrow(dat)) {
  samplename = dat$samplename[i]
  if (sum(purity_mustonen$samplename==samplename)==1) {
    dat$purity_mustonen[i] = purity_mustonen$purity[purity_mustonen$samplename==samplename]
  } else {
    dat$purity_mustonen[i] = NA
  }
}

jeff_adjustments_rounded_20161101 = read.table("purity_ploidy_table/jeff_adjustments_rounded_20161101.txt", header=T, stringsAsFactors=F)
held_back_corrected = read.table("purity_ploidy_table/samplenames_hold_back_corrected.txt", header=T, stringsAsFactors=F)
held_back_corrected = dat[dat$samplename %in% held_back_corrected$samplename,]

held_back_lt30 = read.table("purity_ploidy_table/samplenames_hold_back_lt30.txt", header=T, stringsAsFactors=F)
held_back_lt30 = dat[dat$samplename %in% held_back_lt30$samplename,]

held_back = rbind(held_back_corrected, held_back_lt30)

dat_orig = dat
dat = dat[!dat$samplename %in% held_back$samplename,]


###############################################################################
# Plot agreement levels
###############################################################################
if (make_plots) {
  dat_s = dat[, c("samplename", "agreement_level_a")]
  dat_s = dat_s[with(dat_s, order(-agreement_level_a)),]
  dat_s$samplename = factor(dat_s$samplename, levels=dat_s$samplename)
  
  p = ggplot(dat_s) + aes_string(x="samplename", y="agreement_level_a") + geom_point()
  createPng(p, "agreement_level_a.png", width=700, height=400)
  
  dat_s = dat[, c("samplename", "agreement_level_b")]
  dat_s = dat_s[with(dat_s, order(-agreement_level_b)),]
  dat_s$samplename = factor(dat_s$samplename, levels=dat_s$samplename)
  
  p = ggplot(dat_s) + aes_string(x="samplename", y="agreement_level_b") + geom_point()
  createPng(p, "agreement_level_b.png", width=700, height=400)
  
  dat_s = dat[, c("samplename", "agreement_level_c")]
  dat_s = dat_s[with(dat_s, order(-agreement_level_c)),]
  dat_s$samplename = factor(dat_s$samplename, levels=dat_s$samplename)
  
  p = ggplot(dat_s) + aes_string(x="samplename", y="agreement_level_c") + geom_point()
  createPng(p, "agreement_level_c.png", width=700, height=400)
  
  dat_s = dat[, c("samplename", "agreement_level_d")]
  dat_s = dat_s[with(dat_s, order(-agreement_level_d)),]
  dat_s$samplename = factor(dat_s$samplename, levels=dat_s$samplename)
  
  p = ggplot(dat_s) + aes_string(x="samplename", y="agreement_level_d") + geom_point()
  createPng(p, "agreement_level_d.png", width=700, height=400)
}

###############################################################################
# Get the purity
###############################################################################
# Samples released
purity = dat[, grepl("purity", colnames(dat))]
colnames(purity) = unlist(lapply(colnames(purity), function(x) gsub("purity_", "", x)))
purity = as.matrix(purity[, sort(colnames(purity))])
exclude = dat[, grepl("exclude", colnames(dat))]
colnames(exclude) = unlist(lapply(colnames(exclude), function(x) gsub("exclude_", "", x)))
exclude = as.matrix(exclude[, sort(colnames(exclude))])

purity[exclude] = NA

# Samples with lt30 agreement
purity_held_back_lt30 = held_back_lt30[, grepl("purity", colnames(held_back_lt30))]
colnames(purity_held_back_lt30) = unlist(lapply(colnames(purity_held_back_lt30), function(x) gsub("purity_", "", x)))
purity_held_back_lt30 = purity_held_back_lt30[,colnames(purity)]
purity_held_back_lt30 = as.matrix(purity_held_back_lt30[, sort(colnames(purity_held_back_lt30))])
exclude = held_back_lt30[, grepl("exclude", colnames(held_back_lt30))]
colnames(exclude) = unlist(lapply(colnames(exclude), function(x) gsub("exclude_", "", x)))
exclude = as.matrix(exclude[, sort(colnames(exclude))])

purity_held_back_lt30[exclude] = NA

# Samples corrected
purity_held_back_corrected = held_back_corrected[, grepl("purity", colnames(held_back_corrected))]
colnames(purity_held_back_corrected) = unlist(lapply(colnames(purity_held_back_corrected), function(x) gsub("purity_", "", x)))
purity_held_back_corrected = purity_held_back_corrected[,colnames(purity)]
purity_held_back_corrected = as.matrix(purity_held_back_corrected[, sort(colnames(purity_held_back_corrected))])
exclude = held_back_corrected[, grepl("exclude", colnames(held_back_corrected))]
colnames(exclude) = unlist(lapply(colnames(exclude), function(x) gsub("exclude_", "", x)))
exclude = as.matrix(exclude[, sort(colnames(exclude))])

purity_held_back_corrected[exclude] = NA

###############################################################################
# Make plots for purity
###############################################################################
if (make_plots) {
  has_clonal_aberration = as.matrix(dat[,grepl("clonal_aberration", colnames(dat))])
  has_clonal_aberration[is.na(has_clonal_aberration)] = F
  
  makePurityPlot = function(selection, purity, dat, figure_name, alt=F) {
    large = as.data.frame(purity[selection,])
    large_med = apply(large[, c(1:5)], 1, median, na.rm=T)
    large$samplename = dat$samplename[selection]
    
    large_med = data.frame(samplename=large$samplename, large_med)
    large_med = large_med[with(large_med, order(-large_med)),]
    
    large_med$samplename = factor(large_med$samplename, levels=large_med$samplename)
    large$samplename = factor(large$samplename, levels=large_med$samplename)
    large.m = melt(large)
    
    if (!alt) {
      p = ggplot(large.m) + aes(x=samplename, y=value) + geom_point() + geom_point(data=large_med, mapping=aes(x=samplename, y=large_med, colour="orange"))
    } else {
      p = ggplot(large.m) + aes(x=samplename, y=value, colour=variable) + geom_point() + geom_point(data=large_med, mapping=aes(x=samplename, y=large_med), colour="black", shape=2)
    }
    createPng(p, figure_name, width=700, height=400)
  }
  
  makePurityPlot(has_clonal_aberration[,3], purity, dat, "purities_large_aberration.png")
  makePurityPlot(has_clonal_aberration[,2], purity, dat, "purities_largish_aberration.png")
  makePurityPlot(has_clonal_aberration[,1], purity, dat, "purities_any_aberration.png")
  
  makePurityPlot(rep(T, nrow(held_back_lt30)), purity_held_back_lt30, held_back_lt30, "purities_any_aberration_held_back_lt30.png")
  makePurityPlot(rep(T, nrow(held_back_corrected)), purity_held_back_corrected, held_back_corrected, "purities_any_aberration_held_back_corrected.png")
  
  makePurityPlot(rep(T, nrow(held_back_lt30)), purity_held_back_lt30, held_back_lt30, "purities_any_aberration_held_back_lt30_2.png", alt=T)
  makePurityPlot(rep(T, nrow(held_back_corrected)), purity_held_back_corrected, held_back_corrected, "purities_any_aberration_held_back_corrected_2.png", alt=T)
}
###############################################################################
# Produce the final output file
###############################################################################

purity_med = apply(purity[, c(1:5)], 1, median, na.rm=T)
output = data.frame(samplename=dat$samplename, purity=purity_med, ploidy=dat$ploidy_consensus)
output$purity = format(output$purity,digits = 2,scientific = FALSE)
output$ploidy = format(output$ploidy,digits = 3,scientific = FALSE)

write.table(output, file="purity_ploidy_table/consnsus_purity_ploidy.txt", sep="\t", quote=F, row.names=F)

colnames(purity) = c("absolute", "aceseq", "clonehd", "sclust", "battenberg")
purity = purity[,c("absolute", "aceseq", "battenberg", "clonehd", "sclust")]
write.table(data.frame(output, purity), file="purity_ploidy_table/consnsus_purity_ploidy_annotated.txt", sep="\t", quote=F, row.names=F)
