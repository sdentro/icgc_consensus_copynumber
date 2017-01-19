
kept = c()
output = data.frame()
infiles = list.files("summary_stats/", full.names=T)
for (infile in infiles) {
	d = readr::read_tsv(infile)

	d = d[, c(1, which(grepl("exclude", colnames(d))), which(grepl("purity", colnames(d))), which(grepl("ploidy", colnames(d))), which(grepl("has_", colnames(d))))]
	if (F) {
	if (ncol(d)!=70) {
d$dkfz_bafConfVanloowedge=NA
d$vanloowedge_bafConfVanloowedge=NA
d$peifer_bafConfVanloowedge=NA
d$mustonen_bafConfVanloowedge=NA
d$broad_bafConfVanloowedge=NA
d$jabba_bafConfVanloowedge=NA
d$dkfz_bafConfBroad=NA
d$vanloowedge_bafConfBroad=NA
d$peifer_bafConfBroad=NA
d$mustonen_bafConfBroad=NA
d$broad_bafConfBroad=NA
d$jabba_bafConfBroad=NA
d$numsegments_bafConf=NA
d = d[, !grepl("purities_tested", colnames(d))]
	}
	}
#	if (ncol(d)!=ncol(output)) {
#		kept = c(kept, infile)
#	} else {

		output = rbind(output, d)
#	}
}
write.table(output, file="summary_stats.txt", row.names=F, sep="\t", quote=F)


print(kept)
