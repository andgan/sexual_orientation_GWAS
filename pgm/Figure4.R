library(ggplot2)


col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061")))


### PANEL A

dat <- read.table('data/phewas_heatmap_nhb.tsv', header=T, stringsAsFactors=F, sep="\t", comment.char="@")

dat$shortname[dat$shortname=="# incorrect matches"] <- "# errors in memory test"
dat$shortname[dat$shortname=="Balding pattern 4"] <- "Extensive baldness"

pdf('fig/figure4_panelA.pdf', width=6.5, height=4)
lim = 30
ggplot(dat, aes(rsid, reorder(shortname, minpadj), fill=pmin(lim, -log10(abs(pval_sign)))*sign(pval_sign))) + geom_tile() + xlab('') + ylab('') + theme(axis.text.x=element_text(angle=90), strip.text.x = element_text(angle=90), strip.text.y = element_text(angle=0)) + scale_fill_gradientn(colors=col(50), limits=c(-lim,lim), name='signed -log10 p-value') + ggtitle('Non-heterosexual behavior') + geom_text(aes(rsid, shortname, label=ifelse(big_pval_adj < .05, '*', '')), size=2.7) + facet_grid(category ~ pheno2, scales='free', space='free')
dev.off()



### PANEL B

dat <- read.table('data/phewas_heatmap_nsp.tsv', header=T, stringsAsFactors=F, sep="\t")

dat$shortname[dat$shortname=="A levels"] <- "Secondary education"
dat$shortname[dat$shortname=="Balding pattern 3"] <- "Front/scalp baldness"
dat$shortname[dat$shortname=="Balding pattern 4"] <- "Extensive baldness"
dat$shortname[dat$shortname=="Bone mineral density (heel)"] <- "Bone mineral density (both heels)"
dat$shortname[dat$shortname=="FEV1"] <- "Force respiratory volume (1sec)"

lim = 30
pdf('fig/figure4_panelB.pdf', width=10, height=5)
ggplot(dat, aes(rsid, reorder(shortname, minpadj), fill=pmin(lim, -log10(abs(pval_sign)))*sign(pval_sign))) + geom_tile() + xlab('') + ylab('') + theme(
  axis.text.x=element_text(angle=90, size=7), 
  strip.text.y = element_text(angle=0),
  ) + scale_fill_gradientn(colors=col(50), limits=c(-lim,lim), name='signed -log10 p-value') + ggtitle('Number of heterosexual partners') + geom_text(aes(rsid, shortname, label=ifelse(big_pval_adj < .05, '*', '')), size=2.7) + facet_grid(category ~ pheno2, scales='free', space='free')
dev.off()



