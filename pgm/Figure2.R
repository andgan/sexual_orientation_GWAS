library(gtools)
library(ggplot2)
library(data.table)


mplot <- function(gwas, bp = NA, chrom = NA, pvalue = NA, intervals=NA, ymax=NA, minval=0.01) 
	{
		dfm <- as.data.frame(gwas)
		dfm$chrom <- dfm[,chrom]
		dfm$bp <- as.numeric(as.character(dfm[,bp]))
		dfm$pvalue <- as.numeric(as.character(dfm[,pvalue]))
		dfm$chrom <- as.character(dfm$chrom)
		dfm$marker <- -log10(dfm$pvalue)
		## Remove P-value > 0.001
		dfm <- dfm[dfm$marker > -log10(minval),]
		##add index
		dfm <- dfm[order(dfm$bp),]
		dfm <- dfm[mixedorder(dfm$chrom),]
		dfm$index <- 1:nrow(dfm)
		chrtable <- data.frame(table(dfm$chrom))
		chrtable$Var1 <- as.character(chrtable$Var1)
		chrtable <- chrtable[mixedorder(chrtable$Var1),]
		oddchrom <- as.character(chrtable$Var1[seq(1,nrow(chrtable),2)])
		dfm$chrom_alt <- replace(dfm$chrom, dfm$chrom %in% oddchrom, 0)
		dfm$chrom_alt <- replace(dfm$chrom_alt, dfm$chrom_alt != 0,1)
		dfm$chrom_altA <- ifelse(dfm$chrom_alt=="1",1,2)
		if (is.na(ymax))
		{
		ymax <- max(-log10(dfm$pvalue)) + 1
		}
		dfmsplit <- split(dfm, dfm$chrom)
		xbreaks <- sapply(dfmsplit,function(x) x$index[length(x$index)/2])
		dfm$inint <- 1
		dfm$shape <- 16
		dfm$fill <- 0

		dfm$linepos <- NA
		dfm$linepos[dfm$index[!duplicated(dfm$CHR)]] <- dfm$index[!duplicated(dfm$CHR)]


		for (i in 1:nrow(intervals))
		{
			dfm$inint[dfm$CHR==intervals$chr[i] & dfm$bp<=intervals$end[i] & dfm$bp>=intervals$start[i]] <- 2
			dfm$chrom_altA[dfm$CHR==intervals$chr[i] & dfm$bp<=intervals$end[i] & dfm$bp>=intervals$start[i]] <- 2
			dfm$fill[dfm$CHR==intervals$chr[i] & dfm$bp<=intervals$end[i] & dfm$bp>=intervals$start[i]] <- ifelse(intervals$sex[i]=="M",1,ifelse(intervals$sex[i]=="F",2,ifelse(intervals$sex[i]=="MF",3,NA)))
			dfm$shape[dfm$CHR==intervals$chr[i] & dfm$bp<=intervals$end[i] & dfm$bp>=intervals$start[i]] <- 21
			ind_min = intersect(which(dfm$pvalue == min(dfm$pvalue[dfm$CHR==intervals$chr[i] & dfm$bp<=intervals$end[i] & dfm$bp>=intervals$start[i]])), seq(dfm$pvalue)[dfm$CHR==intervals$chr[i] & dfm$bp<=intervals$end[i] & dfm$bp>=intervals$start[i]])
			dfm$shape[ind_min] <- ifelse(intervals$sex[i]=="M",25,
							   ifelse(intervals$sex[i]=="F",24,
							   ifelse(intervals$sex[i]=="MF",23,NA)))
			dfm$inint[ind_min] <- 3
			dfm <- rbind(dfm,dfm[ind_min,])
		}

		p1 <- ggplot(dfm, aes(x = index,y = marker)) +
	                geom_point(aes(size=as.factor(inint),
	                	shape=shape,
	                	fill=as.factor(fill),
	                	alpha=as.factor(chrom_altA),
	                	colour=as.factor(CHR))) + 
	                scale_shape_identity() +
	                scale_x_continuous(breaks = xbreaks, labels = names(xbreaks), expand=c(0,0)) +
	                scale_y_continuous(expand = c(0,0), limits = c(-log10(minval),ymax), breaks=c(2,4,6,8,10,12),labels=c("2","4","6","8","10","12"))+
	                expand_limits(x = 22.3) +
	                guides(colour = FALSE,alpha=FALSE, size=FALSE, fill=FALSE) +
	                labs(x = "chromosome", y = "-log10(P value)") + 
	                theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "red")) +
	            	geom_hline(aes(yintercept= -log10(0.00000005)),colour = "red", lwd=0.6, linetype = 5)  + 
	            	scale_colour_manual(values = c("red3","red","orange","orange2","orange","yellow","yellow2","yellow","green2","green","green2","green","blue","blue2","blue","blue2","blue","purple","purple3","purple","purple3","purple")) + 
	            	scale_alpha_manual(values=c(1,1)) +
	            	scale_size_manual(values=c(0.4,0.9,3)) + 
	            	scale_fill_manual(values=c("0"="white","1"="blue","2"="pink3","3"="red")) +
	            	geom_segment(aes(x = linepos, y = 2, xend = linepos, yend = 3.5), linetype=1, col="grey34", size=0.3)
}



### !!! Notice that since we are allowed to share only 10000 SNPs from 23andMe the plot for non-heterosexual behaviour is actually not reproducible ####

### Here we show the plot using only UK Biobank data, it will not match the plot in the paper ! ###

### Summary statistics can be downloaded: https://geneticsexorientation.info/summary-statistics/



### PANEL A

nonhetM <- fread("zcat < UKB_nonhetero_cc_M.tsv.gz", header = TRUE, stringsAsFactor=FALSE)
nonhetF <- fread("zcat < UKB_nonhetero_cc_F.tsv.gz", header = TRUE, stringsAsFactor=FALSE)
nonhetMF <- fread("zcat < UKB_nonhetero_cc_MF.tsv.gz", header = TRUE, stringsAsFactor=FALSE)

nonhetFIN <- rbind(nonhetM,nonhetF,nonhetMF)

p_fin <- mplot(gwas=nonhetFIN,pvalue="PVAL",chrom="CHR",bp="BP",intervals=data.frame(chr=c(15,11,12,7),start=c(56958297,59071040,82045652,114940147),end=c(57583301,59247803,82068452,115112776),sex=c("M","M","MF","MF")),ymax=-log10(0.00000000005))
png("fig/figure2_panelA.png", type="cairo", width = 11, height = 5, units = 'in', res = 800)
p_fin
dev.off()




### PANEL B

nhetM <- fread("zcat < UKB_nsp_M.tsv.gz", header = TRUE, stringsAsFactor=FALSE)
nhetF <- fread("zcat < UKB_nsp_F.tsv.gz", header = TRUE, stringsAsFactor=FALSE)
nhetMF <- fread("zcat < UKB_nsp_MF.tsv.gz", header = TRUE, stringsAsFactor=FALSE)

# This file contains the significant intervals
nhet <- read.csv("data/nonhetpartners_sign_beta.csv", stringsAsFactor=F)
nhet$start <- sapply(strsplit(sapply(strsplit(nhet$locus,":"),"[",2),"-"),"[",1)
nhet$end <- sapply(strsplit(sapply(strsplit(nhet$locus,":"),"[",2),"-"),"[",2)

# Depending if sex-specific or not, decide what to plot
nhet$sex <- ifelse(nhet$whichplot=="Combined","MF",ifelse(nhet$whichplot=="Females","F","M"))
df <- data.frame(chr=nhet$chr,start=as.numeric(nhet$start),end=as.numeric(nhet$end),sex=nhet$sex)
df <- df[order(df$sex, decreasing = TRUE),]
df$sex <- as.character(df$sex)

nhetFIN <- rbind(nhetM,nhetF,nhetMF)

p_fin <- mplot(gwas=nhetFIN,pvalue="P_BOLT_LMM_INF",chrom="CHR",bp="BP",intervals=df,ymax=-log10(0.00000000000005))
png("fig/figure2_panelB.png", type="cairo", width = 11, height = 5, units = 'in', res = 800)
p_fin
dev.off()

