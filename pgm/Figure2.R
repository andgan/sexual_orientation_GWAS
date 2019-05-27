library(gtools)
library(ggplot2)
library(data.table)


mplot <- function(gwas, bp = NA, chrom = NA, pvalue = NA, intervals=NA, ymax=NA, minval=0.01) 
	{
		dfm <- as.data.frame(gwas)


		posmin <- tapply(dfm[,bp],dfm[,chrom], min)
    	posmax <- tapply(dfm[,bp],dfm[,chrom], max)
    	posshift <- head(c(0,cumsum(as.numeric(posmax))),-1)

	    names(posshift) <- names(posmin)
	    for (k in unique(dfm[,chrom]))
	    {
	        dfm$pos_new[dfm[,chrom]==k] <-  dfm[,bp][dfm[,chrom]==k] + posshift[names(posshift) == k]
	    }

	  	dfmsplit <- split(dfm, dfm[,chrom])
	    xbreaks <- sapply(dfmsplit,function(x) x$pos_new[length(x$pos_new)/2])

	    dfm$marker <- -log10(dfm[,pvalue])
	    df_manhattan <- dfm[dfm$marker > -log10(minval),]


	    ymax <- ifelse(max(df_manhattan$marker) < -log10(5e-8), -log10(5e-8), max(df_manhattan$marker)) + 0.2


		chrtable <- data.frame(table(df_manhattan[,chrom]))
		chrtable$Var1 <- as.character(chrtable$Var1)
		chrtable <- chrtable[mixedorder(chrtable$Var1),]
		oddchrom <- as.character(chrtable$Var1[seq(1,nrow(chrtable),2)])
		df_manhattan$chrom_alt <- replace(df_manhattan[,chrom], df_manhattan[,chrom] %in% oddchrom, 0)
		df_manhattan$chrom_alt <- replace(df_manhattan$chrom_alt, df_manhattan$chrom_alt != 0,1)
		df_manhattan$chrom_altA <- ifelse(df_manhattan$chrom_alt=="1",1,2)
		if (is.na(ymax)) {ymax <- max(-log10(df_manhattan[,pvalue])) + 1}

		df_manhattan$inint <- 1
		df_manhattan$shape <- 16
		df_manhattan$fill <- 0
		for (i in 1:nrow(intervals))
		{
			df_manhattan$inint[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]] <- 1
			df_manhattan$chrom_altA[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]] <- 1
			df_manhattan$fill[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]] <- ifelse(intervals$sex[i]=="M",1,ifelse(intervals$sex[i]=="F",2,ifelse(intervals$sex[i]=="MF",3,NA)))
			df_manhattan$shape[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]] <- 16

			ind_min = intersect(which(df_manhattan[,pvalue] == min(df_manhattan[,pvalue][df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]])), seq(df_manhattan[,pvalue])[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]])

			df_manhattan$shape[ind_min] <- ifelse(intervals$sex[i]=="M",25,
							   ifelse(intervals$sex[i]=="F",24,
							   ifelse(intervals$sex[i]=="MF",23,NA)))

			print(df_manhattan$shape[ind_min])
			df_manhattan$inint[ind_min] <- 3
			df_manhattan <- rbind(df_manhattan,df_manhattan[ind_min,])

		}
		p1 <- ggplot(df_manhattan, aes(x = pos_new,y = marker)) +
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
	                theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "transparent"),  plot.background = element_rect(fill = "transparent", color = NA),axis.line = element_line(colour = "red")) +
	            	geom_hline(aes(yintercept= -log10(0.00000005)),colour = "red", lwd=0.6, linetype = 5)  + 
	            	scale_colour_manual(values = rep(c("grey67","grey22"),11)) + 
	            	scale_alpha_manual(values=c(1,1)) +
	            	scale_size_manual(values=c(0.4,3)) + 
	            	scale_fill_manual(values=c("0"="white","1"="#41BB8A","2"="steelblue4","3"="red"))
	}


#######################################################

### !!! Notice that since we are allowed to share only 10000 SNPs from 23andMe the plot for non-heterosexual behaviour is actually not reproducible ####

###  We show the plot using only UK Biobank data, it will not match the plot in the paper ! ###

### Summary statistics can be downloaded: https://geneticsexorientation.info/summary-statistics/

###########################################################


nonhetM <- fread("zcat < UKB_nonhetero_cc_M.tsv.gz", header = TRUE, stringsAsFactor=FALSE)
nonhetF <- fread("zcat < UKB_nonhetero_cc_F.tsv.gz", header = TRUE, stringsAsFactor=FALSE)
nonhetMF <- fread("zcat < UKB_nonhetero_cc_MF.tsv.gz", header = TRUE, stringsAsFactor=FALSE)

nonhetFIN <- rbind(nonhetM,nonhetF,nonhetMF)

p_fin <- mplot(gwas=nonhetFIN,pvalue="PVAL",chrom="CHR",bp="BP",intervals=data.frame(chr=c(4,15,11,12,7),start=c(36963942,56999901,59040414,81989337,114940147),end=c(37032454,57583301,59233752,82068452,115314917),sex=c("F","M","M","MF","MF")),ymax=-log10(0.0000000005))

png("fig/figure2.png", type="cairo", width = 11, height = 5, units = 'in', res = 800)
p_fin
dev.off()




