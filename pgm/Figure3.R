## Load libraries ##
library(ggplot2)
library(ggrepel)

## Use this values (calculated from UK Biobank) to convert from observed scale to liability scale
riskliab <- 0.2693821
handliab <- 0.09506235
eversmkliab <- 0.5979527
miliab <- 0.03600239
collegeliab <- 0.3273284
diabliab <- 0.04142383
mainliab <- 0.034

## Load h2 calculated from familial relationships

# Binary
fam <- read.table("data/family_h2.tsv", header=T, stringsAsFactor=F)
fam <- fam[,c("trait","SA","se.SA")]


## Load h2 calculated from SNP data

d <- read.csv("data/snp_h2.csv", header=T, stringsAsFactor=F)
d <- d[d$i %in% fam$trait & d$Category=="univariate",]
d <- d[,c("i","Prop._h2","Prop._h2_std_error")]

## Merge and process
mm <- merge(fam,d,by.x="trait",by.y="i")
mm$liab <- NA
mm$liab[mm$trait=="risk"] <- riskliab
mm$liab[mm$trait=="hand"] <- handliab
mm$liab[mm$trait=="eversmk"] <- eversmkliab
mm$liab[mm$trait=="mi"] <- miliab
mm$liab[mm$trait=="college"] <- collegeliab
mm$liab[mm$trait=="diab"] <- diabliab
mm$liab[mm$trait=="nonhet"] <- mainliab

mm$groupalpha <- factor(ifelse( mm$trait%in% c("nonhet","nhetpartners"),1,0))

mm$yup <- mm$SA + 1.96*mm$se.SA
mm$ydo <- mm$SA - 1.96*mm$se.SA

mm$h2 <- ifelse(!is.na(mm$liab),mm$Prop._h2*mm$liab*(1-mm$liab)/(dnorm(qnorm(mm$liab))^2),mm$Prop._h2)
mm$se.h2 <- ifelse(!is.na(mm$liab),mm$Prop._h2_std_error*mm$liab*(1-mm$liab)/(dnorm(qnorm(mm$liab))^2),mm$Prop._h2_std_error)
mm$xup <- mm$h2 + 1.96*mm$se.h2
mm$xdo <- mm$h2 - 1.96*mm$se.h2

mm$label <- ifelse(mm$trait=="nonhet","Non-heterosexual behaviour",
ifelse(mm$trait=="nhetpartners","N. of opposite-sex partners in heterosexuals",NA))


pdf("fig/figure3.pdf",width=5,height=5)
ggplot(aes(y=SA,x=h2), data=mm) + geom_point(aes(col=groupalpha,alpha=groupalpha),size=3) +
geom_errorbar(aes(ymin=ydo,ymax=yup,col=groupalpha,alpha=groupalpha), width=0) + 
geom_errorbarh(aes(xmin=xdo,xmax=xup,col=groupalpha,alpha=groupalpha),height =0) + 
scale_x_continuous(bquote("SNP-based" ~h^2), limits=c(-0.2,1.1)) +
scale_y_continuous(bquote("Family-based narrow sense"~ h^2), limits=c(-0.2,1.1)) +
geom_smooth( method="lm", se=F, lwd=0.5, col="black", linetype=2) +
geom_abline(aes(intercept=0,slope=1), lwd=0.5) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
scale_colour_manual(guide=F, values=c("black","red")) +
scale_alpha_manual(guide=F, values=c(0.1,1)) +
coord_cartesian(xlim=c(0,0.5), ylim=c(0,1)) +
geom_text_repel(aes(label=label),size=2.8,nudge_y = 0.05,nudge_x = 0.11, alpha=0.7, segment.size=0.2, segment.color="grey48")
dev.off()
