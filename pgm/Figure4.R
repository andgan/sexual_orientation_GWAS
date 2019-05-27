library(ggplot2)


traitmap <- rbind(
c("ADHD_2017","ADHD","Mental health"),
c("AgeFirstBirth_Female","Age at first birth (Females)","Reproductive traits"),
c("AgeFirstBirth_Male","Age at first birth (Males)","Reproductive traits"),
c("openness","Openness to experience","Personality traits"),
c("MDD2018_ex23andMe","Major depressive disorder","Mental health"),
c("Neuroticism_Full","Neuroticism","Personality traits"),
c("NumberChildrenEverBorn_Female","Number of children (Females)","Reproductive traits"),
c("NumberChildrenEverBorn_Male","Number of children (Males)","Reproductive traits"),
c("SCZ2","Schizophrenia","Mental health"),
c("SWB_Full","Subjective well-being","Mental health"),
c("age_at_menarche","Age at menarche","Reproductive traits"),
c("age_at_menopauze","Age at menopause","Reproductive traits"),
c("alcohol_clarke","Alcohol use","Risky behaviours"),
c("anorexia","Anorexia","Mental health"),
c("anxiety","Anxiety","Mental health"),
c("autism_2017","Autism","Mental health"),
c("bipolar","Bipolar disorder","Mental health"),
c("birth_weight","Birth weight","Physical traits"),
c("cannabis_ever_2018","Cannabis use","Risky behaviours"),
c("height_combined","Height","Physical traits"),
c("loneliness","Loneliness","Personality traits"),
c("number_sexual_partners","Number of sex partners","Reproductive traits"),
c("risk_behavior","Risk behaviour","Risky behaviours"),
c("self_rated_health","Self-rated health","Mental health"),
c("smoking_ever_vs_never","Smoking: ever smoking","Risky behaviours"),
c("whr_females","Waist-to-hip ratio (Females)","Physical traits"),
c("whr_males","Waist-to-hip ratio (Males)","Physical traits"),
c("finger2d4d","2D:4D digit ratio","Physical traits"))


d <- read.csv("data/same_sex_sexual_behavior.csv", stringsAsFactor=F)

cat_new <- NULL
for (i in d$p2)
{
  cat_new <- c(cat_new,traitmap[traitmap[,1]==i,3])
}

p2_new <- NULL
for (i in d$p2)
{
  p2_new <- c(p2_new,traitmap[traitmap[,1]==i,2])
}

ds <- d
ds$cat_new <- cat_new
ds$p2_new <- p2_new
ds$ymin <- ds$rg-1.96*ds$se
ds$ymax <- ds$rg+1.96*ds$se

ds$newgroup <- ordered(ds$sex,levels=c("Male","Female"))
ds$cat_newgroup <- ordered(ds$cat_new,levels=c("Risky behaviours","Mental health","Personality traits","Reproductive traits","Physical traits"))

tt <- aggregate(ds$rg,list(ds$p2_new),mean)

ds$p2_newgroup <- ordered(ds$p2_new,levels=tt$Group.1[order(tt$x)])
ds$gsign <- factor(ifelse(ds$p < (0.05/nrow(d)),"*",""))
dsO <- ds[rev(order(ds$cat_newgroup,ds$p2_newgroup,ds$newgroup)),]

x <- 1:84
dsO$pos <- x[x%%3!=0]


pdf("fig/figure4.pdf",width=5,height=8)

p1 <- ggplot(dsO, aes(x=pos, y=rg, ymin=ymin, ymax=ymax)) +
  geom_hline(yintercept=0.05/nrow(d), linetype="dashed") + # se vuoi una linea verticale nel plot, se no togli
  geom_point(aes(color=newgroup), size=3) +
  geom_errorbar(width = .2,aes(color=newgroup)) + # questo usa ymin e ymax
  theme_bw() +
  coord_flip() + # mette tutto orizzontale
  scale_x_continuous('',breaks=seq(1.5,84.5,3), labels=unique(dsO$p2_new)) +
  scale_y_continuous('Genetic correlation', limits=c(-0.5,1), breaks=c(seq(-0.5,1,0.2)), labels=c(as.character(round(seq(-0.5,1,0.2),1)))) +
   geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p, digits=2))), hjust=-.25, vjust=-0.4, size=2, parse=TRUE) + geom_text(aes(label=gsign), col="yellow", vjust=+0.8, size=6) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))  + scale_colour_manual(values=c("#41BB8A","steelblue4"), guide=FALSE) + expand_limits(y = c(-0.5, 1))

p1

dev.off()





