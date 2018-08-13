library(ggplot2)


traitmap <- rbind(
c("ADHD_2017","ADHD"),
c("AgeFirstBirth_Female","Age at first birth (Females)"),
c("AgeFirstBirth_Male","Age at first birth (Males)"),
c("openness","Openness to experience"),
c("MDD2018_ex23andMe","Major depressive disorder"),
c("Neuroticism_Full","Neuroticism"),
c("NumberChildrenEverBorn_Female","Number of children (Females)"),
c("NumberChildrenEverBorn_Male","Number of children (Males)"),
c("SCZ2","Schizophrenia"),
c("SWB_Full","Subjective well-being"),
c("age_at_menarche","Age at menarche"),
c("age_at_menopauze","Age at menopause"),
c("alcohol_clarke","Alcohol use"),
c("anorexia","Anorexia"),
c("anxiety","Anxiety"),
c("autism_2017","Autism"),
c("bipolar","Bipolar disorder"),
c("birth_weight","Birth weight"),
c("cannabis","Cannabis use"),
c("height_combined","Height"),
c("loneliness_333k","Loneliness"),
c("nr_sex_partners_hetero_females","Number of sex partners (Females)"),
c("nr_sex_partners_hetero_males","Number of sex partners (Males)"),
c("risk_behavior","Risk behaviour"),
c("self_rated_health","Self-rated health"),
c("smoking_cigs_per_day","Smoking: cigarettes per day"),
c("smoking_ever_vs_never","Smoking: ever smoking"),
c("whr_females","Waist-to-hip ratio (Females)"),
c("whr_males","Waist-to-hip ratio (Males)"))



### PANEL A ###

d <- read.csv("data/rg_nonheterosexual.csv", stringsAsFactor=F)
d$sex <- ifelse(sapply(strsplit(d$p1,"_|\\."),"[",2)=="female","Female","Male")

d$p2 <- sapply(strsplit(d$p2,"\\."),"[",1)

p2_new <- NULL
for (i in d$p2)
{
  p2_new <- c(p2_new,traitmap[traitmap[,1]==i,2])
}

d$p2_new <- p2_new
d_signt2 <- d$p2_new[d$p < (0.05/nrow(d))]

ds <- d[d$p2_new %in% d_signt2,]
ds$ymin <- ds$rg-1.96*ds$se
ds$ymax <- ds$rg+1.96*ds$se

#ds$groupcat <- factor(ds$Category)

ds$newgroup <- ordered(ds$sex,levels=c("Male","Female"))
ds$gsign <- factor(ifelse(ds$p < (0.05/nrow(d)),"*",""))
dsO <- ds[rev(order(ds$p2_new,ds$newgroup)),]
dsO$finmatch <- paste0(dsO$p2_new,"_",dsO$newgroup)

x <- 1:36
dsO$pos <- x[x%%3!=0]



pdf("fig/figure5_panelA.pdf",width=5,height=6)

p1 <- ggplot(dsO, aes(x=pos, y=rg, ymin=ymin, ymax=ymax)) +
  geom_hline(yintercept=0.05/nrow(d), linetype="dashed") + # se vuoi una linea verticale nel plot, se no togli
  geom_point(aes(color=newgroup), size=3) +
  geom_errorbar(width = .2,aes(color=newgroup)) + # questo usa ymin e ymax
  theme_bw() +
  coord_flip() + # mette tutto orizzontale
  scale_x_continuous('',breaks=seq(1.5,34.5,3), labels=unique(dsO$p2_new)) +
  scale_y_continuous('Genetic correlation', limits=c(-0.5,1), breaks=c(seq(-0.5,1,0.2)), labels=c(as.character(round(seq(-0.5,1,0.2),1)))) +
   geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p, digits=2))), hjust=-.25, vjust=-0.4, size=2.5, parse=TRUE) + geom_text(aes(label=gsign), col="yellow", vjust=+0.8, size=6) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))  + scale_colour_manual(values=c("blue","pink3"), guide=FALSE) + expand_limits(y = c(-0.5, 1))

p1

dev.off()





### PANEL B ###

d <- read.csv("data/rg_nr_sex_partners.csv", stringsAsFactor=F)
d$sex <- ifelse(sapply(strsplit(d$p1,"_|\\."),"[",5)=="females","Female","Male")

d$p2 <- sapply(strsplit(d$p2,"\\."),"[",1)

p2_new <- NULL
for (i in d$p2)
{
  p2_new <- c(p2_new,traitmap[traitmap[,1]==i,2])
}

d$p2_new <- p2_new
d_signt2 <- d$p2_new[d$p < (0.05/nrow(d))]

ds <- d[d$p2_new %in% d_signt2,]
ds$ymin <- ds$rg-1.96*ds$se
ds$ymax <- ds$rg+1.96*ds$se

#ds$groupcat <- factor(ds$Category)

ds$newgroup <- ordered(ds$sex,levels=c("Male","Female"))
ds$gsign <- factor(ifelse(ds$p < (0.05/nrow(d)),"*",""))
dsO <- ds[rev(order(ds$p2_new,ds$newgroup)),]
dsO$finmatch <- paste0(dsO$p2_new,"_",dsO$newgroup)

x <- 1:45
dsO$pos <- x[x%%3!=0]


pdf("fig/figure5_panelB.pdf",width=5,height=7)

p1 <- ggplot(dsO, aes(x=pos, y=rg, ymin=ymin, ymax=ymax)) +
  geom_hline(yintercept=0.05/nrow(d), linetype="dashed") + # se vuoi una linea verticale nel plot, se no togli
  geom_point(aes(color=newgroup), size=3) +
  geom_errorbar(width = .2,aes(color=newgroup)) + # questo usa ymin e ymax
  theme_bw() +
  coord_flip() + # mette tutto orizzontale
  scale_x_continuous('',breaks=seq(1.5,44.5,3), labels=unique(dsO$p2_new)) +
  scale_y_continuous('Genetic correlation', limits=c(-0.5,1), breaks=c(seq(-0.5,1,0.2)), labels=c(as.character(round(seq(-0.5,1,0.2),1)))) +
   geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p, digits=2))), hjust=-.25, vjust=-0.4, size=2.5, parse=TRUE) + geom_text(aes(label=gsign), col="yellow", vjust=+0.8, size=6) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))  + scale_colour_manual(values=c("blue","pink3"), guide=FALSE) + expand_limits(y = c(-0.5, 1))

p1

dev.off()


