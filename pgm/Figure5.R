library(ggplot2)
library(stringi)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(cowplot)


### PANEL A ###

# This is simply a summary of the genetic correlations across the different phenotypes and can be obtained from the supplementary tables in the paper #

### PANEL B ###

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
c("cannabis_ever_2018","Cannabis use"),
c("height_combined","Height"),
c("loneliness","Loneliness"),
c("number_sexual_partners","Number of sex partners"),
c("risk_behavior","Risk behaviour"),
c("self_rated_health","Self-rated health"),
c("smoking_ever_vs_never","Smoking: ever smoking"),
c("whr_females","Waist-to-hip ratio (Females)"),
c("whr_males","Waist-to-hip ratio (Males)"),
c("finger2d4d","2D:4D digit ratio"))


d <- read.csv("data/same_sex_sexual_behavior_vs_same_sex_among_non_heterosexuals.csv", stringsAsFactor=F)

p2_new <- NULL
for (i in d$p2)
{
  p2_new <- c(p2_new,traitmap[traitmap[,1]==i,2])
}

d$p2_new <- p2_new
d$ymin <- d$rg-1.96*d$se
d$ymax <- d$rg+1.96*d$se

d$group <- ifelse(grepl("NHB_partners",d$p1),"All non-heterosexuals","Excluding exclusively same-sex")

a <- d[d$group=="All non-heterosexuals",c("rg","ymin","ymax","p2_new","sex","group")]
a$m_id <- paste0(a$p2_new,"_",a$sex)
b <- d[d$group=="Excluding exclusively same-sex",c("rg","ymin","ymax","p2_new","sex","group")]
b$m_id <- paste0(b$p2_new,"_",b$sex)

ab <- merge(a,b,by="m_id")


pdf("fig/figure5_panelB.pdf",width=6,height=5)
ggplot(ab, aes(x=rg.x, y=rg.y)) +
geom_point(size=2,aes(col=sex.x)) +
geom_errorbarh(aes(xmin = ymin.x,xmax = ymax.x, col=sex.x), alpha=0.5) + 
geom_errorbar(aes(ymin = ymin.y,ymax = ymax.y, col=sex.x), alpha=0.5) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Genetic correlations \n proportion same-sex partners among non-heterosexuals") + xlab("Genetic correlations - same-sex sexual behavior") + scale_colour_manual(values=c("#41BB8A","steelblue4"), guide=FALSE) + xlim(-1,1) + ylim(-1,1) + geom_smooth(method="lm", se=F, col="black", linetype=2) 
dev.off()



### PANEL C ###

## Functions ##
get_upper_tri <- function(cormat){
cormat[lower.tri(cormat)]<- NA
return(cormat)
}

reorder_cormat2 <- function(cormat){
ordnam <- rev(c("sex partners","sexual attr","sexual identity","sexual fantasies","emotional attr","community pref","gender socializing"))
cormat <-cormat[order(match(rownames(cormat),ordnam)), order(match(colnames(cormat),ordnam))]
return(cormat)
}

process_geneti_corr <- function(melted_cormat)
{
  melted_cormat <- rbind(melted_cormat,
                cbind(unique(c(melted_cormat$V1,melted_cormat$V2)),unique(c(melted_cormat$V1,melted_cormat$V2)),NA,1,NA,NA))
  melted_cormat$value <- round(as.numeric(as.character(melted_cormat$V4)),2)
  tt <- dcast(melted_cormat,V1~V2, value.var="value")
  rownames(tt) <- tt$V1
  tt <- tt[,2:ncol(tt)]
  tt <- tt[order(apply(tt,1,function(x){sum(is.na(x))})),order(apply(tt,2,function(x){sum(is.na(x))}),decreasing=T)]
  b <- matrix(0,7,7)
  b[upper.tri(b, diag=FALSE)] <- tt[upper.tri(tt, diag=F)]
  b <- t(b)
  tt[lower.tri(tt, diag=F)] <- b[b!=0]
  cormat <- reorder_cormat2(tt)
  upper_tri <- get_upper_tri(cormat)
  outmatrix <- melt(as.matrix(upper_tri), na.rm = TRUE)
  outmatrix$value <- round(as.numeric(as.character(outmatrix$value)),2)
  return(list(outmatrix,upper_tri))
}


## Read data and process ##
d <- read.csv("data/genetic_correlation_23andme.csv", stringsAsFactor=F)
d <- d[!(grepl("_cc",d$V1) | grepl("_cc",d$V2)),]
d$V4[d$V4 > 1] <- 1 # If genetic correlation > 1 then set to 1


d$V1[grepl("motional_attraction",d$V1)] <- "emotional attr"
d$V2[grepl("motional_attraction",d$V2)] <- "emotional attr"

d$V1[grepl("gender_socialize",d$V1)] <- "gender socializing"
d$V2[grepl("gender_socialize",d$V2)] <- "gender socializing"

d$V1[grepl("sexual_attraction",d$V1)] <- "sexual attr"
d$V2[grepl("sexual_attraction",d$V2)] <- "sexual attr"

d$V1[grepl("sexual_community",d$V1)] <- "community pref"
d$V2[grepl("sexual_community",d$V2)] <- "community pref"

d$V1[grepl("sexual_fantasies",d$V1)] <- "sexual fantasies"
d$V2[grepl("sexual_fantasies",d$V2)] <- "sexual fantasies"

d$V1[grepl("sexual_identity",d$V1)] <- "sexual identity"
d$V2[grepl("sexual_identity",d$V2)] <- "sexual identity"

d$V1[grepl("sexual_partners",d$V1)] <- "sex partners"
d$V2[grepl("sexual_partners",d$V2)] <- "sex partners"


melted_cormat_MF <- d[d$V3=="MF",] # Consider males and females combined

melted_cormat <- process_geneti_corr(melted_cormat_MF)[[1]]
melted_cormat <- melted_cormat[!is.na(melted_cormat$value),]
melted_cormat$Var1 <- ordered(melted_cormat$Var1,levels=colnames(process_geneti_corr(melted_cormat_MF)[[2]]))
melted_cormat$Var2 <- ordered(melted_cormat$Var2,levels=colnames(process_geneti_corr(melted_cormat_MF)[[2]]))

# I set  white the negative correlations
melted_cormat$valuel <- melted_cormat$value
melted_cormat$value[melted_cormat$value < 0] <- 0

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "white", high = "red", limit = c(0,1), space = "Lab", 
    name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1,  hjust = 1))+
 coord_fixed()


gg <- ggheatmap + 
geom_text(aes(Var2, Var1, label = valuel), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.5, 0.7),
  legend.direction = "horizontal")+ guides(fill=FALSE)


pdf("fig/figure5_panelC.pdf", height=5,width=6)
gg
dev.off()
