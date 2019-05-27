## Load library
library(ggplot2)

## Define the main outcome (variables are coded as in UK Biobank as f.xxxx.0.0, where xxxx represent the variable number)
df$nonhet <- ifelse(df$f.2159.0.0=="Yes",1,ifelse(df$f.2159.0.0=="No",0,NA))


##### PANEL A ######
df_plot1 <-  aggregate(df$nonhet,list(df$f.34.0.0,df$f.31.0.0), table)
df_plot1 <- df_plot1[!df_plot1$Group.1 %in% c(1934,1936,1937,1970,1971),] # Exclude because too few individuals
df_plot1$ratio <- sapply(df_plot1$x,"[",2) / (sapply(df_plot1$x,"[",1)+sapply(df_plot1$x,"[",2))

pdf("fig/figure1_panelA.pdf",width=6.5,height=4)
ggplot(aes(Group.1,ratio*100, group=Group.2), data=df_plot1) + geom_point(alpha=0.5,aes(col=Group.2))  + xlab("Year of birth") + ylab("% individuals reporting same-sex sexual behaviour") + geom_smooth(method = "auto", aes(col=Group.2),se=T, lwd=2, alpha=0.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))  + scale_colour_manual(guide=F, values=c("steelblue4","#41BB8A"))
dev.off()




###### PANEL B ######

df_plot2 <- df[!is.na(df$nonhet) & !df$f.34.0.0 %in% c(1934,1936,1970,1971,1937),] # Exclude years with too few individuals
df_plot2 <- df_plot2[df_plot2$f.2149.0.0 <=100,] # Exclude individuals with more than  100 sexual partners

# Define Number of childrens
df_plot2$nchild <- ifelse(df_plot2$f.31.0.0=="Female",df_plot2$f.2734.0.0,df_plot2$f.2405.0.0)
df_plot2$nchild[df_plot2$nchild %in% c(-3,-1)] <- NA
df_plot2$nchild[df_plot2$f.21003.0.0 < 45 & df_plot2$f.31.0.0=="Female"] <- NA
df_plot2$nchild[df_plot2$f.21003.0.0 < 55 & df_plot2$f.31.0.0=="Male"] <- NA

# Define ratio number of same-sex partners over number of lifetime sex partners
homno <- ifelse(df_plot2$f.3669.0.0 %in% c(-1,-3) | df_plot2$f.2129.0.0 == "Skip this section" | df_plot2$f.2139.0.0 == -2 ,NA,ifelse(df_plot2$f.2159.0.0=="No" ,0, df_plot2$f.3669.0.0))
lifeno<- ifelse(df_plot2$f.2149.0.0 %in% c(-1,-3) | df_plot2$f.2129.0.0 == "Skip this section" | df_plot2$f.2139.0.0 == -2,NA,df_plot2$f.2149.0.0)

ratio_for_homo <- homno/lifeno
ratio_for_homo[ratio_for_homo>1] <- 1 # There are few individuals that report more same-sex partners than lifetime partners, we assume is an error and assign ratio=1

df_plot2$nonhet_ratio <- ifelse(ratio_for_homo==0,"0",
          ifelse(ratio_for_homo==1,"1",as.character(cut(ratio_for_homo,4))))
df_plot2$nonhet_ratio <- ordered(df_plot2$nonhet_ratio,levels=c("0","(-0.001,0.25]","(0.25,0.5]","(0.5,0.75]","(0.75,1]","1")) # Cut in quartiles

# Fit model separately for males and females
df_plot2M <- df_plot2[df_plot2$f.31.0.0=="Male",]
df_plot2F <- df_plot2[df_plot2$f.31.0.0=="Female",]

options (contrasts = rep("contr.treatment", 2))

modF <- glm(nchild ~ nonhet_ratio + f.34.0.0, data=df_plot2F, family="poisson")
modM <- glm(nchild ~ nonhet_ratio + f.34.0.0, data=df_plot2M, family="poisson")

dfM <- data.frame(nonhet_ratio=unique(df_plot2M$nonhet_ratio[!is.na(df_plot2M$nonhet_ratio)]),f.34.0.0=1960)
dfM <- dfM[order(dfM$nonhet_ratio),]
dfF <- data.frame(nonhet_ratio=unique(df_plot2F$nonhet_ratio[!is.na(df_plot2F$nonhet_ratio)]),f.34.0.0=1960)
dfF <- dfF[order(dfF$nonhet_ratio),]

pM <- predict(modM,newdata=dfM,type="response", se.fit=T)
pF <- predict(modF,newdata=dfF,type="response", se.fit=T)


df_plot2P <- data.frame(y=c(pF$fit,pM$fit),ymin=c(pF$fit-1.96*pF$se.fit,pM$fit-1.96*pM$se.fit),ymax=c(pF$fit+1.96*pF$se.fit,pM$fit+1.96*pM$se.fit),group=c(rep("Female",6),rep("Male",6)), x=c("0%","(1-25%]","(25-50%]","(50-75%]","(75-100%]","100%"))
df_plot2P$x <- ordered(df_plot2P$x,levels=c("0%","(1-25%]","(25-50%]","(50-75%]","(75-100%]","100%"))


pdf("fig/figure1_panelB.pdf",width=6.5,height=4)
ggplot(aes(x=x,y=y,ymin=ymin,ymax=ymax),data=df_plot2P) + geom_col(aes(fill=group)) + geom_errorbar(colour="black", width=.1) +
  theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background =element_rect(fill="gray95"),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
      )  + xlab("% same-sex partners among 1+ lifetime sexual partners") + facet_wrap(~group) + scale_y_continuous("Average number of children *",breaks=c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8),labels=c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8)) + scale_fill_manual(guide=F, values=c("steelblue4","#41BB8A"))
dev.off()





