## Load library
library(ggplot2)

## Define the main outcome (variables are coded as in UK Biobank as f.xxxx.0.0, where xxxx represent the variable number)
df$nonhet <- ifelse(df$f.2159.0.0=="Yes",1,ifelse(df$f.2159.0.0=="No",0,NA))


##### PANEL A ######
df_plot1 <-  aggregate(df$nonhet,list(df$f.34.0.0,df$f.31.0.0), table)
df_plot1 <- df_plot1[!df_plot1$Group.1 %in% c(1934,1936,1937,1970,1971),] # Exclude because too few individuals
df_plot1$ratio <- sapply(df_plot1$x,"[",2) / (sapply(df_plot1$x,"[",1)+sapply(df_plot1$x,"[",2))

pdf("fig/figure1_panelA.pdf",width=6.5,height=4)
ggplot(aes(Group.1,ratio*100, group=Group.2), data=df_plot1) + geom_point(alpha=0.5,aes(col=Group.2))  + xlab("Year of birth") + ylab("% individuals reporting non-heterosexual behaviour") + geom_smooth(method = "auto", aes(col=Group.2),se=T, lwd=2, alpha=0.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))  + scale_colour_manual(guide=F, values=c("pink3","blue"))
dev.off()




###### PANEL B ######

df_plot2 <- df[df$f.2149.0.0 <=100,] # Exclude individuals with more than  100 sexual partners
homno <- ifelse(df_plot2$f.3669.0.0 %in% c(-1,-3) | df_plot2$f.2129.0.0 == "Skip this section" | df_plot2$f.2139.0.0 == -2 ,NA,ifelse(df_plot2$f.2159.0.0=="No" ,0, df_plot2$f.3669.0.0))
lifeno<- ifelse(df_plot2$f.2149.0.0 %in% c(-1,-3) | df_plot2$f.2129.0.0 == "Skip this section" | df_plot2$f.2139.0.0 == -2,NA,df_plot2$f.2149.0.0)

## Ratio between number of same-sex partners over lifetime partners ##
df_plot2$ratio_for_homo <- homno/lifeno
df_plot2$ratio_for_homo[df_plot2$ratio_for_homo>1] <- 1 # There are few individuals that report more same-sex partners than lifetime partners, we assume is an error and assign ratio=1
df_plot2$ratio_for_homoC <- df_plot2$ratio_for_homoC <- ifelse(df_plot2$ratio_for_homo==0,NA,df_plot2$ratio_for_homo)


pdf("fig/figure1_panelB.pdf",width=6.5,height=4)
ggplot(aes(ratio_for_homoC), data=df_plot2) + geom_density(aes(fill=f.31.0.0,col=f.31.0.0),alpha = 0.5) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Density") + scale_fill_manual(guide=F, values=c("pink3","blue")) + scale_colour_manual(guide=F, values=c("pink3","blue")) + scale_x_continuous("% same-sex partners among 1+ lifetime sexual partners" ,breaks=c(0.01,0.25,0.5,0.75,1),label=c("1% \n (mostly \n heterosexual)","25%","50% \n (bisexual)", "75%", "100% \n (mostly \n non-heterosexual)"))
dev.off()




###### PANEL C ######

options (contrasts = rep("contr.treatment", 2))

df_plot3 <- df[df$f.2149.0.0 <=100 & !is.na(df$nonhet) & !df$f.34.0.0 %in% c(1934,1936,1937,1970,1971),] # Exclude individuals with more than 100 sexual partners and birth years with too few individuals


# Defione number of children but exclude individuals still in reproductive age
df_plot3$nchild <- ifelse(df_plot3$f.31.0.0=="Female",df_plot3$f.2734.0.0,df_plot3$f.2405.0.0)
df_plot3$nchild[df_plot3$nchild %in% c(-3,-1)] <- NA
df_plot3$nchild[df_plot3$f.21003.0.0 < 45 & df_plot3$f.31.0.0=="Female"] <- NA
df_plot3$nchild[df_plot3$f.21003.0.0 < 55 & df_plot3$f.31.0.0=="Male"] <- NA


## Ratio between number of same-sex partners over lifetime partners ##
homno <- ifelse(df_plot3$f.3669.0.0 %in% c(-1,-3) | df_plot3$f.2129.0.0 == "Skip this section" | df_plot3$f.2139.0.0 == -2 ,NA,ifelse(df_plot3$f.2159.0.0=="No" ,0, df_plot3$f.3669.0.0))
lifeno<- ifelse(df_plot3$f.2149.0.0 %in% c(-1,-3) | df_plot3$f.2129.0.0 == "Skip this section" | df_plot3$f.2139.0.0 == -2,NA,df_plot3$f.2149.0.0)
df_plot3$ratio_for_homo <- homno/lifeno
df_plot3$ratio_for_homo[df_plot3$ratio_for_homo>1] <- 1 # There are few individuals that report more same-sex partners than lifetime partners, we assume is an error and assign ratio=1
df_plot3$nonhet_ratio <- ifelse(df_plot3$ratio_for_homo==0,"0",
          ifelse(df_plot3$ratio_for_homo==1,"1",as.character(cut(df_plot3$ratio_for_homo,4))))
df_plot3$nonhet_ratio <- ordered(df_plot3$nonhet_ratio,levels=c("0","(-0.001,0.25]","(0.25,0.5]","(0.5,0.75]","(0.75,1]","1"))

# Separate in males and females and run a poisson model to estimate the fertility ratio
df_plot3M <- df_plot3[df_plot3$f.31.0.0=="Male",]
df_plot3F <- df_plot3[df_plot3$f.31.0.0=="Female",]


modF <- summary(glm(nchild ~ nonhet_ratio + f.34.0.0, data=df_plot3F, family="poisson")) # Adjusting for year at birth
modM <- summary(glm(nchild ~ nonhet_ratio + f.34.0.0, data=df_plot3M, family="poisson"))


df_plot3_final <- data.frame(y=c(coef(modF)[2:6,1],coef(modM)[2:6,1]),ymin=c(coef(modF)[2:6,1]-1.96*coef(modF)[2:6,2],coef(modM)[2:6,1]-1.96*coef(modM)[2:6,2]),ymax=c(coef(modF)[2:6,1]+1.96*coef(modF)[2:6,2],coef(modM)[2:6,1]+1.96*coef(modM)[2:6,2]),group=c(rep("Female",5),rep("Male",5)), x=c("(1-25%]","(25-50%]","(50-75%]","(75-100%]","100%"))


pdf("fig/figure1_panelC.pdf",width=6.5,height=4)
ggplot(aes(x=x,y=exp(y),ymin=exp(ymin),ymax=exp(ymax)),data=df_plot3_final) + geom_point() + geom_errorbar(colour="black", width=.1) +
  theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background =element_rect(fill="gray95"),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
      )  + xlab("% same-sex partners among 1+ lifetime sexual partners") + facet_wrap(~group) + scale_y_continuous("Fertility Ratio",breaks=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),labels=c("0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"), limits=c(0.15,1.1)) 
dev.off()




