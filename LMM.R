library(ggpubr)
library(ggthemes)
library(ggeffects)
library(data.table)
library(dplyr)
library("car")
library("survival")
library("nlme")
library("MASS")
library("plyr")
library("ggplot2")
# library("ggfortify")
library("lme4")
library("lmerTest")


## MEGA data
Subject <- fread("Sample.txt",data.table = F)
Subject <- Subject[!is.na(Subject$STUDY_NAME),]

Longitude <- fread("MEGACombine-LONGITUDINAL-MOCA.txt",data.table = F)

LongSub<-left_join(Longitude,Subject, by=c("ID","STUDY_NAME","STUDY_SITE"))

LongSub$AGEVISIT <- LongSub$AGE + LongSub$YEARS
LongSub$DURVISIT <- LongSub$YEARS + LongSub$DURBASE


## mocha
sample_mocha <- fread("MEGA_moCha_LOY.tsv",
                      data.table = F)
sample_mocha <- sample_mocha[order(sample_mocha$sample_id,sample_mocha$cf,decreasing = T),]
sample_mocha <- sample_mocha[!duplicated(sample_mocha$sample_id),]

PPMI_type <- fread("PPMI_moCha_LOY.tsv",
                   data.table = F)
sample_mocha <- rbind(sample_mocha,PPMI_type)
colnames(sample_mocha)[1] <- "IID"
sample_mocha$cf <- sample_mocha$cf * 100

sample_type <- sample_mocha
LongSub<-left_join(LongSub,sample_type, by=c("IID"))

LongSub$type <- LongSub$event

LongSub <- LongSub[!is.na(LongSub$ID),]
table(LongSub$type)

LongSub$type<-factor(LongSub$type,
                     levels = c("Male_LOY","Female","Male_nonLOY"))

table(LongSub$type,LongSub$SEX)

LongSub <- LongSub[!is.na(LongSub$type),]
length(unique(LongSub$IID))


## 加入GBA apoe phs
GBA <- fread("MEGACombine-SUBJECT-GBAAPOEPHS.txt",data.table = F)

LongSub <- left_join(LongSub,GBA)
library("bootpredictlme4")


Subject.LMMModel<-unique(LongSub[,c("IID","ID","GBA.APOE4.PHS","AGE","YEARSEDUC","DURBASE","SEX","type")])

Subject.Y=subset(Subject.LMMModel,type=='Male_LOY')
Subject.Y.Count=length(unique(Subject.Y$IID))

Subject.H=subset(Subject.LMMModel,type=='Male_nonLOY')
Subject.H.Count=length(unique(Subject.H$IID))

Subject.F=subset(Subject.LMMModel,type=='Female')
Subject.F.Count=length(unique(Subject.F$IID))

predictdf<- function(df, sequence=5, length=11, CLASS, nameid,statclass,stds="-1,0,1"){
  # df <- Subject.Y
  Temp<-df
  Temp<-subset(Temp, select=c(AGE,YEARSEDUC,DURBASE))
  Temps<-data.frame(t(colMeans(Temp,na.rm = T)))
  
  Temps$type<-CLASS
  Temps$ID<-CLASS
  Temps$DURBASE<-mean(df$DURBASE, na.rm=TRUE)+2*stds*sd(df$DURBASE,na.rm=TRUE)  
  Temps$SEX<-nameid
  Temps$GBA.APOE4.PHS<-mean(df$GBA.APOE4.PHS, na.rm=TRUE)
  Temps<-merge(Temps,seq(0,sequence,length=length))
  names(Temps)[names(Temps) == 'y'] <- 'YEARS'
  # Temp[,"LED"]<-mean(df$LED,na.rm=TRUE)
  Temps[,"STAT"]<-statclass
  return(Temps)
}

Population<-NULL

Population<-rbind(Population,predictdf(Subject.Y,8,17,"Male_LOY","M","mean",0))
Population<-rbind(Population,predictdf(Subject.H,8,17,"Male_nonLOY","M","mean",0))
Population<-rbind(Population,predictdf(Subject.F,8,17,"Female","F","mean",0))

## MMSE 
fm1.MMSE<-lmer(MMSE ~ GBA.APOE4.PHS+YEARSEDUC+SEX+AGE+DURBASE+YEARS+type*YEARS+(1+YEARS|ID)+(1|STUDY_NAME), 
               data=LongSub, REML=FALSE)
summary(fm1.MMSE)

p_MMSE <- predict(fm1.MMSE, newdata=Population, re.form=NA, se.fit=TRUE, nsim=200)
Population$MMSE=predict(fm1.MMSE,newdata=Population,re.form=NA)

Population$MMSE_selo <- Population$MMSE-p_MMSE$se.fit
Population$MMSE_sehi <- Population$MMSE+p_MMSE$se.fit



## MDSIII
fm1.MDSIII<-lmer(MDSIII ~ SEX+AGE+DURBASE+YEARS+type*YEARS+(1+YEARS|ID)+(1|STUDY_NAME), data=LongSub, REML=FALSE)
summary(fm1.MDSIII)

p_MDSIII <- predict(fm1.MDSIII, newdata=Population, re.form=NA, se.fit=TRUE, nsim=200)
Population$MDSIII=predict(fm1.MDSIII,newdata=Population,re.form=NA)

Population$MDSIII_selo <- Population$MDSIII-p_MDSIII$se.fit
Population$MDSIII_sehi <- Population$MDSIII+p_MDSIII$se.fit

names <- c("MMSE","MDSIII")
# 
Text.YCarrier <- list()
Text.HCarrier <- list()
Text.FCarrier <- list()
for (i in 1:length(names)) {
  # i=2
  temp.LMMModel<-unique(LongSub[,c("IID","ID","AGE","YEARSEDUC","DURBASE","SEX","type",names[i])])
  temp.LMMModel <- temp.LMMModel[!is.na(temp.LMMModel[,8]),]
  temp.Y=subset(temp.LMMModel,type=='Male_LOY')
  temp.Y.Count=length(unique(temp.Y$IID))
  
  temp.H=subset(temp.LMMModel,type=='Male_nonLOY')
  temp.H.Count=length(unique(temp.H$IID))
  
  temp.F=subset(temp.LMMModel,type=='Female')
  temp.F.Count=length(unique(temp.F$IID))
  
  Text.YCarrier[[i]]<-paste("male_LOY (N = ",temp.Y.Count,")", sep="")
  names(Text.YCarrier)[i] <- names[i]
  Text.HCarrier[[i]]<-paste("male_nonLOY (N = ",temp.H.Count,")", sep="")
  names(Text.HCarrier)[i] <- names[i]
  Text.FCarrier[[i]]<-paste("Female (N = ",temp.F.Count,")", sep="")
  names(Text.FCarrier)[i] <- names[i]
}

## MMSE plot
# summary(fm1.MMSE)


MMSE.LMMplot<-ggplot()+
  geom_ribbon(data=Population, aes(x=YEARS, y=MMSE, ymin=MMSE_selo,
                                   ymax=MMSE_sehi, fill=type), alpha=0.3)+
  geom_line(data=Population, aes(x=YEARS, y=MMSE, color=type))+theme_bw()+
  xlab("Time in study (years)")+ylab("MMSE")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.border = element_blank(),
        axis.line=element_line())+
  labs(title="MMSE")+
  scale_colour_manual(values = c("#E90D8C","#2CABE1","red"))+
  scale_fill_manual(values = c("#E90D8C","#2CABE1","red"))+
  # theme(legend.position="none")+
  scale_y_continuous(breaks = seq(20,32,1),limits = c(24,30),expand = c(0,0))+
  annotate("text",x=1,y=c(29,29.3,29.6),hjust=0,
           label=c(Text.YCarrier[["MMSE"]],Text.HCarrier[["MMSE"]],Text.FCarrier[["MMSE"]]),
           col=c("#2CABE1","red","#E90D8C"),
           size=4)+
  annotate("text", x = 1 , y = 25,parse=F,size = 4,color = "#FC8D62",
           label = c(paste("male_LOY vs male_nonLOY",
                           round(summary(fm1.MMSE)$coefficients["YEARS:typeMale_nonLOY","Pr(>|t|)"],6),
                           collapse = "")),hjust = 0) +
  annotate("text", x = 1 , y = 25.5,parse=F,size = 4,color = "#78AEEF",
           label = c(paste("male_LOY vs Female",
                           round(summary(fm1.MMSE)$coefficients["YEARS:typeFemale","Pr(>|t|)"],6),
                           collapse = "")),hjust = 0) +
  scale_x_continuous(expand=c(0,0))

MMSE.LMMplot


## MDSIII plot
# summary(fm1.MDSIII)

MDSIII.LMMplot<-ggplot()+
  geom_ribbon(data=Population, aes(x=YEARS, y=MDSIII, ymin=MDSIII_selo,
                                   ymax=MDSIII_sehi, fill=type), alpha=0.3)+
  geom_line(data=Population, aes(x=YEARS, y=MDSIII, color=type))+theme_bw()+
  xlab("Time in study (years)")+ylab("MDSIII")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.border = element_blank(),
        axis.line=element_line())+
  labs(title="MDSIII")+
  scale_colour_manual(values = c("#E90D8C","#2CABE1","red"))+
  scale_fill_manual(values = c("#E90D8C","#2CABE1","red"))+
  # theme(legend.position="none")+
  scale_y_continuous(breaks = seq(25,50,5),limits = c(25,50),expand = c(0,0))+
  annotate("text",x=1,y=c(47,48,49),hjust=0,
           label=c(Text.YCarrier[["MDSIII"]],Text.HCarrier[["MDSIII"]],Text.FCarrier[["MDSIII"]]),
           col=c("#2CABE1","red","#E90D8C"),
           size=4)+
  annotate("text", x = 1 , y = 44,parse=F,size = 4,color = "#FC8D62",
           label = c(paste("male_LOY vs male_nonLOY",
                           round(summary(fm1.MDSIII)$coefficients["YEARS:typeMale_nonLOY","Pr(>|t|)"],6),
                           collapse = "")),hjust = 0) +
  annotate("text", x = 1 , y = 45,parse=F,size = 4,color = "#78AEEF",
           label = c(paste("male_LOY vs Female",
                           round(summary(fm1.MDSIII)$coefficients["YEARS:typeFemale","Pr(>|t|)"],6),
                           collapse = "")),hjust = 0) +
  scale_x_continuous(expand=c(0,0))

MDSIII.LMMplot


p1 <- MMSE.LMMplot +theme(legend.position = "none")
p2 <- MDSIII.LMMplot +theme(legend.position = "none")
legend <- get_legend(MMSE.LMMplot)
cowplot::plot_grid(p1, p2,labels = c('A.', 'B.'))

