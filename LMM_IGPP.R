
##########  MEGA only male  ######
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
library("bootpredictlme4")


## MEGA data
Subject <- fread("c:/Users/bioin/Desktop/ppn/Project/Sample.txt",data.table = F)
Subject <- Subject[!is.na(Subject$STUDY_NAME),]
# table(Subject$STUDY_NAME)

## mocha
sample_mocha <- fread("c:/Users/bioin/Desktop/ppn/Project/MEGA_moCha_LOY.tsv",
                      data.table = F)
sample_mocha <- sample_mocha[order(sample_mocha$sample_id,sample_mocha$cf,decreasing = T),]
sample_mocha <- sample_mocha[!duplicated(sample_mocha$sample_id),]

PPMI_type <- fread("c:/Users/bioin/Desktop/ppn/Project/PPMI_moCha_LOY.tsv",
                   data.table = F)
sample_mocha <- rbind(sample_mocha,PPMI_type)
colnames(sample_mocha)[1] <- "IID"
sample_mocha$cf <- sample_mocha$cf * 100

length(intersect(sample_mocha$IID,Subject$IID))



##   mLRRY
sample_mLRRY <- fread("c:/Users/bioin/Desktop/ppn/Project/MEGA.madloy_230316_ALL.txt",
                      data.table = F)
sample_mLRRY <- sample_mLRRY[!is.na(sample_mLRRY$mLRR),]
colnames(sample_mLRRY)[1] <- "IID"


## merge
sample_type <- full_join(sample_mLRRY,sample_mocha,by = c("IID"))

Subject <- inner_join(Subject,sample_type,by = c("IID"))
length(unique(Subject$IID))


Longitude <- fread("c:/Users/bioin/Desktop/ppn/Project/MEGACombine-LONGITUDINAL-MOCA.txt",data.table = F)
Longitude <- Longitude[Longitude$STUDY_NAME != "BannerHealth",]

LongSub<-inner_join(Longitude,Subject, by=c("ID","STUDY_NAME","STUDY_SITE"))

LongSub$AGEVISIT <- LongSub$AGE + LongSub$YEARS
LongSub$DURVISIT <- LongSub$YEARS + LongSub$DURBASE

LongSub$type <- "Male_nonLOY"
LongSub$type[LongSub$Fosberg_new99 == "Male_LOY" |
               LongSub$event == "Male_LOY" ]  = "Male_LOY"
LongSub$type[LongSub$SEX == "F"] <- "Female"

LongSub <- LongSub[!is.na(LongSub$ID),]
table(LongSub$type)




### only male
LongSub <- LongSub[LongSub$SEX != "F",]

LongSub$type<-factor(LongSub$type,
                     levels = c("Male_LOY","Male_nonLOY"))

table(LongSub$type,LongSub$SEX)

LongSub <- LongSub[!is.na(LongSub$type),]
length(unique(LongSub$ID))
length(unique(LongSub$IID))

## join GBA apoe prs
GBA <- fread("c:/Users/bioin/Desktop/ppn/Project/MEGACombine-SUBJECT-GBAAPOEPHS.txt",data.table = F)

LongSub <- left_join(LongSub,GBA)


Subject.LMMModel<-unique(LongSub[,c("IID","ID","GBA.APOE4.PHS","AGE","YEARSEDUC","DURBASE","SEX","type")])

Subject.Y=subset(Subject.LMMModel,type=='Male_LOY')
Subject.Y.Count=length(unique(Subject.Y$IID))

Subject.H=subset(Subject.LMMModel,type=='Male_nonLOY')
Subject.H.Count=length(unique(Subject.H$IID))


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


## MMSE
fm1.MMSE<-lmer(MMSE ~ GBA.APOE4.PHS+YEARSEDUC+AGE+DURBASE+YEARS+type*YEARS+(1+YEARS|ID)+(1|STUDY_NAME),
               data=LongSub, REML=FALSE)
summary(fm1.MMSE)


p_MMSE <- predict(fm1.MMSE, newdata=Population, re.form=NA, se.fit=TRUE, nsim=200)
Population$MMSE=predict(fm1.MMSE,newdata=Population,re.form=NA)

Population$MMSE_selo <- Population$MMSE-p_MMSE$se.fit
Population$MMSE_sehi <- Population$MMSE+p_MMSE$se.fit


## MDSII
fm1.MDSII<-lmer(MDSII ~ AGE+DURBASE+YEARS+type*YEARS+(1+YEARS|ID)+(1|STUDY_NAME),
               data=LongSub, REML=FALSE)

summary(fm1.MDSII)

p_MDSII <- predict(fm1.MDSII, newdata=Population, re.form=NA, se.fit=TRUE, nsim=200)
Population$MDSII=predict(fm1.MDSII,newdata=Population,re.form=NA)

Population$MDSII_selo <- Population$MDSII-p_MDSII$se.fit
Population$MDSII_sehi <- Population$MDSII+p_MDSII$se.fit


## MDSIII
fm1.MDSIII<-lmer(MDSIII ~ AGE+DURBASE+YEARS+type*YEARS+(1+YEARS|ID)+(1|STUDY_NAME), data=LongSub, REML=FALSE)
summary(fm1.MDSIII)

p_MDSIII <- predict(fm1.MDSIII, newdata=Population, re.form=NA, se.fit=TRUE, nsim=200)
Population$MDSIII=predict(fm1.MDSIII,newdata=Population,re.form=NA)

Population$MDSIII_selo <- Population$MDSIII-p_MDSIII$se.fit
Population$MDSIII_sehi <- Population$MDSIII+p_MDSIII$se.fit

names <- c("MMSE","MDSII","MDSIII")
Text.YCarrier <- list()
Text.HCarrier <- list()
for (i in 1:length(names)) {
  temp.LMMModel<-unique(LongSub[,c("IID","ID","AGE","YEARSEDUC","DURBASE","SEX","type",names[i])])
  temp.LMMModel <- temp.LMMModel[!is.na(temp.LMMModel[,8]),]
  temp.Y=subset(temp.LMMModel,type=='Male_LOY')
  temp.Y.Count=length(unique(temp.Y$IID))

  temp.H=subset(temp.LMMModel,type=='Male_nonLOY')
  temp.H.Count=length(unique(temp.H$IID))

  
  Text.YCarrier[[i]]<-paste("male_LOY (N = ",temp.Y.Count,")", sep="")
  names(Text.YCarrier)[i] <- names[i]
  Text.HCarrier[[i]]<-paste("male_nonLOY (N = ",temp.H.Count,")", sep="")
  names(Text.HCarrier)[i] <- names[i]
 }

## MMSE plot

MMSE.LMMplot<-ggplot()+
  geom_ribbon(data=Population, aes(x=YEARS, y=MMSE, ymin=MMSE_selo,
                                   ymax=MMSE_sehi, fill=type), alpha=0.3)+
  geom_line(data=Population, aes(x=YEARS, y=MMSE, color=type))+theme_bw()+
  xlab("Time in study (years)")+ylab("MMSE")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.border = element_blank(),
        axis.line=element_line())+
  labs(title="MMSE")+
  scale_colour_manual(values = c("#2CABE1","red"))+
  scale_fill_manual(values = c("#2CABE1","red"))+
  # theme(legend.position="none")+
  scale_y_continuous(breaks = seq(20,32,1),limits = c(24,30),expand = c(0,0))+
  annotate("text",x=1,y=c(29,29.3),hjust=0,
           label=c(Text.YCarrier[["MMSE"]],Text.HCarrier[["MMSE"]]),
           col=c("#2CABE1","red"),
           size=4)+
  annotate("text", x = 1 , y = 25,parse=F,size = 4,color = "#FC8D62",
           label = c(paste("male_LOY vs male_nonLOY",
                           round(summary(fm1.MMSE)$coefficients["YEARS:typeMale_nonLOY","Pr(>|t|)"],6),
                           collapse = "")),hjust = 0) +
  scale_x_continuous(expand=c(0,0))

MMSE.LMMplot


## MDSII plot

MDSII.LMMplot<-ggplot()+
  geom_line(data=Population, aes(x=YEARS, y=MDSII, color=type))+theme_bw()+
  xlab("Time in study (years)")+ylab("MDSII")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.border = element_blank(),
        axis.line=element_line())+
  labs(title="MDSII")+
  scale_colour_manual(values = c("#2CABE1","red"))+
  scale_fill_manual(values = c("#2CABE1","red"))+
  scale_y_continuous(breaks = seq(10,20,2),limits = c(10,20),expand = c(0,0))+
  annotate("text",x=4,y=c(11,12),hjust=0,
           label=c(Text.YCarrier[["MDSII"]],Text.HCarrier[["MDSII"]]),
           col=c("#2CABE1","red"),
           size=4)+
  annotate("text", x = 1 , y = 18,parse=F,size = 4,color = "#FC8D62",
           label = c(paste("male_LOY vs male_nonLOY",
                           round(summary(fm1.MDSII)$coefficients["YEARS:typeMale_nonLOY","Pr(>|t|)"],6),
                           collapse = "")),hjust = 0) +
  scale_x_continuous(expand=c(0,0))

MDSII.LMMplot

## MDSIII plot

MDSIII.LMMplot<-ggplot()+
  geom_ribbon(data=Population, aes(x=YEARS, y=MDSIII, ymin=MDSIII_selo,
                                   ymax=MDSIII_sehi, fill=type), alpha=0.3)+
  geom_line(data=Population, aes(x=YEARS, y=MDSIII, color=type))+theme_bw()+
  xlab("Time in study (years)")+ylab("MDSIII")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.border = element_blank(),
        axis.line=element_line())+
  labs(title="MDSIII")+
  scale_colour_manual(values = c("#2CABE1","red"))+
  scale_fill_manual(values = c("#2CABE1","red"))+
  # theme(legend.position="none")+
  scale_y_continuous(breaks = seq(25,60,5),limits = c(25,55),expand = c(0,0))+
  annotate("text",x=1,y=c(47,48),hjust=0,
           label=c(Text.YCarrier[["MDSIII"]],Text.HCarrier[["MDSIII"]]),
           col=c("#2CABE1","red"),
           size=4)+
  annotate("text", x = 1 , y = 44,parse=F,size = 4,color = "#FC8D62",
           label = c(paste("male_LOY vs male_nonLOY",
                           round(summary(fm1.MDSIII)$coefficients["YEARS:typeMale_nonLOY","Pr(>|t|)"],6),
                           collapse = "")),hjust = 0) +
  scale_x_continuous(expand=c(0,0))

MDSIII.LMMplot


table(Subject$MADthres)

p1 <- MMSE.LMMplot +theme(legend.position = "none")
p3 <- MDSIII.LMMplot +theme(legend.position = "none")
legend <- get_legend(MMSE.LMMplot)

pdf("c:/Users/bioin/Desktop/ppn/Project/LMM.pdf")
print(p1)
print(p3)
dev.off()
cowplot::plot_grid(p1, p3,
                   ncol=2,labels = c('A.', 'B.'))


