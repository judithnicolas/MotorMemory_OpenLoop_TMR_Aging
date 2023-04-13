
setwd("D:/Documents/Research/TMR/openLoop/bimanual_open_loop/data/OL_CA/group")
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)
library(ungeviz)
library(BayesFactor)
library(readr)
# load(file = 'MSLSummary.RData')

questionnaires_Summarize_CA <- read_delim("questionnaires Summarize CA.csv", ";", escape_double = FALSE, trim_ws = TRUE)
questionnaires_Summarize_CA$Sub = as.factor(questionnaires_Summarize_CA$Sub)
questionnaires_Summarize_CA$Sexe = as.factor(questionnaires_Summarize_CA$Sexe)

questionnaires_Summarize_CA = questionnaires_Summarize_CA[questionnaires_Summarize_CA$Sub!='OL_CA_16',]
offLineGain=matrix(,length(allSub)*length(allSequence)*3,8)
colnames(offLineGain)=c(colnames(MSLSummary)[c(1,4,5)],"Confond",'Time','Gain_RT','Gain_Acc','Gain_PI')
offLineGain=as.data.frame(offLineGain)
startBlock = c(18,21,25)
endBlock = c(20,24,28)
# startBlock = c(13,21,25)
# endBlock = c(17,24,28)

counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_seq in 1:length(allSequence))
  {
    
    meanPreRT   = mean(MSLSummary$Mean[MSLSummary$Sub==allSub[idx_sub] & 
                                         MSLSummary$Session==' TestPreNap' &
                                         MSLSummary$Block %in% as.factor(c(startBlock[1]:endBlock[1])) &
                                         MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    meanEarlyRT = mean(MSLSummary$Mean[MSLSummary$Sub==allSub[idx_sub] & 
                                         MSLSummary$Session==' TestPostNap' &
                                         MSLSummary$Block %in% as.factor(c(startBlock[2]:endBlock[2])) &
                                         
                                         MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    meanLateRT  = mean(MSLSummary$Mean[MSLSummary$Sub==allSub[idx_sub] & 
                                         MSLSummary$Session==' TestPostNight' &
                                         MSLSummary$Block %in% as.factor(c(startBlock[3]:endBlock[3])) &
                                         MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    
    meanPreAcc  = mean(MSLSummary$Acc[MSLSummary$Sub==allSub[idx_sub] & 
                                        MSLSummary$Session==' TestPreNap' &
                                        MSLSummary$Block %in% as.factor(c(startBlock[1]:endBlock[1])) &
                                        
                                        MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    meanEarlyAcc= mean(MSLSummary$Acc[MSLSummary$Sub==allSub[idx_sub] & 
                                        MSLSummary$Session==' TestPostNap' &
                                        MSLSummary$Block %in% as.factor(c(startBlock[2]:endBlock[2])) &
                                        
                                        MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    meanLateAcc = mean(MSLSummary$Acc[MSLSummary$Sub==allSub[idx_sub] & 
                                        MSLSummary$Session==' TestPostNight' &
                                        MSLSummary$Block %in% as.factor(c(startBlock[3]:endBlock[3])) &
                                        MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    
    meanPrePI   = mean(MSLSummary$PI[MSLSummary$Sub==allSub[idx_sub] & 
                                       MSLSummary$Session==' TestPreNap' &
                                       MSLSummary$Block %in% as.factor(c(startBlock[1]:endBlock[1])) &
                                       MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    meanEarlyPI = mean(MSLSummary$PI[MSLSummary$Sub==allSub[idx_sub] & 
                                       MSLSummary$Session==' TestPostNap' &
                                       MSLSummary$Block %in% as.factor(c(startBlock[2]:endBlock[2])) &
                                       
                                       MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    meanLatePI  = mean(MSLSummary$PI[MSLSummary$Sub==allSub[idx_sub] & 
                                       MSLSummary$Session==' TestPostNight' &
                                       MSLSummary$Block %in% as.factor(c(startBlock[3]:endBlock[3])) &
                                       MSLSummary$Sequence==allSequence[idx_seq]],na.rm = T)
    
    offLineGain$Sub[c(counter,counter+1,counter+2)]       = allSub[idx_sub]
    offLineGain$Sequence[c(counter,counter+1,counter+2)]  = allSequence[idx_seq]
    offLineGain$Condition[c(counter,counter+1,counter+2)] = as.character(MSLSummary$Condition[MSLSummary$Sub==allSub[idx_sub] & 
                                                                                                MSLSummary$Sequence==allSequence[idx_seq]][1])
    offLineGain$Confond[c(counter,counter+1,counter+2)] = as.character(MSLSummary$Confond[MSLSummary$Sub==allSub[idx_sub] & 
                                                                                            MSLSummary$Sequence==allSequence[idx_seq]][1])
    
    offLineGain$Time[counter]      = 'early'
    offLineGain$Gain_RT[counter]   = ((meanPreRT-meanEarlyRT)/meanPreRT)*100
    offLineGain$Gain_Acc[counter]  = ((meanEarlyAcc-meanPreAcc)/meanPreAcc)*100
    offLineGain$Gain_PI[counter]   = (meanEarlyPI-meanPrePI)/meanPrePI
    
    offLineGain$Time[counter+1] = 'late'
    offLineGain$Gain_RT[counter+1]   = ((meanPreRT-meanLateRT)/meanPreRT)*100
    offLineGain$Gain_Acc[counter+1]  = ((meanLateAcc-meanPreAcc)/meanPreAcc)*100
    offLineGain$Gain_PI[counter+1]   = (meanLatePI-meanPrePI)/meanPrePI
    
    offLineGain$Time[counter+2] = 'both'
    offLineGain$Gain_RT[counter+2]   = mean(c(offLineGain$Gain_RT[counter+1],offLineGain$Gain_RT[counter]))
    offLineGain$Gain_Acc[counter+2]  = mean(c(offLineGain$Gain_Acc[counter+1],offLineGain$Gain_Acc[counter]))
    offLineGain$Gain_PI[counter+2]   = mean(c(offLineGain$Gain_PI[counter+1],offLineGain$Gain_PI[counter]))
    

    offLineGain$Age[counter]         = questionnaires_Summarize_CA$Age[questionnaires_Summarize_CA$Sub==allSub[idx_sub]]
    offLineGain$Age[counter+1]       = questionnaires_Summarize_CA$Age[questionnaires_Summarize_CA$Sub==allSub[idx_sub]]
    offLineGain$Age[counter+2]       = questionnaires_Summarize_CA$Age[questionnaires_Summarize_CA$Sub==allSub[idx_sub]]
    
    
    counter = counter+3
  }
}



offLineGain$Sub       = as.factor(offLineGain$Sub)
offLineGain$Sequence  = as.factor(offLineGain$Sequence)
offLineGain$Condition = as.factor(offLineGain$Condition)
offLineGain$Time      = as.factor(offLineGain$Time)
allTime = levels(offLineGain$Time)
limInfEarly = mean(offLineGain$Gain_RT[offLineGain$Time=='early'])-3*sd(offLineGain$Gain_RT[offLineGain$Time=='early'])
limSupEarly = mean(offLineGain$Gain_RT[offLineGain$Time=='early'])+3*sd(offLineGain$Gain_RT[offLineGain$Time=='early'])

limInfLate = mean(offLineGain$Gain_RT[offLineGain$Time=='late'])-3*sd(offLineGain$Gain_RT[offLineGain$Time=='late'])
limSupLate = mean(offLineGain$Gain_RT[offLineGain$Time=='late'])+3*sd(offLineGain$Gain_RT[offLineGain$Time=='late'])

limInfBoth = mean(offLineGain$Gain_RT[offLineGain$Time=='both'])-3*sd(offLineGain$Gain_RT[offLineGain$Time=='both'])
limSupBoth = mean(offLineGain$Gain_RT[offLineGain$Time=='both'])+3*sd(offLineGain$Gain_RT[offLineGain$Time=='both'])

offLineGain$Gain_RT[offLineGain$Time=='early']<limSupEarly & offLineGain$Gain_RT[offLineGain$Time=='early']>limInfEarly
offLineGain$Gain_RT[offLineGain$Time=='late']<limSupLate & offLineGain$Gain_RT[offLineGain$Time=='late']>limInfLate
offLineGain$Gain_RT[offLineGain$Time=='both']<limSupBoth & offLineGain$Gain_RT[offLineGain$Time=='both']>limInfBoth



ezANOVA(offLineGain[offLineGain$Time!='both' ,] , 
                  dv = .(Gain_RT),
                  wid=.(Sub), 
                  within = .(Time,Condition), 
                  detailed=T)

#Power analysis Condition Condition
tmp = summarySE(offLineGain[offLineGain$Time!='both',], measurevar="Gain_RT", 
                groupvars=c("Sub","Condition"))
#etaSquare = 0.151265
cor.test(tmp$Gain_RT[tmp$Condition!=' react' ],tmp$Gain_RT[tmp$Condition!=' notReact' ])
#cor = 0.67838 

#=> Power = 0.9947708

#Power analysis  Time
tmp = summarySE(offLineGain[offLineGain$Time!='both',], measurevar="Gain_RT", 
                groupvars=c("Sub","Time"))
#etaSquare = 0.1730364
cor.test(tmp$Gain_RT[tmp$Time=='early' ],tmp$Gain_RT[tmp$Time=='late' ])
#cor = 0.72097

#=> Power = 0.9989834

#Power analysis Condition Time
tmp = summarySE(offLineGain[offLineGain$Time!='both',], measurevar="Gain_RT", 
                groupvars=c("Sub","Time","Condition"))
#etaSquare = 0.1730364
cor.test(tmp$Gain_RT[tmp$Time=='early' ],tmp$Gain_RT[tmp$Time=='late' ])
#cor = 0.72097

#=> Power = 0.9999998


tmp = c()
tmpEarly = rowMeans(cbind(offLineGain$Gain_RT[offLineGain$Time=='early' &offLineGain$Condition==' react'],
            offLineGain$Gain_RT[offLineGain$Time=='early' &offLineGain$Condition==' notReact']))
tmpLate = rowMeans(cbind(offLineGain$Gain_RT[offLineGain$Time=='late' &offLineGain$Condition==' react'],
                          offLineGain$Gain_RT[offLineGain$Time=='late' &offLineGain$Condition==' notReact']))
t.test(tmpEarly,tmpLate,paired = T)

tmp = c()
tmpReact = rowMeans(cbind(offLineGain$Gain_RT[offLineGain$Time=='early' &offLineGain$Condition==' react'],
                          offLineGain$Gain_RT[offLineGain$Time=='late' &offLineGain$Condition==' react']))
tmpNotReac = rowMeans(cbind(offLineGain$Gain_RT[offLineGain$Time=='early' &offLineGain$Condition==' notReact'],
                         offLineGain$Gain_RT[offLineGain$Time=='late' &offLineGain$Condition==' notReact']))
t.test(tmpReact,tmpNotReac,paired = T)

tmp = summarySE(offLineGainCA[offLineGainCA$Time!='both',], measurevar="Gain_RT", 
                groupvars=c("Sub"))
t.test(tmp$Gain_RT,alternative = 'greater')

tmp = summarySE(offLineGainYA[offLineGainYA$Time=='both',], measurevar="Gain_RT", 
                groupvars=c("Sub"))
t.test(tmp$Gain_RT,alternative = 'greater')

####### Plot offline gain by block and sequence condition. 
#RT
# Individual data points + line separatly for condition
ggplot(summarySE(offLineGain[offLineGain$Time!='both',], measurevar="Gain_RT", 
                 groupvars=c("Sub","Condition","Time")), aes( y=Gain_RT, fill = Condition,x = Condition)) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,color=c("darkblue","darkviolet","darkblue","darkviolet"))+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),color=c("darkblue","darkviolet","darkblue","darkviolet") ) +
  scale_fill_manual(values=c("#3b3bf5","magenta"))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
  geom_point(aes(fill=Condition,group=Sub),size=0.5,shape=21, position = position_dodge(0.2)) +
  coord_cartesian(ylim=c(-50,50))+
  facet_grid(. ~ Time)+
  theme_classic() 


#Acc
 ezANOVA(offLineGain[offLineGain$Time!='both',], 
                      dv = .(Gain_Acc),
                      wid=.(Sub), 
                      within = .(Time,Condition), 
                      detailed=T)

 ggplot(summarySE(offLineGain[offLineGain$Time!='both',], measurevar="Gain_Acc", 
                  groupvars=c("Sub","Condition","Time")), aes( y=Gain_Acc, fill = Condition,x = Condition)) +
   geom_boxplot(position=position_dodge(1),outlier.shape = NA,color=c("darkblue","darkviolet","darkblue","darkviolet"))+
   stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),color=c("darkblue","darkviolet","darkblue","darkviolet") ) +
   scale_fill_manual(values=c("#AAAAEF","magenta"))+
   geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
   geom_point(aes(fill=Condition,group=Sub),size=0.5,shape=21, position = position_dodge(0.2)) +
   # coord_cartesian(ylim=c(-50,50))+
   facet_grid(. ~ Time)+
   theme_classic() 
 
 #PI
 ezANOVA(offLineGain[offLineGain$Time!='both' & offLineGain$Sub!='OL_CA_21',], 
         dv = .(Gain_PI),
         wid=.(Sub), 
         within = .(Time,Condition), 
         detailed=T)
 
 ggplot(summarySE(offLineGain[offLineGain$Time!='both' & offLineGain$Sub!='OL_CA_21',], measurevar="Gain_PI", 
                  groupvars=c("Sub","Condition","Time")), aes( y=Gain_PI, fill = Condition,x = Condition)) +
   geom_boxplot(position=position_dodge(1),outlier.shape = NA,color=c("darkblue","darkviolet","darkblue","darkviolet"))+
   stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),color=c("darkblue","darkviolet","darkblue","darkviolet") ) +
   scale_fill_manual(values=c("#AAAAEF","magenta"))+
   geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
   geom_point(aes(fill=Condition,group=Sub),size=0.5,shape=21, position = position_dodge(0.2)) +
   coord_cartesian(ylim=c(-1,3))+
   facet_grid(. ~ Time)+
   theme_classic() 

 
 t.test(offLineGain$Gain_Acc[offLineGain$Time=='early' & offLineGain$Condition==' react'], 
       offLineGain$Gain_Acc[offLineGain$Time=='early' & offLineGain$Condition==' notReact'],
       paired = T,alternative = "greater")
t.test(offLineGain$Gain_RT[offLineGain$Time=='late' & offLineGain$Condition==' react'], 
       offLineGain$Gain_RT[offLineGain$Time=='late' & offLineGain$Condition==' notReact'],
       paired = T,alternative = "greater")
p.adjust(c(0.08143, 0.02984),method = 'fdr')


#PI
aovGain_PI = ezANOVA(offLineGain[offLineGain$Time!='both',], 
                      dv = .(Gain_PI),
                      wid=.(Sub), 
                      within = .(Time,Condition), 
                      detailed=T)

ggplot(summarySE(offLineGain[offLineGain$Time!='both',], measurevar="Gain_PI", 
                 groupvars=c("Sub","Condition","Time")), aes(x=Time, y=Gain_PI, fill = Condition)) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,color=c("darkblue","darkviolet","darkblue","darkviolet"))+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),color=c("darkblue","darkblue","darkviolet","darkviolet") ) +
  
  scale_fill_manual(values=c("#AAAAEF","magenta"))+
  ylim(-2.5,2.55)+
  theme_classic() 





## Analysis with young adults
load(file ="D:/Documents/Research/TMR/openLoop/bimanual_open_loop/data/group/OL/offLineGain.RData")


offLineGainCA = offLineGain

offLineGainCA$Age = 'old'
offLineGainYA$Age = 'young'
offLineGainAll = rbind(offLineGainCA,offLineGainYA)

sorted_Age   = c('young',  'old')
offLineGainAll$Age = factor(offLineGainAll$Age, levels = sorted_Age)

ezANOVA(offLineGainAll[offLineGainAll$Time!='both' ,] , 
        dv = .(Gain_RT),
        wid=.(Sub), 
        within = .(Time,Condition), 
        between = .(Age), 
        detailed=T)

ezANOVA(offLineGainAll[offLineGainAll$Time!='both' &  offLineGainAll$Condition!=' react',] , 
        dv = .(Gain_RT),
        wid=.(Sub), 
        within = .(Time), 
        between = .(Age), 
        detailed=T)

ploplot = ezPlot(offLineGainAll[offLineGainAll$Time!='both' ,] , 
                 dv = .(Gain_RT   ),
                 wid=.(Sub), 
                 within = .(Time,Condition), 
                 between = .(Age),
                 x = .(),
                 split = .(Age),
                 col = .(Time))
print(ploplot)

ploplot = ezPlot(offLineGainAll[offLineGainAll$Time!='both' &  offLineGainAll$Condition!=' notReact',] , 
                 dv = .(Gain_RT   ),
                 wid=.(Sub), 
                 within = .(Time), 
                 between = .(Age),
                 x = .(),
                 split = .(Age))
print(ploplot)

#Time * Age Effect
tmp = c()
tmpEarly = rowMeans(cbind(offLineGainCA$Gain_RT[offLineGainCA$Time=='early' & offLineGainCA$Condition==' react'],
                          offLineGainCA$Gain_RT[offLineGainCA$Time=='early' & offLineGainCA$Condition==' notReact']))
tmpLate = rowMeans(cbind(offLineGainCA$Gain_RT[offLineGainCA$Time=='late' & offLineGainCA$Condition==' react'],
                            offLineGainCA$Gain_RT[offLineGainCA$Time=='late' & offLineGainCA$Condition==' notReact']))
timeOld = t.test(tmpEarly,tmpLate,paired = T)
ttestBF(tmpEarly,tmpLate,paired = T)

mean(tmpEarly);sd(tmpEarly)
mean(tmpLate);sd(tmpLate)
cor.test(tmpEarly,tmpLate,paired = T)
# d = 0.4437733

tmp = c()
tmpEarly = rowMeans(cbind(offLineGainYA$Gain_RT[offLineGainYA$Time=='early' & offLineGainYA$Condition==' react'],
                          offLineGainYA$Gain_RT[offLineGainYA$Time=='early' & offLineGainYA$Condition==' notReact']))
tmpLate = rowMeans(cbind(offLineGainYA$Gain_RT[offLineGainYA$Time=='late' & offLineGainYA$Condition==' react'],
                            offLineGainYA$Gain_RT[offLineGainYA$Time=='late' & offLineGainYA$Condition==' notReact']))
timeYoung = t.test(tmpEarly,tmpLate,paired = T)
ttestBF(tmpEarly,tmpLate,paired = T)

mean(tmpEarly);sd(tmpEarly)
mean(tmpLate);sd(tmpLate)
cor.test(tmpEarly,tmpLate,paired = T)
# d = 1.392462

#Condition * Age Effect
tmp = c()
tmpReact = rowMeans(cbind(offLineGainCA$Gain_RT[offLineGainCA$Time=='early' & offLineGainCA$Condition==' react'],
                          offLineGainCA$Gain_RT[offLineGainCA$Time=='late' & offLineGainCA$Condition==' react']))
tmpNotReac = rowMeans(cbind(offLineGainCA$Gain_RT[offLineGainCA$Time=='early' & offLineGainCA$Condition==' notReact'],
                            offLineGainCA$Gain_RT[offLineGainCA$Time=='late' & offLineGainCA$Condition==' notReact']))
condOld = t.test(tmpReact,tmpNotReac,paired = T)
ttestBF(tmpReact,tmpNotReac,paired = T)
mean(tmpReact);sd(tmpReact)
mean(tmpNotReac);sd(tmpNotReac)
cor.test(tmpReact,tmpNotReac,paired = T)
# d = 0.4095611

tmp = c()
tmpReact = rowMeans(cbind(offLineGainYA$Gain_RT[offLineGainYA$Time=='early' & offLineGainYA$Condition==' react'],
                          offLineGainYA$Gain_RT[offLineGainYA$Time=='late' & offLineGainYA$Condition==' react']))
tmpNotReac = rowMeans(cbind(offLineGainYA$Gain_RT[offLineGainYA$Time=='early' & offLineGainYA$Condition==' notReact'],
                            offLineGainYA$Gain_RT[offLineGainYA$Time=='late' & offLineGainYA$Condition==' notReact']))
condYoung = t.test(tmpReact,tmpNotReac,paired = T)
ttestBF(tmpReact,tmpNotReac,paired = T)
mean(tmpReact);sd(tmpReact)
mean(tmpNotReac);sd(tmpNotReac)
cor.test(tmpReact,tmpNotReac,paired = T)
# d = 0.4450764


p.adjust(c(timeOld$p.value,timeYoung$p.value,condOld$p.value,condYoung$p.value),n=4,'fdr')

ggplot(summarySE(offLineGainAll[offLineGainAll$Time!='both',], measurevar="Gain_RT", 
                 groupvars=c("Sub","Condition","Age")), aes( y=Gain_RT, fill = Condition,x = Condition)) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,
               color=c("darkblue","darkviolet","darkblue","darkviolet"), alpha=c(0.2,0.2,1,1))+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),
               color=c("darkblue","darkviolet","darkblue","darkviolet"), alpha=c(0.5,0.5,1,1)) +
  scale_fill_manual(values=c("#AAAAEF","magenta"))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
  geom_point(aes(fill=Condition,group=Sub),size=0.5,shape=21, position = position_dodge(0.2)) +
  coord_cartesian(ylim=c(-100,50))+
  facet_grid(. ~ Age)+
  theme_classic() 



ggplot(summarySE(offLineGainAll[offLineGainAll$Time!='both',], measurevar="Gain_RT", 
                 groupvars=c("Sub","Time","Age")), aes( y=Gain_RT, fill = Time,x = Time)) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,
               color=c("olivedrab4","orangered4","olivedrab4","orangered4"), alpha=c(0.2,0.2,1,1))+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Time), position=position_dodge(1),
               color=c("olivedrab4","orangered4","olivedrab4","orangered4"), alpha=c(0.5,0.5,1,1)) +
  scale_fill_manual(values=c("olivedrab2","orangered3"))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
  geom_point(aes(fill=Time,group=Sub),size=0.5,shape=21, position = position_dodge(0.2)) +
  coord_cartesian(ylim=c(-100,50))+
  facet_grid(. ~ Age)+
  theme_classic() 



### TMR index
allTime = levels(offLineGain$Time)

TMRIndex=matrix(,length(allSub)*length(allSequence)*3,4)
colnames(TMRIndex)=c(colnames(MSLSummary)[c(1)],'Age','Time','index_RT')
TMRIndex=as.data.frame(TMRIndex)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_time in 1:length(allTime))
  {
    
    offlineGainReact    = offLineGain$Gain_RT[offLineGain$Sub==allSub[idx_sub] & offLineGain$Condition==' react' & offLineGain$Time==allTime[idx_time]]
    offlineGainNotReact = offLineGain$Gain_RT[offLineGain$Sub==allSub[idx_sub] & offLineGain$Condition==' notReact' & offLineGain$Time==allTime[idx_time]]
    
    TMRIndex$Sub[counter]       = allSub[idx_sub]
    TMRIndex$Time[counter]      = allTime[idx_time]
    TMRIndex$index_RT[counter]  = (offlineGainReact-offlineGainNotReact)
    
    TMRIndex$Age[counter]       = questionnaires_Summarize_CA$Age[questionnaires_Summarize_CA$Sub==allSub[idx_sub]]
    counter = counter+1
  }
  
  
}


TMRIndex$Sub       = as.factor(TMRIndex$Sub)
TMRIndex$Time      = as.factor(TMRIndex$Time)

TMRIndex=na.omit(TMRIndex)



limInfEarly = mean(TMRIndex$index_RT[TMRIndex$Time=='early'])-3*sd(TMRIndex$index_RT[TMRIndex$Time=='early'])
limSupEarly = mean(TMRIndex$index_RT[TMRIndex$Time=='early'])+3*sd(TMRIndex$index_RT[TMRIndex$Time=='early'])

limInfLate = mean(TMRIndex$index_RT[TMRIndex$Time=='late'])-3*sd(TMRIndex$index_RT[TMRIndex$Time=='late'])
limSupLate = mean(TMRIndex$index_RT[TMRIndex$Time=='late'])+3*sd(TMRIndex$index_RT[TMRIndex$Time=='late'])

limInfBoth = mean(TMRIndex$index_RT[TMRIndex$Time=='both'])-3*sd(TMRIndex$index_RT[TMRIndex$Time=='both'])
limSupBoth = mean(TMRIndex$index_RT[TMRIndex$Time=='both'])+3*sd(TMRIndex$index_RT[TMRIndex$Time=='both'])

TMRIndex$index_RT[TMRIndex$Time=='early']<limSupEarly & TMRIndex$index_RT[TMRIndex$Time=='early']>limInfEarly
TMRIndex$index_RT[TMRIndex$Time=='late']<limSupLate & TMRIndex$index_RT[TMRIndex$Time=='late']>limInfLate
TMRIndex$index_RT[TMRIndex$Time=='both']<limSupBoth & TMRIndex$index_RT[TMRIndex$Time=='both']>limInfBoth




ggplot(TMRIndex[TMRIndex$Time=='both',], aes( x=Age,y=index_RT)) +
  geom_point(color=c("black")) +
  geom_smooth(method = lm,color="magenta",alpha=0.2,fill='magenta') +
  # coord_cartesian(ylim = c(-0.4, 0.4 ))+
  stat_cor(method = "pearson")+
  theme_light()


shapiro.test(TMRIndex$index_RT[TMRIndex$Time=='both'])
shapiro.test(TMRIndex$Age[TMRIndex$Time=='both'])


t.test(TMRIndex$index_RT[TMRIndex$Time=='early' ],alternative = 'greater')
t.test(TMRIndex$index_RT[TMRIndex$Time=='late'],alternative = 'greater')
t.test(TMRIndex$index_RT[TMRIndex$Time=='both'],alternative = 'greater')


tmp = summarySE(generationSummary[generationSummary$Session== ' GenerationPre',],
                measurevar="Accuracy", groupvars=c("Sub"))
cor.test(tmp$Accuracy,
         TMRIndex$index_RT[TMRIndex$Time=='both'],paired =TRUE)

## Analysis with young adults
load(file ="D:/Documents/Research/TMR/openLoop/bimanual_open_loop/data/group/OL/TMRIndex.RData")


TMRIndexCA = TMRIndex

TMRIndexCA$Age = 'old'
TMRIndexYA$Age = 'young'
TMRIndexAll = rbind(TMRIndexCA,TMRIndexYA)

sorted_Age   = c('young',  'old')
TMRIndexAll$Age = factor(TMRIndexAll$Age, levels = sorted_Age)

ezANOVA(TMRIndexAll[TMRIndexAll$Time!='both' ,] , 
        dv = .(index_RT   ),
        wid=.(Sub), 
        within = .(Time), 
        between = .(Age), 
        detailed=T)

ploplot = ezPlot(TMRIndexAll[TMRIndexAll$Time!='both' ,]  , 
                 dv = .(index_RT   ),
                 wid=.(Sub), 
                 within = .(Time), 
                 between = .(Age),
                 x = .(Age),
                 split = .(Time))
print(ploplot)

ploplot = ezPlot(offLineGainAll[offLineGainAll$Time!='both' ,]  , 
                 dv = .(Gain_RT   ),
                 wid=.(Sub), 
                 within = .(Time,Condition), 
                 between = .(Age),
                 x = .(Time),
                 split = .(),
                 col = .(Age))
print(ploplot)

#Relation Age ofline gain / TMR index

t.test(offLineGain$Age[offLineGain$Gain_RT>median(offLineGain$Gain_RT[offLineGain$Time=='both' & offLineGain$Condition==' notReact']) & offLineGain$Time=='both' & offLineGain$Condition==' notReact'],
       offLineGain$Age[offLineGain$Gain_RT<=median(offLineGain$Gain_RT[offLineGain$Time=='both' & offLineGain$Condition==' notReact']) & offLineGain$Time=='both' & offLineGain$Condition==' notReact'])
cor.test(offLineGain$Age[offLineGain$Time=='both' & offLineGain$Condition==' notReact'],offLineGain$Gain_RT[offLineGain$Time=='both' & offLineGain$Condition==' notReact'])

t.test(offLineGain$Age[offLineGain$Gain_RT>=median(offLineGain$Gain_RT[offLineGain$Time=='both' & offLineGain$Condition==' react']) & offLineGain$Time=='both' & offLineGain$Condition==' react'],
       offLineGain$Age[offLineGain$Gain_RT<median(offLineGain$Gain_RT[offLineGain$Time=='both' & offLineGain$Condition==' react']) & offLineGain$Time=='both' & offLineGain$Condition==' react'])
cor.test(offLineGain$Age[offLineGain$Time=='both' & offLineGain$Condition==' react'],offLineGain$Gain_RT[offLineGain$Time=='both' & offLineGain$Condition==' react'])


t.test(TMRIndex$Age[TMRIndex$index_RT<median(TMRIndex$index_RT[TMRIndex$Time=='both']) & TMRIndex$Time=='both'],
       TMRIndex$Age[TMRIndex$index_RT>=median(TMRIndex$index_RT[TMRIndex$Time=='both']) & TMRIndex$Time=='both'])

cor.test(TMRIndex$Age[TMRIndex$Time=='both'] , TMRIndex$index_RT[TMRIndex$Time=='both'] )
ggplot(TMRIndex[TMRIndex$Time=='both',], aes( x=Age,y=index_RT)) +
  geom_point(color=c("black")) +
  geom_smooth(method = lm,color="grey",alpha=0.2,fill='grey') +
  stat_cor(method = "pearson")+
  theme_light()

ggplot(offLineGain[offLineGain$Time=='both' & offLineGain$Condition==' react',], aes( x=Age,y=Gain_RT)) +
  geom_point(color=c("magenta")) +
  geom_smooth(method = lm,color="darkviolet",alpha=0.2,fill='darkviolet') +
  stat_cor(method = "pearson")+
  coord_cartesian(ylim = c(-30,40))+
  theme_light()

ggplot(offLineGain[offLineGain$Time=='both' & offLineGain$Condition==' notReact',], aes( x=Age,y=Gain_RT)) +
  geom_point(color=c("darkblue")) +
  geom_smooth(method = lm,color="darkblue",alpha=0.2,fill='darkblue') +
  stat_cor(method = "pearson")+
  coord_cartesian(ylim = c(-30,40))+
  theme_light()
