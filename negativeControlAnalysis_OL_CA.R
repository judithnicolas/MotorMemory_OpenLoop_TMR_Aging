setwd("D:/Documents/Research/TMR/openLoop/bimanual_open_loop/data/OL_CA/group")
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)
library(readr)

### Equivalent baseline performance between the two movement sequences
# load(file = 'MSLSummary.RData')

#RT

ezANOVA (MSLSummary[MSLSummary$Session ==" TrainingPreNap",], dv = .(Mean), wid = .(Sub),
         within= .(Sequence,Block), detailed=T)

ezANOVA (MSLSummary[MSLSummary$Session ==" TestPreNap",], dv = .(Mean), wid = .(Sub),
         within= .(Sequence,Block), detailed=T)

ezANOVA (MSLSummary[MSLSummary$Block %in% as.factor(c(18:20)),], dv = .(Mean), wid = .(Sub),
         within= .(Sequence,Block), detailed=T)


ggplot(summarySE(MSLSummary[MSLSummary$Block %in% as.factor(c(1:20)) ,], measurevar="Mean", groupvars=c(  "Block", "Sequence "),na.rm=T),
       aes(x=Block, y=Mean, group = Sequence )       ) + 
  geom_point(aes(color=Sequence))+
  geom_line(aes(color=Sequence))+
  geom_ribbon(aes(ymin = Mean-se,
                  ymax = Mean+se,fill=factor(Sequence)), alpha = 0.3)+
  ylab("Reaction time (ms)")+
  xlab("Block of training")+
  # coord_cartesian(xlim=c(1,20))+
  theme_classic() 


# Accuracy
ezANOVA (MSLSummary[MSLSummary$Session ==" TrainingPreNap",], dv = .(Acc), wid = .(Sub),
         within= .(Sequence,Block), detailed=T)

ezANOVA (MSLSummary[MSLSummary$Session ==" TestPreNap",], dv = .(Acc), wid = .(Sub),
         within= .(Sequence,Block), detailed=T)

ezANOVA (MSLSummary[MSLSummary$Block %in% as.factor(c(18:20)),], dv = .(Acc), wid = .(Sub),
         within= .(Sequence,Block), detailed=T)

ggplot(summarySE(MSLSummary[MSLSummary$Block %in% as.factor(c(1:20)) ,], measurevar="Acc", groupvars=c(  "Block", "Sequence "),na.rm=T),
       aes(x=Block, y=Acc, group = Sequence )       ) + 
  geom_point(aes(color=Sequence,shape=Sequence,size=Sequence))+
  geom_line(aes(color=Sequence))+
  scale_shape_manual(values=c(15, 16))+
  scale_size_manual(values=c(2,4))+
  geom_ribbon(aes(ymin = Acc-se,
                  ymax = Acc+se,fill=factor(Sequence)), alpha = 0.3)+
  ylab("Accuracy")+
  xlab("Block of training")+
  coord_cartesian(ylim=c(0.9,1))+
  theme_classic() 


### Equivalent baseline performance between the two conditions 
#RT

ezANOVA (MSLSummary[MSLSummary$Session ==" TrainingPreNap",], dv = .(Mean), wid = .(Sub),
         within= .(Condition,Block), detailed=T)

ezANOVA (MSLSummary[MSLSummary$Session ==" TestPreNap" ,], dv = .(Mean), wid = .(Sub),
         within= .(Condition,Block), detailed=T)

ezANOVA (MSLSummary[MSLSummary$Block %in% as.factor(c(18:20)),], dv = .(Mean), wid = .(Sub),
         within= .(Condition,Block), detailed=T)


ggplot(summarySE(MSLSummary[MSLSummary$Block %in% as.factor(c(1:20)) ,], measurevar="Mean", groupvars=c(  "Block", "Condition "),na.rm=T),
       aes(x=Block, y=Median, group = Condition )       ) + 
  geom_point(aes(color=Condition,shape=Condition,size=Condition))+
  geom_line(aes(color=Condition))+
  scale_shape_manual(values=c(15, 16))+
  scale_size_manual(values=c(2,4))+
  geom_ribbon(aes(ymin = Median-se,
                  ymax = Median+se,fill=factor(Condition)), alpha = 0.3)+
  scale_color_manual(values=c("darkblue","magenta")) +
  scale_fill_manual(values=c("darkblue","magenta")) +
  ylab("Reaction time (ms)")+
  xlab("Block of training")+
  # coord_cartesian(xlim=c(1,20))+
  theme_classic() 


# Accuracy
ezANOVA (MSLSummary[MSLSummary$Session ==" TrainingPreNap",], dv = .(Acc), wid = .(Sub),
         within= .(Condition,Block), detailed=T)
ezANOVA (MSLSummary[MSLSummary$Session ==" TestPreNap",], dv = .(Acc), wid = .(Sub),
         within= .(Condition,Block), detailed=T)

ezANOVA (MSLSummary[MSLSummary$Block %in% as.factor(c(18:20)),], dv = .(Acc), wid = .(Sub),
         within= .(Condition,Block), detailed=T)

ggplot(summarySE(MSLSummary[MSLSummary$Block %in% as.factor(c(1:20)) ,], measurevar="Acc", groupvars=c(  "Block", "Condition "),na.rm=T),
       aes(x=Block, y=Acc, group = Condition )       ) + 
  geom_point(aes(color=Condition,shape=Condition,size=Condition))+
  geom_line(aes(color=Condition))+
  scale_shape_manual(values=c(15, 16))+
  scale_size_manual(values=c(2,4))+
  geom_ribbon(aes(ymin = Acc-se,
                  ymax = Acc+se,fill=factor(Condition)), alpha = 0.3)+
  scale_color_manual(values=c("darkblue","magenta")) +
  scale_fill_manual(values=c("darkblue","magenta")) +
  ylab("Accuracy")+
  xlab("Block of training")+
  # coord_cartesian(xlim=c(1,20))+
  theme_classic() 

### Motor exectution performances
# load(file = 'randomSummary.RData')
learningRate=matrix(,length(allSub)*6,3)
colnames(learningRate)=c('Sub','Task',"DV")
learningRate=as.data.frame(learningRate)
counter = 1
for (idx_sub in 1:length(allSub))
{
  
  tmpMSL = MSLSummary[MSLSummary$Sub==allSub[idx_sub] & 
                        MSLSummary$Block %in% as.factor(c(1:4,17:20, 37:40)),]
  
  meanMSLpre = mean(tmpMSL$Mean[tmpMSL$Session==" TrainingPreNap"])
  meanMSLpost = mean(tmpMSL$Mean[tmpMSL$Session==" TestPreNap"])
  meanMSLpostLate = mean(tmpMSL$Mean[tmpMSL$Session==" TestPostNight"])
  
  
  tmpRandom = randomSummary[randomSummary$Sub==allSub[idx_sub] ,]
  meanRandomPre = mean(tmpRandom$Mean[tmpRandom$Session==" preNap"])
  meanRandomPost = mean(tmpRandom$Mean[tmpRandom$Session==" postNight"])
  
  
  learningRate$Sub[counter]   = allSub[idx_sub]
  learningRate$Task[counter]  = 'MSL_rate'
  learningRate$DV[counter]  =  (meanMSLpre-meanMSLpost)/meanMSLpre*100
  
  learningRate$Sub[counter+1]   = allSub[idx_sub]
  learningRate$Task[counter+1]  = 'MSL_early'
  learningRate$DV[counter+1]  =  meanMSLpre
  
  learningRate$Sub[counter+2]   = allSub[idx_sub]
  learningRate$Task[counter+2]  = 'MSL_late'
  learningRate$DV[counter+2]  =  meanMSLpost
  
  
  learningRate$Sub[counter+3]   = allSub[idx_sub]
  learningRate$Task[counter+3]  = 'Random_rate'
  learningRate$DV[counter+3]  =  (meanRandomPre-meanRandomPost)/meanRandomPre*100
  
  learningRate$Sub[counter+4]   = allSub[idx_sub]
  learningRate$Task[counter+4]  = 'Random_early'
  learningRate$DV[counter+4]  =  meanRandomPre
  
  learningRate$Sub[counter+5]   = allSub[idx_sub]
  learningRate$Task[counter+5]  = 'Radom_late'
  learningRate$DV[counter+5]  =  meanRandomPost
  
  counter = counter +6
  
  
}


learningRate$Sub       = as.factor(learningRate$Sub)
learningRate$Task  = as.factor(learningRate$Task)


t.test(learningRate$DV[learningRate$Task=='MSL_rate'],learningRate$DV[learningRate$Task=='Random'],paired = T)
mean(learningRate$DV[learningRate$Task=='Random']);sd(learningRate$DV[learningRate$Task=='Random'])
mean(learningRate$DV[learningRate$Task=='MSL_late']);sd(learningRate$DV[learningRate$Task=='MSL_late'])
cor.test(learningRate$DV[learningRate$Task=='Random'],learningRate$DV[learningRate$Task=='MSL_late'],paired = T)


## Analysis with young adults
load(file ="D:/Documents/Research/TMR/openLoop/bimanual_open_loop/data/group/OL/learningRate.RData")


learningRateCA = learningRate

learningRateCA$Age = 'old'
learningRateYA$Age = 'young'
learningRateAll = rbind(learningRateCA,learningRateYA)

sorted_Age   = c('young',  'old')
learningRateAll$Age = factor(learningRateAll$Age, levels = sorted_Age)

t.test(learningRateAll$DV [learningRateAll$Task =='MSL_rate' & learningRateAll$Age=='old'],
       learningRateAll$DV [learningRateAll$Task =='MSL_rate' & learningRateAll$Age=='young'])
mean(learningRateAll$Rate [learningRateAll$Task =='MSL_early' & learningRateAll$Age=='young']);sd(learningRateAll$Rate [learningRateAll$Task =='MSL_early' & learningRateAll$Age=='young'])
mean(learningRateAll$Rate [learningRateAll$Task =='MSL_early' & learningRateAll$Age=='old']);sd(learningRateAll$Rate [learningRateAll$Task =='MSL_early' & learningRateAll$Age=='old'])


ggplot(summarySE(learningRateAll[learningRateAll$Task=='MSL_early' | learningRateAll$Task=='MSL_late',], measurevar="DV", 
                 groupvars=c("Sub","Task","Age")), aes( y=DV, fill = Task,x = Task)) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,
               color=c("deepskyblue","indianred","deepskyblue","indianred"), alpha=c(0.7,0.7,1,1))+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Task), position=position_dodge(1),
               color=c("deepskyblue","indianred","deepskyblue","indianred"), alpha=c(0.7,0.7,1,1)) +
  scale_fill_manual(values=c("deepskyblue4","indianred4"))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
  geom_point(aes(fill=Task,group=Sub),size=0.5,shape=21, position = position_dodge(0.2)) +
  # coord_cartesian(ylim=c(-100,50))+
  facet_grid(. ~ Age)+
  theme_classic() 

ggplot(summarySE(learningRateAll[learningRateAll$Task=='MSL_rate' ,],
                 measurevar="DV", groupvars=c("Sub","Age"),na.rm=T),
       aes(x=Age, y=DV,fill=Age )) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,color=c("hotpink","hotpink"),alpha=c(0.5,1))+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Age), 
               position=position_dodge(1),color=c("hotpink","hotpink"),alpha=c(0.5,1) ) +
  scale_fill_manual(values=c("hotpink4",'hotpink4'))+
  theme_classic() 



### PVT
PVT = read.csv(file = 'PVT.csv',header = T,sep = ';');

PVT$DV=PVT$DV*1000
PVT$Session= as.factor(PVT$Session)
tapply(PVT$DV,list(PVT$Session),CI)

ezANOVA (PVT, dv = .(DV), wid=.(Sub), within = .(Session), detailed=T)

### SSS
questionnaires_Summarize_CA <- read_delim("questionnaires Summarize CA.csv", ";", escape_double = FALSE, trim_ws = TRUE)
questionnaires_Summarize_CA$Sub = as.factor(questionnaires_Summarize_CA$Sub)
questionnaires_Summarize_CA$Sexe = as.factor(questionnaires_Summarize_CA$Sexe)
allSub = levels(questionnaires_Summarize_CA$Sub)

SSS=matrix(,length(allSub)*3,3)
colnames(SSS)=c('Sub','Session',"DV")
SSS=as.data.frame(SSS)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_sess in 1: 3)
  {
    SSS$Sub[counter]     = allSub[idx_sub]
    SSS$Session[counter] = idx_sess
    SSS$DV[counter]      = as.numeric(questionnaires_Summarize_CA[questionnaires_Summarize_CA$Sub==allSub[idx_sub],8+idx_sess])

    counter = counter +1
  }
}


SSS$Session=as.factor(SSS$Session)
SSS$Sub=as.factor(SSS$Sub)

ezANOVA (SSS[SSS$Sub!='OL_CA_16',], dv = .(DV), wid=.(Sub), within = .(Session), detailed=T)

tapply(SSS$DV[SSS$Sub!='OL_CA_16'],SSS$Session[SSS$Sub!='OL_CA_16'],CI)



### Learning between sequences

learningRate=matrix(,length(allSub)*length(allSequence),5)
colnames(learningRate)=c('Sub','Sequence', "Condition","Confond","Rate")
learningRate=as.data.frame(learningRate)
counter = 1
for (idx_sub in 1:length(allSub))
{
  
  for (idx_seq in 1:length(allSequence))
  {
    tmp = MSLSummary$Mean[MSLSummary$Sub==allSub[idx_sub] & 
                            MSLSummary$Block %in% as.factor(c(1,20)) &
                            MSLSummary$Sequence==allSequence[idx_seq]]
    
    learningRate$Condition[counter]  = as.character(MSLSummary$Condition[MSLSummary$Sub==allSub[idx_sub] & 
                                                                           MSLSummary$Block %in% as.factor(c(1,16)) &
                                                                           MSLSummary$Sequence==allSequence[idx_seq]][1])
    learningRate$Sub[counter]       = allSub[idx_sub]
    learningRate$Sequence[counter]  = allSequence[idx_seq]
    learningRate$Confond[counter]   = as.character(MSLSummary$Confond[MSLSummary$Sub==allSub[idx_sub] & 
                                                                        MSLSummary$Block %in% as.factor(c(1,20)) &
                                                                        MSLSummary$Sequence==allSequence[idx_seq]][1])
    
    
    #% of change
    learningRate$Rate[counter]      = (tmp[1]-tmp[2])/tmp[1]*100
    
    
    counter = counter +1
    
  }
}


learningRate$Sub       = as.factor(learningRate$Sub)
learningRate$Sequence  = as.factor(learningRate$Sequence)
learningRate$Confond   = as.factor(learningRate$Confond)

boxplot(learningRate$Rate~learningRate$Sequence)
t.test(learningRate$Rate[learningRate$Sequence=='1'],learningRate$Rate[learningRate$Sequence=='2'],paired = T)

boxplot(learningRate$Rate~learningRate$Condition)
t.test(learningRate$Rate[learningRate$Condition==' react'],learningRate$Rate[learningRate$Condition==' notReact'],paired = T)


ezANOVA(learningRate, dv = .(Rate),wid=.(Sub),within = .(Condition),
        between_covariates  = .(Sequence),detailed=T)


limInfSeq1 = mean(learningRate$Rate[learningRate$Sequence=='1'])-3*sd(learningRate$Rate[learningRate$Sequence=='1'])
limSupSeq1 = mean(learningRate$Rate[learningRate$Sequence=='1'])+3*sd(learningRate$Rate[learningRate$Sequence=='1'])

limInfSeq2 = mean(learningRate$Rate[learningRate$Sequence=='2'])-3*sd(learningRate$Rate[learningRate$Sequence=='2'])
limSupSeq2 = mean(learningRate$Rate[learningRate$Sequence=='2'])+3*sd(learningRate$Rate[learningRate$Sequence=='2'])


learningRate$Rate[learningRate$Sequence=='2']<limInfSeq2 & learningRate$Rate[learningRate$Sequence=='2']>limSupSeq2
learningRate$Rate[learningRate$Sequence=='1']<limInfSeq1 & learningRate$Rate[learningRate$Sequence=='1']>limSupSeq1


#misc

