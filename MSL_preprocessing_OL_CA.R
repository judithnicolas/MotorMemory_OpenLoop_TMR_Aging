setwd("D:/Documents/Research/TMR/openLoop/bimanual_open_loop/data/OL_CA/group")
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)
MSLData = read.csv(file = 'sequentialSRTT.csv',header = T,sep = ';')

MSLData$Acc=MSLData$Cue==MSLData$Rep
MSLData$Acc = as.numeric(MSLData$Acc)
MSLData$Sequence   = as.factor(MSLData$Sequence)
MSLData$Block      = as.factor(MSLData$Block)
MSLData$OrdinalPos = as.factor(MSLData$OrdinalPos)
MSLData$Cue        = as.factor(MSLData$Cue)
MSLData$Rep        = as.factor(MSLData$Rep)
MSLData$Repetition = as.factor(MSLData$Repetition)

sorted_Block  = paste(sort(as.integer(levels(MSLData$Block))))
MSLData$Block = factor(MSLData$Block, levels = sorted_Block)

sorted_Session  = c(' TrainingPreNap',  ' TestPreNap',  ' TestPostNap',  ' TestPostNight')
MSLData$Session = factor(MSLData$Session, levels = sorted_Session)

allSess       = levels(MSLData$Session)
allCond       = levels(MSLData$Condition)
allSequence   = levels(MSLData$Sequence)
allBlock      = levels(MSLData$Block)
allOrdinalPos = levels(MSLData$OrdinalPos)
allCue        = levels(MSLData$Cue)
allRep        = levels(MSLData$Rep)
allSub        = levels(MSLData$Sub)

for (idx in 1:length(MSLData$Condition))
{
  if (MSLData$Condition[idx]==' react' & MSLData$Sequence[idx]=='1')
  {
    MSLData$Confond[idx] = 'reactivated_1'
  }
  if (MSLData$Condition[idx]==' react' & MSLData$Sequence[idx]=='2')
  {
    MSLData$Confond[idx] = 'reactivated_2'
  }
  
  if (MSLData$Condition[idx]==' notReact' & MSLData$Sequence[idx]=='1')
  {
    MSLData$Confond[idx] = 'reactivated_2'
  }
  if (MSLData$Condition[idx]==' notReact' & MSLData$Sequence[idx]=='2')
  {
    MSLData$Confond[idx] = 'reactivated_1'
  }
  
}


Seq1 = as.factor(c(1,	6,	3,	5,	4,	8,	2,	7))
Seq2 = as.factor(c(7,	2,	6,	4,	5,	1,	8,	3))
blockOfInterst = as.factor(c(17:28))


####### Compute Accuracy by Block, mean rt per key press PI 

MSLSummary=matrix(,length(allSub)*length(allBlock)*length(allSequence),10)
colnames(MSLSummary)=c(colnames(MSLData)[c(1:5,8)],"Mean", "Median", "Acc","Confond")
MSLSummary=as.data.frame(MSLSummary)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_block in 1 : length(allBlock))
  {
    
    for (idx_seq in 1:length(allSequence))
    {
      
     
      
       tmp = MSLData[MSLData$Sub==allSub[idx_sub] &
                      MSLData$Block==allBlock[idx_block] &
                      MSLData$Sequence==allSequence[idx_seq] ,]
      limInf = mean(tmp$RT[tmp$Acc=='1' ])-3*sd(tmp$RT[tmp$Acc=='1' ])
      limSup = mean(tmp$RT[tmp$Acc=='1' ])+3*sd(tmp$RT[tmp$Acc=='1' ])
      
      
      MSLSummary$Sub[counter]       = allSub[idx_sub]
      MSLSummary$Block[counter]     = allBlock[idx_block]
      MSLSummary$Sequence[counter]  = allSequence[idx_seq]
      MSLSummary$Session[counter]   = as.character(tmp$Session[1])
      MSLSummary$Condition[counter] = as.character(tmp$Condition[1])
      MSLSummary$Confond[counter] = as.character(tmp$Confond[1])
      MSLSummary$SeqinBlock[counter] = as.character(tmp$SeqinBlock[1])
      
      
      # % Accuracy per block
      MSLSummary$Acc[counter]       = sum(tmp$Acc)/length(tmp$Acc)
      MSLData$PerCorr[MSLData$Sub==allSub[idx_sub] & 
                        MSLData$Block==allBlock[idx_block] &
                        MSLData$Sequence==allSequence[idx_seq]] = repmat(sum(tmp$Acc)/length(tmp$Acc),length(tmp$Sequence),1)
      
      #Mean RT per key presses
      MSLSummary$Mean[counter]      = mean(tmp$RT[tmp$Acc=='1' & tmp$RT>limInf & tmp$RT<limSup ])
      MSLSummary$Median[counter]    = median(tmp$RT[tmp$Acc=='1' & tmp$RT>limInf & tmp$RT<limSup ])
      MSLData$Outlier[MSLData$Sub==allSub[idx_sub] & 
                        MSLData$Block==allBlock[idx_block] &
                        MSLData$Sequence==allSequence[idx_seq]] = tmp$Acc=='1' & tmp$RT>limInf & tmp$RT<limSup 
      
      

      seqDuration = mean(tmp$timeRep[tmp$OrdinalPos==8] - tmp$timeRep[tmp$OrdinalPos==1] )
      counterOrdPos = 1
      nbCorrectSeq = 0
      for (idx_rep in 1: (length(tmp$Acc)/length(Seq2)))
      {
        if (sum(tmp$Acc[ seq(counterOrdPos,counterOrdPos +length(Seq2)-1,1)])==8)
        {
          nbCorrectSeq = nbCorrectSeq+1
        }
        counterOrdPos= counterOrdPos+length(Seq2)
        
      }
      
      B = ((length(tmp$Acc)/length(Seq2))-nbCorrectSeq)/(length(tmp$Acc)/length(Seq2))
      MSLSummary$PI[counter]       = exp(-seqDuration)*exp(-B)*100
      counter = counter+1
      

    }
  }
}


MSLSummary$Sub       = as.factor(MSLSummary$Sub)
MSLSummary$Sequence  = as.factor(MSLSummary$Sequence)
MSLSummary$Session   = as.factor(MSLSummary$Session)
MSLSummary$Condition = as.factor(MSLSummary$Condition)
MSLSummary$Block     = as.factor(MSLSummary$Block)
MSLSummary$Session   = as.factor(MSLSummary$Session)
MSLSummary$Confond   = as.factor(MSLSummary$Confond)
MSLSummary$SeqinBlock   = as.factor(MSLSummary$SeqinBlock)
sorted_Block  = paste(sort(as.integer(levels(MSLSummary$Block))))
MSLSummary$Block = factor(MSLSummary$Block, levels = sorted_Block)
sorted_Session  = c(' TrainingPreNap',  ' TestPreNap',  ' TestPostNap',  ' TestPostNight')
MSLSummary$Session = factor(MSLSummary$Session, levels = sorted_Session)


100 -length(MSLData$Acc[MSLData$Acc=='1'])/length(MSLData$Acc)*100

100 -length(MSLData$Acc[MSLData$Acc=='1' & MSLData$Outlier==TRUE])/length(MSLData$Acc[MSLData$Acc=='1'])*100

ezANOVA(MSLSummary[MSLSummary$Block %in% as.factor(c(1:20,21:24)),],
        dv = .(Mean),
        wid=.(Sub), 
        within = .(Condition,Session), detailed=T)

ezANOVA(MSLSummary[MSLSummary$Block %in% as.factor(c(18:20)) & MSLSummary$Condition==' notReact',],
        dv = .(Median),
        wid=.(Sub), 
        within = .(Block), detailed=T)

ezANOVA(MSLSummary[MSLSummary$Session == ' TestPostNight',],
        dv = .(Mean),
        wid=.(Sub), 
        within = .(Condition,Block), detailed=T)

####### Plot raw by block and sequence condition. 
#RT
ggplot(summarySE(MSLSummary, measurevar="Mean", groupvars=c(  "Block", "Condition"),na.rm=T),
       aes(x=Block, y=Mean, group = Condition )       ) + 
  geom_point(aes(color=Condition))+
  geom_line(aes(color=Condition))+
  geom_ribbon(aes(ymin = Mean-se,
                  ymax = Mean+se,fill=factor(Condition)), alpha = 0.4)+
  scale_color_manual(values=c("#AAAAEF","magenta")) +
  scale_fill_manual(values=c("#AAAAEF","magenta")) +
  ylab("REaction Time (ms)")+
  xlab("Block of training")+
  theme_classic() 


ggplot(summarySE(MSLSummary[MSLSummary$Block %in% as.factor(c(18:20,21:24)),], measurevar="Median", 
                 groupvars=c("Sub","Condition","Session")), aes(x=Session, y=Median, fill = Condition)) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,color=c("darkblue","darkviolet","darkblue","darkviolet"))+
  scale_fill_manual(values=c("#AAAAEF","magenta"))+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),color=c("darkblue","darkblue","darkviolet","darkviolet") ) +
  theme_classic() 

#RT with Random SRTT
ggplot(rbind(summarySE(MSLSummary, measurevar="Mean", groupvars=c(  "Block", "Condition"),na.rm=T),summarySE(randomSummary, measurevar="Mean", groupvars=c(  "Block","Condition"),na.rm=T)),
       aes(x=Block, y=Mean, group = Condition )) + 
  geom_point(aes(color=Condition),size =3)+
  geom_line(aes(color=Condition))+
  geom_ribbon(aes(ymin = Mean-se,
                  ymax = Mean+se,fill=factor(Condition)),alpha = c(0.3))+
  scale_color_manual(values=c("#3b3bf5","magenta","Black","Black")) +
  scale_fill_manual(values=c("#3b3bf5","magenta","Black","Black")) +
  ylab("Reaction time (ms)")+
  xlab("Block of training")+
  coord_cartesian(xlim = c(17,24))+
  theme_classic() 

  # ylim(150,650)

#Accuracy with Random SRTT
ggplot(rbind(summarySE(MSLSummary, measurevar="Acc", groupvars=c(  "Block", "Condition"),na.rm=T),summarySE(randomSummary, measurevar="Acc", groupvars=c(  "Block","Condition"),na.rm=T)),
       aes(x=Block, y=Acc, group = Condition )       ) + 
  geom_point(aes(color=Condition))+
  geom_line(aes(color=Condition))+
  geom_ribbon(aes(ymin = Acc-se,
                  ymax = Acc+se,fill=factor(Condition)), alpha = 0.3)+
  scale_color_manual(values=c("#3b3bf5","magenta","Black","Black")) +
  scale_fill_manual(values=c("#3b3bf5","magenta","Black","Black")) +
  ylab("Accuracy (%)")+
  xlab("Block of training")+
  coord_cartesian(xlim = c(17,24))+
  theme_classic() +
  ylim(0.9,1)



#PI
ggplot(summarySE(MSLSummary[MSLSummary$Sub!='OL_CA_21',], measurevar="PI", groupvars=c(  "Block", "Condition"),na.rm=T),
       aes(x=Block, y=PI, group = Condition )       ) + 
  geom_point(aes(color=Condition))+
  geom_ribbon(aes(ymin = PI-se,
                  ymax = PI+se,fill=factor(Condition)), alpha = 0.4)+
  scale_color_manual(values=c("#AAAAEF","magenta")) +
  scale_fill_manual(values=c("#AAAAEF","magenta")) +
  ylab("REaction Time (ms)")+
  xlab("Block of training")+
  coord_cartesian(ylim = c(0,3))+
  theme_classic() 

####### Plot raw by block and sequence  
#RT
ggplot(summarySE(MSLSummary, measurevar="Median", groupvars=c(  "Block", "Sequence"),na.rm=T),
       aes(x=Block, y=Median, group = Sequence )       ) + 
  geom_point(aes(color=Sequence))+
  geom_line(aes(color=Sequence))+
  geom_ribbon(aes(ymin = Median-se,
                  ymax = Median+se,fill=factor(Sequence)), alpha = 0.4)+
  ylab("REaction Time (ms)")+
  xlab("Block of training")+
  theme_classic() 

#Acc
ggplot(summarySE(MSLSummary[MSLSummary$Block %in% as.factor(c(1:20)),], measurevar="Acc", groupvars=c(  "Block", "Sequence"),na.rm=T),
       aes(x=Block, y=Acc, group = Sequence )       ) + 
  geom_point(aes(color=Sequence))+
  geom_line(aes(color=Sequence))+
  geom_ribbon(aes(ymin = Acc-se,
                  ymax = Acc+se,fill=factor(Sequence)), alpha = 0.4)+
  ylab("REaction Time (ms)")+
  xlab("Block of training")+
  coord_cartesian(xlim = c(0,20))+
  coord_cartesian(ylim = c(0.90,1))+
  theme_classic() 

### Learning magnitude
# load(file = 'randomSummary.RData')

LearnMag=matrix(,length(allSub)*3,4)
colnames(LearnMag)=c('Sub',"Age",'Condition',"DV")
LearnMag=as.data.frame(LearnMag)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_cond in 1:length(allCond))
  {
  
  tmpMSL = MSLSummary[MSLSummary$Sub==allSub[idx_sub] & 
                        MSLSummary$Block %in% as.factor(c(1:4,17:20)) &
                        MSLSummary$Condition ==allCond[idx_cond],]
  
  meanMSLpre = mean(tmpMSL$Mean[tmpMSL$Session==" TrainingPreNap"])
  meanMSLpost = mean(tmpMSL$Mean[tmpMSL$Session==" TestPreNap"])

  
  LearnMag$Sub[counter]       = allSub[idx_sub]
  LearnMag$Condition[counter] = allCond[idx_cond]
  LearnMag$DV[counter]        = (meanMSLpre-meanMSLpost)/meanMSLpre*100
  LearnMag$Age[counter]    = questionnaires_Summarize_CA$Age[questionnaires_Summarize_CA$Sub==allSub[idx_sub]]
  
  counter = counter +1
  }
  tmpMSL = MSLSummary[MSLSummary$Sub==allSub[idx_sub] & 
                        MSLSummary$Block %in% as.factor(c(1:4,17:20))  ,]
  
  meanMSLpre = mean(tmpMSL$Mean[tmpMSL$Session==" TrainingPreNap"])
  meanMSLpost = mean(tmpMSL$Mean[tmpMSL$Session==" TestPreNap"])
  
  
  LearnMag$Sub[counter]       = allSub[idx_sub]
  LearnMag$Condition[counter] = "all"
  LearnMag$DV[counter]        = (meanMSLpre-meanMSLpost)/meanMSLpre*100
  LearnMag$Age[counter]    = questionnaires_Summarize_CA$Age[questionnaires_Summarize_CA$Sub==allSub[idx_sub]]
  
  counter = counter +1
  
}


LearnMag$Sub       = as.factor(LearnMag$Sub)
LearnMag$Condition = as.factor(LearnMag$Condition)


shapiro.test(LearnMag$DV[LearnMag$Condition=='all'])

ggplot(LearnMag[LearnMag$Condition=='all',], aes( x=Age,y=DV)) +
  geom_point(color=c("black")) +
  geom_smooth(method = lm,color="grey",alpha=0.2,fill='grey') +
  stat_cor(method = "pearson")+
  # coord_cartesian(ylim = c(80, 120 ))+
  theme_light()

cor.test(LearnMag$Age[LearnMag$Condition=='all'],LearnMag$DV[LearnMag$Condition=='all'])

ggplot(LearnMag[LearnMag$Condition==' react',], aes( x=Age,y=DV)) +
  geom_point(color=c("magenta")) +
  geom_smooth(method = lm,color="darkviolet",alpha=0.2,fill='darkviolet') +
  # stat_cor(method = "pearson")+
  coord_cartesian(ylim = c(-11, 42 ))+
  theme_light()

ggplot(LearnMag[LearnMag$Condition==' notReact',], aes( x=Age,y=DV)) +
  geom_point(color=c("darkblue")) +
  geom_smooth(method = lm,color="darkblue",alpha=0.2,fill='darkblue') +
  # stat_cor(method = "pearson")+
  coord_cartesian(ylim = c(-11, 40 ))+
  theme_light()
