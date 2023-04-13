setwd(dir = '/Users/u0129763/Documents/Research/TMR/openLoop/bimanual_open_loop/data/group//OL_CA/')
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)

### Preprocessing Generation
generation  = read.csv(file = 'generation_allKeys.csv',header = T,sep = ';')
generation$attempts = as.factor(generation$attempts)
generation$Sequence = as.factor(generation$Sequence)
generation$OrdinalPosition = as.factor(generation$OrdinalPosition)

allSequence = levels(generation$Sequence)
allSub      = levels(generation$Sub)
allSession  = levels(generation$Session)
allPos  = levels(generation$OrdinalPosition)

sorted_Session  = c(' GenerationPre',  ' GenerationPost')
generation$Session = factor(generation$Session, levels = sorted_Session)


generationSummary=matrix(,length(allSub)*length(allSequence)*length(allSession)*9,6)
colnames(generationSummary)=c('Sub','Session','Sequence','Condition','OrdinalPosition','Accuracy')
generationSummary=as.data.frame(generationSummary)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_sess in 1:length(allSession))
  {
    for (idx_seq in 1:length(allSequence))
    {    
      for (idx_pos in 1:length(allPos))
      {
      
        tmp = generation[generation$Sub== allSub[idx_sub] 
                                  & generation$Session== allSession[idx_sess] 
                                  & generation$Sequence == allSequence[idx_seq]
                                  & generation$OrdinalPosition   == allPos[idx_pos]
                         ,]
        generationSummary$Sub[counter]       = allSub[idx_sub]
        generationSummary$Session[counter]   = allSession[idx_sess]
        generationSummary$Sequence[counter]  = allSequence[idx_seq]
        generationSummary$Condition[counter] = as.character(tmp$Condition[1])
        generationSummary$OrdinalPosition[counter] = allPos[idx_pos]
        
        generationSummary$Accuracy[counter] = sum(tmp$Accuracy)/length(tmp$Accuracy)*100
        
        counter = counter+1
      }

      generationSummary$Accuracy[counter] = mean(generationSummary$Accuracy[generationSummary$Sub== allSub[idx_sub] 
                                                                            & generationSummary$Session== allSession[idx_sess] 
                                                                            & generationSummary$Sequence == allSequence[idx_seq]],na.rm=T)
      
      generationSummary$Sub[counter]       = allSub[idx_sub]
      generationSummary$Session[counter]   = allSession[idx_sess]
      generationSummary$Sequence[counter]  = allSequence[idx_seq]
      generationSummary$Condition[counter] = as.character(tmp$Condition[1])
      generationSummary$OrdinalPosition[counter] = 'all'
      
      
      counter = counter+1
      
    }
  }
}
generationSummary$Sub= as.factor(generationSummary$Sub)
generationSummary$Sequence= as.factor(generationSummary$Sequence)
generationSummary$Condition= as.factor(generationSummary$Condition)
sorted_Session  = c(' GenerationPre',  ' GenerationPost')
generationSummary$Session = factor(generationSummary$Session, levels = sorted_Session)


shapiro.test(generationSummary$Accuracy[generationSummary$OrdinalPosition=='all' 
                                        & generationSummary$Session==' GenerationPre'  
                                        & generationSummary$Condition==' react' ])
shapiro.test(offLineGain$Gain_RT[offLineGain$Time=='both' & offLineGain$Condition==' react' ])

cor.test(generationSummary$Accuracy[generationSummary$OrdinalPosition=='all' 
                                    & generationSummary$Session==' GenerationPre'  
                                    & generationSummary$Condition==' react' ],
         TMRIndex$index_RT[TMRIndex$Time=='both' ], paired = T)


plot(generationSummary$Accuracy[generationSummary$OrdinalPosition=='all' 
                                    & generationSummary$Session==' GenerationPre'  
                                    & generationSummary$Condition==' react' ],
         offLineGain$Gain_RT[offLineGain$Time=='both' & offLineGain$Condition==' react' ])

ezANOVA (generationSummary[generationSummary$OrdinalPosition==1,], dv = .(Accuracy), wid = .(Sub),
         within= .(Condition,Session), detailed=T)

ggplot(generationSummary,
       aes(x=Session, y=Accuracy, fill = Condition)) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,color=c("darkblue","darkviolet","darkblue","darkviolet"))+
  scale_fill_manual(values=c("#AAAAEF","magenta"))+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),color=c("darkblue","darkblue","darkviolet","darkviolet") ) +
  coord_cartesian(ylim=c(0,100))+
  theme_classic() 

ezANOVA (generationSummary, dv = .(Accuracy), wid = .(Sub),
         within= .(Sequence,Session), detailed=T)

ggplot(summarySE(generationSummary, measurevar="Accuracy", groupvars=c("OrdinalPosition","Sub")),
       aes(x=OrdinalPosition, y=Accuracy)) +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA)+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=OrdinalPosition), position=position_dodge(1) ) +
  ylim(0,100)+
  theme_classic() 

pairwise.t.test(generationSummary$Accuracy,generationSummary$OrdinalPosition,p.adjust.method = "fdr")

summarySequence = summarySE(generationSummary, measurevar="Accuracy", groupvars=c("Condition","Session","Sub"))
cor.test(summarySequence$Accuracy[summarySequence$Condition== ' react' & summarySequence$Session== ' GenerationPre'],
           TMRIndex$index_RT[TMRIndex$Time=='both'],paired=T) 

plot(summarySequence$Accuracy[summarySequence$Condition== ' react' & summarySequence$Session== ' GenerationPre'],
   TMRIndex$index_RT[TMRIndex$Time=='both'])

