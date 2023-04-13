setwd(dir = '/Users/u0129763/Documents/Research/TMR/openLoop/bimanual_open_loop/data/group//')
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)

### Preprocessing Generation
# load('generationSummary.RData')
# load('offLineGain.RData')


session = " GenerationPre"
condition = " react"


#pre registered
cor.test(generationSummary$Accuracy[generationSummary$Session==session & generationSummary$Condition==condition],TMRIndex$index_RT [TMRIndex$Time=='both'],paired = T)



cor.test(generationSummary$Accuracy[generationSummary$Session==session & generationSummary$Condition==condition],offLineGain$Gain_RT [offLineGain$Time=='both'  & offLineGain$Condition==condition],paired = T)




library("ggpubr")
ggplot(corMatrix, aes(x=genrationPre, y=tmrIndexEarly)) +
  geom_point(color="magenta") +
  geom_smooth(method = lm,color="magenta", fill="blue") +
  stat_cor(method = "pearson", label.x = 20)

ggplot(corMatrix, aes(x=genrationPre, y=TMRIndex$index_RT [TMRIndex$Time=='late'  ])) +
  geom_point(color="blue") +
  geom_smooth(method = lm,color="blue", fill="magenta") +
  stat_cor(method = "pearson", label.x = 20)


ggplot(corMatrix, aes(x=genrationPre, y=TMRIndex$index_RT [TMRIndex$Time=='both'  ])) +
  geom_point(color="blue") +
  geom_smooth(method = lm,color="blue", fill="magenta") +
  stat_cor(method = "pearson", label.x = 20)



