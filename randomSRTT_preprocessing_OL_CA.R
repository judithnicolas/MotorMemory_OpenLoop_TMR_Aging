
setwd("D:/Documents/Research/TMR/openLoop/bimanual_open_loop/data/OL_CA/group")
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)

randomSRTT = read.csv(file = 'RandomSRTT.csv',header = T,sep = ';');



randomSRTT$Block   = as.factor(randomSRTT$Block)
randomSRTT$Cue     = as.factor(randomSRTT$Cue)
randomSRTT$Rep     = as.factor(randomSRTT$Rep)
sorted_Session     = c(" preNap",  " postNight")
randomSRTT$Session = factor(randomSRTT$Session, levels = sorted_Session)

allSess       = levels(randomSRTT$Session)
allBlock      = levels(randomSRTT$Block)
allCue        = levels(randomSRTT$Cue)
allRep        = levels(randomSRTT$Rep)
allSub        = levels(randomSRTT$Sub)

####### OUtlier

for (idx_sub in 1:length(allSub))
{
  for (idx_block in 1 : length(allBlock))
  {
    
    # RT outlyers
    tmpLim = randomSRTT[randomSRTT$Sub==allSub[idx_sub] & randomSRTT$Acc=='1' & 
                          randomSRTT$Block==allBlock[idx_block],]
    
    
    IQR = quantile(tmpLim$RT)[4]-quantile(tmpLim$RT)[2]
    limInf = quantile(tmpLim$RT)[2]-(1.5*IQR)
    limSup = quantile(tmpLim$RT)[4]+(1.5*IQR)
    
    
    randomSRTT$Outlier[randomSRTT$Sub==allSub[idx_sub] & 
                         randomSRTT$Block==allBlock[idx_block]] = randomSRTT$Acc[randomSRTT$Sub==allSub[idx_sub] & 
                                                                                   randomSRTT$Block==allBlock[idx_block]]=='1' & 
      randomSRTT$RT[randomSRTT$Sub==allSub[idx_sub] & 
                      randomSRTT$Block==allBlock[idx_block]]>limInf     & randomSRTT$RT[randomSRTT$Sub==allSub[idx_sub] & 
                                                                                          randomSRTT$Block==allBlock[idx_block]]<limSup 
    
  }
}

####### Compute Accuracy by Block, mean rt per key press PI 

randomSummary=matrix(,length(allSub)*length(allBlock)*length(allSess),6)
colnames(randomSummary)=c(colnames(randomSRTT)[c(1:3)],"Mean","Median", "Acc")
randomSummary=as.data.frame(randomSummary)
counter = 1

for (idx_sub in 1:length(allSub))
{
  for (idx_sess in 1 : length(allSess))
  {
    
    for (idx_block in 1 : length(allBlock))
    {
      
      # RT outlyers
      tmp = randomSRTT[randomSRTT$Sub==allSub[idx_sub] & 
                         randomSRTT$Session==allSess[idx_sess] & 
                         randomSRTT$Block==allBlock[idx_block] ,]
      limInf = mean(tmp$RT[tmp$Acc=='1' ])-3*sd(tmp$RT[tmp$Acc=='1' ])
      limSup = mean(tmp$RT[tmp$Acc=='1' ])+3*sd(tmp$RT[tmp$Acc=='1' ])
      
      
      randomSummary$Sub[counter]       = allSub[idx_sub]
      if (allSess[idx_sess] ==  " postNight")
      {
        randomSummary$Block[counter]     = as.numeric(allBlock[idx_block])+36  
      }
      if (allSess[idx_sess] ==  " preNap")
      {
        randomSummary$Block[counter]     = allBlock[idx_block]
      }
      
      randomSummary$Session[counter]   = as.character(tmp$Session[1])
      
      
      # % Accuracy per block
      randomSummary$Acc[counter]       = sum(tmp$Acc)/length(tmp$Acc)
      
      #Mean RT per key presses
      randomSummary$Median[counter]      = median(tmp$RT[tmp$Acc=='1' & tmp$RT>limInf & tmp$RT<limSup ])
      randomSummary$Mean[counter]      = mean(tmp$RT[tmp$Acc=='1' & tmp$RT>limInf & tmp$RT<limSup ])
      
      counter = counter+1 
    }
    
  }
}


randomSummary$Sub       = as.factor(randomSummary$Sub)
randomSummary$Session   = as.factor(randomSummary$Session)
randomSummary$Block     = as.factor(randomSummary$Block)
sorted_Block     = c('1', '2', '3', '4', '37', '38', '39', '40')
randomSummary$Block = factor(randomSummary$Block, levels = sorted_Block)
randomSummary$Session   = as.factor(randomSummary$Session)
sorted_Session     = c(' preNap',  ' postNight')
randomSummary$Session = factor(randomSummary$Session, levels = sorted_Session)

randomSummary$Condition = randomSummary$Session

100 -length(randomSRTT$Acc[randomSRTT$Acc=='1'])/length(randomSRTT$Acc)*100

100 -length(randomSRTT$Acc[randomSRTT$Acc=='1' & randomSRTT$Outlier==TRUE])/length(randomSRTT$Acc[randomSRTT$Acc=='1'])*100

