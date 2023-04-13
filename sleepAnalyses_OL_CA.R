
setwd("D:/Documents/Research/TMR/openLoop/bimanual_open_loop/data/OL_CA/group")
library('ez')
library('Rmisc')


ci_perso = function (x, ci = 0.95) 
{
  a <- mean(x)
  s <- sd(x)
  n <- length(x)
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(c( mean = round(a,2), lower = round(a - error,2),upper =round( a + error,2)))
}

############ SLEEP Duration
#sleepDuration.csv is included in the data repository 

duration = read.csv(file = 'sleepDuration.csv',header = T,sep = ';');

tapply(duration$DV/60,list(duration$Stage),CI)
CI(duration$Percentage[duration$Stage == ' sleep' ])

CI(duration$Percentage[duration$Stage == ' S2']  + 
     duration$Percentage[ duration$Stage == ' S3' ] + 
     duration$Percentage[ duration$Stage == ' REM'] )

all_sub= levels(duration$Sub)
all_stage = levels(duration$Stage)

############ SLEEP Latendy
latency = read.csv(file = 'sleepLatency.csv',header = T,sep = ';');
 CI(latency$S1/60)

############ SLEEP arousal
arousal = read.csv(file = 'sleepArousal.csv',header = T,sep = ';');


############ STIM REPARTION
cuesCount = read.csv(file = 'stimEfficiency.csv',header = T,sep = ';');

tapply(cuesCount$Percentage[cuesCount$Stage==' NREM'],cuesCount$type[cuesCount$Stage==' NREM'],CI)

tapply(cuesCount$DV[cuesCount$type==' all'],list(cuesCount$Stage[cuesCount$type==' all']),CI)
tapply(cuesCount$DV[cuesCount$type==' stim'],list(cuesCount$Stage[cuesCount$type==' stim']),CI)
tapply(cuesCount$DV[cuesCount$type==' random'],list(cuesCount$Stage[cuesCount$type==' all']),CI)

ezANOVA (cuesCount[cuesCount$type!=' all' & (cuesCount$Stage!=' total' & cuesCount$Stage!=' S2' & cuesCount$Stage!=' S3'),],
         dv = .(DV),wid=.(Sub), within = .(type,Stage), detailed=T)

############ EEG trial count
trlCount <- read.csv("trlCount.csv", header=T,sep = '\t')
trlCount$PercentageRel = (trlCount$NREM_Rel /trlCount$Sent_Rel)*100
trlCount$PercentageIrrel = (trlCount$NREM_Irrel/trlCount$Sent_Irrel)*100
trlCount$PercentageAll = (trlCount$NREM_All/trlCount$Sent_All)*100


t.test(trlCount$Sent_Rel,trlCount$Sent_Irrel,paired = T)
ci_perso(trlCount$Sent_Rel)
ci_perso(trlCount$Sent_Irrel)

t.test(trlCount$NREM_Rel,trlCount$NREM_Irrel,paired = T)
ci_perso(trlCount$NREM_Rel)
ci_perso(trlCount$NREM_Irrel)

cor.test(trlCount$NREMRel,TMRIndex$index_RT[TMRIndex$Time =='late'])

############ SLEEP instrucitons
questionnaires_Summarize_CA <- read_delim("questionnaires Summarize CA.csv", ";", escape_double = FALSE, trim_ws = TRUE)
questionnaires_Summarize_CA$Sub = as.factor(questionnaires_Summarize_CA$Sub)
questionnaires_Summarize_CA = 
  questionnaires_Summarize_CA[questionnaires_Summarize_CA$Sub != 'OL_CA_16',] #Participants removed due to experimental error
questionnaires_Summarize_CA$Sexe = as.factor(questionnaires_Summarize_CA$Sexe)


CI(questionnaires_Summarize_CA$sleep_duration_estimate_StMary_Night3)/60/60
CI(questionnaires_Summarize_CA$sleep_quality_St_Mary_Night3)
CI(questionnaires_Summarize_CA$sleep_duration_estimate_StMary_Night4)
CI(questionnaires_Summarize_CA$sleep_quality_St_Mary_Night4)


t.test(questionnaires_Summarize_CA$sleep_duration_estimate_StMary_Night3,questionnaires_Summarize_CA$sleep_duration_estimate_StMary_Night4,paired = T)
t.test(questionnaires_Summarize_CA$sleep_quality_St_Mary_Night3,questionnaires_Summarize_CA$sleep_quality_St_Mary_Night4,paired = T)


### Acti + sleep journal
CI(c(questionnaires_Summarize_CA$`Actiwatch night 1`,questionnaires_Summarize_CA$Sleepjournal_night_1))/60/60
CI(c(questionnaires_Summarize_CA$`Actiwatch night 2`,questionnaires_Summarize_CA$Sleepjournal_night_2))/60/60
CI(c(questionnaires_Summarize_CA$`Actiwatch night 3`,questionnaires_Summarize_CA$Sleepjournal_night_3))/60/60
CI(c(questionnaires_Summarize_CA$`Actiwatch night 4`,questionnaires_Summarize_CA$Sleepjournal_night_4))/60/60


Nights=matrix(,length(allSub)*4,3)
colnames(Nights)=c('Sub','Session',"DV")
Nights=as.data.frame(Nights)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_sess in 1: 4)
  {
    if (idx_sess==1)
    {
      tmp = mean(c(questionnaires_Summarize_CA$`Actiwatch night 1`[questionnaires_Summarize_CA$Sub == allSub[idx_sub]],questionnaires_Summarize_CA$Sleepjournal_night_1[questionnaires_Summarize_CA$Sub == allSub[idx_sub]]))/60/60  
    }
    if (idx_sess==2)
    {
      tmp = mean(c(questionnaires_Summarize_CA$`Actiwatch night 2`[questionnaires_Summarize_CA$Sub == allSub[idx_sub]],questionnaires_Summarize_CA$Sleepjournal_night_2[questionnaires_Summarize_CA$Sub == allSub[idx_sub]]))/60/60  
    }
    if (idx_sess==3)
    {
      tmp = mean(c(questionnaires_Summarize_CA$`Actiwatch night 3`[questionnaires_Summarize_CA$Sub == allSub[idx_sub]],questionnaires_Summarize_CA$Sleepjournal_night_3[questionnaires_Summarize_CA$Sub == allSub[idx_sub]]))/60/60  
    }
    if (idx_sess==4)
    {
      tmp = mean(c(questionnaires_Summarize_CA$`Actiwatch night 4`[questionnaires_Summarize_CA$Sub == allSub[idx_sub]],questionnaires_Summarize_CA$Sleepjournal_night_4[questionnaires_Summarize_CA$Sub == allSub[idx_sub]]))/60/60  
    }
    
    
    Nights$Sub[counter]     = allSub[idx_sub]
    Nights$Session[counter] = idx_sess
    Nights$DV[counter]      = as.numeric(tmp)
    
    counter = counter +1
  }
}


Nights$Session=as.factor(Nights$Session)
Nights$Sub=as.factor(Nights$Sub)

ezANOVA (Nights, dv = .(DV), wid=.(Sub), within = .(Session), detailed=T)



AOVplot = ezPlot(
  data = Nights,
  , dv = .(DV)
  , wid = .(Sub)
  , within = .(Session)
  , x = .(Session)
  , x_lab =    'Block'
  , y_lab =    'RT' )
print(AOVplot)

############ Participant characteristics

CI(questionnaires_Summarize_CA$Edinburgh_handeness)
CI(questionnaires_Summarize_CA$Daytime_sleepiness)

CI(questionnaires_Summarize_CA$Beck_depression)
CI(questionnaires_Summarize_CA$Beck_anxiety)
CI(questionnaires_Summarize_CA$PQSI)
CI(questionnaires_Summarize_CA$Chronotype)



