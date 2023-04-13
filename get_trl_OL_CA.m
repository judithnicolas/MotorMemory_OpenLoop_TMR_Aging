

listSub= getScoredDatasets_CA;

nbTrialsSent = [];

for idx_sub =1 : length(listSub)
    
    sub = listSub{idx_sub};

    FqSampling = 1000;
    
    % Getting the good trial matrix (/!\ depends on subject /!\)
    preStim  = 1; %sec
    postStim = 3; %sec
    
    %Via the vmkr
    eventfile     = [initPath.Exp '\data\OL_CA\' sub '\exp\' sub '.vmrk'];
    event = ft_read_event(eventfile);
    sel = strcmp({event.value}, {'T479'});
    smp = [event(sel).sample];
    
    selCDL  = strcmp({event.value}, {'cdl'});
    smpCDL  = [event(selCDL).sample];
    selODL  = strcmp({event.value}, {'odl'});
    smpODL  = [event(selODL).sample];
    
    %Via the vmkr
    eventfilePsychToolBox     = [initPath.Exp '\data\OL_CA\' sub '\exp\MSL_BilateralClosedLoop_' sub '_Pauselog.txt'];
    eventPsychToolBox= tdfread(eventfilePsychToolBox);
    
    if ismember(sub, 'OL_CA_02')
        eventPsychToolBox.Time(1:108,2)= eventPsychToolBox.Time(1:108,1)-eventPsychToolBox.Time(1,1);
        eventPsychToolBox.Time(109,2)= eventPsychToolBox.Time(108,2)+(event(219).sample-event(217).sample)/FqSampling;
        eventPsychToolBox.Time(110:end,2)=eventPsychToolBox.Time(110:end,1)-eventPsychToolBox.Time(109,1)+eventPsychToolBox.Time(109,2);
    elseif ismember(sub, 'OL_CA_05')
        eventPsychToolBox.Time(1:92,2)= eventPsychToolBox.Time(1:92,1)-eventPsychToolBox.Time(1,1);
        eventPsychToolBox.Time(93,2)= eventPsychToolBox.Time(92,2)+(event(195).sample-event(193).sample)/FqSampling;
        eventPsychToolBox.Time(94:end,2)=eventPsychToolBox.Time(94:end,1)-eventPsychToolBox.Time(93,1)+eventPsychToolBox.Time(93,2);
    else 
        eventPsychToolBox.Time(:,2)=eventPsychToolBox.Time(:,1)-eventPsychToolBox.Time(1,1);
    end
    
    if ismember(sub, 'OL_CA_19')
            eventPsychToolBox.Sample = eventPsychToolBox.Time(:,2);%*FqSampling; Because pauselog from vmkr
            
    else
            eventPsychToolBox.Sample = eventPsychToolBox.Time(:,2)*FqSampling;

    end
    
    counterPsychToolBox = length(eventPsychToolBox.Sample);
    coef= [];
    
    if ismember(sub, 'OL_CA_14')
        trig = 'R 14';
    else
        trig = 'T479';
    end
    
    for i=1:length(event)
        if strcmp(event(i).value, trig)
            if isempty(coef)
                coef= event(i).sample;
            end
        end
    end
    eventPsychToolBox.SampleEEG = eventPsychToolBox.Sample+coef;
    
    trl = [];
    trl(:,1) = round(eventPsychToolBox.SampleEEG-preStim*FqSampling);
    trl(:,2) = round(eventPsychToolBox.SampleEEG+postStim*FqSampling);
    trl(:,3) = repmat(-FqSampling*preStim,length(eventPsychToolBox.SampleEEG(:,1)),1);
    trl(:,4) = (eventPsychToolBox.Block(:,1)=='s'); % stim =1 ; random = 0;
    

    
    if trl(end,1)>smpODL
        ft_warning(['check trial matrix: events after end of recording ' sub])
        trl(trl(:,1)>smpODL,:)=[];
    end
    
    nbTrialsSent = vertcat(nbTrialsSent, [sum(trl(:,4)) (length(trl(:,4))-sum(trl(:,4)))] );
    save ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trlRaw.mat'], 'trl')

    
end