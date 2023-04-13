listSub= getScoredDatasets_CA;

avg =0;

for idx_sub =  2% : length(listSub)
    
    sub = listSub{idx_sub};
    fprintf(['\nload ' sub ' data\n'])
    
    load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_preprocessed_continuous.mat'])
    load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trl_epoch.mat'])
    
    trl_corrected = trl;%(find(trl(:,5) == 2 | trl(:,5) == 3),:);
    
    stim = [trl_corrected(:,1)+1*data.fsample trl_corrected(:,4)];
    
    stim(stim(:,1)> data.sampleinfo(end,2),:) = [];
    
    if avg ==1
        cfg = [];
        cfg.channel        = [1:6];
        cfg.avgoverchan    = 'yes';
        data = ft_selectdata(cfg, data);
        eegIdx = 5;
        nameFile = '_withRest_yasa_avgOverChan';
    else
        eegIdx = [5:10];
        nameFile = '_withRest_yasa';
        
    end
    
    
    cfg = [];
    cfg.sampleindex     = 'no';
    cfg.resamplefs      = 500;
    cfg.detrend         = 'no';
    data_resample = ft_resampledata(cfg, data);
    
    if data_resample.fsample==500
        step = 2;
    else
        step = 1;
    end
    if data.sampleinfo(1,1)==1
        data.sampleinfo(1,1)=0;
        data_resample.sampleinfo(1,1)=0;
        data.sampleinfo(1,2)=29999;
        data_resample.sampleinfo(1,2)=29999;  
        data.trial{1,1}(:,29999:30000)=nan(6,2);
        data_resample.trial{1,1}(:,end+1)=nan(6,1);
    end

    allSamples =horzcat([data.sampleinfo(1,1):step:data.sampleinfo(end,2)]',...
        nan(length([data.sampleinfo(1,1):step:data.sampleinfo(end,2)]),4));

    for idx_epoch = 1 : length(data.sampleinfo)
        
        t1= find(allSamples(:,1)==data.sampleinfo(idx_epoch,1));
        t2= find(allSamples(:,1) - data.sampleinfo(idx_epoch,2)==-1);
        
        if ~isempty(t1)
            
            allSamples(t1:t2,2)=...
                repmat(data.trialinfo(idx_epoch,1), 1,length([data.sampleinfo(idx_epoch,1):step:data.sampleinfo(idx_epoch,2)]));
            
            allSamples(t1:t2,3)=...
                repmat(2,1,length([data.sampleinfo(idx_epoch,1):step:data.sampleinfo(idx_epoch,2)]));
            
            allSamples(t1:t2,4)=...
                repmat(idx_epoch,1,length([data.sampleinfo(idx_epoch,1):step:data.sampleinfo(idx_epoch,2)]));
            
            allSamples(t1:t2,eegIdx)=...
                data_resample.trial{idx_epoch}';
        end
        
    end
    
    tmp = find(isnan(allSamples(:,5)));
    
    allSamples(tmp, 4)=0;
    allSamples(tmp, 2:3)=5;
    allSamples(tmp, eegIdx)=1;
    allSamples = horzcat(allSamples, nan(size(allSamples,1),1));

    for idx_stim = 1 : length(stim)
        
        tmpStim = find(allSamples(:,1)-stim(idx_stim,1)==1 | allSamples(:,1)-stim(idx_stim,1)==0);
        if ~isempty(tmpStim)
            if idx_stim~= length(stim)
                tmpStimNext = find(allSamples(:,1)-stim(idx_stim+1,1)==1 | allSamples(:,1)-stim(idx_stim+1,1)==0);% find(allSamples(:,1)== stim(idx_stim+1,1));
                if stim(idx_stim+1,1)-stim(idx_stim,1) > 10 * data.fsample
                    allSamples([tmpStim:1:tmpStim+5*data_resample.fsample],3) = stim(idx_stim,2);
                else
                    allSamples([tmpStim:1:tmpStimNext],3) = stim(idx_stim,2);
                end
            else
                if tmpStim+5*data.fsample< length(allSamples)
                    allSamples([tmpStim:1:tmpStim+5*data_resample.fsample],3) = stim(idx_stim,2);
                else
                    allSamples([tmpStim:1:end],3) = stim(idx_stim,2);
                end
            end
            allSamples(tmpStim,end)= stim(idx_stim,2);

        end
    end
    
    allSamples(end+1,1) = data_resample.fsample;
    

    
    fprintf(['\nsaving ' sub ' data\n'])
    dlmwrite([initPath.Exp '\data\OL_CA\' sub '\exp\' sub nameFile '.txt'],allSamples)
    
end