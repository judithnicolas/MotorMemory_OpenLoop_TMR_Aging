listSub = getScoredDatasets_CA;


for idx_sub = 2%1 : length(listSub)
    
    sub = listSub{idx_sub};
    load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_preprocessed_continuous.mat'])
    load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trl_epoch.mat'])
    
    trl_corrected = trl;
    if idx_sub<8 
        trl_corrected(:,1:2) = trl(:,1:2) + 0.21*data.fsample;
    else
        trl_corrected(:,1:2) = trl(:,1:2) + 0.52*data.fsample;
    end

    cfg = [];
    cfg.trl      = trl_corrected(find(trl_corrected(:,5) == 2 | trl_corrected(:,5) == 3),:);
    data_epoched = ft_redefinetrial(cfg, data);

    
    disp (['loading ' sub ' dataset'])
    
    cfg = [];
    cfg.channel        = [1:6];
    cfg.baselinewindow = [-0.3 -0.1];
    cfg.demean         = 'yes';
    data_epoched = ft_preprocessing(cfg, data_epoched);
    
    cfg            = [];
    cfg.resamplefs = 500;
    data_epoched = ft_resampledata(cfg, data_epoched);
    
    cfg = [];
    cfg.channel        = [1:6];
    cfg.avgoverchan    = 'yes';
    data_epoched = ft_selectdata(cfg, data_epoched);
    
    data_epoched.trialinfo= trl_corrected(trl_corrected(:,5)== 2 | trl_corrected(:,5)== 3,4:5);

    data_pac=[];
%     for idx_chan = 1 : length(data.label)
        
        for idx_trl = 1 :length(data_epoched.trial)
            
            data_pac = vertcat(data_pac,[data_epoched.trial{idx_trl}(1,:) data_epoched.trialinfo(idx_trl,:) idx_trl]);

        end
        
%     end
    

    dlmwrite([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_ERpac_avgOverChan.txt'],data_pac)
    
end