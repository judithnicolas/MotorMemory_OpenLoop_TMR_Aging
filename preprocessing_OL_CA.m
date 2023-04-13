
%% Continuous data pre-processing
listSub=   getCompleteDatasets_CA;


for idx_sub =1 : length(listSub)

    sub = listSub{idx_sub};
        
        load([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trl_epoch.mat'])
        load([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '.mat'])
        load([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_scored_epoch.mat'])
        
        % re-ref
        cfg = [];
        cfg.trl = scoreSleep(scoreSleep(:,4)==2 |scoreSleep(:,4)==3,: );
        cfg.dataset = [initPath.Exp '\data\OL_CA\' sub '\exp\' sub '.edf']; % Take the edf to be the same as the scored files
        cfg.reref       = 'yes';
        cfg.channel     = 'all';
        cfg.implicitref = 'REF1';         % the implicit (non-recorded) reference channel is added to the data representation
        cfg.refchannel  = {'REF1', 'REF2'}; % the average of these two is used as the new reference, channel '53' corresponds to the right mastoid (M2)
        data_reref= ft_preprocessing(cfg);
        
        % Filtering
        cfg = [];
        cfg.channel       = 1:6;
        cfg.hpfilter      = 'yes'; % highpass filter (default = 'no')
        cfg.hpfreq        = 0.1;   % highpass frequency in Hz
        cfg.hpfiltord     = 4;   % filter order
        cfg.lpfilter      = 'yes'; % lowpass filter (default = 'no')
        cfg.lpfreq        = 30;    % lowpass  frequency in Hz
        cfg.continuous    = 'yes';
        data_EEG_filtered = ft_preprocessing(cfg, data_reref);
        
        cfg = [];
        cfg.channel       = 7:10;
        cfg.hpfilter      = 'yes'; % highpass filter (default = 'no')
        cfg.hpfreq        = 2;   % highpass frequency in Hz
        cfg.lpfilter      = 'yes'; % lowpass filter (default = 'no')
        cfg.lpfreq        = 15;    % lowpass  frequency in Hz
        cfg.continuous    = 'yes';
        data_EOG_filtered = ft_preprocessing(cfg, data_reref);
        
        cfg = [];
        cfg.channel       = 11:12;
        cfg.hpfilter      = 'yes'; % highpass filter (default = 'no')
        cfg.hpfreq        = 50;   % highpass frequency in Hz
        cfg.lpfilter      = 'yes'; % lowpass filter (default = 'no')
        cfg.lpfreq        = 150;    % lowpass  frequency in Hz
        cfg.continuous    = 'yes';
        data_EMG_filtered = ft_preprocessing(cfg, data_reref);
        
        cfg = [];
        data_filtered= ft_appenddata(cfg,data_EEG_filtered,data_EOG_filtered,data_EMG_filtered);
        
        cfg=[];
        cfg.artfctdef.reject = 'nan'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
        cfg.artfctdef.muscle.artifact = D.other.CRC.score{6,1}*D.Fsample;
        data_arousal_artifacts = ft_rejectartifact(cfg,data_filtered);
        
        
        % Artefact inspection Data browser
        cfg = [];
        cfg.viewmode   = 'vertical';
        cfg.ylim       = [-150 150];
        cfg.continuous = 'yes';
        cfg.blocksize  = 30;
        artf=ft_databrowser(cfg,data_arousal_artifacts);
        
        
        cfg=[];
        cfg.artfctdef.reject = 'nan'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
        cfg.artfctdef.visual.artifact =  artf.artfctdef.visual.artifact;
        cfg.artfctdef.muscle.artifact = D.other.CRC.score{6,1}*D.Fsample;
        data_artifacts = ft_rejectartifact(cfg,data_arousal_artifacts);
        
        artifacts = {};
        artifacts.muscle = D.other.CRC.score{6,1}*D.Fsample;
        artifacts.visual =  artf.artfctdef.visual.artifact;
        
        cfg = [];
        ic_data = ft_componentanalysis(cfg,data_arousal_artifacts);
        
        cfg          = [];
        cfg.viewmode = 'component';
        cfg.layout =  [initPath.FieldTrip '\template\layout\EEG1020.lay'];
        ft_databrowser(cfg, ic_data)
        
        [components] = input('Value(s) for component ([x1 x2...] if more than one values): ');
        
        cfg = [];
        cfg.component = components; % to be removed component(s)
        data = ft_rejectcomponent(cfg, ic_data, data_artifacts);
        
        cfg=  [];
        cfg.channel = 1:6;
        data = ft_selectdata(cfg,data);
        
        save ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_preprocessed_continuous.mat'], 'data')
    
end
