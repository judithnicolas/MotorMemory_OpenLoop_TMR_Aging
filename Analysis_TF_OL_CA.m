%% 

% Results presented in 3.2.1 

%%
clc
listSub = getScoredDatasets_CA;
%% Load data & baseline

allSubTF = {};


for idx_sub = 1 : length(listSub)
    
    sub = listSub{idx_sub};
    
    load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_preprocessed_continuous.mat'])
    load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trl_epoch.mat'])
    data = rmfield(data, 'trialinfo');

    trl_corrected = trl;

    if idx_sub<8 
        trl_corrected(:,1:2) = trl(:,1:2) + 0.21*data.fsample;% correction for trigger lag
    else
        trl_corrected(:,1:2) = trl(:,1:2) + 0.52*data.fsample;% correction for trigger lag
    end
    
    disp (['loading ' sub ' dataset'])
    
    cfg = [];
    cfg.trl      = trl_corrected(find(trl_corrected(:,5) == 2 | trl_corrected(:,5) == 3),1:3);
    data_epoched = ft_redefinetrial(cfg, data);
    
    cfg = [];
    cfg.channel        = [1:6];
    cfg.resamplefs      = 100;
    cfg.demean          = 'yes';
    data_epoched = ft_resampledata(cfg, data_epoched);
    
    data_epoched.trialinfo = trl_corrected(find(trl_corrected(:,5) == 2 | trl_corrected(:,5) == 3),4);
     
    % NREM up vs down and all

    foi = 5:0.5:30; % 0.1 Hz steps
    toi = -1:0.01:3; % 0.1 s steps
    
    %asso
    cfg = [];
    
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi;
    cfg.toi        = toi;
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;   
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)==1);
    timelock_data  = ft_freqanalysis(cfg, data_epoched);
    allSubTF.nbTrl(idx_sub,1) =  length(cfg.trials);
    
    cfg = [];
    cfg.baseline     = [-0.3 -0.1] ;
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);
    
    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSubTF.react{idx_sub}= timelock_data;
    
    %not reqct
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi; % 0.1 Hz steps
    cfg.toi        = toi; % 0.1 s steps
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;    
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)==0 );
    timelock_data  = ft_freqanalysis(cfg, data_epoched);
    
    allSubTF.nbTrl(idx_sub,2) =  length(cfg.trials);
     
    cfg = [];
    cfg.baseline     = [-0.3 -0.1] ;
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);
    
    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSubTF.notReact{idx_sub}= timelock_data;
    
    %All
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi; % 0.1 Hz steps
    cfg.toi        = toi; % 0.1 s steps
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;   
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)==1 | data_epoched.trialinfo(:,1)==0);
    timelock_data  = ft_freqanalysis(cfg, data_epoched);
    
    allSubTF.nbTrl(idx_sub,3) =  length(cfg.trials);
    
    cfg = [];
    cfg.baseline     = [-0.3 -0.1] ;
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);

    timelock_data.cfg        = [];
    timelock_data.freq       = foi; 
    allSubTF.all{idx_sub}= timelock_data;

    allSubTF.zero{idx_sub} = allSubTF.all{idx_sub};
    allSubTF.zero{idx_sub}.powspctrm = zeros(size(allSubTF.all{idx_sub}.powspctrm,1),size(allSubTF.all{idx_sub}.powspctrm,2),size(allSubTF.all{idx_sub}.powspctrm,3));
 

    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x1-x2';
    allSubTF.subtracted{idx_sub}= ft_math(cfg,allSubTF.react{idx_sub},allSubTF.notReact{idx_sub});
end


cfg = [];
grdAvgTF.react = ft_freqgrandaverage(cfg, allSubTF.react{:});
grdAvgTF.notReact = ft_freqgrandaverage(cfg, allSubTF.notReact{:});
grdAvgTF.all = ft_freqgrandaverage(cfg, allSubTF.all{:});
grdAvgTF.subtracted = ft_freqgrandaverage(cfg, allSubTF.subtracted{:});

% clearvars -except initPath grdAvg allSub allSubERP  grdAvgERP listSubEEG 
%% Plots
%From ERPs
cfg = [];
cfg.avgoverchan         = 'yes';
cfg.latency             = [0 2.5];
tmpup = ft_selectdata(cfg,grdAvgERP.react);
tmpdown = ft_selectdata(cfg,grdAvgERP.notReact);
tmpAll= ft_selectdata(cfg,grdAvgERP.all);


cfg = [];
cfg.layout      = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.xlim     = [-0.10 2.5];
cfg.zlim     = [-0.1 0.25];
% cfg.ylim     = [12 16];

figure;ft_multiplotTFR(cfg, grdAvgTF.react)
figure;ft_multiplotTFR(cfg, grdAvgTF.notReact)
figure;ft_multiplotTFR(cfg, grdAvgTF.subtracted)

cfg.avgoverchan = 'yes';
cfg.channel = 'Fz';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlegend  = 'yes';

figure;ft_singleplotTFR(cfg, grdAvgTF.react)
figure;ft_singleplotTFR(cfg, grdAvgTF.notReact)
figure;ft_singleplotTFR(cfg, grdAvgTF.subtracted)

hold on 
yyaxis right
plot(tmpAll.time,squeeze(mean(tmpAll.individual)),'-k')
ylim([-2.1 2])
%% Stat

% Comparison

load neighboursPerso.mat
cfg                     = [];
cfg.design(1,1:2*length(listSub))  = [ones(1,length(listSub)) 2*ones(1,length(listSub))];
cfg.design(2,1:2*length(listSub))  = [1:length(listSub) 1:length(listSub)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
cfg.method              = 'montecarlo';       
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.minnbchan           = 0;              
cfg.avgoverchan         = 'yes';
cfg.neighbours          = neighbours_perso; 
cfg.alpha               = 0.025; 
cfg.clusteralpha        = 0.05;     
cfg.numrandomization    = 500;      % number of draws from the permutation distribution
cfg.latency             = [0 2.5];
% cfg.frequency           = [12 16];

[statTF] = ft_freqstatistics(cfg,  allSubTF.react{:}, allSubTF.notReact{:});



%% Correlation % Results presented in 3.3

% compute statistics with correlationT

% TMR index both early and late with TFR relevant - TFR irrelevant
% Correlation between SO PAC and behaviour
% TMR index both early and late with TFR relevant - TFR irrelevant

[TMRindex , reactGains , notReactGains] = getBehav_OL_CA();

design=[];
design(1,1:length(listSub))          = TMRindex;

cfg = [];
cfg.statistic           = 'ft_statfun_correlationT';
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.neighbours          = neighbours_perso; 
cfg.minnbchan           = 0;              
cfg.alpha               = 0.025; 
cfg.clusteralpha        = 0.05;     
cfg.latency             = [0 2.5];
cfg.frequency           = [5 30];
cfg.numrandomization    = 500;
cfg.design              = design;
cfg.ivar                = 1;

statCorAll= ft_freqstatistics(cfg, allSubTF.subtracted{:});
