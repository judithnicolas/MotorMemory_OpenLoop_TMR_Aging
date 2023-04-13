%% 

% Results presented in 3.2.1 and Figure S4


%%
clc

listSub = getScoredDatasets_CA;

%% Load data & baseline

allSubERP = {};

for idx_sub = 1 : length(listSub)
    
    sub = listSub{idx_sub};
    load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_preprocessed_continuous.mat'])
    load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trl_epoch.mat'])

    trl_corrected = trl;
    if idx_sub<8 
        trl_corrected(:,1:2) = trl(:,1:2) + 0.21*data.fsample;% correction for trigger lag
    else
        trl_corrected(:,1:2) = trl(:,1:2) + 0.52*data.fsample;% correction for trigger lag
    end
    
    cfg = [];
    cfg.trl      = trl_corrected(find(trl_corrected(:,5) == 2 | trl_corrected(:,5) == 3),1:4);
    data_epoched = ft_redefinetrial(cfg, data);
    
    cfg = [];
    cfg.resamplefs      = 100;
    data_epoched = ft_resampledata(cfg, data_epoched);   

    
    %NREM react vs not react and all
    cfg = [];
    cfg.baselinewindow = [-0.3 -0.1];
    cfg.demean         = 'yes';
    data_epoched = ft_preprocessing(cfg, data_epoched);

    %react
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1) == 1);
    timelock_data  = ft_timelockanalysis(cfg, data_epoched);
    allSubERP.nbTrl(idx_sub,1) =  length(cfg.trials);

    timelock_data.cfg=[];
    allSubERP.react{idx_sub}= timelock_data;
    
    % not react
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)== 0 );
    timelock_data  = ft_timelockanalysis(cfg, data_epoched);
    allSubERP.nbTrl(idx_sub,2) =  length(cfg.trials);

    timelock_data.cfg=[];
    allSubERP.notReact{idx_sub}= timelock_data;
    
    % all
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     =find(data_epoched.trialinfo(:,1)== 0  |  data_epoched.trialinfo(:,1)== 1  );
    timelock_data  = ft_timelockanalysis(cfg, data_epoched);
    allSubERP.nbTrl(idx_sub,3) =  length(cfg.trials);
 
    timelock_data.cfg=[];
    allSubERP.all{idx_sub} = timelock_data;

    allSubERP.zero{idx_sub} = allSubERP.all{idx_sub};
    allSubERP.zero{idx_sub}.avg = zeros(size(allSubERP.all{idx_sub}.avg,1),size(allSubERP.all{idx_sub}.avg,2));
    cfg = [];
    cfg.parameter = 'avg';
    cfg.operation = 'x1-x2';
    allSubERP.subtracted{idx_sub}= ft_math(cfg,allSubERP.react{idx_sub},allSubERP.notReact{idx_sub});

    
end

cfg = [];
cfg.keepindividual = 'yes';
grdAvgERP.react = ft_timelockgrandaverage(cfg, allSubERP.react{:});
grdAvgERP.notReact = ft_timelockgrandaverage(cfg, allSubERP.notReact{:});
grdAvgERP.all = ft_timelockgrandaverage(cfg, allSubERP.all{:});
grdAvgERP.zero = ft_timelockgrandaverage(cfg, allSubERP.zero{:});
grdAvgERP.subtracted = ft_timelockgrandaverage(cfg, allSubERP.subtracted{:});


%% Plot all ERP

cfg = [];
cfg.channel= 'all';
cfg.xlim = [-0.31 2.5];
cfg.ylim = [-5.1 5.1];
figure
cfg.color = 'k';
h_plot_erf(cfg,allSubERP.all');

% From cluster based analysis
xline(0.37,'Color',[227, 227, 227]/255);xline(0.53,'Color',[227, 227,227]/255)
xline(1.03,'Color',[227, 227, 227]/255);xline(1.35,'Color',[227, 227,227]/255)
xline(0,'g','Color',[227, 227, 227]/255);yline(5,'Color',[227, 227, 227]/255);yline(-5,'Color',[227, 227, 227]/255)

set(gca,'TickDir','out');

%% Plot two conditions ERP


cfg = [];
cfg.channel= 'all';
cfg.xlim = [-0.31 2.5];
cfg.ylim = [-5.1 5.1];
figure

cfg.color = 'm';
h_plot_erf(cfg,allSubERP.react')
cfg.color = 'y';
h_plot_erf(cfg,allSubERP.notReact');hold on

% xline(0.37,'Color',[227, 227, 227]/255);xline(0.53,'Color',[227, 227,227]/255)
% xline(1.03,'Color',[227, 227, 227]/255);xline(1.35,'Color',[227, 227,227]/255)
xline(0,'g','Color',[227, 227, 227]/255);yline(5,'Color',[227, 227, 227]/255);yline(-5,'Color',[227, 227, 227]/255);yline(0,'Color',[227, 227, 227]/255)
set(gca,'TickDir','out');

%% Plots per channel

for idx_channel = 5:6% length(allSubERP.all{1}.label)
    
    cfg = [];
    cfg.channel= allSubERP.all{1}.label{idx_channel};%'all'
    cfg.xlim = [-0.31 2.5];
    cfg.ylim = [-5 7];
    figure
%     cfg.color = 'k';
%     cfg.title = allSubERP.all{1}.label{idx_channel};
%     h_plot_erf(cfg,allSubERP.all');hold on
%     
%     xline(0.37,'Color',[227, 227, 227]/255);xline(0.53,'Color',[227, 227,227]/255)
%     xline(1.03,'Color',[227, 227, 227]/255);xline(1.35,'Color',[227, 227,227]/255)
%     xline(0,'g','Color',[227, 227, 227]/255);yline(5,'Color',[227, 227, 227]/255);yline(-5,'Color',[227, 227, 227]/255)
%     set(gca,'TickDir','out');
%     
%     hold off
% 
%     
    figure
    cfg.color = 'm';
    h_plot_erf(cfg,allSubERP.react')
    cfg.color = 'y';
    h_plot_erf(cfg,allSubERP.notReact');hold on
    
%     xline(0.37,'Color',[227, 227, 227]/255);xline(0.53,'Color',[227, 227,227]/255)
%     xline(1.03,'Color',[227, 227, 227]/255);xline(1.35,'Color',[227, 227,227]/255)
    xline(0,'g','Color',[227, 227, 227]/255);yline(0,'Color',[227, 227, 227]/255)
    yline(5,'Color',[227, 227, 227]/255);yline(-5,'Color',[227, 227, 227]/255)
    
    hold off
end


%% Stat

load neighboursPerso.mat

cfg                     = [];
cfg.design(1,1:2*(length(listSub)))  = [ones(1,(length(listSub))) 2*ones(1,(length(listSub)))];
cfg.design(2,1:2*(length(listSub)))  = [1:(length(listSub)) 1:(length(listSub))];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
cfg.method              = 'montecarlo';       
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.minnbchan           = 0;              
cfg.neighbours          = neighbours_perso; 
cfg.tail                = 0;                    
cfg.clustertail         = 0;
cfg.alpha               = 0.05; 
cfg.clusteralpha        = 0.05;  
cfg.avgoverchan         = 'yes';
cfg.numrandomization    = 500;      % number of draws from the permutation distribution
cfg.latency             = [0 2.5];

[statERP] = ft_timelockstatistics(cfg,  allSubERP.react{:}, allSubERP.notReact{:});
[statAll] = ft_timelockstatistics(cfg,  allSubERP.all{:}, allSubERP.zero{:});


cfg = [];
cfg.avgoverchan         = 'yes';
cfg.latency             = [0 2.5];
tmpup = ft_selectdata(cfg,grdAvgERP.react);
tmpdown = ft_selectdata(cfg,grdAvgERP.notReact);
tmpAll= ft_selectdata(cfg,grdAvgERP.all);
tmpSubtracted= ft_selectdata(cfg,grdAvgERP.subtracted);
tmpAll.mask = statAll.mask;

cfg = [];
cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.maskparameter = 'mask';
cfg.maskfacealpha = 0.5;
cfg.ylim = [-2.5 2.5];
figure;ft_singleplotER(cfg, tmpAll)


erpNegPeak = [];

for idx_sub = 1 : length(allSubERP.all)

    cfg = [];
    cfg.avgoverchan         = 'yes';
    cfg.latency             = [0.36 0.56];
    cfg.avgovertime         = 'yes';
    tmp = ft_selectdata(cfg,allSubERP.react{idx_sub});

    erpNegPeak(idx_sub,1) = tmp.avg;

    tmp = ft_selectdata(cfg,allSubERP.notReact{idx_sub});
    erpNegPeak(idx_sub,2) = tmp.avg;
  
    
    cfg.latency             = [1.02 1.42];
    tmp = ft_selectdata(cfg,allSubERP.react{idx_sub});

    erpNegPeak(idx_sub,3) = tmp.avg;

    tmp = ft_selectdata(cfg,allSubERP.notReact{idx_sub});
    erpNegPeak(idx_sub,4) = tmp.avg;

end

csvwrite([initPath.Exp '\data\OL_CA\group\erpNegPeak.csv'],erpNegPeak)

% Computes effect size (needs grand Avg option cfg.keepindividual = 'yes';)

cfg = [];
cfg.channel = 'all';
cfg.latency = [0.36 0.56];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
all = ft_selectdata(cfg, grdAvgERP.all);
zero  = ft_selectdata(cfg, grdAvgERP.zero);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSub)))  = [ones(1,(length(listSub))) 2*ones(1,(length(listSub)))];
cfg.design(2,1:2*(length(listSub)))  = [1:(length(listSub)) 1:(length(listSub))];
effect_roiTrough = ft_timelockstatistics(cfg, all, zero);
disp(effect_roiTrough)
%0.7059

cfg = [];
cfg.latency = [1.02 1.42];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
all = ft_selectdata(cfg, grdAvgERP.all);
zero  = ft_selectdata(cfg, grdAvgERP.zero);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSub)))  = [ones(1,(length(listSub))) 2*ones(1,(length(listSub)))];
cfg.design(2,1:2*(length(listSub)))  = [1:(length(listSub)) 1:(length(listSub))];

effect_roiPeak = ft_timelockstatistics(cfg, all, zero);
disp(effect_roiPeak)
%0.7005

