%% 

% Results presented in 3.2.3 and Figure S6 and S7 (supp)

%% Utilities
parulaHalfClear = [0.120903225806452,0.754387096774194,0.697670967741935;0.184429032258065,0.771745161290323,0.639109677419355;0.232335483870968,0.788816129032258,0.571925806451613;0.321229032258064,0.799632258064516,0.494625806451613;0.425522580645161,0.802896774193548,0.406574193548387;0.543361290322581,0.796035483870968,0.318687096774194;0.656267741935484,0.781867741935484,0.233203225806452;0.761322580645161,0.762416129032258,0.170558064516129;0.853903225806452,0.742577419354839,0.157364516129032;0.932683870967742,0.729283870967742,0.202977419354839;0.994006451612903,0.740158064516129,0.239851612903226;0.995622580645161,0.786170967741936,0.204903225806452;0.979764516129032,0.836200000000000,0.177732258064516;0.961296774193548,0.887400000000000,0.154329032258065;0.962651612903226,0.936567741935484,0.127070967741935;0.976900000000000,0.983900000000000,0.0805000000000000];


%% Cue-Locked preferred phase
% Import data from text file
% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "VarName1";
opts.VariableTypes = "double";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
PP_ER_irrel= readtable([initPath.Exp '\data\OL_CA\group\PP_ER_irrel_avgOverChan.txt'], opts);

PP_ERS0_irrel= readtable([initPath.Exp '\data\OL_CA\group\PP_irrel_avgOverChan_detectionFz.txt'], opts);
PP_ERS0_rel= readtable([initPath.Exp '\data\OL_CA\group\PP_rel_avgOverChan_detectionFz.txt'], opts);
PP_ERS0_rest= readtable([initPath.Exp '\data\OL_CA\group\PP_rest_avgOverChan_detectionFz.txt'], opts);

% Convert to output type
PP_ERS0_irrel = table2array(PP_ERS0_irrel);
PP_ERS0_rel = table2array(PP_ERS0_rel);
PP_ERS0_rest = table2array(PP_ERS0_rest);

clear opts



figure; circ_plot(PP_ERS0_rel,'hist',[],15,true,true,'linewidth',2,'color','r'); 
 [pval z] =circ_rtest(PP_ERS0_rel)
figure; circ_plot(PP_ERS0_irrel,'hist',[],15,true,true,'linewidth',2,'color','r'); 
[pval z] =circ_rtest(PP_ERS0_irrel)
figure; circ_plot(PP_ERS0_rest,'hist',[],15,true,true,'linewidth',2,'color','r'); 
[pval z] =circ_rtest(PP_ERS0_rest)
[mu_ERSO ul_ERSO ll_ERSO] = circ_mean([PP_ERS0_rel;PP_ERS0_irrel;PP_ERS0_rest]);

rad2deg([mu_ERSO ul_ERSO ll_ERSO])+360


 [pval table] = circ_wwtest(PP_ERS0_rel, PP_ERS0_irrel)
 [pval table] = circ_wwtest(PP_ERS0_rel, PP_ERS0_rest)
 [pval table] = circ_wwtest(PP_ERS0_irrel, PP_ERS0_rest)


%%

listSub = getScoredDatasets_CA;
erPAC_SO_allSub_avg = {};


counterSub = 1;
for idx_sub = 1 : length(listSub)
    sub = listSub{idx_sub};
            

    if idx_sub == 1
        load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_preprocessed_continuous.mat'])
        load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trl_epoch.mat'])
        
        cfg = [];
        cfg.trl      = trl(find(trl(:,5) == 2 | trl(:,5) == 3),1:4);
        data_epoched = ft_redefinetrial(cfg, data);

        cfg         = [];
        cfg.channel = {'Fz'};
        data = ft_selectdata(cfg,data);
        idx=find(ismember(trl(:,1),data.sampleinfo(:,1)));
        data.trialinfo= trl(idx,4:5);

        cfg = [];
        cfg.method     = 'mtmconvol';
        cfg.pad        = 'nextpow2';
        cfg.taper      = 'hanning';
        cfg.keeptrials = 'no';
        cfg.trials     = 1;
        
        cfg.foi = 5.25:0.5:29.25; 
        cfg.toi = -3:0.002:3; 
        cfg.t_ftimwin  = 5./cfg.foi;
        cfg.tapsmofrq  = 0.4 *cfg.foi;

        % er PAC
        erPAC_template  = ft_freqanalysis(cfg, data);
        erPAC_template.freq  = 7.25:0.5:29.25; 
        erPAC_template.time= -1:0.002:1.999; 
        erPAC_template.powspctrm= [];

    end
%    
    if ~ismember(idx_sub ,[4 12 13])
        data = textread([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_pacERSO_rel_detectionFz_PACall_new.txt'],'','delimiter',',');
        erPAC_SO_allSub_avg{counterSub,1} = erPAC_template;
        erPAC_SO_allSub_avg{counterSub,1}.powspctrm(1,:,1:length(erPAC_template.time)) = data(:,1:length(erPAC_template.time));

        data = textread([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_pacERSO_irrel_detectionFz_PACall_new.txt'],'','delimiter',',');
        erPAC_SO_allSub_avg{counterSub,2} = erPAC_template;
        erPAC_SO_allSub_avg{counterSub,2}.powspctrm(1,:,1:length(erPAC_template.time)) = data(:,1:length(erPAC_template.time));

        data = textread([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_pacERSO_rest_detectionFz_PACall_new.txt'],'','delimiter',',');
        erPAC_SO_allSub_avg{counterSub,3} = erPAC_template;
        erPAC_SO_allSub_avg{counterSub,3}.powspctrm(1,:,1:length(erPAC_template.time)) = data(:,1:length(erPAC_template.time));
        counterSub = counterSub +1;
    end

end



for idx_sub = 1 : length(erPAC_SO_allSub_avg)
    
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x1-x2';
    erPAC_SO_allSub_avg{idx_sub,4} = ft_math(cfg, erPAC_SO_allSub_avg{idx_sub,1},erPAC_SO_allSub_avg{idx_sub,2});
    erPAC_SO_allSub_avg{idx_sub,5} = ft_math(cfg, erPAC_SO_allSub_avg{idx_sub,1},erPAC_SO_allSub_avg{idx_sub,3});
    erPAC_SO_allSub_avg{idx_sub,6} = ft_math(cfg, erPAC_SO_allSub_avg{idx_sub,2},erPAC_SO_allSub_avg{idx_sub,3});
    erPAC_SO_allSub_avg{idx_sub,7} = erPAC_SO_allSub_avg{idx_sub,3};
    erPAC_SO_allSub_avg{idx_sub,7}.powspctrm = zeros(size(erPAC_SO_allSub_avg{idx_sub,3}.powspctrm,1),size(erPAC_SO_allSub_avg{idx_sub,3}.powspctrm,2),size(erPAC_SO_allSub_avg{idx_sub,3}.powspctrm,3));

    
end

cfg = [];
cfg.keepindividual = 'yes';

grdAvg_avg.erPAC.stim                   = ft_freqgrandaverage(cfg, erPAC_SO_allSub_avg{:,1});
grdAvg_avg.erPAC.random                 = ft_freqgrandaverage(cfg, erPAC_SO_allSub_avg{:,2});
grdAvg_avg.erPAC.rest                   = ft_freqgrandaverage(cfg, erPAC_SO_allSub_avg{:,3});
grdAvg_avg.erPAC.subtracted             = ft_freqgrandaverage(cfg, erPAC_SO_allSub_avg{:,4});
grdAvg_avg.erPAC.subtractedStimRest     = ft_freqgrandaverage(cfg, erPAC_SO_allSub_avg{:,5});
grdAvg_avg.erPAC.subtractedRandoRest    = ft_freqgrandaverage(cfg, erPAC_SO_allSub_avg{:,6});

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x1-x2';
grdAvg_avg.erPAC.subtractedPlot = ft_math(cfg, grdAvg_avg.erPAC.random,grdAvg_avg.erPAC.stim);

cfg = [];
cfg.layout      = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlegend  = 'yes';
cfg.xlim        = [-1 2];
cfg.zlim        = [0 0.1];
cfg.colormap    = colormap(parulaHalfClear);
ft_singleplotTFR(cfg, grdAvg_avg.erPAC.stim);set(gca,'TickDir','out');
ft_singleplotTFR(cfg, grdAvg_avg.erPAC.random);set(gca,'TickDir','out');
ft_singleplotTFR(cfg, grdAvg_avg.erPAC.rest);set(gca,'TickDir','out');


cfg.zlim        = [-0.05 0.05];
cfg.colormap    = colormap(parula);

ft_singleplotTFR(cfg, grdAvg_avg.erPAC.subtracted)
ft_singleplotTFR(cfg, grdAvg_avg.erPAC.subtractedRandoRest)
ft_singleplotTFR(cfg, grdAvg_avg.erPAC.subtractedStimRest)



load neighboursPerso.mat


cfg                     = [];
cfg.design(1,1:2*(length(erPAC_SO_allSub_avg)))  = [ones(1,(length(erPAC_SO_allSub_avg))) 2*ones(1,(length(erPAC_SO_allSub_avg)))];
cfg.design(2,1:2*(length(erPAC_SO_allSub_avg)))  = [1:(length(erPAC_SO_allSub_avg)) 1:(length(erPAC_SO_allSub_avg))];
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
cfg.alpha               = 0.025; 
cfg.clusteralpha        = 0.05;     
cfg.numrandomization    = 500;      % number of draws from the permutation distribution
cfg.latency             = [-1 2];

% 1 = associated; 2 = unassociated; 3 = rest;

[stat12] = ft_freqstatistics(cfg,  erPAC_SO_allSub_avg{:,1}, erPAC_SO_allSub_avg{:,2});
[stat13] = ft_freqstatistics(cfg,  erPAC_SO_allSub_avg{:,1}, erPAC_SO_allSub_avg{:,3});%0.025/0.05
[stat23] = ft_freqstatistics(cfg,  erPAC_SO_allSub_avg{:,2}, erPAC_SO_allSub_avg{:,3});


cfg = [];
cfg.latency             = [-1 2];
tmp = ft_selectdata(cfg,grdAvg_avg.erPAC.subtractedPlot);
tmp.mask   = stat13.prob<0.05;

cfg = [];
cfg.layout      = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlegend  = 'yes';
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.4;
cfg.maskparameter = 'mask';
cfg.zlim        = [-0.05 0.05];
figure;ft_singleplotTFR(cfg, tmp)
tmpSO = csvread([initPath.Exp 'data\group\OL\grdAvg_SO_allstage.csv']);tmpSO= tmpSO(1:end-1,2);

hold on
yyaxis right
plot(tmp.time,tmpSO','-k')
ylim([-90 30])
hold off
set(gca,'TickDir','out');

cfg = [];
cfg.channel = 'all';
cfg.latency = [-0.43 0.22];
cfg.frequency = [20.5 22];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq = 'yes';  % this "squeezes" the time dimension out of the data
unasso = ft_selectdata(cfg, grdAvg_avg.erPAC.stim);
rest  = ft_selectdata(cfg, grdAvg_avg.erPAC.rest);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(erPAC_SO_allSub_avg)))  = [ones(1,(length(erPAC_SO_allSub_avg))) 2*ones(1,(length(erPAC_SO_allSub_avg)))];
cfg.design(2,1:2*(length(erPAC_SO_allSub_avg)))  = [1:(length(erPAC_SO_allSub_avg)) 1:(length(erPAC_SO_allSub_avg))];
effect_roiAssovsRest= ft_freqstatistics(cfg, unasso, rest);
disp(effect_roiAssovsRest)


PACAverage = [];
for idx_sub = 1 : length(erPAC_SO_allSub_avg)
    
        cfg = [];
        cfg.avgovertime         = 'yes';
        cfg.avgoverfreq         = 'yes';
        cfg.avgoverchan         = 'yes';

        cfg.latency = [-0.43 0.22];
        cfg.frequency = [20.5 22];


        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,1});
        PACAverage(idx_sub,1)   = tmp.powspctrm;
        
        
end

csvwrite([initPath.Exp 'data\OL_CA\group\extract_PACSO_CBP_AssvsRest_OL_CA.csv'],PACAverage)




cfg = [];
cfg.latency             = [-1 2];
tmp = ft_selectdata(cfg,grdAvg_avg.erPAC.subtractedRandoRest);
tmp.mask   = stat23.mask;

cfg = [];
cfg.layout      = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlegend  = 'yes';
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.5;
cfg.maskparameter = 'mask';
cfg.zlim        = [-0.05 0.05];
figure;ft_singleplotTFR(cfg, tmp)
tmpSO = csvread([initPath.Exp 'data\group\OL\grdAvg_SO_allstage.csv']);tmpSO= tmpSO(1:end-1,2);

hold on
yyaxis right
plot(tmp.time,tmpSO','-k')
ylim([-90 30])
set(gca,'TickDir','out');

cfg = [];
cfg.channel = 'all';
cfg.latency = [-1 -0.02];
cfg.frequency = [9 12];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq = 'yes';  % this "squeezes" the time dimension out of the data
unasso = ft_selectdata(cfg, grdAvg_avg.erPAC.random);
rest  = ft_selectdata(cfg, grdAvg_avg.erPAC.rest);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(erPAC_SO_allSub_avg)))  = [ones(1,(length(erPAC_SO_allSub_avg))) 2*ones(1,(length(erPAC_SO_allSub_avg)))];
cfg.design(2,1:2*(length(erPAC_SO_allSub_avg)))  = [1:(length(erPAC_SO_allSub_avg)) 1:(length(erPAC_SO_allSub_avg))];
effect_roiUnassovsRest1= ft_freqstatistics(cfg, unasso, rest);
disp(effect_roiUnassovsRest1)


cfg = [];
cfg.channel = 'all';
cfg.latency = [-0.58 0.09];
cfg.frequency = [17.5 20];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq = 'yes';  % this "squeezes" the time dimension out of the data
unasso = ft_selectdata(cfg, grdAvg_avg.erPAC.random);
rest  = ft_selectdata(cfg, grdAvg_avg.erPAC.rest);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(erPAC_SO_allSub_avg)))  = [ones(1,(length(erPAC_SO_allSub_avg))) 2*ones(1,(length(erPAC_SO_allSub_avg)))];
cfg.design(2,1:2*(length(erPAC_SO_allSub_avg)))  = [1:(length(erPAC_SO_allSub_avg)) 1:(length(erPAC_SO_allSub_avg))];
effect_roiUnassovsRest2= ft_freqstatistics(cfg, unasso, rest);
disp(effect_roiUnassovsRest2)



cfg = [];
cfg.channel = 'all';
cfg.latency = [-0.49 0.09];
cfg.frequency = [26 30];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq = 'yes';  % this "squeezes" the time dimension out of the data
unasso = ft_selectdata(cfg, grdAvg_avg.erPAC.random);
rest  = ft_selectdata(cfg, grdAvg_avg.erPAC.rest);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(erPAC_SO_allSub_avg)))  = [ones(1,(length(erPAC_SO_allSub_avg))) 2*ones(1,(length(erPAC_SO_allSub_avg)))];
cfg.design(2,1:2*(length(erPAC_SO_allSub_avg)))  = [1:(length(erPAC_SO_allSub_avg)) 1:(length(erPAC_SO_allSub_avg))];
effect_roiUnassovsRest3= ft_freqstatistics(cfg, unasso, rest);
disp(effect_roiUnassovsRest3)


PACAverage = [];
for idx_sub = 1 : length(erPAC_SO_allSub_avg)
    
        cfg = [];
        cfg.avgovertime         = 'yes';
        cfg.avgoverfreq         = 'yes';
        cfg.avgoverchan         = 'yes';
        
        cfg.latency = [-1 -0.02];
        cfg.frequency = [9 12];
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,2});
        PACAverage(idx_sub,1) = tmp.powspctrm;
        
        cfg.latency = [-0.58 0.09];
        cfg.frequency = [17.5 20];
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,2});
        PACAverage(idx_sub,2) = tmp.powspctrm;
        
        cfg.latency = [-0.49 0.09];
        cfg.frequency = [26 30];
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,2});
        PACAverage(idx_sub,3) = tmp.powspctrm;
end

csvwrite([initPath.Exp 'data\OL_CA\group\extract_PACSO_CBP_UnassvsRest_OL_CA.csv'],PACAverage)



PACAverage = [];
for idx_sub = 1 : length(erPAC_SO_allSub_avg)
    
        cfg = [];
        cfg.avgovertime         = 'yes';
        cfg.avgoverfreq         = 'yes';
        cfg.avgoverchan         = 'yes';
        
        
        cfg.latency = [-0.58 0.22];
        cfg.frequency = [17.5 22];
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,1});
        PACAverage(idx_sub,1) = tmp.powspctrm;
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,2});
        PACAverage(idx_sub,2) = tmp.powspctrm;
     
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,4});
        PACAverage(idx_sub,3) = tmp.powspctrm;
end

csvwrite([initPath.Exp 'data\OL_CA\group\extract_PACSO_forAgeCorr_OL_CA.csv'],PACAverage)


%% Correlation % Results presented in 3.3
[TMRindex , reactGains , notReactGains] = getBehav_OL_CA();
TMRindex = TMRindex([1:3 5:11 14:17]);% remove participants with no SW
reactGains = reactGains([1:3 5:11 14:17]);
notReactGains = notReactGains([1:3 5:11 14:17]);

% Correlation between SO PAC and behaviour
% TMR index both early and late with TFR relevant - TFR irrelevant
design=[];
design(1,1:length(erPAC_SO_allSub_avg))       = TMRindex;

cfg=[];
cfg.statistic           = 'ft_statfun_correlationT';
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.neighbours          = neighbours_perso; 
cfg.minnbchan           = 0;              
cfg.alpha               = 0.025; %0.25
cfg.clusteralpha        = 0.025;     %0.25
cfg.latency             = [-1 2];
% cfg.frequency           = [10 20];
cfg.numrandomization    = 500;
cfg.design              = design;
cfg.ivar                = 1;

statCorAssovsUnasso = ft_freqstatistics(cfg, erPAC_SO_allSub_avg{:,4});
statCorAssovsRest = ft_freqstatistics(cfg, erPAC_SO_allSub_avg{:,5});
statCorUnassovsRest = ft_freqstatistics(cfg, erPAC_SO_allSub_avg{:,6});

cfg = [];
cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive   = 'yes';
cfg.showoutline   = 'yes';  
cfg.showlegend    = 'yes';
cfg.avgoverchan   = 'no';
cfg.xlim          = [-1 2];
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.5;
cfg.maskparameter = 'mask';
cfg.parameter = 'rho';
cfg.zlim = [-1 1];
figure;ft_singleplotTFR(cfg, statCorAssovsRest)
figure;ft_singleplotTFR(cfg, statCorUnassovsRest)
tmpSO = csvread([initPath.Exp 'data\group\OL\grdAvg_SO_allstage.csv']);tmpSO= tmpSO(1:end-1,2);

hold on
yyaxis right
plot(tmp.time,tmpSO','-k')
ylim([-90 30])
set(gca,'TickDir','out');







PACAverage = [];
for idx_sub = 1 : length(erPAC_SO_allSub_avg)
    
        cfg = [];
        cfg.avgovertime         = 'yes';
        cfg.avgoverfreq         = 'yes';
        cfg.avgoverchan         = 'yes';
        %For Cor with TMR
        cfg.latency             = [0.47 1.7];
        cfg.frequency           = [23 28];
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,4});
        PACAverage(idx_sub,1) = tmp.powspctrm;
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,1});
        PACAverage(idx_sub,2) = tmp.powspctrm;
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,2});
        PACAverage(idx_sub,3) = tmp.powspctrm;
        
        cfg.latency             = [1.14 2];
        cfg.frequency           = [13 17];
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,4});
        PACAverage(idx_sub,4) = tmp.powspctrm;
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,1});
        PACAverage(idx_sub,5) = tmp.powspctrm;
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub_avg{idx_sub,2});
        PACAverage(idx_sub,6) = tmp.powspctrm;

end

csvwrite([initPath.Exp 'data\OL_CA\group\extract_PACSO_CBPcor_TMRvsPACdiff_OL_CA.csv'],PACAverage)



% Offline gain in react both early and late with TFR relevant 
design=[];
design(1,1:length(erPAC_SO_allSub_avg))       = reactGains;

cfg=[];
cfg.statistic           = 'ft_statfun_correlationT';
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.neighbours          = neighbours_perso; 
cfg.minnbchan           = 0;              
cfg.alpha               = 0.025; 
cfg.clusteralpha        = 0.05;     
cfg.latency             = [-1 2];
cfg.numrandomization    = 500;
cfg.design              = design;
cfg.ivar                = 1;

statCorReact = ft_freqstatistics(cfg, erPAC_SO_allSub_avg{:,1});

cfg = [];
cfg.latency             = [-1 2];
tmp = ft_selectdata(cfg,grdAvg_avg.erPAC.stim);
tmp.mask   = statCorReact.prob<0.1;

cfg = [];
cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive   = 'yes';
cfg.showoutline   = 'yes';  
cfg.showlegend    = 'yes';
cfg.avgoverchan   = 'no';
cfg.xlim          = [-1 2];
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.5;
cfg.maskparameter = 'mask';
% cfg.parameter = 'rho';
% cfg.zlim = [-1 1];
figure;ft_singleplotTFR(cfg, tmp)
tmpSO = csvread([initPath.Exp 'data\group\OL\grdAvg_SO_allstage.csv']);tmpSO= tmpSO(1:end-1,2);

hold on
yyaxis right
plot(tmp.time,tmpSO','-k')
ylim([-90 30])
set(gca,'TickDir','out');

% Offline gain in not react both early and late with TFR irrelevant 
design=[];
design(1,1:length(erPAC_SO_allSub_avg))       = notReactGains;

cfg=[];
cfg.statistic           = 'ft_statfun_correlationT';
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.neighbours          = neighbours_perso; 
cfg.minnbchan           = 0;              
cfg.alpha               = 0.05; 
cfg.clusteralpha        = 0.025;     
cfg.latency             = [-1 2];
cfg.numrandomization    = 500;
cfg.design              = design;
cfg.ivar                = 1;

statCornotReact = ft_freqstatistics(cfg, erPAC_SO_allSub_avg{:,2});



cfg = [];
cfg.latency             = [-1 2];
tmp = ft_selectdata(cfg,grdAvg_avg.erPAC.random);
tmp.mask   = statCornotReact.prob<0.05;

cfg = [];
cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive   = 'yes';
cfg.showoutline   = 'yes';  
cfg.showlegend    = 'yes';
cfg.avgoverchan   = 'no';
cfg.xlim          = [-1 2];
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.5;
cfg.maskparameter = 'mask';
% cfg.parameter = 'rho';
% cfg.zlim = [-1 1];
figure;ft_singleplotTFR(cfg, tmp)
tmpSO = csvread([initPath.Exp 'data\group\OL\grdAvg_SO_allstage.csv']);tmpSO= tmpSO(1:end-1,2);

hold on
yyaxis right
plot(tmp.time,tmpSO','-k')
ylim([-90 30])
