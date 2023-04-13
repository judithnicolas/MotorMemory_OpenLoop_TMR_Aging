%Utilities
parulaHalfClear = [0.0963952380952380,0.750000000000000,0.712038095238095;0.140771428571429,0.758400000000000,0.684157142857143;0.171700000000000,0.766961904761905,0.655442857142857;0.193766666666667,0.775766666666667,0.625100000000000;0.216085714285714,0.784300000000000,0.592300000000000;0.246957142857143,0.791795238095238,0.556742857142857;0.290614285714286,0.797290476190476,0.518828571428572;0.340642857142857,0.800800000000000,0.478857142857143;0.390900000000000,0.802871428571429,0.435447619047619;0.445628571428572,0.802419047619048,0.390919047619048;0.504400000000000,0.799300000000000,0.348000000000000;0.561561904761905,0.794233333333333,0.304480952380953;0.617395238095238,0.787619047619048,0.261238095238095;0.671985714285714,0.779271428571429,0.222700000000000;0.724200000000000,0.769842857142857,0.191028571428571;0.773833333333333,0.759804761904762,0.164609523809524;0.820314285714286,0.749814285714286,0.153528571428571;0.863433333333333,0.740600000000000,0.159633333333333;0.903542857142857,0.733028571428571,0.177414285714286;0.939257142857143,0.728785714285714,0.209957142857143;0.972757142857143,0.729771428571429,0.239442857142857;0.995647619047619,0.743371428571429,0.237147619047619;0.996985714285714,0.765857142857143,0.219942857142857;0.995204761904762,0.789252380952381,0.202761904761905;0.989200000000000,0.813566666666667,0.188533333333333;0.978628571428571,0.838628571428572,0.176557142857143;0.967647619047619,0.863900000000000,0.164290476190476;0.961009523809524,0.889019047619048,0.153676190476191;0.959671428571429,0.913457142857143,0.142257142857143;0.962795238095238,0.937338095238095,0.126509523809524;0.969114285714286,0.960628571428571,0.106361904761905;0.976900000000000,0.983900000000000,0.0805000000000000];
%%
listSub = getScoredDatasets_CA;
erPAC_SO_allSub = {};
Cond = {'rel' 'irrel' 'rest'};

counterSub = 1;

for idx_sub = 1: length(listSub)
    sub = listSub{idx_sub};
    if ~ismember(idx_sub ,[4  12:13 ])
        
        if idx_sub == 1
            load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_preprocessed_continuous.mat'])
            load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trl_epoch.mat'])
            
            cfg         = [];
            cfg.channel = [1:6];%{'Fz'};
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
            cfg.toi = -3:0.002:2.999;
            cfg.t_ftimwin  = 5./cfg.foi;
            cfg.tapsmofrq  = 0.4 *cfg.foi;
            
            % er PAC
            erPAC_template  = ft_freqanalysis(cfg, data);
            erPAC_template.freq  = 7.25:0.5:29.25;
            erPAC_template.time= -1:0.002:1.999;
            
            erPAC_template.powspctrm= [];
            
        end
        
        for idx_condition =  1 : length(Cond)
            erPAC_SO_allSub{counterSub,idx_condition} = erPAC_template;
            for idx_channel = 1 : length(erPAC_template.label)
                
                if ~isempty(dir([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_pacERSO_' Cond{idx_condition} '_detectionFz_' erPAC_template.label{idx_channel} '.txt']))
                    data = textread([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_pacERSO_' Cond{idx_condition} '_detectionFz_' erPAC_template.label{idx_channel} '.txt'],'','delimiter',',');
                    erPAC_SO_allSub{counterSub,idx_condition}.powspctrm(idx_channel,:,1:length(erPAC_template.time)) = data(:,1:length(erPAC_template.time));
                else
                end
            end
            indx_4_avg(idx_sub,idx_condition) = 1;
            
        end
    end
    counterSub = counterSub +1;

end

for idx_sub = 1 : length(erPAC_SO_allSub)
    
        
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.operation = 'x1-x2';
        erPAC_SO_allSub{idx_sub,4} = ft_math(cfg, erPAC_SO_allSub{idx_sub,1},erPAC_SO_allSub{idx_sub,2});
        erPAC_SO_allSub{idx_sub,5} = ft_math(cfg, erPAC_SO_allSub{idx_sub,1},erPAC_SO_allSub{idx_sub,3});
        erPAC_SO_allSub{idx_sub,6} = ft_math(cfg, erPAC_SO_allSub{idx_sub,2},erPAC_SO_allSub{idx_sub,3});

end

cfg = [];
% cfg.keepindividual = 'yes';

grdAvg.erPAC.stim                   = ft_freqgrandaverage(cfg, erPAC_SO_allSub{:,1});
grdAvg.erPAC.random                 = ft_freqgrandaverage(cfg, erPAC_SO_allSub{:,2});
grdAvg.erPAC.rest                   = ft_freqgrandaverage(cfg, erPAC_SO_allSub{:,3});

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x1-x2';

grdAvg.erPAC.subtracted             = ft_math(cfg, grdAvg.erPAC.stim,grdAvg.erPAC.random);
grdAvg.erPAC.subtractedStimRest     = ft_math(cfg, grdAvg.erPAC.stim,grdAvg.erPAC.random);
grdAvg.erPAC.subtractedRandoRest    = ft_math(cfg, grdAvg.erPAC.random,grdAvg.erPAC.rest);



%%
cfg = [];
cfg.layout      = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive = 'yes';
% cfg.showoutline = 'yes';
cfg.showlegend  = 'yes';
cfg.xlim        = [-1 2];
cfg.zlim        = [0 0.1];
cfg.colormap    = colormap(parulaHalfClear);


figure;
cfg.title = 'stim';
ft_multiplotTFR(cfg, grdAvg.erPAC.stim)
figure;
cfg.title = 'random';
ft_multiplotTFR(cfg, grdAvg.erPAC.random)
figure;
cfg.title = 'rest';
ft_multiplotTFR(cfg, grdAvg.erPAC.rest)


for idx_channel = 1 : length(grdAvg.erPAC.rest.label)
    
    cfg = [];
    cfg.channel= grdAvg.erPAC.rest.label{idx_channel};%'all'
    cfg.layout      = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
    cfg.interactive = 'yes';
    % cfg.showoutline = 'yes';
    cfg.showlegend  = 'yes';
    cfg.xlim        = [-1 2];
    cfg.zlim        = [0 0.1];
    cfg.colormap    = colormap(parulaHalfClear);
    
    figure;
    cfg.title = ['stim ' grdAvg.erPAC.rest.label{idx_channel};];
    ft_singleplotTFR(cfg, grdAvg.erPAC.stim)
    set(gca,'TickDir','out');
    
    figure;
    cfg.title = ['random '  grdAvg.erPAC.rest.label{idx_channel};];
    ft_singleplotTFR(cfg, grdAvg.erPAC.random)
    set(gca,'TickDir','out');

    figure;
    cfg.title = ['rest '  grdAvg.erPAC.rest.label{idx_channel};];
    ft_singleplotTFR(cfg, grdAvg.erPAC.rest)
    set(gca,'TickDir','out');

end

cfg.zlim        = [-0.2 0.2];
figure;
ft_multiplotTFR(cfg, grdAvg.erPAC.subtracted)
figure;
ft_multiplotTFR(cfg, grdAvg.erPAC.subtractedRandoRest)
figure;
ft_multiplotTFR(cfg, grdAvg.erPAC.subtractedStimRest)


% cfg.channel  = 'yes';
% cfg.avgoverchan  = 'yes';
cfg.zlim        = [0 0.05];
figure;
cfg.title = 'stim';
ft_singleplotTFR(cfg, grdAvg.erPAC.stim)
figure;
cfg.title = 'random';
ft_singleplotTFR(cfg, grdAvg.erPAC.random)
figure;
cfg.title = 'rest';
ft_singleplotTFR(cfg, grdAvg.erPAC.rest)

load neighboursPerso.mat


cfg                     = [];
cfg.design(1,1:2*(length(erPAC_SO_allSub)))  = [ones(1,(length(erPAC_SO_allSub))) 2*ones(1,(length(erPAC_SO_allSub)))];
cfg.design(2,1:2*(length(erPAC_SO_allSub)))  = [1:(length(erPAC_SO_allSub)) 1:(length(erPAC_SO_allSub))];
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



[stat12] = ft_freqstatistics(cfg,  erPAC_SO_allSub{:,1}, erPAC_SO_allSub{:,2});
[stat13] = ft_freqstatistics(cfg,  erPAC_SO_allSub{:,1}, erPAC_SO_allSub{:,3});
[stat23] = ft_freqstatistics(cfg,  erPAC_SO_allSub{:,2}, erPAC_SO_allSub{:,3});

cfg = [];
cfg.latency             = [-1 2];
cfg.channel             = [1:3 5:6];    
tmp = ft_selectdata(cfg,grdAvg.erPAC.subtractedRandoRest);
tmp.mask   = stat23.mask;

cfg = [];
cfg.layout      = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive = 'yes';
cfg.showlegend  = 'yes';
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.3;
cfg.maskparameter = 'mask';
figure;ft_multiplotTFR(cfg, tmp)



cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSub)-2))  = [ones(1,length(listSub)-2) 2*ones(1,length(listSub)-2)];
cfg.design(2,1:2*(length(listSub)-2))  = [1:length(listSub)-2 1:length(listSub)-2];
effect_roiStimvsRamdom = ft_freqstatistics(cfg, stim, random);
disp(effect_roiStimvsRamdom)

cfg = [];
cfg.channel = 'all';
cfg.latency = [-1 0.5];
cfg.frequency = [13.5 20];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq = 'yes';  % this "squeezes" the time dimension out of the data
random = ft_selectdata(cfg, grdAvg.erPAC.random);
rest  = ft_selectdata(cfg, grdAvg.erPAC.rest);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSub)-2))  = [ones(1,length(listSub)-2) 2*ones(1,length(listSub)-2)];
cfg.design(2,1:2*(length(listSub)-2))  = [1:length(listSub)-2 1:length(listSub)-2];
effect_roiRandomvsRest= ft_freqstatistics(cfg, random, rest);
disp(effect_roiRandomvsRest)

%%
[TMRindex , reactGains , notReactGains] = getBehav_OL_CA();
TMRindex = TMRindex([1:3 5:11 14:17]);
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
cfg.alpha               = 0.025; 
cfg.clusteralpha        = 0.025;     
cfg.latency             = [-1 2];
cfg.numrandomization    = 500;
cfg.design              = design;
cfg.ivar                = 1;

statCor = ft_freqstatistics(cfg, erPAC_SO_allSub{:,4});


cfg = [];
cfg.latency             = [-1 2];
tmp = ft_selectdata(cfg,grdAvg.erPAC.subtracted);
statCor.mask   = statCor.prob<0.026;

cfg = [];
cfg.layout      = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlegend  = 'yes';
cfg.maskstyle     = 'opacity';
cfg.maskalpha        = 0.4;
cfg.maskparameter = 'mask';
cfg.parameter = 'rho';
cfg.zlim        = [-1 1];

figure;ft_multiplotTFR(cfg, statCor)


figure;ft_singleplotTFR(cfg, statCor)
tmpSO = csvread('data\group\OL\grdAvg_SO_allstage.csv');tmpSO= tmpSO(:,2);

hold on
yyaxis right
plot(tmp.time,tmpSO','-k')
ylim([-90 30])

cfg = [];
cfg.latency             = [0.5 1];
cfg.frequency           = [14 17];
cfg.avgoverfreq         = 'yes';
cfg.avgovertime         = 'yes';
tmp = ft_selectdata(cfg,statCor);


PACAverage = [];
for idx_sub = 1 : length(listSub)
    
        cfg = [];
        cfg.avgovertime         = 'yes';
        cfg.avgoverfreq         = 'yes';
        cfg.avgoverchan         = 'yes';
        
        cfg.latency             = [0.68 0.82];
        cfg.frequency           = [14.5 17];
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub{idx_sub,4});
        PACAverage(idx_sub,1) = tmp.powspctrm;
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub{idx_sub,1});
        PACAverage(idx_sub,2) = tmp.powspctrm;
        
        tmp                     = ft_selectdata(cfg,erPAC_SO_allSub{idx_sub,2});
        PACAverage(idx_sub,3) = tmp.powspctrm;
end

csvwrite([initPath.Exp 'data\group\OL\TMRvsPAC.csv'],PACAverage)


% Offline gain in react both early and late with TFR relevant 
design(1,1:length(listSub)-2)       = [  0.2116431234  0.3405858412  0.3214817219 ...
    0.0004167532  0.1442537994  0.2172953504 -0.0225442508  0.3362848209 ...
   -0.0844974713  0.1020246232  0.2200232825  ...
 0.0359798153  0.1039230677  0.2935166032   0.1323151041  0.1695819974  ...
 0.1177051794  0.0511158980  0.0908001297  0.1864793652  0.1650714297 -0.0067169097]; 

cfg.design           = design;

statCor = ft_freqstatistics(cfg, erPAC_SO_allSub{:,1});

% Offline gain in not react both early and late with TFR irrelevant 
design(1,1:length(listSub)-2)       =  [ -0.36279555  0.19205582  0.28019562 ...
    -0.01665761  0.16103859  0.01546686 -0.26225657  0.40485138 ...
  0.06620860  0.09779284  0.11531575  -0.02538434  0.09021207  ...
  0.16310537 -0.11078920  0.05585621 0.25243604  0.10025795  0.01879572  ...
  0.23137372  0.02086235  0.13766224]; 

cfg.design           = design;

statCor = ft_freqstatistics(cfg, erPAC_SO_allSub{:,2});



PACAverage = [];
for idx_sub = 1 : length(listSub)

    cfg = [];
    cfg.avgovertime         = 'yes';
    cfg.avgoverfreq         = 'yes';
    cfg.avgoverchan         = 'yes';
    
    cfg.latency             = [0.5 1];
    cfg.frequency           = [14 18];
    tmp                     = ft_selectdata(cfg,erPAC_SO_allSub{idx_sub,1});
    PACAverage(idx_sub,1) = tmp.powspctrm;

    tmp                     = ft_selectdata(cfg,erPAC_SO_allSub{idx_sub,2});
    PACAverage(idx_sub,2) = tmp.powspctrm;

    tmp                     = ft_selectdata(cfg,erPAC_SO_allSub{idx_sub,4});
    PACAverage(idx_sub,3) = tmp.powspctrm;
    
end