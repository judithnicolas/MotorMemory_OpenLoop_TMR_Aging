%% 

% Results presented in 3.2.3 and Figure S6 and S7 (supp)

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
PP_ER_rel= readtable([initPath.Exp '\data\OL_CA\group\PP_ER_rel_avgOverChan.txt'], opts);

% Convert to output type
PP_ER_rel = table2array(PP_ER_rel);
PP_ER_irrel = table2array(PP_ER_irrel);

clear opts

figure; circ_plot(PP_ER_rel,'hist',[],15,true,true,'linewidth',2,'color','r'); 
[pval z] =circ_rtest(PP_ER_rel)
figure; circ_plot(PP_ER_irrel,'hist',[],15,true,true,'linewidth',2,'color','r'); 
[pval z] =circ_rtest(PP_ER_irrel)

[pval table] = circ_wwtest(PP_ER_rel, PP_ER_irrel)
[mu_ER ul_ER ll_ER] = circ_mean([PP_ER_rel;PP_ER_irrel]);
rad2deg([mu_ER ul_ER ll_ER])+360



%% Cue-Locked related PAC 


listSub = getScoredDatasets_CA;
erPAC_allSub = {};

for idx_sub = 1 : length(listSub)
    sub = listSub{idx_sub};

    if idx_sub == 1
        
        load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_preprocessed_continuous.mat'])
        load ([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_trl_epoch.mat'])
      
        trl_corrected = trl;
        trl_corrected(:,1:2) = trl(:,1:2) + 0.44*data.fsample;
        
        
        cfg = [];
        cfg.trl      = trl_corrected(find(trl_corrected(:,5) == 2 | trl_corrected(:,5) == 3),1:3);
        data = ft_redefinetrial(cfg, data);

        cfg         = [];
        cfg.channel = {'Fz'};
        data = ft_selectdata(cfg,data);
        
        data.trialinfo= trl(:,4:5);

        cfg = [];
        cfg.method     = 'mtmconvol';
        cfg.pad        = 'nextpow2';
        cfg.taper      = 'hanning';
        cfg.keeptrials = 'no';
        cfg.foi = 7.25:0.5:29.25; 
        cfg.toi = -1:0.002:3;
        cfg.t_ftimwin  = 5./cfg.foi;
        cfg.tapsmofrq  = 0.4 *cfg.foi;

        cfg.trials     = 1;
       
        % er PAC
        erPAC_template  = ft_freqanalysis(cfg, data);
        erPAC_template.freq  = 7.25:0.5:29.25; 

        erPAC_template.powspctrm= [];


    end

    
    
    data = textread([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_pacER_rel_avgOverChan.txt'],'','delimiter',',');
    erPAC_allSub{idx_sub,1} = erPAC_template;
    erPAC_allSub{idx_sub,1}.powspctrm(1,:,:) = data;
    
    data = textread([initPath.Exp '\data\OL_CA\' sub '\exp\' sub '_pacER_irrel_avgOverChan.txt'],'','delimiter',',');
    erPAC_allSub{idx_sub,2} = erPAC_template;
    erPAC_allSub{idx_sub,2}.powspctrm(1,:,:) = data;
 

end



for idx_sub = 1 : length(erPAC_allSub)

        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.operation = 'x1-x2';
        erPAC_allSub{idx_sub,3} = ft_math(cfg, erPAC_allSub{idx_sub,1},erPAC_allSub{idx_sub,2});
        erPAC_allSub{idx_sub,4}= erPAC_template;
        erPAC_allSub{idx_sub,4}.powspctrm = zeros(size(erPAC_allSub{idx_sub,3}.powspctrm,1),size(erPAC_allSub{idx_sub,3}.powspctrm,2),size(erPAC_allSub{idx_sub,3}.powspctrm,3));

end


cfg = [];
grdAvg.erPAC.stim = ft_freqgrandaverage(cfg, erPAC_allSub{:,1});
grdAvg.erPAC.random = ft_freqgrandaverage(cfg, erPAC_allSub{:,2});

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x1-x2';

grdAvg.erPAC.subtracted = ft_math(cfg, grdAvg.erPAC.stim,grdAvg.erPAC.random);


load neighboursPerso.mat


cfg                     = [];
cfg.design(1,1:2*length(listSub))  = [ones(1,length(listSub)) 2*ones(1,length(listSub))];
cfg.design(2,1:2*length(listSub))  = [1:length(listSub) 1:length(listSub)];
cfg.ivar                = 1; 
cfg.uvar                = 2;
cfg.method              = 'montecarlo';       
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.minnbchan           = 0;              
cfg.avgoverchan         = 'yes';
cfg.neighbours          = neighbours_perso; 
cfg.tail                = 0;                    
cfg.clustertail         = 0;
cfg.alpha               = 0.025; 
cfg.clusteralpha        = 0.05;     
cfg.numrandomization    = 500;      

cfg.latency             = [0 2.5];
% 1 = associated; 2 = unassociated; 
[statERPAC] = ft_freqstatistics(cfg,  erPAC_allSub{:,1}, erPAC_allSub{:,2});


% Restults presented in 3.3

[TMRindex , reactGains , notReactGains] = getBehav_OL_CA();

% Correlation between ER PAC and behaviour
% TMR index both early and late with PAC associated - PAC unassociated 
design=[];
design(1,1:length(listSub))       = TMRindex;
cfg                     = [];
cfg.statistic           = 'ft_statfun_correlationT';
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.neighbours          = neighbours_perso; 
cfg.minnbchan           = 0;              
cfg.alpha               = 0.025; 
cfg.clusteralpha        = 0.05;     
cfg.latency             = [0 2.5];
cfg.numrandomization    = 500;
cfg.avgoverchan         = 'yes';
cfg.design              = design;
cfg.ivar                = 1;

statERPAC_cor = ft_freqstatistics(cfg, erPAC_allSub{:,3});
