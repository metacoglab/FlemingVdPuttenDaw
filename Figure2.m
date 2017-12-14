% Plots performance and confidence data together with posterior predictive
% simulations from STAN model fits
% Reproduces Figure 2 in Fleming, van der Putten & Daw
%
% Steve Fleming stephen.fleming@ucl.ac.uk 2017
clear all
close all
fs = filesep;

dataset = input('Which dataset? 2=behav, 3=fmri ');
model = 'ideal_dt';

Ndraws = 1000;    % how many draws from subject-level posteriors / simulated trial sequences to average over
conf_sigma = 0.025; % Fixed parameter relating model to observed confidence

cwd = pwd;
baseDir = '~/Dropbox/Research/Metacognition/stateactionexpt/FlemingVdPuttenDaw/stan/modelfits/indiv';
dirData = '~/Dropbox/Research/Metacognition/stateactionexpt/FlemingVdPuttenDaw/data';
fitsDir = [baseDir fs model fs];

if dataset == 2
    filename = 'fMRI_pilotData_sub_';    % use for fMRI
    suffix = ''; % use for fMRI data out-of-scanner
    subjects = [12:28 30:37];
elseif dataset == 3
    filename = 'fMRI_pilotData_sub_';    % use for fMRI
    suffix = '_fMRI';
    subjects = [12:19 23:28 30:37];
end

%% Sort subject data
allConf_cor = cell(1,9);   % stores aggregate confidence data, 9 x (Ntrials*Nsubjects) matrix; slowest changing factor is post
allConf_err = cell(1,9);
allModelConf_cor = cell(1,9);
allModelConf_err = cell(1,9);

for s = 1:length(subjects)

    %% Load data for this subject
    datafile = [filename num2str(subjects(s)) suffix '_2.mat'];
    cd(dirData);
    load(datafile);
    cd(cwd);
    
    precoh_index = [];
    postcoh_index = [];
    
    precoh = locDATA.dots_coherence';
    postcoh = locDATA.post_coherence';
    dir = locDATA.dots_direction/360;
    dir(dir==0.5) = -1;
    action = locDATA.button_response - 1;
    conf = locDATA.mouse_response;
    logrt = log(locDATA.reaction_time_button);
    conf_rt = log(locDATA.reaction_time_mouse);
    transformed_action = action;
    transformed_action(action == 0) = -1;
    acc = dir == transformed_action;
    
    % Add indicator variables for pre/post confidence
    coherence = unique(precoh);
    
    for i = 1:3
        precoh_index(locDATA.dots_coherence==coherence(i))=i;
    end
    for i = 1:3
        postcoh_index(locDATA.post_coherence==coherence(i))=i;
    end
    
    %% check for missed responses
    err_type1(s) = sum(isnan(action));
    err_type2(s) = sum(isnan(conf) & ~isnan(action));
    err_tot(s) = sum(isnan(action) | isnan(conf));
    
    %% Confidence/performance
    % 3x3
    j=1;
    for post = 1:3
        for pre = 1:3
            allConf_cor{j} = [allConf_cor{j} conf(precoh_index == pre & postcoh_index == post & acc == 1)];
            allConf_err{j} = [allConf_err{j} conf(precoh_index == pre & postcoh_index == post & acc == 0)];
            perf(s,j) = nanmean(acc(precoh_index == pre & postcoh_index == post))*100;
            conf_cor(s,j) = nanmean(conf(precoh_index == pre & postcoh_index == post & acc == 1));
            conf_err(s,j) = nanmean(conf(precoh_index == pre & postcoh_index == post & acc == 0));
            j=j+1;
        end
    end
    
    % Marginals
    for p = 1:3
        mean_conf_pre_cor(p, s) = nanmean(conf(precoh_index == p & acc == 1));
        mean_conf_pre_err(p, s) = nanmean(conf(precoh_index == p & acc == 0));
        mean_conf_post_cor(p, s) = nanmean(conf(postcoh_index == p & acc == 1));
        mean_conf_post_err(p, s) = nanmean(conf(postcoh_index == p & acc == 0));
    end
    
    %% Load model parameters
    modelfile = [model '_sub' num2str(subjects(s)) '_dataset' num2str(dataset) '_fit.mat'];
    cd(fitsDir)
    load(modelfile)
    cd(cwd)
    
    model_conf_vec = [];
    model_perf_vec = [];
    model_precoh_vec = [];
    model_postcoh_vec = [];
    %% Generate model data and point log-likelihoods from subject-level parameters using actual trial sequence
    for rep = 1:Ndraws
        % draw parameters for this replication
        k1_rep = normrnd(sub_param.k1, sub_param.k1_sd);
        m_rep = normrnd(sub_param.m, sub_param.m_sd);
        switch model
            case {'ideal', 'random'} % this is a hack, parameters have no influence for random
                sim = metaModelFit_generate_trialSequence(dir, action, conf, coherence, precoh_index, postcoh_index, conf_sigma, model, 'k1', k1_rep, 'm', m_rep);
            case {'choicebias'}
                b_rep = normrnd(sub_param.w, sub_param.w_sd);
                sim = metaModelFit_generate_trialSequence(dir, action, conf, coherence, precoh_index, postcoh_index, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'b', b_rep);
            case {'weighted', 'choiceweighted'}
                w1_rep = normrnd(sub_param.w1, sub_param.w1_sd);
                w2_rep = normrnd(sub_param.w2, sub_param.w2_sd);
                sim = metaModelFit_generate_trialSequence(dir, action, conf, coherence, precoh_index, postcoh_index, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'w1', w1_rep, 'w2', w2_rep);
            case 'mapping'
                gamma_rep = normrnd(sub_param.gamma, sub_param.gamma_sd);
                sim = metaModelFit_generate_trialSequence(dir, action, conf, coherence, precoh_index, postcoh_index, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'gamma', gamma_rep);
            case 'ideal_dt'
                brt_rep = normrnd(sub_param.brt, sub_param.brt_sd);
                sim = metaModelFit_generate_trialSequence(dir, action, conf, coherence, precoh_index, postcoh_index, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'brt', brt_rep, 'logrt', logrt);
            case {'choicebias_dt'}
                b_rep = normrnd(sub_param.w, sub_param.w_sd);
                brt_rep = normrnd(sub_param.brt, sub_param.brt_sd);
                sim = metaModelFit_generate_trialSequence(dir, action, conf, coherence, precoh_index, postcoh_index, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'b', b_rep, 'brt', brt_rep, 'logrt', logrt);
            case {'weighted_dt', 'choiceweighted_dt'}
                w1_rep = normrnd(sub_param.w1, sub_param.w1_sd);
                w2_rep = normrnd(sub_param.w2, sub_param.w2_sd);
                brt_rep = normrnd(sub_param.brt, sub_param.brt_sd);
                sim = metaModelFit_generate_trialSequence(dir, action, conf, coherence, precoh_index, postcoh_index, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'w1', w1_rep, 'w2', w2_rep, 'brt', brt_rep, 'logrt', logrt);
            case 'mapping_dt'
                gamma_rep = normrnd(sub_param.gamma, sub_param.gamma_sd);
                brt_rep = normrnd(sub_param.brt, sub_param.brt_sd);
                sim = metaModelFit_generate_trialSequence(dir, action, conf, coherence, precoh_index, postcoh_index, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'gamma', gamma_rep, 'brt', brt_rep, 'logrt', logrt);
        end
        % Store matrix of model output with rows = samples, columns =
        % trials
        model_conf_vec = [model_conf_vec sim.modelConf];
        model_perf_vec = [model_perf_vec sim.modelAcc];
        model_precoh_vec = [model_precoh_vec precoh_index];
        model_postcoh_vec = [model_postcoh_vec postcoh_index];
        
    end
    %% Average over simulations
    disp(['Finished simulation for model ' model ', subject ' num2str(s)])
    
    % Get model predictions by condition per subject
    j=1;
    for post = 1:3
        for pre = 1:3
            model_conf_cor(s, j) = nanmean(model_conf_vec(model_perf_vec == 1 & model_precoh_vec == pre & model_postcoh_vec == post));
            model_conf_err(s, j) = nanmean(model_conf_vec(model_perf_vec == 0 & model_precoh_vec == pre & model_postcoh_vec == post));
            model_perf(s, j) = nanmean(model_perf_vec(model_precoh_vec == pre & model_postcoh_vec == post)).*100;
            allModelConf_cor{j} = [allModelConf_cor{j} model_conf_vec(model_perf_vec == 1 & model_precoh_vec == pre & model_postcoh_vec == post)];
            allModelConf_err{j} = [allModelConf_err{j} model_conf_vec(model_perf_vec == 0 & model_precoh_vec == pre & model_postcoh_vec == post)];
            j=j+1;
        end
    end
    
end

% Ensure that NaN confidence entries in behav due to ceiling performance
% are matched in model output
temp_nan = isnan(conf_err);
model_conf_err(temp_nan) = NaN;

%% PLOTS
%
%% 3 x 3 performance
h1a = figure;
xpos = [1.2 2 2.8 4.2 5 5.8 7.2 8 8.8];
set(gcf, 'Position', [200 200 500 300])
mean_perf = mean(perf);
sem_perf = std(perf)./sqrt(length(subjects));
ci_perf = sem_perf.*tinv(0.975,length(subjects)-1);
mean_model_perf = mean(model_perf);
sem_model_perf = std(model_perf)./sqrt(length(subjects));
ci_model_perf = sem_perf.*tinv(0.975,length(subjects)-1);
% Subjects
boxplot(perf, 'positions', xpos, 'outliersize', 8, 'symbol', 'o', 'widths', 0.3, 'jitter', 0, 'colors', 'k')
hold on
% Model
for group = 1:3
    % Add error bars
    x = xpos((group*3)-2:group*3);
    mu = mean_model_perf((group*3)-2:group*3);
    se = ci_model_perf((group*3)-2:group*3);
    for p = 1:length(x)-1
        pX = [x(p) x(p+1) x(p+1) x(p)];
        pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
        h = fill(pX,pY,[0.5 0.5 0.5]);
        set(h,'edgecolor',[0.5 0.5 0.5]);
    end
    plot(x, mu, 'k-', 'LineWidth', 1.5)
end
boxplot(perf, 'positions', xpos, 'outliersize', 8, 'symbol', 'o', 'widths', 0.3, 'jitter', 0, 'colors', 'k')
set(gca, 'YLim', [50 100], 'XLim', [0.5 9.5], 'XTick', xpos, 'XTickLabel', {'PreL', 'PreM', 'PreH'}, 'FontSize', 14);
text(xpos(2)-0.4, 45, 'PostL', 'FontSize', 14)
text(xpos(5)-0.4, 45, 'PostM', 'FontSize', 14)
text(xpos(8)-0.4, 45, 'PostH', 'FontSize', 14)
ylabel('Performance (% correct)', 'FontSize', 18);
box off

%% 3 x 3 confidence
h1b = figure;
xpos = [1.2 2 2.8 4.2 5 5.8 7.2 8 8.8];
set(gcf, 'Position', [200 200 500 300])
hold on

mean_model_conf_cor = nanmean(model_conf_cor);
sem_model_conf_cor = nanstd(model_conf_cor)./sqrt(length(subjects));
ci_model_conf_cor = sem_model_conf_cor.*tinv(0.975,length(subjects)-1);
mean_model_conf_err = nanmean(model_conf_err);
sem_model_conf_err = nanstd(model_conf_err)./sqrt(length(subjects));
ci_model_conf_err = sem_model_conf_err.*tinv(0.975,length(subjects)-1);

% Corrects - subjects
boxplot(conf_cor, 'positions', xpos, 'outliersize', 8, 'symbol', 'o', 'widths', 0.3, 'jitter', 0, 'colors', 'g')
% Corrects - model
for group = 1:3
    x = xpos((group*3)-2:group*3);
    mu = mean_model_conf_cor((group*3)-2:group*3);
    se = ci_model_conf_cor((group*3)-2:group*3);
    for p = 1:length(x)-1
        pX = [x(p) x(p+1) x(p+1) x(p)];
        pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
        h = fill(pX,pY,[0.5 1 0.5]);
        set(h,'edgecolor',[0.5 1 0.5]);
    end
    plot(x, mu, 'g-', 'LineWidth', 1.5)
end
% Corrects - subjects (2nd time around to ensure is displayed on top)
boxplot(conf_cor, 'positions', xpos, 'outliersize', 8, 'symbol', 'o', 'widths', 0.3, 'jitter', 0, 'colors', 'g')
% Errors - model
for group = 1:3
    x = xpos((group*3)-2:group*3) + 0.4;
    mu = mean_model_conf_err((group*3)-2:group*3);
    se = ci_model_conf_err((group*3)-2:group*3);
    for p = 1:length(x)-1
        pX = [x(p) x(p+1) x(p+1) x(p)];
        pY = [mu(p)-se(p) mu(p+1)-se(p+1) mu(p+1)+se(p+1) mu(p)+se(p)];
        h = fill(pX,pY,[1 0.5 0.5]);
        set(h,'edgecolor',[1 0.5 0.5]);
    end
    plot(x, mu, 'r-', 'LineWidth', 1.5)
end
boxplot(conf_err, 'positions', xpos+0.4, 'outliersize', 8, 'symbol', 'o', 'widths', 0.3, 'jitter', 0, 'colors', 'r')

set(gca, 'YLim', [0 1], 'XLim', [0.5 9.5], 'XTick', xpos+0.2, 'XTickLabel', {'PreL', 'PreM', 'PreH'}, 'FontSize', 14);
text(xpos(2), -0.1, 'PostL', 'FontSize', 14)
text(xpos(5), -0.1, 'PostM', 'FontSize', 14)
text(xpos(8), -0.1, 'PostH', 'FontSize', 14)
ylabel('Confidence', 'FontSize', 16);
box off

%% Densities
h3 = figure;
j=1;
bincenters = linspace(0.05,0.95,10);
for post = 1:3
    for pre = 1:3
        subplot(3,3,j);
        [n x] = hist(allConf_cor{j},bincenters);
        n = n./length(allConf_cor{j});  % normalise by number of corrects/incorrects
        bar(x+0.015, n, 0.35, 'g');
        hold on
        
        [n x] = hist(allConf_err{j},bincenters);
        n = n./length(allConf_err{j});
        bar(x-0.015, n, 0.35, 'r')
        set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 14);
        
        if pre == 1 & post == 2
            ylabel('Probability', 'FontSize', 20)
        end
        if pre == 2 & post == 3
            xlabel('Confidence', 'FontSize', 20)
        end
        
        j=j+1;
    end
end

h4 = figure;
j=1;
bincenters = linspace(0.05,0.95,10);
for post = 1:3
    for pre = 1:3
        subplot(3,3,j);
        [n x] = hist(allModelConf_cor{j},bincenters);
        n = n./length(allModelConf_cor{j});
        bar(x+0.015, n, 0.35, 'g');
        hold on
        [n x] = hist(allModelConf_err{j},bincenters);
        n = n./length(allModelConf_err{j});
        bar(x-0.015, n, 0.35, 'r')
        set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 14);
        
        if pre == 1 & post == 2
            ylabel('Probability', 'FontSize', 20)
        end
        if pre == 2 & post == 3
            xlabel('Confidence', 'FontSize', 20)
        end
        j=j+1;
    end
end

