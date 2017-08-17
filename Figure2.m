% Plots performance and confidence data together with posterior predictive
% simulations from STAN model fits
% Reproduces Figure 2 in Fleming, van der Putten & Daw
%
% Steve Fleming stephen.fleming@ucl.ac.uk 2017

clear all
close all
fs = filesep;

saveRegressors = 0;
dataset = input('Which dataset? 2=behav, 3=fmri ');
model = 'ideal'; % which model to simulate predictions from

simTrials = 5000;    % how many trials to simulate per parameter draw per condition (9 conditions)
Ndraws = 100;    % how many draws from subject-level posteriors
conf_sigma = 0.025; % fixed parameter relating model to observed confidence

cwd = pwd;
baseDir = '~/Dropbox/Research/Metacognition/stateactionexpt/github/stan/modelfits'; %% path to stan model fits

%% Load fitted parameters and data
fitsDir = [baseDir fs model fs];

switch model
    case {'weighted', 'accweighted'}
        cd(fitsDir)
        w = csvread(['subParams_w_dataset' num2str(dataset) '.csv'],1,1);
        w_sd = csvread(['subParams_w_sd_dataset' num2str(dataset) '.csv'],1,1);
    case 'mapping'
        cd(fitsDir)
        gamma = csvread(['subParams_gamma_dataset' num2str(dataset) '.csv'],1,1);
        gamma_sd = csvread(['subParams_gamma_sd_dataset' num2str(dataset) '.csv'],1,1);
end
cd(fitsDir)
k1 = csvread(['subParams_k1_dataset' num2str(dataset) '.csv'],1,1);
k1_sd = csvread(['subParams_k1_sd_dataset' num2str(dataset) '.csv'],1,1);
m = csvread(['subParams_m_dataset' num2str(dataset) '.csv'],1,1);
m_sd = csvread(['subParams_m_sd_dataset' num2str(dataset) '.csv'],1,1);
cd(cwd);

if dataset == 2
    dirData = '~/Dropbox/Research/Metacognition/stateactionexpt/github/data/'; % fmri experiment data out of scanner
    filename = 'fMRI_pilotData_sub_';    % use for fMRI
    suffix = ''; % use for fMRI data out-of-scanner
    subjects = [12:28 30:37];
elseif dataset == 3
    dirData = '~/Dropbox/Research/Metacognition/stateactionexpt/github/data/'; % fmri experiment data
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
    
    datafile = [filename num2str(subjects(s)) suffix '_2.mat'];
    cd(dirData);
    load(datafile);
    cd(cwd);
    
    conf = [];
    acc = [];
    precoh_index = [];
    postcoh_index = [];
    
    precoh = locDATA.dots_coherence';
    postcoh = locDATA.post_coherence';
    dir = locDATA.dots_direction/360;
    dir(dir==0.5) = -1;
    action = locDATA.button_response - 1;
    conf = locDATA.mouse_response;
    rt = log(locDATA.reaction_time_button);
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
    
    %% confidence/performance 3x3 for each subject
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
    
    %% generate model data from subject-level parameter draws
    for draw = 1:Ndraws
        conf = [];
        acc = [];
        precoh_index = [];
        postcoh_index = [];
        j=1;
        % draw parameters for this replication
        k1_rep = normrnd(k1(s), k1_sd(s));
        m_rep = normrnd(m(s), m_sd(s));
        switch model
            case {'weighted', 'accweighted'}
                w_rep = normrnd(w(s), w_sd(s));
                if w_rep < 0
                    w_rep = 0;
                elseif w_rep > 1;
                    w_rep = 1;
                end
            case 'mapping'
                gamma_rep = normrnd(gamma(s), gamma_sd(s));
        end
        
        for post = 1:3
            for pre = 1:3
                
                switch model
                    case 'ideal'
                        simdata = metaModelFit_generate_data(simTrials, coherence, pre, post, conf_sigma, model, 'k1', k1(s), 'm', m(s));
                    case {'weighted', 'accweighted'}
                        simdata = metaModelFit_generate_data(simTrials, coherence, pre, post, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'w', w_rep);
                    case 'mapping'
                        simdata = metaModelFit_generate_data(simTrials, coherence, pre, post, conf_sigma, model, 'k1', k1_rep, 'm', m_rep, 'gamma', gamma_rep);
                end
                
                % get big vectors concatenated over conditions, conditional
                % on accuracy from the model's choices
                conf = [conf simdata.conf];
                acc = [acc simdata.acc];
                precoh_index = [precoh_index ones(1,length(simdata.conf)).*pre];
                postcoh_index = [postcoh_index ones(1,length(simdata.conf)).*post];
                allModelConf_cor{j} = [allModelConf_cor{j} simdata.conf(simdata.acc == 1)];
                allModelConf_err{j} = [allModelConf_err{j} simdata.conf(simdata.acc == 0)];
                
                % take averages per condition
                model_perf_bydraw(s,j,draw) = nanmean(simdata.acc)*100;
                model_loglikPost_cor_bydraw(pre,post,draw) = nanmean(simdata.loglik_post(simdata.acc == 1));
                model_loglikPost_err_bydraw(pre,post,draw) = nanmean(simdata.loglik_post(simdata.acc == 0));
                model_loglikPre_cor_bydraw(pre,post,draw) = nanmean(simdata.loglik_pre(simdata.acc == 1));
                model_loglikPre_err_bydraw(pre,post,draw) = nanmean(simdata.loglik_pre(simdata.acc == 0));
                model_conf_cor_bydraw(s,j,draw) = nanmean(simdata.conf(simdata.acc == 1));
                model_conf_err_bydraw(s,j,draw) = nanmean(simdata.conf(simdata.acc == 0));
                
                j=j+1;
            end
        end
        
        disp(['Parameter draw ' num2str(draw) ' out of ' num2str(Ndraws)])
    end
    
    % marginalise over parameter draws
    model_perf = nanmean(model_perf_bydraw, 3);
    model_loglikPre_cor = nanmean(model_loglikPre_cor_bydraw, 3);
    model_loglikPre_err = nanmean(model_loglikPre_err_bydraw, 3);
    model_loglikPost_cor = nanmean(model_loglikPost_cor_bydraw, 3);
    model_loglikPost_err = nanmean(model_loglikPost_err_bydraw, 3);
    model_conf_cor = nanmean(model_conf_cor_bydraw, 3);
    model_conf_err = nanmean(model_conf_err_bydraw, 3);
    
    if sum(any(isnan(model_loglikPost_err))) | sum(any(isnan(model_loglikPost_cor)))    % not enough trials to get estimate
        warning('Not enough draws to get a good estimate, increase simTrials')
    end
        
    if saveRegressors && dataset == 3
        % store regressors for use in mediation modeling of fMRI data
        cd regressors
        save([model num2str(subjects(s)) '_loglikPost.mat'], 'model_loglikPre_cor', 'model_loglikPre_err', 'model_loglikPost_cor', 'model_loglikPost_err');
        cd(cwd);
    end
    
    disp(['Finished simulation for subject ' num2str(s)])
    
    
end

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

mean_conf_cor = nanmean(conf_cor);
sem_conf_cor = nanstd(conf_cor)./sqrt(length(subjects));
ci_conf_cor = sem_conf_cor.*tinv(0.975,length(subjects)-1);
mean_model_conf_cor = nanmean(model_conf_cor);
sem_model_conf_cor = nanstd(model_conf_cor)./sqrt(length(subjects));
ci_model_conf_cor = sem_model_conf_cor.*tinv(0.975,length(subjects)-1);
mean_conf_err = nanmean(conf_err);
sem_conf_err = nanstd(conf_err)./sqrt(length(subjects));
ci_conf_err = sem_conf_err.*tinv(0.975,length(subjects)-1);
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

