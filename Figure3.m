% Plots data from pMFC ROI
% Reproduces panels from Figure 3 in Fleming, van der Putten & Daw
%
% Steve Fleming stephen.fleming@ucl.ac.uk 2017

clear all
close all

cwd = pwd;

fs = filesep;
addpath(genpath('~/Dropbox/Utils/fmri/MetaLabCore/'));
dir_roiData = '~/Dropbox/Published/FlemingNN2018/Code&analysis/FlemingVdPuttenDaw/fmri/roidata_post';
dir_roiData_upsampled = '~/Dropbox/Published/FlemingNN2018/Code&analysis/FlemingVdPuttenDaw/fmri/roidata_upsample';
dir_datafiles = '~/Dropbox/Published/FlemingNN2018/Code&analysis/FlemingVdPuttenDaw/data';
dir_modelfiles = '~/Dropbox/Published/FlemingNN2018/Code&analysis/FlemingVdPuttenDaw/regressors';

model = 'ideal_dt';
name_subj = {'sub12','sub13','sub14','sub15','sub16','sub17','sub18','sub19','sub23',...
    'sub24','sub25','sub26','sub27','sub28','sub30','sub31','sub32','sub33','sub34','sub35','sub36','sub37'};
rois = {'pMFC'};

c.cor =  [0 0.45 0.74];
c.err =  [0.85 0.33 0.10];
xcenters = [1 2 3];
allBetas = [];
allSE = [];

% get exclusion vectors for all subjects
load singleTrial_exclusions_post.mat

for r = 1:length(rois)
    
    bigData = [];
    
    for s = 1:length(name_subj)
        
        cd(dir_roiData)
        datafile = ['ROI_' rois{r} '_' name_subj{s}];
        load(datafile)
        
        %% Get behavioural data for X, Y and covariates
        sub = name_subj{s}(4:5);
        % Define behavioural data path
        datafile = ['fMRI_pilotData_sub_' num2str(sub) '_fMRI_2.mat'];
        cd(dir_datafiles);
        load(datafile)
        cd(cwd);
        
        % reshape some of the locDATA variables
        rt =  locDATA.reaction_time_button;
        conf_rt =  locDATA.reaction_time_mouse;
        precoh = locDATA.dots_coherence;
        postcoh = locDATA.post_coherence;
        accuracy = locDATA.accuracy;
        accuracy(accuracy == 0) = -1;
        conf = locDATA.mouse_response;
        err = isnan(rt) | isnan(accuracy) | isnan(conf);
        
        % Extract out ordinal variables for coherence
        index = unique(sort(locDATA.dots_coherence));
        for i = 1:3
            precoh(precoh == index(i)) = i;
            postcoh(postcoh == index(i)) = i;
        end
        
        % Note that trial exclusions have already been completed by
        % ROI_analysis script, all behavioural variables are indexed relative to this
        % structure within ROI
        Ynew = zscore(ROI.Y);
        qsr = max(1-(1 - ROI.conf).^2, 1-(0 - ROI.conf).^2);
        
        %% Get model log-odds update for this subject
        modelfile = [model num2str(sub) '_loglikPost.mat'];
        cd(dir_modelfiles);
        load(modelfile)
        cd(cwd);
        
        % Assign model log odds to trials
        logOddsPre = [];
        logOddsPost = [];
        for i = 1:length(precoh)
            if accuracy(i) == 1
                logOddsPre(i) = model_loglikPre_cor(precoh(i), postcoh(i));
                logOddsPost(i) = model_loglikPost_cor(precoh(i), postcoh(i));
            else
                logOddsPre(i) = model_loglikPre_err(precoh(i), postcoh(i));
                logOddsPost(i) = model_loglikPost_err(precoh(i), postcoh(i));
            end
        end
        
        %% Compute relevant means as function of post-decision coherence/accuracy
        for i = 1:3
            bold_correct(s,i) = mean(Ynew(ROI.acc == 1 & ROI.postcoh == i));
            bold_error(s,i) = mean(Ynew(ROI.acc ~= 1 & ROI.postcoh == i));
        end
        
        % Get linear fit
        p_correct = polyfit(xcenters, bold_correct(s,:), 1);
        p_error = polyfit(xcenters, bold_error(s,:), 1);
        fit_correct(s,:) = polyval(p_correct, xcenters);
        fit_error(s,:) = polyval(p_error, xcenters);
        
        %% Compute relevant means per bin of post-decision evidence
        subjectIndex = find(strcmp(name_subj{s}, output.subj));  % exclude BOLD outliers from model variables
        postDE = logOddsPost(~err & ~output.exclusion{subjectIndex});
        postDEBins = linspace(min(postDE), max(postDE), 5);
        for b = 1:length(postDEBins)-1
            bold_postDE(s,b) = mean(Ynew(postDE > postDEBins(b) & postDE <= postDEBins(b+1)));
            postDEcenter(s,b) = (postDEBins(b+1) + postDEBins(b))/2;
        end
        
        preDE = logOddsPre(~err & ~output.exclusion{subjectIndex});
        preDEBins = linspace(min(preDE), max(preDE), 5);
        for b = 1:length(preDEBins)-1
            bold_preDE(s,b) = mean(Ynew(preDE > preDEBins(b) & preDE <= preDEBins(b+1)));
            preDEcenter(s,b) = (preDEBins(b+1) + preDEBins(b))/2;
        end
             
        % Get linear fits
        missingpoints = isnan(bold_postDE(s,:));
        p_postDE = polyfit(postDEcenter(s,~missingpoints), bold_postDE(s,~missingpoints), 1);
        fit_postDE(s,:) = polyval(p_postDE, postDEcenter(s,:));
        
        %% Upsampled data
        % Get upsampled ROI data
        cd(dir_roiData_upsampled)
        datafile = ['ROI_' rois{r} '_' name_subj{s}];
        load(datafile)
        cd(cwd)
        
        subjectMean(s,:) = nanmean(ROI.time_series);
        subjectXtime(s,:) = ([1:length(ROI.time_series(1,:))].*ROI.spec.samp_reso) - ROI.spec.samp_reso;  % this is the same for all subjects, subtract one bin to start at 0
        
        % Extract out condition means as function of accuracy and
        % post-decision motion coherence
        for i = 1:3
            bold_correct_ups{i}(s,:) = nanmean(ROI.time_series(~err & accuracy == 1 & postcoh == i,:));
            bold_error_ups{i}(s,:) = nanmean(ROI.time_series(~err & accuracy ~= 1 & postcoh == i,:));
        end
        
        
    end
    
   % Post-decision motion strength figure
    h = figure;
    set(gcf, 'Position', [300 300 250 275])
    corcolor = [0.7 1 0.7; 0.35 1 0.35; 0 1 0];
    errcolor = [1 0.7 0.7; 1 0.35 0.35; 1 0 0];
    mean_bold_correct = mean(bold_correct);
    sem_bold_correct = std(bold_correct)./sqrt(length(name_subj));
    mean_bold_error = mean(bold_error);
    sem_bold_error = std(bold_error)./sqrt(length(name_subj));
    hold on
    plot(xcenters-0.2, mean(fit_correct), 'Color', c.cor, 'LineWidth', 2)
    plot(xcenters+0.2, mean(fit_error), 'Color', c.err, 'LineWidth', 2)
    hb = boxplot(bold_correct, 'positions', xcenters-0.2, 'outliersize', 8, 'symbol', 'o', 'widths', 0.3, 'jitter', 0, 'colors', c.cor);
    set(hb,'LineWidth', 2)
    hb = boxplot(bold_error, 'positions', xcenters+0.2, 'outliersize', 8, 'symbol', 'o', 'widths', 0.3, 'jitter', 0, 'colors', c.err);
    set(hb,'LineWidth', 2)
    set(gca, 'XLim', [0.2 3.8], 'XTick', xcenters, 'XTickLabel', {'Low', 'Med', 'High'}, 'FontSize', 16);
    ylabel('BOLD a.u.')
    box off
    
    % Post-decision evidence figure
    h = figure;
    set(gcf, 'Position', [300 300 400 275])
    mean_bold_postDE = nanmean(bold_postDE);
    sem_bold_postDE = nanstd(bold_postDE)./sqrt(length(name_subj));
    hold on
    hb = boxplot(bold_postDE, 'positions', 1:length(mean_bold_postDE), 'outliersize', 8, 'symbol', 'o', 'widths', 0.3, 'jitter', 0, 'colors', 'k');
    set(hb, 'LineWidth', 2)
    plot(1:length(mean_bold_postDE), nanmean(fit_postDE), 'k', 'LineWidth', 2)
    set(gca, 'FontSize', 16, 'XTick', 1:length(mean_bold_postDE));
    ylabel('BOLD a.u.')
    xlabel('PDE bins (\Delta LO_{correct})')
    box off
    
    % Upsampled post-decision motion figure
    h = figure;
    set(gcf, 'Position', [600 300 400 275])
    corcolor = [0.52 0.65 0.74; 0.26 0.55 0.74; 0 0.45 0.74];
    errcolor = [0.85 0.69 0.63; 0.85 0.51 0.36; 0.85 0.33 0.1];
    
    xMean = mean(subjectXtime);
    for i = 1:3
        boldMean = mean(bold_correct_ups{i});
        boldSE = std(bold_correct_ups{i})./sqrt(length(name_subj));
        plot(xMean, boldMean, 'Color', corcolor(i,:), 'LineWidth', i+0.5);
        hold on
        boldMean = mean(bold_error_ups{i});
        boldSE = std(bold_error_ups{i})./sqrt(length(name_subj));
        plot(xMean, boldMean, 'Color', errcolor(i,:), 'LineWidth', i+0.5)
        hold on
    end
    ylabel('BOLD a.u.')
    xlabel('Time (s)')
    legend({'Cor low', 'Err low', 'Cor med', 'Err med', 'Cor high', 'Err high'}, 'Location', 'NorthEast')
    legend boxoff
    set(gca, 'FontSize', 16, 'XLim', [0 max(xMean)])
    box off
    
end

cd(cwd)
