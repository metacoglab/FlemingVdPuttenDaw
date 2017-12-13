% Implements multi-level mediation models on ROI data
% Reproduces statistics used to create Figure 5 in Fleming, van der Putten & Daw
% Stored in stats_forward{r}, stats_reverse{r}, where r indexes ROI
%
% Requires Tor Wager's Mediation Toolbox and functions from CanLabCore available here:
% https://canlabweb.colorado.edu/wiki/doku.php/help/mediation/m3_mediation_fmri_toolbox
%
% Steve Fleming stephen.fleming@ucl.ac.uk 2017

clear all
close all

cwd = pwd;

addpath(genpath('~/Dropbox/Utils/fmri/MetaLabCore/'));
dir_roiData_upsampled = '~/Dropbox/Research/Metacognition/stateactionexpt/FlemingVdPuttenDaw/fmri/roidata_upsample';
dir_datafiles = '~/Dropbox/Research/Metacognition/stateactionexpt/FlemingVdPuttenDaw/data';
dir_stats = '~/Dropbox/Research/Metacognition/stateactionexpt/FlemingVdPuttenDaw/stats';

name_subj = {'sub12','sub13','sub14','sub15','sub16','sub17','sub18','sub19','sub23',...
    'sub24','sub25','sub26','sub27','sub28','sub30','sub31','sub32','sub33','sub34','sub35','sub36','sub37'};
rois = {'pMFC','union_FPL','union_46','union_FPm'};

permutationtest = 0;    % If 0, load previously computed permutation test stats from .mat file; if 1, compute again
npermute = 1000;

xcenters = [1 2 3];
allBetas = [];
allSE = [];
for r = 1:length(rois)
    
    % Load up regression outputs from R hierarchical model of confidence/QSR, sems stored as
    % second part of vector
    cd(dir_stats)
    betaFile = ['regression_conf_split_' rois{r} '.csv'];
    betas = csvread(betaFile,1,1);
    segmented_ffx = betas([2 3])';
    segmented_se = betas([6 7])';
    
    % Load up regression outputs from R hierarchical segmented regression model, sems stored as
    % second part of vector
    h = figure;
    set(gcf, 'Position', [600 300 220 250])
    bar(1, segmented_ffx(1), 'FaceColor', [0.6 0.6 1])
    hold on
    bar(2, segmented_ffx(2), 'FaceColor', [0.8 0.8 0.8])
    errorbar(1, segmented_ffx(1), segmented_se(1), 'o ', 'Color', 'k', 'LineWidth', 2)
    errorbar(2, segmented_ffx(2), segmented_se(2), 'o ', 'Color', 'k', 'LineWidth', 2)
    set(gca, 'XTickLabel', {'',''}, 'YLim', [-0.7 0.5], 'FontSize', 16)
    ylabel('Regression coefficient')
    box off
    title(rois{r})
    
    % Create upsampled regressions
    for s = 1:length(name_subj)
        
        % Get behavioural data
        sub = name_subj{s}(4:5);
        datafile = ['fMRI_pilotData_sub_' num2str(sub) '_fMRI_2.mat'];        
        cd(dir_datafiles);
        load(datafile)
        cd(cwd);
        
        % extract locDATA variables
        rt =  locDATA.reaction_time_button;
        conf_rt =  locDATA.reaction_time_mouse;
        precoh = locDATA.dots_coherence;
        postcoh = locDATA.post_coherence;
        accuracy = locDATA.accuracy;
        accuracy(accuracy == 0) = -1;
        conf = locDATA.mouse_response;
        qsr = max(1-(1 - conf).^2, 1-(0 - conf).^2);
        err = isnan(rt) | isnan(accuracy) | isnan(conf);
        
        % Extract out ordinal variables for coherence
        index = unique(sort(locDATA.dots_coherence));
        for i = 1:3
            precoh(precoh == index(i)) = i;
            postcoh(postcoh == index(i)) = i;
        end
        
        % Get upsampled ROI data
        cd(dir_roiData_upsampled)
        datafile = ['ROI_' rois{r} '_' name_subj{s}];
        load(datafile)
        cd(cwd)
        
        subjectXtime(s,:) = ([1:length(ROI.time_series(1,:))].*ROI.spec.samp_reso) - ROI.spec.samp_reso;  % same for all subjects, subtract one bin to start at 0

        % Regression by timepoint on confidence / conf + qsr
        for t = 1:length(ROI.time_series(1,:))
            
            % Extract out signal over trials for this timepoint
            bold_t = ROI.time_series(:,t);
            
            % include qsr, store for permutation test
            regress(t).X{s} = [conf' qsr' log(rt)'];
            regress(t).X{s} = scale(regress(t).X{s});
            regress(t).Y{s} = bold_t;    % zscore per timepoint to get standardised coefficients
            b = glmfit(regress(t).X{s},regress(t).Y{s});
            confValBeta1(s,t) = b(2);
            confValBeta2(s,t) = b(3);
            
        end
        % End subject loop
    end
    
    %% Group-level analyses
    % Create null distribution of betas for group-level test by permuting
    % Y trial ordering within-subject, don't break neighbourhood of Y over
    % time
    if permutationtest
        for n = 1:npermute
            permute_confValBeta1 = [];
            permute_confValBeta2 = [];
            for s = 1:length(name_subj)
                permuteY = randperm(length(ROI.time_series(:,1)));
                for t = 1:length(ROI.time_series(1,:))
                    X = regress(t).X{s};
                    Y = regress(t).Y{s};
                    bold_permute = Y(permuteY);
                    b = glmfit(X,bold_permute);
                    permute_confValBeta1(s,t) = b(2);
                    permute_confValBeta2(s,t) = b(3);
                end
            end
            % Store t-values against zero (at group level) for this iteration
            [h p ci stats] = ttest(permute_confValBeta1);
            permute_tBeta1(n,:) = stats.tstat;
            [h p ci stats] = ttest(permute_confValBeta2);
            permute_tBeta2(n,:) = stats.tstat;
        end
        % compute significance bounds
        upperCI_beta1 = quantile(permute_tBeta1, 0.975);
        lowerCI_beta1 = quantile(permute_tBeta1, 0.025);
        upperCI_beta2 = quantile(permute_tBeta2, 0.975);
        lowerCI_beta2 = quantile(permute_tBeta2, 0.025);
        
    else
        cd(dir_stats)
        load([rois{r} '_confval_permutationTest.mat']);
        cd(cwd)
    end
    [h p ci stats] = ttest(confValBeta1);
    actual_tBeta1 = stats.tstat;
    [h p ci stats] = ttest(confValBeta2);
    actual_tBeta2 = stats.tstat;
    sigBeta1 = (actual_tBeta1 > upperCI_beta1) | (actual_tBeta1 < lowerCI_beta1);
    sigBeta2 = (actual_tBeta2 > upperCI_beta2) | (actual_tBeta2 < lowerCI_beta2);
    
    %% PLOT FIGURES
    xMean = mean(subjectXtime);
    
    %% Confidence / value beta
    h = figure;
    set(gcf, 'Position', [600 300 300 250])
    betaMean = mean(confValBeta2);
    betaSE = std(confValBeta2)./sqrt(length(name_subj));
    errorarea(xMean, betaMean, betaMean-betaSE, betaMean+betaSE, 'y', [0.8 0.8 0.4]);
    minmaxBeta1 = [min(betaMean) max(betaMean)];
    hold on
    betaMean = mean(confValBeta1);  % plot confidence second on top
    betaSE = (std(confValBeta1)./sqrt(length(name_subj)));
    errorarea(xMean, betaMean, betaMean-betaSE, betaMean+betaSE, [0.91 0.30 0.06], [0.91 0.44 0.20]);
    minmaxBeta2 = [min(betaMean) max(betaMean)];
    minmaxAll = [min([minmaxBeta1(1) minmaxBeta2(1)]) max([minmaxBeta1(2) minmaxBeta2(2)])];
    ymax = minmaxAll(2) + 0.15;
    ymin = minmaxAll(1) - 0.15;
    
    set(gca, 'FontSize', 16, 'YLim', [ymin ymax], 'XLim', [0 max(xMean)])
    line([0 max(xMean)], [0 0], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
    line([1.9 1.9], [ymin ymin+0.05], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1.5) % add line at onset of PDE
    
    % Add significance below plot
    hold on
    minBeta = min(betaMean);
    totSig1 = sum(sigBeta1);
    plot(xMean(sigBeta1), ones(1,totSig1).*(ymin+0.03), 'o ', 'Color', [0.91 0.44 0.20])
    totSig2 = sum(sigBeta2);
    plot(xMean(sigBeta2), ones(1,totSig2).*(ymin+0.015), 'yo ')
    ylabel('Regression coefficient')
    xlabel('Time (s)')
    box off
    title(rois{r})
end
cd(cwd)