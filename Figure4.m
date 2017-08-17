% Implements multi-level mediation models on ROI data
% Reproduces statistics used to create Figure 4 in Fleming, van der Putten & Daw
% Stored in stats_forward{r}, stats_reverse{r}, where r indexes ROI
%
% Requires Tor Wager's Mediation Toolbox and functions from CanLabCore available here:
% https://canlabweb.colorado.edu/wiki/doku.php/help/mediation/m3_mediation_fmri_toolbox
%
% Steve Fleming stephen.fleming@ucl.ac.uk 2017

addpath(genpath('~/Dropbox/Utils/fmri/MediationToolbox/'));
addpath(genpath('~/Dropbox/Utils/fmri/CanlabCore/'));

dir_datafiles = '~/Dropbox/Research/Metacognition/stateactionexpt/github/data';
dir_modelfiles = '~/Dropbox/Research/Metacognition/stateactionexpt/github/regressors';
dir_roiData = '~/Dropbox/Research/Metacognition/stateactionexpt/github/fmri/roidata_conf';

allRoiName = {'union_46', 'union_FPL', 'union_FPm'};    % which ROIs are Y variables
model = 'ideal'; % which model to extract out PDE predictions from

name_subj = {'sub12','sub13','sub14','sub15','sub16','sub17','sub18','sub19','sub23',...
    'sub24','sub25','sub26','sub27','sub28','sub30','sub31','sub32','sub33','sub34','sub35','sub36','sub37'};
load singleTrial_exclusions_conf.mat

for r = 1:length(allRoiName)
    
    for s = 1:length(name_subj)
        
        % get behavioural data
        sub = name_subj{s}(4:5);
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
        
        % get timeseries in ROI
        
        cd(dir_roiData)
        datafile = ['ROI_' allRoiName{r} '_' name_subj{s}];
        load(datafile)
        ts = zscore(ROI.Y);
        
        % Get model log-odds update for this subject
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
        
        % make X/Y/M and covariate vectors and center/scale
        subjectIndex = find(strcmp(name_subj{s}, output.subj));
        X{s} = zscore(logOddsPost(~err & ~output.exclusion{subjectIndex}))';
        Y{s} = zscore(conf(~err & ~output.exclusion{subjectIndex}))';
        C{s} = [zscore(log(rt(~err & ~output.exclusion{subjectIndex})))' zscore(logOddsPre(~err & ~output.exclusion{subjectIndex}))'];
        M{s} = ts;
        
        disp(['Extracting timecourses for subject ' name_subj{s}])
    end
    
    %% Run multi-level mediation
    [paths_forward{r}, stats_forward{r}] = mediation(X, Y, M, 'covs', C, 'boot', 'plots', 'verbose', 'names', {'PDE', 'Confidence', allRoiName{r}}, 'bootsamples', 10000);
    close all
    %% Run multi-level mediation in reverse (swap X and Y)
    [paths_reverse{r}, stats_reverse{r}] = mediation(Y, X, M, 'covs', C, 'boot', 'plots', 'verbose', 'names', {'Confidence', 'PDE', allRoiName{r}}, 'bootsamples', 10000);
    close all
end