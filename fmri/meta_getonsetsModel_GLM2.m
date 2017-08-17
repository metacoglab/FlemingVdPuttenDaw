% Script to get design matrix onsets
%
% Loads into spm 1st level design matrix and estimates job
%
% Steve Fleming 2008-2017 stephen.fleming@ucl.ac.uk

cwd = pwd;

% Options for processing
compute_conditions = 1;
construct_design = 1;
estimate = 1;

spmdir =  '/Users/sfleming/Dropbox/Utils/spm12';
addpath(spmdir);

%% Change according to your directory strucutre and scan parameters
dir_behav = '~/Dropbox/Research/Metacognition/stateactionexpt/task/locFullData';
dir_base    = '/Users/sfleming/Documents/Data/Postdecision/Data';
dir_epi     = 'Functional';
dir_stats = 'StatsDir_GLM2';
sess_prfx   = 'sess';
name_subj = {'sub12','sub13','sub14','sub15','sub16','sub17','sub18','sub19','sub23',...
    'sub24','sub25','sub26','sub27','sub28','sub30','sub31','sub32','sub33','sub34','sub35','sub36','sub37'};

fs = filesep;
blockNo = 4;
nslices = 42;
TR = 2.3352;
hpcutoff = 128;

for n = 1:length(name_subj)
    
    
    % Go to subject directory
    cd([dir_base fs name_subj{n}]);
    % if the intended stats directory does not exits, make it and go into it
    if exist(dir_stats) ~= 7
        mkdir(dir_stats);
        cd(dir_stats);
    else
        % if it does exist, go into it and delete existing contrast files
        % as we are re-estimating
        cd(dir_stats);
        delete('SPM.mat','*.nii','*.hdr');
    end
    
    outputDir = [dir_base fs name_subj{n} fs dir_stats];
    
    % Define behavioural data path
    datafile = ['fMRI_pilotData_sub_' name_subj{n}(4:5) '_fMRI_2.mat'];
    cd(dir_behav);
    load(datafile)
    cd(cwd);
    
    %=========================================================================
    %% Get onsets in scans from behavioural datafiles and define images
    %======================================================================
    % Reshape data into N blocks
    coherence = sort(unique(locDATA.dots_coherence));
    ntrials = length(locDATA.dots_direction);
    d = reshape(locDATA.dots_direction, ntrials/blockNo, blockNo)';
    rt =  reshape(locDATA.reaction_time_button, ntrials/blockNo, blockNo)';
    precoh = reshape(locDATA.dots_coherence, ntrials/blockNo, blockNo)';
    postcoh = reshape(locDATA.post_coherence, ntrials/blockNo, blockNo)';
    conf = reshape(locDATA.mouse_response, ntrials/blockNo, blockNo)';
    conf_rt =  reshape(locDATA.reaction_time_mouse, ntrials/blockNo, blockNo)';
    accuracy = reshape(locDATA.accuracy, ntrials/blockNo, blockNo)';
    err = isnan(rt) | isnan(accuracy) | isnan(conf);
    startTime = reshape(locDATA.timing.blockStart, ntrials/blockNo, blockNo)';
    
    for k = 1:blockNo;
        
        disp(['Computing event onsets for subject ' name_subj{n} ', session ' num2str(k)]);
        
        qsr = max(1-(1 - conf(k,:)).^2, 1-(0 - conf(k,:)).^2);
        
        scanStart = startTime(k,1);
        preOnsets = startTime(k,:) - scanStart;
        postOnsets = startTime(k,:) - scanStart + 1.9;
        confOnsets = startTime(k,:) - scanStart + 2.4;
        
        names = [];
        onsets = [];
        durations = [];
        pmod = [];
        orth = [];
        
        % Change names to allow marsbar analysis
        names{1} = 'Pre';
        onsets{1} = preOnsets(~err(k,:));
        durations{1} = 0.3;
        pmod(1).name{1} = 'RT';
        pmod(1).param{1} = zscore(log(rt(k, ~err(k,:))));
        pmod(1).poly{1} = 1;
        
        names{2} = 'Post';
        onsets{2} = postOnsets(~err(k,:));
        durations{2} = 0.3;
        
        names{3} = 'Conf';
        onsets{3} = confOnsets(~err(k,:));
        durations{3} = 0;
        pmod(3).name{1} = 'Confidence';
        pmod(3).param{1} = zscore(conf(k, ~err(k,:)));
        pmod(3).poly{1} = 1;
        
        cd(outputDir);
        conditionFile = sprintf('conditions%d.mat',k);
        
        save(conditionFile, 'onsets', 'names', 'durations', 'pmod', 'orth');
        
        %==========================================================================
        %% Construct design matrix
        %==========================================================================
        
        % Load files we have just created
        epiDir = [dir_base fs name_subj{n} fs dir_epi fs sess_prfx num2str(k)];
        conditionFile = sprintf('conditions%d.mat',k);
        
        % Get text file with movement regressors and concatenate with
        % first derivative
        mCorrFile = spm_select('List',epiDir,'^rp_af.*\.txt$');
        M = textread([epiDir fs mCorrFile]);
        R = [M [zeros(1,6); diff(M)]];
        
        cd(outputDir);
        multiFile = sprintf('multireg%d.mat',k);
        save(multiFile, 'R');
        
        % Assign .mat file with onsets/names/pmods in to path
        conditionPath = [outputDir fs conditionFile];
        multiregPath = [outputDir fs multiFile];
        
        % get epi files for this session
        epiDir = [dir_base fs name_subj{n} fs dir_epi fs sess_prfx num2str(k)];
        % select scans and concatenate
        f   = spm_select('List', epiDir, '^swuaf.*\.nii$');     % Select smoothed normalised images
        files  = cellstr([repmat([epiDir fs],size(f,1),1) f]);
        % clear temporary variables for next run
        jobs{1}.stats{1}.fmri_spec.sess(k).scans = files;
        f = []; files = [];
        
        jobs{1}.stats{1}.fmri_spec.sess(k).multi = {conditionPath};
        jobs{1}.stats{1}.fmri_spec.sess(k).multi_reg = {multiregPath};
        % high pass filter
        jobs{1}.stats{1}.fmri_spec.sess(k).hpf = hpcutoff;
        jobs{1}.stats{1}.fmri_spec.sess(k).regress = struct([]);
    end
    
    %==========================================================================
    %======================================================================
    jobs{1}.stats{1}.fmri_spec.dir = {outputDir};
    % timing variables
    jobs{1}.stats{1}.fmri_spec.timing.units     = 'secs';
    jobs{1}.stats{1}.fmri_spec.timing.RT        = TR;
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t    = nslices;
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t0   = nslices/2;
    
    % basis functions
    jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs         = [0 0];
    % model interactions (Volterra) OPTIONS: 1|2 = order of convolution
    jobs{1}.stats{1}.fmri_spec.volt                     = 1;
    % global normalisation
    jobs{1}.stats{1}.fmri_spec.global                   = 'None';
    % explicit masking
    jobs{1}.stats{1}.fmri_spec.mask                     = {[spmdir fs 'tpm/mask_ICV.nii']};
    % serial correlations
    jobs{1}.stats{1}.fmri_spec.cvi                      = 'AR(1)';
    % no factorial design
    jobs{1}.stats{1}.fmri_spec.fact = struct('name', {}, 'levels', {});
    
    
    %==========================================================================
    %% run model specification
    %==========================================================================
    if construct_design
        cd(outputDir);
        % save and run job
        save specify.mat jobs
        disp(['RUNNING model specification for subject ' name_subj{n}]);
        spm_jobman('run','specify.mat');
        clear jobs
    end
    
    % Ensure orthogonalisation and implicit masking is switched off
    load SPM
    SPM.xM.TH = repmat(-Inf, length(SPM.xM.TH), 1);
    SPM.xM.I = 0;
    save SPM SPM
    
    %% Estimate
    % setup job structure for model estimation and estimate
    % ---------------------------------------------------------------------
    if estimate
        jobs{1}.stats{1}.fmri_est.spmmat = {[outputDir fs 'SPM.mat']};
        save estimate.mat jobs
        disp(['RUNNING model estimation for subject ' name_subj{n}])
        spm_jobman('run','estimate.mat');
        clear jobs
    end
end
cd(cwd);
