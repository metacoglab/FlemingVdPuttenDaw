% Script to set up conrasts for 1st level analysis
%
% Requires that model estimation for all subjects being analysed be
% completed first by getonsetsModel.m
%
% Steve Fleming 2008-2017 stephen.fleming@ucl.ac.uk

% =========================== Set up contrasts ============================
% =========================================================================

cwd = pwd;

%% Change according to your scan parameters and directory structure
dir_base    = '/Users/sfleming/Documents/Data/Postdecision/Data';
dir_epi     = 'Functional';
dir_stats = 'StatsDir_GLM1';
sess_prfx   = 'sess';
name_subj = {'sub12','sub13','sub14','sub15','sub16','sub17','sub18','sub19','sub23',...
    'sub24','sub25','sub26','sub27','sub28','sub30','sub31','sub32','sub33','sub34','sub35','sub36','sub37'};
blockNo = 4;

cwd = pwd;

fs = filesep;
cd(cwd);
k = 12; % number of motion regs

for s0 = 1 : length(name_subj)
    % Contrast names
    clear T
    T.contrasts = {'Post-acc int'};
    
    % You need as many rows here as entries in the cell array above
    T.contrastVectors(1,:) =    [repmat([0 0 0 0 1 0 -1 0 zeros(1,k)],1,blockNo) zeros(1,blockNo)];
    j = 1;
    jobs{1}.stats{1}.con.spmmat    = {[dir_base fs name_subj{s0} fs dir_stats fs 'SPM.mat']};
    
    %     Specify tcontrasts
    for cont_nr = 1:length(T.contrasts)
        contrasttype = T.contrasts{cont_nr};
        contr_input = T.contrastVectors(cont_nr,:);
        
        % setup job structure for contrasts
        jobs{1}.stats{1}.con.consess{j}.tcon.name           = contrasttype;
        jobs{1}.stats{1}.con.consess{j}.tcon.convec         = contr_input;
        jobs{1}.stats{1}.con.consess{j}.tcon.sessrep        = 'none';
        j=j+1;
    end
    
    % if 1 then all existing contrasts are deleted
    jobs{1}.stats{1}.con.delete                         = 1;
    
    outputDir = [dir_base fs name_subj{s0} fs dir_stats];
    cd (outputDir);
    % save and run job
    save contrasts.mat jobs
    disp(['RUNNING contrast specification for subject number  ' name_subj{s0}]);
    spm_jobman('run','contrasts.mat');
    disp(['Contrasts created for ' num2str(s0) ' subjects']);
    clear jobs
    
end   % end of subject loop

cd(cwd);
