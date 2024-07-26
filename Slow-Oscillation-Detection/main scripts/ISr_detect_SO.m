%% Slow oscillation detection
% Arnaud Boutin (CIAMS, Université Paris-Saclay)
% Adrien Conessa (CIAMS, Université Paris-Saclay)

%%%%%%
% Create "raw" and "output" folders
% The raw folder includes eeg files (.dat or .eeg, .vmrk, .vhdr) and scoring files (modify the suffix of the files)
% EEG files: Binary format and multiplexed 
%%%%%%

clc
clearvars
close all

subject_to_process=2:31;
method='Sta';  % Alb  Sta
%% Directories

redo=1; % 1 if you want to ecrase previous sps files

MainDir = 'your dir'; % your main dir

% Load eeglab and the spindle-detection package with all required functions
addpath('C:\toolbox\eeglab2021.1\');
addpath('your dir\dectection functions\');

EEGData =[MainDir 'EEG_sps\'];%EEG data DIR
OutputDir = '\'; % output folder slow oscillations
ScoringDir='\'; % folder with scoring files

eeglab
close(gcf)
clc

load('ISr_channel_infos.mat');


chan2process={'C4'};
% chan2process={'Fp1','Fp2',...
%     'F7','F3','Fz','F4','F8',...
%     'FT9','FC5','FC1','FC2','FC6','FT10',...
%     'T7','C3','Cz','C4','T8',...
%     'CP5','CP1','CP2','CP6',...
%     'P7','P3','Pz','P4','P8',...
%     'O1','Oz','O2'};


pos_chan2pro=find(ismember({Channel.Name},chan2process)); % postition of channels to process

all_files = dir([EEGData '*.eeg']);
baseline_or_conso={'baseline','conso'};
night_to_process='conso';

%% Spindle detection 

for num_subject=subject_to_process % loop over each subject
    tic
    subject=num2str(num_subject);

    if num_subject<10
        subject_str=['00' subject];
    else
        subject_str=['0' subject];
    end


    if ~strcmp(night_to_process,'both')
        files_to_process=contains({all_files.name},[subject_str '_' night_to_process]);
    else
        files_to_process=contains({all_files.name},subject_str);
    end
    
    file_name={all_files(files_to_process).name};

    fprintf(['\nprocessing subject' subject_str '\nfiles : ' file_name{1} '\n'])
    % name of the eeg files
    vhdr_name = strrep(file_name,'eeg','vhdr'); 
    out_name = strrep(file_name,'.eeg','_raw_SO.mat');


    for i_night=1:numel(file_name) % loop over each night

        if exist([OutputDir out_name{i_night}],'file') && redo==0
            disp(['Already done : ' out_name{i_night}]);
        else
            
            disp(['Analysis : ' file_name{i_night}]);
   
            if strcmp(night_to_process,'both')
                night_type=baseline_or_conso{i_night};
            else
                night_type=night_to_process;
            end

            try % Read Scoring file
                i_scoringFile = [ScoringDir 'hypno_S' subject_str '_night_' night_type '.mat']; % scoring file suffix %  'score_test' 
                load (i_scoringFile)
            catch
                disp(['Error no scoring file found : hypno_S' subject_str '_night_' night_type '.mat'])
            end
      
            % load EEG data
            EEG = pop_loadbv(EEGData, vhdr_name{i_night}, [], pos_chan2pro); % specify channels to load 
             
            for nCh=1:length(EEG.chanlocs)
                disp(['Channel : ' EEG.chanlocs(nCh).labels]);
            end
    
            % Info initialization
            Info.Recording.dataDim = size(EEG.data);
            Info.Recording.sRate = EEG.srate;
            % get the default settings for SO detection
            try
                Info = ISr_SO_getInfoDefaults(Info, ['SO_' method]);
            catch
                disp(['Error ISr_SO_getInfoDefaults function: ' file_name{i_night}])
            end
     

            % Check here if trouble with multiple ch detection
            for nChan =length(EEG.chanlocs):-1:1  % number of channels
                Ind = find(strcmp({Channel.Name},EEG.chanlocs(nChan).labels));
                if isempty(Ind)
                    disp(['Error no ' EEG.chanlocs(nChan).labels ' ch found in channel infos'])
                else
                    currentChanInfos(1,nChan) = Channel(Ind);
                end
            end
    
    
            Info.Electrodes = currentChanInfos;
            if size(Info.Electrodes,2)>1
                [~,sortElec] = sort(convertStringsToChars({Info.Electrodes(:).Name}));
                Info.Electrodes = Info.Electrodes(sortElec);
            else
                sortElec=1;
            end
                
            Data.SSRef = EEG.data(sortElec,:); % channel's data
            
            if strcmp(method,'Sta') % detection
                SO_core = ISr_SO_FindSORef_Sta(Data, Info, hypno_night);
            end

            save([OutputDir 'somhtd_' subject_str '_' night_type '_SO_raw_' method '.mat'],'SO_core','Info')

            clear Data EEG SO_Core Info currentChanInfos hypno_night
        end
    
   
    end
    toc
end
