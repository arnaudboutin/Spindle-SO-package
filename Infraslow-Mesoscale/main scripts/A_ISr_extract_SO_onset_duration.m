%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adrien Conessa (CIAMS, Université Paris-Saclay)
% Arnaud Boutin (CIAMS, Université Paris-Saclay)


% Information extracted for each SO:
%   onset
%   duration

clearvars 
clc
subject2process=2:31;

channels={'Fp1','Fp2',...
    'F7','F3','Fz','F4','F8',...
    'FT9','FC5','FC1','FC2','FC6','FT10',...
    'T7','C3','Cz','C4','T8',...
    'CP5','CP1','CP2','CP6',...
    'P7','P3','Pz','P4','P8',...
    'O1','Oz','O2'};


MainDir         = 'your dir\';   % full path to the main directory
SODataDir      = '\IS_output_SO\';                              % directory with SO data; should be placed under <main_dir_path>
SOOutputDir    = '\onset_duration_extracted\';                 % output directory to save the exctracted SOs ; will be placed under <main_dir_path>

baseline_or_conso={'baseline','conso'};
method='Sta'; % Alb
%% some vars

sleepStage          = 'NREM23'; %'NREM2' - for NREM2/ 'NREM3' - for NREM3/ 'NREM23' - for NREM2 & NREM3 together
% Sleep Stage 
disp(['Sleep stage: ' sleepStage]);
if strcmp(sleepStage,'NREM2')
    indSleepStage = 2;
elseif strcmp(sleepStage,'NREM3')
    indSleepStage = 3;
elseif strcmp(sleepStage,'NREM23')
    indSleepStage = [2 3];
end


SODirPath = fullfile(MainDir, SODataDir);   % full path to the directory that contains SO data
if ~exist(SODirPath, 'dir') 
    warning(['Looking for the SO data directory: ' SODirPath ' ...']);
    error('The direcory with SO data does not exist. CHECK!!!');
end

outputDirPath = fullfile(MainDir, SOOutputDir);  % full path to the output directory
% Create output folder if it do not exist before running the analysis 
if ~exist(outputDirPath, 'dir')
    mkdir(outputDirPath);
end

% get all SO files for all subjects
allFiles = dir(fullfile(SODirPath, '*.mat'));

% store recap in txt file
SO_recap=fopen( [MainDir '\SO_recap.txt'], 'a' );
fprintf(SO_recap,['\n' datestr(datetime) ' - '  method '\n']);
fclose(SO_recap); 
   
for num_subject=subject2process % loop over each subject

    subject=num2str(num_subject);

    if num_subject<10
        subject_str=['00' subject];
    else
        subject_str=['0' subject];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for inight=1:2
        night_type=baseline_or_conso{inight};
        filename=['somhtd_' subject_str '_' night_type '_SO_raw_' method '.mat'];

        if ~any(strcmp({allFiles.name},filename))
            fprintf(['\nno file for ' night_type ' night\n'])
            continue
        end

        output_subj_dir_path    = fullfile(outputDirPath, ['somhtd_' subject_str]);

        if ~exist(output_subj_dir_path, 'dir')
            mkdir(output_subj_dir_path);
        end

        % Display file name
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp(['Filename: ', filename]);

         %load SO Core
        % each line represents one detected SO over all recorded electrodes
        load(fullfile(SODirPath, filename), 'SO_core');
        sampling_rate=SO_core.sampling_rate;
        
        for chan=1:numel(channels) % loop over each channel
            source_EEG=channels{chan};

            SO=SO_core.(source_EEG);
            
            if isempty(SO)
                SO_extracted=[];


                % store in file text
                fprintf(['found 0 SO on electrode ' source_EEG ' during ' sleepStage '\n'])
                SO_recap=fopen([MainDir '\SO_recap.txt'], 'a' );
                fprintf(SO_recap,['somhtd_' subject_str '_' night_type ' : 0 SO on ' source_EEG ' for ' sleepStage '\n']);
                fclose(SO_recap); 

                outputFilePath = fullfile(output_subj_dir_path, ['somhtd_' subject_str '_' night_type '_SO_' source_EEG '_' sleepStage '_' method '.mat']);

                save(outputFilePath,"SO_extracted")
                continue
            end

            % find SOs during the stage of interest
            SO_to_extract=ismember(SO.stage,indSleepStage);
            nSO=sum(SO_to_extract);
            fprintf(['found ' num2str(nSO) ' SO on electrode ' source_EEG ' during ' sleepStage '\n'])

            onset=SO.onset(SO_to_extract)./sampling_rate; % get the onset in secs
            duration=SO.duration(SO_to_extract); % duration in sec
            amplitude=SO.amplitude(SO_to_extract); % amplitude
            stage=SO.stage(SO_to_extract); % sleep stage

            SO_extracted=[onset' duration' amplitude' stage']; % put it in array

            % save extracted data
            outputFilePath = fullfile(output_subj_dir_path, ['somhtd_' subject_str '_' night_type '_SO_' source_EEG '_' sleepStage '_' method '.mat']);
            save(outputFilePath,"SO_extracted")
            
            % Overwrite text file with these infos
            SO_recap=fopen([MainDir '\SO_recap.txt'], 'a' );
            fprintf(SO_recap,['somhtd_' subject_str '_' night_type ' : ' num2str(nSO) ' SO on ' source_EEG ' for ' sleepStage '\n']);
            fclose(SO_recap); 


            clear SO_extracted
        end

        clear SO_core
    end

    
end

