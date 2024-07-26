%% fft power sigma
% Adrien Conessa (CIAMS, Université Paris-Saclay)
% Arnaud Boutin (CIAMS, Université Paris-Saclay)

% clc
clearvars
close all

subject_to_process=2:31;
%% Directories
MainDir = '\'; % your main directory

% eeglab path
addpath('C:\toolbox\eeglab2021.1\');

EEGData =[MainDir 'EEG\']; % eeg data dir
ScoringDir=[MainDir 'Scoring\']; % scoring dir
OutpuDir=[MainDir 'data\']; % output dir

%% initiate eeg lab
eeglab
close(gcf)
clc

%% variables 
freq_ms=0.01:0.01:0.6; % mesoscale frequency band

load([OutpuDir 'ISr_channel_infos.mat']);

all_files = dir([EEGData '*.eeg']);
baseline_or_conso={'baseline','conso'};
night_to_process='baseline';
stage='N2'; % sleep stage to process
method='morlet';
time_specific=0; % if you want to add time constraint (60 to compute on the first 60 minutes of sleep)

all_eeg_channels={Channel.Name};
all_eeg_channels(contains(all_eeg_channels,{'Chin','EoG'}))=[];

chan2process={'C4'}; % channel to process
pos_chan2pro=find(ismember(all_eeg_channels,chan2process));% find position in the channel file


idx_subject=1;

if time_specific==0
    idx_time='';
else
    idx_time=['_time_' num2str(time_specific)];
end


%% load data

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

    vhdr_name = strrep(file_name,'eeg','vhdr'); 
    out_name = strrep(file_name,'.eeg','_raw_SO.mat');


    for i_night=1:numel(file_name)  % loop over each night
        
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
        EEG = pop_loadbv(EEGData, vhdr_name{i_night},[],pos_chan2pro); % specify channels to load 
         

        for ichan=1:numel(chan2process)

            current_chan=EEG.chanlocs(ichan).labels;

            Data=EEG.data(ichan,:);
            % process
            fprintf(['processing channel ' current_chan '\n'])
            [spectrum_mesoscale.(current_chan)(idx_subject,:),spectrum_mesoscale_weighted.(current_chan)(idx_subject,:)]=ISr_compute_fft_sigma_mesoscale(Data,hypno_night,stage,method,freq_ms,time_specific);

        end





        clear EEG Data
   
    end
    toc

    idx_subject=idx_subject+1;


end

% save data
save([OutpuDir 'spectrum_mesoscale_weighted_' night_type '_' stage '_' method idx_time '.mat'],'spectrum_mesoscale_weighted','freq_ms')
save([OutpuDir 'spectrum_mesoscale_' night_type '_' stage '_' method idx_time '.mat'],'spectrum_mesoscale','freq_ms')
