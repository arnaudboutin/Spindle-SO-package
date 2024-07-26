%% compute SO-sps coupling
% Adrien Conessa (CIAMS, Université Paris-Saclay)
% Arnaud Boutin (CIAMS, Université Paris-Saclay)
% script to initialize the SO-sps coupling computation

clearvars
close all
% clc

maxNumCompThreads(6)

%% Initialize constant
subject_to_process=2:31;
plotresult=0; % if you want to plot the polarplot for each subject
sps_band=[11 16]; % sps frequency band
SO_band=[0.5 1.25]; % SO frequency band
SO_check=0; % 1 if you want to plot example of SO-sps coupling

method_SO='Sta';
night_to_process='conso'; % both, baseline or conso
stage='NREM23'; % stage to process 
type=''; % '', 'Grp', 'Iso'
time_specific=0; % put 20 if you want to compute the concomitant for the first 20 minutes of sleep for example (necessitate the corresponding conco SO sps)

chan2process={'Fp1','Fp2',...
    'F7','F3','Fz','F4','F8',...
    'FT9','FC5','FC1','FC2','FC6','FT10',...
    'T7','C3','Cz','C4','T8',...
    'CP5','CP1','CP2','CP6',...
    'P7','P3','Pz','P4','P8',...
    'O1','Oz','O2'};


sampling_rate=250;

time_window=[-1 2.5]; % extracted time window around SO onset

%% Dir

MainDir='\';
DataDir=[MainDir '\data\'];
ConcoDir=[DataDir 'conco_SO_SPS']; % folder with concomitant events stored
EEGDir='\EEG\'; % eeg dir

addpath('C:\toolbox\eeglab2021.1\'); % need eeglab to load the data
eeglab
close(gcf)

%% load EEG
if time_specific==0
    indexation='';
else
    indexation=['_time_' num2str(time_specific)];
end


load ([DataDir 'ISr_channel_infos_sps.mat']) % load channels infos
pos_chan2pro=find(ismember({Channel.Name},chan2process)); % find position of channels to process

all_files = dir([EEGDir '*.eeg']);

% load the output struct to append it, if it already exist
if exist([DataDir 'PAC_conco_SO_SPS' type '_' stage '_' method_SO indexation '.mat'],'file')
    load([DataDir 'PAC_conco_SO_SPS' type '_' stage '_' method_SO indexation '.mat'],'PAC')
end
if exist([DataDir 'all_phase_SO_SPS' type '_' stage '_' method_SO indexation '.mat'],'file')
    load([DataDir 'all_phase_SO_SPS' type '_' stage '_' method_SO indexation '.mat'],'all_phase_SO_SPS')
end


%% process
for num_subject=subject_to_process % loop over each subject

    subject=num2str(num_subject);

    if num_subject<10
        subject_str=['00' subject];
    else
        subject_str=['0' subject];
    end
    
    SubjectName=['somhtd_' subject_str];

    if ~strcmp(night_to_process,'both')
        files_to_process=contains({all_files.name},[subject_str '_' night_to_process]);
    else
        files_to_process=contains({all_files.name},subject_str);
    end
    
    
    file_name={all_files(files_to_process).name};

    fprintf(['\nprocessing subject' subject_str '\nfiles : ' file_name{1} '\n'])
    if numel(file_name)>1
        fprintf(['------- '  file_name{2:end} '\n'])
    end

    % load concomitant events
    if time_specific==0

        load([DataDir 'conco_SO_SPS\' SubjectName '_conco_SO_SPS_' stage '_' type '_' method_SO '.mat'])
    else
        load([DataDir 'conco_SO_SPS\' SubjectName '_conco_SO_SPS_' stage '_' type '_' method_SO indexation '.mat'])
    end

    for ifile=1:numel(file_name) % loop over each night
        night=split(file_name{ifile},'_');
        night=night{3};

        fprintf(['processing ' night ' night\n'])
        tic

        % load EEG data
        vhdr_name = strrep(file_name{ifile},'eeg','vhdr');   
        EEG = pop_loadbv(EEGDir, vhdr_name, [], pos_chan2pro); % specify channels to load 
        channels_order={EEG.chanlocs.labels};
        EEG_channel_f1=cell(1);
        EEG_channel_f2=cell(1);
        EEG_channel=cell(1);
        
        
        for ichan=1:numel(channels_order) % loop over each channel, bandpass according to SO and sps frequency band
            fprintf(['filtering channel ' channels_order{ichan} '\n'])

            if all(isnan(EEG.data(ichan,:)))
                continue
            end

            % filter data
            [EEG_channel_f1{ichan,1},d3]=bandpass(EEG.data(ichan,:),SO_band,sampling_rate,'ImpulseResponse','fir');
            EEG_channel_f1{ichan,2}=channels_order{ichan};

            [EEG_channel_f2{ichan,1},d4]=bandpass(EEG.data(ichan,:),sps_band,sampling_rate,'ImpulseResponse','fir');
            EEG_channel_f2{ichan,2}=channels_order{ichan};

        end


        clear EEG

        % initiate PAC struct for the current subject
        PAC.(night).length(num_subject).name=SubjectName;
        PAC.(night).direction(num_subject).name=SubjectName;
        PAC.(night).sDirection(num_subject).name=SubjectName;

        PAC.(night).fast.direction(num_subject).name=SubjectName;
        PAC.(night).fast.length(num_subject).name=SubjectName;

        for ichan=1:numel(channels_order) % loop again over each channel to compute PAC
            
            source_EEG=channels_order{ichan}; % current channel
            fprintf(['computing PAC for channel ' source_EEG '\n'])
            fprintf(['processing ' stage '\n'])

            if isempty(conco_SO_SPS.(night).(source_EEG)) % if no concomitant events, pass
                fprintf('No conco SO sps\n')
                PAC.(night).length(num_subject).(source_EEG)=nan;
                PAC.(night).direction(num_subject).(source_EEG)=nan;
                PAC.(night).sDirection(num_subject).(source_EEG)=nan;

                PAC.(night).Nevents(num_subject).(source_EEG)=0;
                continue            
            end

            start_SO=[conco_SO_SPS.(night).(source_EEG).onsetSO]; % get the start of each SO
            SO_duration = [conco_SO_SPS.(night).(source_EEG).SODuration]; % get the duration of each SO

            if all(isnan(start_SO))
                fprintf('No conco SO sps\n')
                PAC.(night).length(num_subject).(source_EEG)=nan;
                PAC.(night).direction(num_subject).(source_EEG)=nan;
                PAC.(night).sDirection(num_subject).(source_EEG)=nan;

                PAC.(night).Nevents(num_subject).(source_EEG)=0;
                continue
            end

            start_epoch=(start_SO+time_window(1))*sampling_rate; % start and end of the extracted time window around each SO onset
            end_epoch=(start_SO+time_window(2))*sampling_rate;
            n_epoch=numel(start_epoch);
   
            % preallocale var
            EEG_epoch_f1=zeros(n_epoch,round(end_epoch(1))-round(start_epoch(1))+1);
            EEG_epoch_f2=zeros(n_epoch,round(end_epoch(1))-round(start_epoch(1))+1);

            % extract time window around SO onsets
            for i_epoch=1:n_epoch
                EEG_epoch_f1(i_epoch,:)=EEG_channel_f1{ichan}(round(start_epoch(i_epoch)):round(end_epoch(i_epoch)));
                EEG_epoch_f2(i_epoch,:)=EEG_channel_f2{ichan}(round(start_epoch(i_epoch)):round(end_epoch(i_epoch)));
            end


            % compute PAC
            [PAC.(night).length(num_subject).(source_EEG),...
                PAC.(night).direction(num_subject).(source_EEG),...
                PAC.(night).Nevents(num_subject).(source_EEG),...
                PAC.(night).sDirection(num_subject).(source_EEG),...
                all_phase_SO_SPS.(night).(source_EEG).(['S' num2str(num_subject)])(1,:),...
                all_phase_SO_SPS.(night).(source_EEG).(['S' num2str(num_subject)])(2,:),...
                PAC.(night).fast.length(num_subject).(source_EEG),...
                PAC.(night).fast.direction(num_subject).(source_EEG),...
                PAC.(night).fast.Nevents(num_subject).(source_EEG)]=...
                ISr_compute_PAC_hilbert_conco_SO_SPS(EEG_epoch_f1,EEG_epoch_f2,time_window,sampling_rate,SO_duration,plotresult,SO_check);


        end %end for ichan
    
        toc
    end % end for ifile (==inight)


end % end subject

% save data
save([DataDir 'PAC_conco_SO_SPS' type '_' stage '_' method_SO indexation '.mat'],'PAC')
save([DataDir 'all_phase_SO_SPS' type '_' stage '_' method_SO indexation '.mat'],'all_phase_SO_SPS')

