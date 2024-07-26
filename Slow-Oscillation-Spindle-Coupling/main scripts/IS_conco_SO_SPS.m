%% identify concomitant SO and SPS detected

clearvars
close all
clc


%% some vars
subject_to_process=2:31;

channels={'Fp1','Fp2',...
    'F7','F3','Fz','F4','F8',...
    'FT9','FC5','FC1','FC2','FC6','FT10',...
    'T7','C3','Cz','C4','T8',...
    'CP5','CP1','CP2','CP6',...
    'P7','P3','Pz','P4','P8',...
    'O1','Oz','O2'};

method='Sta';
baseline_or_conso={'baseline','conso'};
type=''; % '' for all, 'Grp' for grouped spindles, 'Iso' for isolated ones
stage='NREM2';  % sleep stage to process
time_specific=0; % put 20 if you want to compute the concomitant for the first 20 minutes of sleep for example (necessitate to extract the spindles before)

switch type
    case ''
        type_sps_str='spsNoBi';
    case 'Grp'
        type_sps_str='spsGrpNoBi';
    case 'Iso'
        type_sps_str='spsIsoNoBi';
end



if time_specific==0 % if we take spindles over the whole night
    sps_str='byStage';
    sps_str_dir='sps_by_stage';
    indexation='';
else  % if we put time constraints
    sps_str=['time_' num2str(time_specific)];
    sps_str_dir='sps_time_specific';
    indexation=['_time_' num2str(time_specific)];
end



%% Directories


MainDir='\SO detection\'; % main path

SODir=[MainDir '\onset_duration_extracted\'];  % SO data extracted 
SPSDir= [MainDir 'IS_spindle_processing\' sps_str_dir '\']; % extracted and clustered spindle data directory
OutputDataDir=[MainDir '\data\']; % output dir
OutputDir=[OutputDataDir 'conco_SO_SPS\'];


if ~exist(OutputDir,'dir')
   mkdir(OutputDir)
end

for num_subject=subject_to_process
    subject=num2str(num_subject);

    if num_subject<10
        subject_str=['00' subject];
    else
        subject_str=['0' subject];
    end

    SubjectName=['somhtd_' subject_str];
    fprintf(['\nprocessing subject' subject_str '\n'])

    % check if the structure already exist, and load if it is the case (to append another night for
    % example)
    try 
        load([OutputDir SubjectName '_conco_SO_SPS_' stage '_' type '_' method indexation])
    catch
        conco_SO_SPS=struct;

    end

    for inight=1:2 % loop over the nights
        night_type=baseline_or_conso{inight};


        for chan=1:numel(channels) % loop for each channel
    
            source_EEG=channels{chan};

            % load sps and SO
            sps_file=[SPSDir SubjectName '\somhtd_' subject_str '_' night_type '_spsGrpIso_' sps_str '_' source_EEG '.mat'];
            SO_file=[SODir SubjectName '\somhtd_' subject_str '_' night_type '_SO_' source_EEG '_NREM23_' method '.mat'];

            if ~exist(sps_file,'file') || ~exist(SO_file,'file') % if one of these file doesn't exist, continue
                fprintf(['SO or SPS file not found for channel ' source_EEG ' or for ' night_type ' night\n'])
                conco_SO_SPS.(night_type).(source_EEG)=[];

                continue
            else
                fprintf(['processing channel ' source_EEG ' for ' night_type ' night\n'])
            end
            
            load(sps_file)
            load(SO_file)

            if time_specific~=0
                spsByStage=spsTimeSpecific;
            end

            if strcmp(type,'Grp')  % be careful with grouped spindles, onsets and onset field are not the same         
                 sps_onset=[spsByStage.(stage).(type_sps_str).onsets];
            else
                sps_onset=[spsByStage.(stage).(type_sps_str).onset];
            end

            if isempty(SO_extracted)
                conco_SO_SPS.(night_type).(source_EEG)=[];
                continue

            end


            SO_intervals=[unique(SO_extracted(:,1)');unique(SO_extracted(:,1)' + SO_extracted(:,2)')]; % start and end of the SOs
            % process
            idx=0;
            for isps=1:numel(sps_onset) % loop over each spindles
                conc_SO=find((sps_onset(isps) > SO_intervals(1,:)) & (sps_onset(isps)< SO_intervals(2,:))); % check if the onset is within the interval of an SO
                if ~isempty(conc_SO) % if we found one, store the data
                    idx=idx+1;
                    conco_SO_SPS.(night_type).(source_EEG)(idx).iSO=conc_SO;
                    conco_SO_SPS.(night_type).(source_EEG)(idx).onsetSO=SO_extracted(conc_SO,1);
                    conco_SO_SPS.(night_type).(source_EEG)(idx).isps=isps;
                    conco_SO_SPS.(night_type).(source_EEG)(idx).onsetsps=sps_onset(isps);
                    conco_SO_SPS.(night_type).(source_EEG)(idx).stage=SO_extracted(conc_SO,4);
                    conco_SO_SPS.(night_type).(source_EEG)(idx).SODuration=SO_extracted(conc_SO,2);

                end
                
            end % for isps
        

            if idx==0 % if no concomitant events
                fprintf(['No conco SO SPS for subject' subject_str ' channel ' source_EEG ' night ' night_type ' ' indexation '\n'])
                conco_SO_SPS.(night_type).(source_EEG)(1).iSO=nan;
                conco_SO_SPS.(night_type).(source_EEG)(1).onsetSO=nan;
                conco_SO_SPS.(night_type).(source_EEG)(1).isps=nan;
                conco_SO_SPS.(night_type).(source_EEG)(1).onsetsps=nan;
                conco_SO_SPS.(night_type).(source_EEG)(1).stage=nan;
                conco_SO_SPS.(night_type).(source_EEG)(1).SODuration=nan;

            
            end


        end % for ichan
        
    end % for inight



   % save data
   save([OutputDir SubjectName '_conco_SO_SPS_' stage '_' type '_' method indexation],'conco_SO_SPS')
   clear conco_SO_SPS SO_extracted spsByStage spsTimeSpecific
end


