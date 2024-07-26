function SO_core = ISr_SO_FindSORef_Sta(Data, Info, hypno_night)
% Adrien Conessa (CIAMS, Université Paris-Saclay)
% Arnaud Boutin (CIAMS, Université Paris-Saclay)


EEG=Data.SSRef;

n_channels=size(EEG,1); % get the numbes of channels

freq_bp=Info.Parameters.Bandpass_freq; % SO frequency band
min_SO_duration= Info.Parameters.min_SO_duration; % minimum SO duration
max_SO_duration=Info.Parameters.max_SO_duration; % maximum SO duration
amplitude_prctile=Info.Parameters.amplitude_prctile; % amplitude criteria (>75% of all SO candidates here)
sampling_rate=Info.Recording.sRate; % sampling rate

channels={Info.Electrodes.Name}; % all channels to processed

n_epoch=numel(hypno_night{1}); % number of sleep epochs
hypno_sampling=(0:n_epoch).*hypno_night{2}.*sampling_rate; % convert seconds into pnts

for ichan=1:n_channels % loop over each channel

    if all(isnan(EEG(ichan,:)))
        SO_core.(channels{ichan})=[];
        continue
    end


    EEG_bp=bandpass(EEG(ichan,:),freq_bp,sampling_rate,'ImpulseResponse','fir'); % filter data


    % keep only NREM23 samples
    EEG_bp_N23= ISr_SO_ReplaceSampleByNan(EEG_bp,hypno_night); 

    % find SO candidates with zero crossing
    ZC_EEG=diff(sign(EEG_bp_N23));
    
    SO_candidate=find(ZC_EEG==-2); % positive to negative (start and end of SO. In the SO the voltage go from negative to positive)

    % check duration of SO candidate
    duration_SO_candidate=diff(SO_candidate);
    duration_SO_confirmed=duration_SO_candidate>min_SO_duration*sampling_rate & duration_SO_candidate<max_SO_duration*sampling_rate;

    % take the start and end of each SO
    start_SO=SO_candidate(duration_SO_confirmed);
    end_SO=SO_candidate([false duration_SO_confirmed]);

    % collect the peak to peak amplitude
    p2p_amplitude=zeros(1,numel(start_SO));
    for SO=1:numel(start_SO)
        
        data_SO=EEG_bp_N23(start_SO(SO):end_SO(SO));

        p2p_amplitude(SO)=max(data_SO)-min(data_SO);
    end

    p2p_percentile=prctile(p2p_amplitude,amplitude_prctile); % amplitude threshold
    SO_confirmed=p2p_amplitude>=p2p_percentile; % keep only SO above the amplitude threshold

    start_SO_confirmed=start_SO(SO_confirmed); % collect start and end of the SO detected
    end_SO_confirmed=end_SO(SO_confirmed);
    amplitude_SO_confirmed=p2p_amplitude(SO_confirmed);

    n_SO=sum(SO_confirmed);

    % loop through all SOs
    for SO=1:n_SO

        SO_stage=hypno_night{1}(find(hypno_sampling<=start_SO_confirmed(SO),1,'last')); % find the sleep stage of each SO
    

        SO_core.(channels{ichan}).onset(SO)=start_SO_confirmed(SO);
        SO_core.(channels{ichan}).end(SO)=end_SO_confirmed(SO);
        SO_core.(channels{ichan}).duration(SO)=(end_SO_confirmed(SO)-start_SO_confirmed(SO))/sampling_rate; % duration in seconds
        SO_core.(channels{ichan}).amplitude(SO)=amplitude_SO_confirmed(SO);
        SO_core.(channels{ichan}).stage(SO)=SO_stage; % SO stage

    end

    fprintf([num2str(n_SO) ' slow oscillations found in channel ' channels{ichan} '\n'])

end
SO_core.sampling_rate=sampling_rate;

