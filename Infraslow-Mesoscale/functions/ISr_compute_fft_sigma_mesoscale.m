function [spectrum_mesoscale,spectrum_mesoscale_weighted]=ISr_compute_fft_sigma_mesoscale(Data, hypno_night,stage,method,freq_ms,time_specific)

% Data : EEG data
% hypno_night : scoring infos
% stage : sleep stage to process
% method : hilbert or morlet (hilbert is WIP)
% freq_ms : mesoscale frequencies
% time_specific: if any constraint on the temporal scale

%OUTPUT
%  spectrum_mesoscale : basic spectrum of the infraslow rhythm
% spectrum_mesoscale_weighted : weighted by the bouts duration

sps_band=[11 16];  % sps frequency band
sampling_rate=250;

n_freqs=numel(freq_ms);

time_step=0.5*sampling_rate;% 0.5 sec here

max_sleep_duration=time_specific; %minutes

epoch_threshold=0.5; % the spectrum will be computed on continous epoch of 30 secondes minimum
epoch_threshold_sec=epoch_threshold*60;


% preallocate

epoch_data=cell(1);
nchan=size(Data,1);

spectrum_mesoscale=zeros(numel(freq_ms),nchan);
spectrum_mesoscale_weighted=zeros(numel(freq_ms),nchan);


for ichan=1:nchan
    %% divise the signal into continuous N23 epoch


%     end_window=210*60*sampling_rate; % first  3.5 hours ?
    if strcmp(method,'hilbert') % if you want to do it with hilber transform
        Data_bp=bandpass(Data(ichan,:),sps_band,sampling_rate,'ImpulseResponse','fir');
        power_data=hilbert(Data_bp);
    elseif strcmp(method,'morlet')
        % compute power time course
        conv_result_fft=ISr_morlet_transform(Data(ichan,:),sps_band(1):0.2:sps_band(2),4,-6:1/250:6,0,250);
        power_data=mean(abs(conv_result_fft).^2);  % take the power, if you want to use the magnitude instead of power, juste remove the .^2
    else
        beep
        fprintf('Unrecognised method\n')
        pause
    end



    % remove all non - NREM2 smaples
    Data_stage=ISr_Infraslow_ReplaceSampleByNan(power_data,hypno_night,stage); 


    start_sleep=find(~isnan(Data_stage),1,'first');
    Data_window=Data_stage(start_sleep:end);


    %% gap identification

    % identification of all NREM2 bouts
    gap=diff(isnan(Data_window));
    start_gap=find(gap==1);
    end_gap=find(gap==-1);

    if numel(end_gap)<numel(start_gap) % if last gap has no end (ie if sleep doesn't end with N2)
        start_gap(end)=[];
    end


    start_end_gap=[start_gap;end_gap];
    onset_epoch=nan(1,size(start_end_gap,2));

    for igap=1:size(start_end_gap,2)+1 % loop over each gap to extract the data (between the gaps)
        
        
        if igap==1 % for the first gap, take the start of data window, until the start of the first gap
            epoch_data{igap}=Data_window(1:start_end_gap(1,1));
            onset_epoch(igap)=1;
        elseif igap==size(start_end_gap,2)+1 % for the last epoch
            epoch_data{igap}=Data_window(start_end_gap(2,igap-1)+1:end);
            onset_epoch(igap)=start_end_gap(2,igap-1)+1;

        else
            epoch_data{igap}=Data_window(start_end_gap(2,igap-1)+1:start_end_gap(1,igap));            
            onset_epoch(igap)=start_end_gap(2,igap-1)+1;

        end


    end

    %% keep only the first x minutes of sleep
    % if you want to put time constraint

    cumul_sleep_duration=cumsum(cellfun(@numel,epoch_data));
    if max_sleep_duration~=0
        thres_sample=max_sleep_duration*60*sampling_rate;
        last_epoch=find(onset_epoch<=thres_sample,1,'last');

    else
        last_epoch=numel(cumul_sleep_duration);
        thres_sample=cumul_sleep_duration(end);
    end

    if isempty(last_epoch)
        last_epoch=1;
    end


    if last_epoch>numel(cumul_sleep_duration) % if not enough sleep for the time constraint
        fprintf('WARNING !! Not enough sleep\n')
        last_epoch=numel(cumul_sleep_duration);
    else
        surplus_last_epoch=cumul_sleep_duration(last_epoch)-thres_sample;
        epoch_data(last_epoch+1:end)=[];
        if surplus_last_epoch>0
            epoch_data{last_epoch}=epoch_data{last_epoch}(1:end-surplus_last_epoch);
        end
    end

    if any(isnan(epoch_data{last_epoch}))
        last_sample=find(isnan(epoch_data{last_epoch}),1,'first')-1;
        epoch_data{last_epoch}= epoch_data{last_epoch}(1:last_sample);
    end

    %% remove too short epoch
    
    size_epoch=cell2mat(cellfun(@numel,epoch_data,'UniformOutput',false));
    epoch_data(~(size_epoch>epoch_threshold*60*sampling_rate))=[];

    n_epoch=numel(epoch_data);


    if isempty(epoch_data)
        spectrum_mesoscale(:,ichan)=nan(n_freqs,1);
        spectrum_mesoscale_weighted(:,ichan)=nan(n_freqs,1);

        fprintf(['mesoscale could not be computed, not large enough ' stage ' epoch\n'])

    else
    
        sleepduration_computation=cumsum(cellfun(@numel,epoch_data));
        fprintf(['mesoscale computed on ' num2str(sleepduration_computation(end)/250/60) ' min\n'])
    
    
        %% now performe fft on all epoch
        spectrum_epoch=[];
        for iepoch=n_epoch:-1:1
    
            tmp_epoch=epoch_data{iepoch};
            power_epoch=tmp_epoch;
            norm_power=(power_epoch-mean(power_epoch))./mean(power_epoch);  % standardize (highlight better the oscillation)
    
            % second time freq decomposition
            conv_result_fft=ISr_morlet_transform(norm_power,freq_ms,4,-epoch_threshold_sec/2:1/250:epoch_threshold_sec/2,0,250); % convolution
    
            % take the mean power for each frequency at each time step (ie, the spectrum density at each time step)
            spectrum_epoch(:,iepoch)=mean(abs(conv_result_fft(:,1:time_step:end)).^2,2);
    
    
        end
        
        
    

        
        bout_duration=(cellfun(@numel,epoch_data)/sampling_rate/60); % take each period duration
        weight_by_duration=bout_duration./(sum(bout_duration)); % compute the ratio
        spectrum_epoch_weighted=spectrum_epoch.*repmat(weight_by_duration,n_freqs,1);  % multiply the spectrum to weight by the duration of each period
    
    
        spectrum_mesoscale(:,ichan)=mean(spectrum_epoch,2); % average 
        spectrum_mesoscale_weighted(:,ichan)=sum(spectrum_epoch_weighted,2);% sum because the spectrum have been weighted by a ratio (sum of all ratio = 1)

    end

end








