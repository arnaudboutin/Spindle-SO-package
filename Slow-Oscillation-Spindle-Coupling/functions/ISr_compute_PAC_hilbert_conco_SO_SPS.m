function [PAC, mean_angle,n_event_conco,sAngle,...
    phase_pic_sps,peak_freq_sps,...
    PAC_fast_sps,mean_angle_fast_sps,n_fast_sps_conco]= ISr_compute_PAC_hilbert_conco_SO_SPS(EEG_f1,EEG_f2,time_window,sampling_rate,plotresult,SO_check)
% Adrien Conessa (CIAMS, Université Paris-Saclay)
% Arnaud Boutin (CIAMS, Université Paris-Saclay)

% used by ISr_init_SO_SPS_hilbert_conco_SO_SPS.m
% OUTPUT
% PAC : mean length of the resultat vector
% mean angle : preferred direction
% n_event_conco : number of concomitant events
% circular standard deviation
% phase_pic_sps : SO phase for each spindle peak power
% peak_freq_sps : peak frequency of each spindle
% PAC_fast_sps, etc. : same thing but for fast spindle (peak > 14 Hz)

% INPUT
% EEG_f1 : EEG epochs, bandpass in the SO frequency band
% EEG_f2 : same but in the sps frequency band
% time_window : extracted time window
% samping_rate of the data
% plotresult = 1 if you want to plot the polarplot
% SO_check =1 for some example of SO-spindle coupling.


if nargin<6
    SO_check=0;
end
if nargin<5
    plotresult=0;
end


i_start_SO=(-time_window(1))*sampling_rate; % get the onset of SO (depend only on the time window here)
time_tot=time_window(2)-time_window(1); % total duration of the extracted time window
n_sample=time_tot*sampling_rate +1; % number of pnts
SO_duration_sample=SO_duration*sampling_rate; % duration of each SO, in number of pnts
end_SO=i_start_SO+SO_duration_sample; % end of each detected SO

if size(EEG_f1,1)~=n_sample % the matrix orientation is important
    EEG_f1=EEG_f1';
    EEG_f2=EEG_f2';
end

n_trial=size(EEG_f1,1);
EEG_bp_SO=EEG_f1;
EEG_bp_SPS=EEG_f2;

% uncomment if you want to put an hamming window
% W=hamming(size(EEG_bp_SO,1));
% W=repmat(W,1,size(EEG_bp_SO,2));
% hilbert_SO=hilbert(EEG_bp_SO.*W);
% hilbert_SPS=hilbert(EEG_bp_SPS.*W);


% hilbert transform
hilbert_SO=hilbert(EEG_bp_SO);
hilbert_SPS=hilbert(EEG_bp_SPS);

% extract phase and magnitude (or power)
power_sps=abs(hilbert_SPS);
phase_SO=angle(hilbert_SO);

% get the peak power with a fft, if you want to dissociated according to
% slow or fast spindles for example, optionnal.
[peak_freq_sps]=ISr_fft_SPS(EEG_bp_SPS,sampling_rate,size(EEG_bp_SPS,1));


[~,i_max_p_sps]=max(power_sps(i_start_SO+1:end,:)); %get the peak power after the SO onset
i_max_p_sps=i_max_p_sps+i_start_SO;

ind = sub2ind(size(phase_SO),i_max_p_sps,1:n_trial);
phase_pic_sps=phase_SO(ind); % corresponding phase

phase_pic_sps(i_max_p_sps>end_SO)=nan; % remove if peak sps is beyond the end of the SO

mean_angle=angle(mean(exp(1i*(phase_pic_sps)),'omitnan')); % get the preferred direction
PAC=abs(mean(exp(1i*(phase_pic_sps)),'omitnan')); % get the mean resultat vector
n_event_conco=sum(~isnan(phase_pic_sps)); % get the number of effective concomitant events
[~,sAngle]=circ_std(phase_pic_sps'); % get the circular standard deviation, necessitate circ_stat toolbox


if PAC==1 % if only one concomitant event, remove. The resultatn vector length would not be reliable
    PAC=nan;
    mean_angle=nan;
    n_event_conco=1;
end


% only for fast spindle
mean_angle_fast_sps=angle(mean(exp(1i*(phase_pic_sps(peak_freq_sps>14))),'omitnan'));
PAC_fast_sps=abs(mean(exp(1i*(phase_pic_sps(peak_freq_sps>14))),'omitnan'));
n_fast_sps_conco=sum(~isnan(phase_pic_sps(peak_freq_sps>14)));


%% check if necessary here
if SO_check
    if n_SO<10
        n_SO_check=n_SO;
    else
        n_SO_check=10;
    end

    for i=1:n_SO_check
        SO_to_check=i;
        
        figure
        plot(EEG_f1(:,SO_to_check));
        hold on
        plot(EEG_f2(:,SO_to_check),'color',[1 0 0 0.5]);
        phase_plot=plot(phase_SO(:,SO_to_check)*5);
        yL=ylim;
        xL=xlim;
        start_SO_plot=plot([i_start_SO i_start_SO],[yL(1) yL(2)],'k');
        pic_sps_plot=plot([i_max_p_sps(SO_to_check) i_max_p_sps(SO_to_check)],[yL(1) yL(2)],'r');
        plot(power_sps(:,SO_to_check),'color',[1 0 0])

%         legend([phase_plot,start_SO_plot,pic_sps_plot],{'phase SO','start SO','pic power sps'})
%         
%         text(xL(1)+abs(0.05*xL(1)),0.99*yL(2),['phase pic = ' num2str(phase_pic_sps(SO_to_check))],...
%             'HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold')
    
    end
end
%% polar plot

if plotresult
    figure
    set(gca, 'FontName', 'Calibri')   
    set(gca,'FontSize',11)
    
    phase_sps=phase_pic_sps;
    phase_sps(isnan(phase_sps))=[];

    PAC_plot=polarplot([mean_angle mean_angle],[0 PAC],'Color',[0 0 0],'LineWidth',2);
    hold on
    pp_sps=polarplot([phase_pic_sps' phase_pic_sps'],[0 1],'Color',[0 0 0 0.1],'LineWidth',1); 
%     pp_sps=polarhistogram(phase_sps',6,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',.3,'EdgeAlpha',.5,'Normalization','probability');

    rlim([0 1])
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
%     h=legend([PAC_plot,pp_sps(1)],{'Phase préférée','Phase pic sps'});
%     set(h,'position',[0.7176 0.8687 0.2482 0.0869]);
    hold off


end









end