function conv_result_fft=ISr_morlet_transform(EEG,freqs,cycle,time,plot_wavelet_fft,sampling_rate)

% compute the convolution of morlet wavelet on the EEG data


if nargin<5
    plot_wavelet_fft=0;
end

n=cycle; % variable cycles
s=n./(2*pi*freqs)';

A=1./((s*sqrt(pi)).^0.5);
gaussian=A.*exp(-time.^2./(2*s.^2)); % gaussian func

complex_sin_wave=zeros(numel(freqs),numel(time)); 
cmw=zeros(numel(freqs),numel(time)); % complex morlet wavelet

idx_f=1;
for f=freqs
    complex_sin_wave(idx_f,:)=exp(2*pi*1i*f.*time);
    cmw(idx_f,:)=gaussian(idx_f,:).*complex_sin_wave(idx_f,:);
    idx_f=idx_f+1;
end

%% compute conv

data=squeeze(EEG(1,:,:));
n_trial=size(EEG,3);

n_wavelet=size(cmw,2);
n_data=size(EEG,2);
n_convolution=n_wavelet+n_data-1;
half_of_wavelet_size = (n_wavelet-1)/2;

conv_result_fft=zeros(numel(freqs),n_data,n_trial);

idx_f=1;

for f=freqs
    fft_wavelet=fft(cmw(idx_f,:),n_convolution);
    fft_data=fft(data,n_convolution);
    
    fft_wavelet=repmat(fft_wavelet,[1 n_trial]);

    ifft_data = ifft(fft_wavelet.*fft_data,n_convolution);

    conv_result_fft(idx_f,:,:) = ifft_data(1,half_of_wavelet_size+1:end-half_of_wavelet_size)*1/sampling_rate;
    idx_f=idx_f+1;

end


if plot_wavelet_fft % if you want to check the wavelet spec


        figure
        subplot(2,1,1)
        plot3(time,real(cmw(1,:)),imag(cmw(1,:)))
        ylim(max(abs(ylim)).*[-1 1])
        view(0,90)
        xlabel('Time (s)')
        ylabel('Real part')
        zlabel('Imag part')
        title(['Freq : ' num2str(freqs(1)) ' Hz'])

        subplot(2,1,2)
        plot3(time,real(cmw(end,:)),imag(cmw(end,:)))
        ylim(max(abs(ylim)).*[-1 1])
        view(0,90)
        xlabel('Time (s)')
        ylabel('Real part')
        zlabel('Imag part')
        title(['Freq : ' num2str(freqs(end)) ' Hz'])
   


end
