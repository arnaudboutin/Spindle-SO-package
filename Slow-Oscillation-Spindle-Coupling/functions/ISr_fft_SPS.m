function [max_freq,P1,f] = ISr_fft_SPS(EEG,Fs,L)
%
%

% EEG : data
% Fs : Sampling frequency                    
% L : Length of signal
f = Fs*(0:(L/2))/L;
T = 1/Fs;             % Sampling period       


Y = fft(EEG);

P2 = abs(Y/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);

[~,idx_max_freq]=max(P1);
max_freq=f(idx_max_freq);


% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of S(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')


end

