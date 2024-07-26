function Data = ISr_Infraslow_ReplaceSampleByNan(Data, hypno_night,stage)

% Replace non NREM23 sample by NaN value to later compute  only on these
% sample


Data_N23 = zeros(1,numel(Data));
Data_N2 = zeros(1,numel(Data));
Data_N3 = zeros(1,numel(Data));


NREM2_samples = floor([hypno_night{3}.NREM2.onset*250 ; (hypno_night{3}.NREM2.onset + hypno_night{3}.NREM2.duration) *250]);% get the pnts during NREM2 sleep stage
NREM3_samples = floor([hypno_night{3}.NREM3.onset*250 ; (hypno_night{3}.NREM3.onset + hypno_night{3}.NREM3.duration) *250]);% get the pnts during NREM3 sleep stage
bad_samples = floor([hypno_night{3}.BadIntervals.onset*250 ; (hypno_night{3}.BadIntervals.onset + hypno_night{3}.BadIntervals.duration) *250]);% get the pnts during bad intervals
bad_samples(bad_samples==0)=1;  % if the bad intervals begin at 0 sec
NREM2_samples(NREM2_samples==0)=1;  % if the bad intervals begin at 0 sec


for NREM2_wd=1:size(NREM2_samples,2)
    
    Data_N23(NREM2_samples(1,NREM2_wd):NREM2_samples(2,NREM2_wd)-1)=1;
    Data_N2(NREM2_samples(1,NREM2_wd):NREM2_samples(2,NREM2_wd)-1)=1;

end

for NREM3_wd=1:size(NREM3_samples,2)
    
    Data_N23(NREM3_samples(1,NREM3_wd):NREM3_samples(2,NREM3_wd)-1)=1;
    Data_N3(NREM3_samples(1,NREM3_wd):NREM3_samples(2,NREM3_wd)-1)=1;

end

for bad_wd=1:size(bad_samples,2)
   
        
    Data_N23(bad_samples(1,bad_wd):bad_samples(2,bad_wd)-1)=0;
    Data_N2(bad_samples(1,bad_wd):bad_samples(2,bad_wd)-1)=0;
    Data_N3(bad_samples(1,bad_wd):bad_samples(2,bad_wd)-1)=0;

end


switch stage
    case 'N23'
    
    % if sum(Data_N23)~=0
        Data(Data_N23==0) = nan;
    % end

    case 'N2'
        Data(Data_N2==0) = nan;

    case 'N3'
        Data(Data_N3==0) = nan;
end


