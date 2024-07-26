function Data = ISr_SO_ReplaceSampleByNan(Data, hypno_night)

% Replace non NREM23 sample by NaN value to later compute only on these
% sample


NREM23_samples = zeros(1,numel(Data));


NREM2_samples = floor([hypno_night{3}.NREM2.onset*250 ; (hypno_night{3}.NREM2.onset + hypno_night{3}.NREM2.duration) *250]); % get the pnts during NREM2 sleep stage
NREM3_samples = floor([hypno_night{3}.NREM3.onset*250 ; (hypno_night{3}.NREM3.onset + hypno_night{3}.NREM3.duration) *250]);% get the pnts during NREM3 sleep stage
bad_samples = floor([hypno_night{3}.BadIntervals.onset*250 ; (hypno_night{3}.BadIntervals.onset + hypno_night{3}.BadIntervals.duration) *250]); % get the pnts during bad intervals
bad_samples(bad_samples==0)=1;  % if the bad intervals begin at 0 sec
NREM2_samples(NREM2_samples==0)=1;  % if the first N2 begin at 0 sec (recording issue, typicaly if the recording start when the subject already sleep)


for NREM2_wd=1:size(NREM2_samples,2) % loop over each N2 periods
    
    NREM23_samples(NREM2_samples(1,NREM2_wd):NREM2_samples(2,NREM2_wd))=1; % keep the corresponding pnts
    
end

for NREM3_wd=1:size(NREM3_samples,2)  % loop over each N3 periods
    
    NREM23_samples(NREM3_samples(1,NREM3_wd):NREM3_samples(2,NREM3_wd))=1; % keep the corresponding pnts
    
end

for bad_wd=1:size(bad_samples,2) % loop over each bad intervals
   
        
    NREM23_samples(bad_samples(1,bad_wd):bad_samples(2,bad_wd))=0;%  remove the corresponding pnts
    
end


if sum(NREM23_samples)~=0
    Data(NREM23_samples==0) = nan; % replace 0 (removed pnts) with nan, easier to work with after.
end
