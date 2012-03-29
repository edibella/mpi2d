function  [newTimeCurve,heartRate,delta_t]=InterpTimeCurve2(timeStampFile, origTimeCurve)
% newTimeCurve will have length totaltime/delta_t    (defined below)
% example call:      bldcurve=InterpTimeCurves(timeStampFileAIF, delta_t, bldcurve);


load(timeStampFile);   % assumes gives variable timeStamp
timeStamps=squeeze(timeStamp(:,1)-timeStamp(1,1));   % choosing timeStamp from slice 1 assumedly
d=timeStamps(2:end)-timeStamps(1:end-1);
%%% to keep the same delta_t as used with 2-comp and Fermi modeling, I'll
%%% keep the same delta-t below...
delta_t=mean(d); % could divide by 2 or other...this is arbitrary...maybe I should use the minimum sample spacing -----ie. delta_t=min(d)
totalTime=timeStamp(length(origTimeCurve),1)-timeStamp(1,1);
nTimes=round(totalTime/delta_t);


t=find(origTimeCurve(1,:)==max(origTimeCurve(1,:))); % to find the peak of the blood (or tissue) curve.

try
    %HR=60/mean(d(t-9:t+10)); % to choose a range of t to average the heart rate during the peak bolus of the study
    HR=60/delta_t; % i'll instead try the mean value across the entire scan
catch
end
% units are in seconds after subtract off first time
origSampleTimes=timeStamps;
origSampleTimes=origSampleTimes(2:length(origTimeCurve)+1);


extrapval=0;
newTimeCurve = interp1(origSampleTimes, origTimeCurve, delta_t*(0:nTimes-1),'cubic',extrapval);
newTimeCurve=newTimeCurve';  %%% NOTE: "interp1.m" seems different in matlab6! Using matlab7 here!
try
    heartRate=HR;
catch
    heartRate=60;
end

return