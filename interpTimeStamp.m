% This routine is to interpolate known deltaSIcurves vs. timeStamp values
% into equally spaced points on the AIF and tissue curves

Slice=3; % first choose what slice you're working with
origSampleTimes=timeStamp(2:end,Slice);
delta_t=1; % this can be changed depending on the time between heart beats
nT=round(max(origSampleTimes./delta_t));
bldcurve=interp1(origSampleTimes,deltaSIcurves(1,:),delta_t*(0:nT-1)','linear',0); % can be nearest, linear, spline, or cubic 
bldcurve=bldcurve';

for i=1:8 % or # of tissue curves
tisscurve(i,:)=interp1(origSampleTimes,deltaSIcurves(i+1,:),delta_t*(0:nT-1)','linear',0);
end

deltaSIcurves1=[bldcurve;tisscurve];
 % then save deltaSIcurves deltaSIcurves
 % and replace a specified deltaSIcurves file with this new one