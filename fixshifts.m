function shifts_out = fixshifts(shifts, timestamps, sigmamultiplierThreshold)
savedShifts = shifts;
fatnesses = findWeakestLink(shifts);
threshold = sigmamultiplierThreshold*std(fatnesses);
shifts(fatnesses > threshold,:) = [];
if(length(shifts) < 1)
    shifts = savedShifts;
else
    timestamps(fatnesses > threshold) = [];
end

[b, stats] = robustfit(shifts(:,1), shifts(:,2));
myx = min(shifts(:,1)):max(shifts(:,1));
theta1 = -atan2(b(2),1);
rotated1 = zeros(length(shifts),2);
for t=1:length(shifts)
    rotated1(t,:) = shifts(t,:)*[cos(theta1) sin(theta1);-sin(theta1) cos(theta1)];
end
[bnew, stats] = robustfit(rotated1(:,1), rotated1(:,2));
theta2 = -atan2(bnew(2),1);
rotated2 = zeros(length(shifts),2);
for t=1:length(shifts)
    rotated2(t,:) = rotated1(t,:)*[cos(theta2) sin(theta2);-sin(theta2) cos(theta2)];
end
[bnew, stats] = robustfit(rotated2(:,1), rotated2(:,2));
tallnesses = rotated2(:,1);
if(any(isnan(tallnesses)))
    x = 0;
end
timestamps = timestamps - timestamps(1);
newtallnesses = interp1(timestamps,tallnesses,0:.1:max(timestamps),'cubic');
myfft = fft(newtallnesses); 
halfLength = round(length(myfft)/2);
myfilter = 1./(1+((1:halfLength)/25).^6);
myfilter = [myfilter fliplr(myfilter(1:(end-1)))];
if(length(myfft)-1 == length(myfilter))
    myfilter = 1./(1+((1:halfLength)/25).^6);
    myfilter = [myfilter fliplr(myfilter)];
end
myfft = myfft.*myfilter;
mynewcurve = real(ifft(myfft));
mynewcurve = interp1(0:.1:max(timestamps), mynewcurve, timestamps);
rotated2(:,1) = mynewcurve';
rotated2(:,2) = 0;
untheta = -theta1-theta2;
mynewderotatedcurve = zeros(length(shifts),2);
shifts_out = savedShifts;
unrotated = zeros(size(rotated2));
for t=1:length(shifts)
    unrotated(t,:) = rotated2(t,:)*[cos(untheta) sin(untheta);-sin(untheta) cos(untheta)];
    newnorms(t) = norm(shifts_out(t,:) - unrotated(t,:));
end
newthreshold = sigmamultiplierThreshold*std(newnorms(~isnan(newnorms)));
for t=1:length(shifts)
    if((newnorms(t) > newthreshold) && ~any(isnan(rotated2(t,:))))
        shifts_out(t,:) = round(unrotated(t,:));
    end
end
