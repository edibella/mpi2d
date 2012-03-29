function [fatnesses,b] = findWeakestLink(shifts)

DC = mean(shifts,1);
for t=1:length(shifts)
    shifts(t,:) = shifts(t,:) - DC;
end
[b, stats] = robustfit(shifts(:,1), shifts(:,2));
myx = min(shifts(:,1)):max(shifts(:,1));
theta = -atan2(b(2),1);
for t=1:length(shifts)
    rotated(t,:) = shifts(t,:)*[cos(theta) sin(theta);-sin(theta) cos(theta)];
end
[bnew, stats] = robustfit(rotated(:,1), rotated(:,2));
theta = -atan2(bnew(2),1);
for t=1:length(shifts)
    rotated(t,:) = rotated(t,:)*[cos(theta) sin(theta);-sin(theta) cos(theta)];
end
[bnew, stats] = robustfit(rotated(:,1), rotated(:,2));
fatnesses = rotated(:,2);
