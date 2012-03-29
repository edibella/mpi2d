function [parValues,delays,smoothed] = model(AIF,TimeDelayWindow,cinemri1,mask)
[sx sy st] = size(cinemri1);
if(~exist('mask'))
    mask = ones(sx,sy);
end
parValues = zeros([sx sy 3]);
delays = zeros([sx sy]);
[r,c] = find(mask);

AIF = AIF - mean(AIF(1:5));
deltaSIcurves = zeros([st sx sy]);
offset = mean(cinemri1(:,:,1:5),3);
disp('extracting deltaSIcurves');
for t=1:st
    deltaSIcurves(t,:,:) = cinemri1(:,:,t) - offset(:,:);
end

bloodCurves = zeros(st);
start = min(TimeDelayWindow);
for t=TimeDelayWindow
    bloodCurves((t-start+1):st,t) = AIF(1:(st-(t-start)));
    integral(:,t) = cumsum(bloodCurves(t,:));
end
smoothed = cinemri1;
for x=1:sx
    for y=1:sy
        if(mask(x,y) ~= 1) continue; end;
        tissueCurve = deltaSIcurves(:,x,y);
        %findTimeDelay
        bestDelay = 0;
        lowestError = 9999999;
        A(:,2) = -cumsum(tissueCurve);
        for t=TimeDelayWindow
            %find the parameters
            A(:,1) = integral(:,t);
            A(:,3) = bloodCurves(:,t);
            par = A\tissueCurve;
            fit = A*par;
            err = sum((fit-tissueCurve).^2);
            if(err < lowestError)
                lowestError = err;
                bestDelay = t;
            end
        end
        if(bestDelay == 0)
            bestDelay = 1;
            disp('Tried to model a curve, but couldn''t find a good delay');
        end

        %find the parameters
        A(:,1) = integral(:,bestDelay);
        A(:,3) = bloodCurves(:,bestDelay);
        par = A\tissueCurve;
        parValues(x,y,:) = par;
        delays(x,y) = bestDelay;
        %create movie based on those parameters
        %cinemri should be overwritten in this phase
        smoothed(x,y,:) = A*par+offset(x,y);
    end
end