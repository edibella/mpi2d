function [shifts,myparsaves,mydelaysaves,pixelwisevariance,cinemri1] = modelRegistration(UpsampledFullCine, rangex, rangey,show,Nitter)
global waith
if(~exist('waith')),waith = 0;end
if(ishandle(waith))
    waith = waitbar(0,waith,'Good Morning Dave');
else
    waith = waitbar(0,'Good Morning Dave');
end

if(~exist('Nitter'))
    Nitter = 4;
end
if(~exist('show'))
    show = 0;
end
if(show)
    subplot(2,2,2), cla,hold on;
end
% h = fspecial('gaussian',[size(UpsampledFullCine,1) size(UpsampledFullCine,2)],1);
% for t=1:size(UpsampledFullCine,3)
%     UpsampledFullCine(:,:,t) = filter2(h,UpsampledFullCine(:,:,t));
% end
%extract out cinemri and make a smooth version
cinemri1 = UpsampledFullCine(rangex, rangey, :);  %Set the default smooth movie to a smooth movie
DC = mean(cinemri1(:));
cinemri1 = cinemri1 - DC;
[X, Y] = find(cinemri1(:,:,1));
UpsampledFullCine = UpsampledFullCine - DC;
[sx sy st] = size(cinemri1);
smoothed= zeros(size(cinemri1));
h = fspecial('gaussian',[1 7],2);
h = reshape(h,[1 1 7]);
h = repmat(h,[sx sy 1]);
disp('Smoothing movie');

for t=1:st
    %smoothed(:,:,t) = sum(cinemri1(:,:,max(1,t-3):min(st,t+3)).*h(:,:,max(1,4-t):max(1,1)),3);
    smoothed(:,:,t) = mean(cinemri1(:,:,max(1,t-3):min(st,t+3)),3);
end

% if(show)
%     close(figure(1));
%     figure(1);
%     title('Smoothed Movie');
%     for t=1:st
%         imagesc(smoothed(:,:,t)), colormap gray, pause(.1);
%     end
% end

%automatically find the RV and reference frame
disp('Finding LV and RV');
[RV, LV] = FindLVRV(cinemri1,0);

%calculate when the RV onsets
RVcurve = squeeze(mean(mean(cinemri1((RV(1)-4):(RV(1)+4), (RV(2)-4):(RV(2)+4),:),1),2));
LVcurve = squeeze(mean(mean(cinemri1((LV(1)-4):(LV(1)+4), (LV(2)-4):(LV(2)+4),:),1),2));
TimeDelay = [0 0];
maxSlope = [0 0];
for t=1:(st-5)
    b = robustfit(1:5,RVcurve(t:(t+4)));
    if(maxSlope(1) < b(2))
        maxSlope(1) = b(2);
        TimeDelay(1) = t;
    end
    b = robustfit(1:5,LVcurve(t:(t+4)));
    if(maxSlope(2) < b(2))
        maxSlope(2) = b(2);
        TimeDelay(2) = t;
    end
end
if(TimeDelay(1) < TimeDelay(2))
    RVupslope = TimeDelay(1);
    LV = RV;   %the RVLV extraction got them backwards
else
    RVupslope = TimeDelay(2);
    RV = LV;   %the RVLV extraction got them backwards
end
%find the range of delays that are understandable
TimeDelayWindow = max(1,round(.5*min(TimeDelay))):min(round(1.5*max(TimeDelay)),st);

if(show)
    subplot(2,2,2);
    plot(RVcurve);
    hold on;
    plot([RVupslope RVupslope],[0 max(RVcurve)],'k');
    plot([max(TimeDelayWindow) max(TimeDelayWindow)],[0 max(RVcurve)],'k');
    hold off
    legend('RV curve through time','Reference Frame');
    pause(.1);
end
%initalize the shift vector, Nitter is used so it can keep a history and
%sum it up at the end
shifts = zeros(Nitter,st,2);

disp('Modeling and registering to smooth movie')
for itter = 1:(Nitter+1)
    disp(['Pass ' num2str(itter)]);
    %do cross correlation with a smoothed movie(smoothed)
    disp('Xcoor')
    cine = UpsampledFullCine(rangex,rangey,:);
    for t=(st-1):-1:1
        c = normxcorr2(cine(:,:,t),smoothed(:,:,t));
        [max_c, imax] = max(abs(c(:)));
        [xpeak, ypeak] = ind2sub(size(c),imax(1));
        shifts(itter,t,:) = [(xpeak-sx) (ypeak-sy)];
         %shift the whole frame and reextract
        UpsampledFullCine(:,:,t) = circshift(UpsampledFullCine(:,:,t), squeeze(shifts(itter,t,:)));
    end
    cinemri1 = UpsampledFullCine(rangex, rangey, :);
    
    if(itter == Nitter+1)
        break;
    end
    if(show)
        title('New Smoothed Movie');
        for t=1:st
            subplot(2,2,1)
            imagesc(smoothed(:,:,t)), colormap gray;
            title('Smoothed');
            subplot(2,2,3)
            imagesc(cinemri1(:,:,t)), colormap gray;  title('Registered');
            pause(.1);
        end
    end
    
    %extract the deltaSIcurves for the movie
    RVcurve = squeeze(mean(mean(cinemri1((RV(1)-4):(RV(1)+4), (RV(2)-4):(RV(2)+4),:),1),2));
    RVcurve = RVcurve - mean(RVcurve(1:5));
    deltaSIcurves = zeros([st sx sy]);
    offset = mean(cinemri1(:,:,1:5),3);
    disp('extracting deltaSIcurves');
    for t=1:st
        deltaSIcurves(t,:,:) = cinemri1(:,:,t) - offset(:,:);
    end

    
    
    %% find parameters for all curves
    mylength = length(deltaSIcurves);
    time=(0:(st-1))/60;
    %make every shifted AIF possible (only possible shifts are forward in
    %time
    bloodCurves = zeros(st);
    start = min(TimeDelayWindow);
    for t=TimeDelayWindow
        bloodCurves((t-start+1):st,t) = RVcurve(1:(st-(t-start)));
        %integral(t,:) = cumtrapz(time,bloodCurves(t,:));
        integral(:,t) = cumsum(bloodCurves(t,:));
    end
    warning off;
    disp('Modeling and creating smooth movie');
    myparsaves = zeros([sx sy 3]);
    mydelaysaves = zeros([sx sy]);
    pixelwisevariance = zeros([sx sy]);
    for x=1:sx
        if(mod(x,4)==0)
            waitbar(x/sx,waith,['Registration - Pass ' num2str(itter) '/' num2str(Nitter)]);
        end
        for y=1:sy
            tissueCurve = deltaSIcurves(:,x,y);
            %findTimeDelay
            bestDelay = 0;
            lowestError = 9999999;
            %A(:,2) = -cumtrapz(time,tissueCurve);
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

            %find the parameters
            A(:,1) = integral(:,bestDelay);
            A(:,3) = bloodCurves(:,bestDelay);
            par = A\tissueCurve;
            myparsaves(x,y,:) = par;
            mydelaysaves(x,y) = bestDelay;
            %create movie based on those parameters
            %cinemri should be overwritten in this phase
            smoothed(x,y,:) = A*par+offset(x,y);
            pixelwisevariance(x,y) = var(abs(squeeze(smoothed(x,y,:)) - tissueCurve));
        end
    end
    %subplot(2,2,4), imagesc(myparsaves(:,:,3)/mean(mean(myparsaves((RV(1)-4):(RV(1)+4), (RV(2)-4):(RV(2)+4),3),1),2)),colormap(1-hsv), colorbar, title('k-trans');
    figure(10);
    %subplot(2,2,4)
    %imagesc(myparsaves(:,:,3) - myparsaves(:,:,1).*myparsaves(:,:,2),[0 .4]),colormap(jet); h = colorbar, title('k-trans');cbfreeze(h);
    imagesc(myparsaves(:,:,3) - myparsaves(:,:,1).*myparsaves(:,:,2),[0 .4]),colormap(jet); colorbar, title('k-trans');
    warning on;
end

shifts = squeeze(sum(shifts,1));