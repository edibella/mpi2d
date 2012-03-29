function [shifts,pixelwisevariance] = model_newer(upsampled, rangex, rangey, Nitter)
%save away the original so we can apply our shifts to this
originalUpsampled = upsampled;

%pixelwise average of the first 9 frames
InitialTissue = mean(upsampled(:,:,1:9),3);
spatialAverage = mean(InitialTissue(:));
%compensate for coils
for t=1:size(upsampled,3)
    upsampled(:,:,t) = upsampled(:,:,t) ./ InitialTissue * spatialAverage;
end
%sometimes there are pixels that have a 0 InitialTissue value
upsampled(isnan(upsampled)) =spatialAverage; 
%subtract out so that the first 9 frams have 0 mean
InitialTissue = mean(upsampled(:,:,1:9),3);
for t=1:size(upsampled,3)
    upsampled(:,:,t) = upsampled(:,:,t) - InitialTissue;
end
cinemri1 = upsampled(rangex, rangey, :);
[sx sy st] = size(cinemri1);


globalmin = min(cinemri1(:));
globalmax = max(cinemri1(:));

%get at the AIF
RVpoint = [50,30];
debugg = 0;
randomPoints = rand(1000,2);
randomPoints(:,1) = ceil(sx*randomPoints(:,1));
randomPoints(:,2) = ceil(sy*randomPoints(:,2));

%construct a filter so that only a window of 20 pixels is considered
c = normxcorr2(cinemri1(:,:,3),cinemri1(:,:,4));
[csx, csy] = size(c);
diskFilter = zeros(csx,csy);
disk = fspecial('disk',20);
disk = disk /(disk(10,10));
startx = (round(csx/2)-round(size(disk,1)/2));
starty = (round(csy/2)-round(size(disk,2)/2));
diskFilter(startx:(startx+size(disk,1)-1),starty:(starty+size(disk,2)-1)) = disk;


shifts = zeros(Nitter+1,st,2);
gaussianFilter = fspecial('gaussian',[csx,csy],10);
%preliminary registration

%I've found that this does most of the work and that model based
%registration actually makes things worse
myvar = var(cinemri1,0,3);
ccScores = zeros(st,1);
for t=1:st
    c = normxcorr2(cinemri1(:,:,t),myvar);
    [vmax, ~] = max(abs(c(:)));
    ccScores(t) = vmax;
end
[~,ReferenceFrame] = max(ccScores);
for t=1:st
    c = normxcorr2(cinemri1(:,:,t),cinemri1(:,:,ReferenceFrame));
    %subplot(1,2,1), surf(c)
    c = c .* gaussianFilter;
    %subplot(1,2,2), surf(c);pause();
    [~, imax] = max(abs(c(:)));
    [xpeak, ypeak] = ind2sub(size(c),imax(1));
    if(norm([(xpeak-sx) (ypeak-sy)]) < 20)
        shifts(1,t,:) = [(xpeak-sx) (ypeak-sy)];
        %shift the whole frame and reextract
        upsampled(:,:,t) = circshift(upsampled(:,:,t), squeeze(shifts(1,t,:)));
    else
        x = 0;
    end
end
old = squeeze(shifts(1,:,:));
cinemri1 = upsampled(rangex, rangey, :);

%what shifts are considered acceptable for the AIF
TimeDelayWindow = 5:20;
start = min(TimeDelayWindow);
last = max(TimeDelayWindow);
pixelwisevariance = zeros(sx,sy);
h = fspecial('gaussian',[st,1],1);
%figure(1)
for itter = 1:Nitter
    
    tInitialTissue = mean(cinemri1(:,:,1:9),3);
    for t=1:size(cinemri1,3)
        cinemri1(:,:,t) = cinemri1(:,:,t) - tInitialTissue;
    end
    %make deltaSIcurves
    deltaSIcurves = cinemri1;
    
    %get a nice smooth AIF
    blood = deltaSIcurves(RVpoint(1),RVpoint(2),:);
    blood = blood + deltaSIcurves(RVpoint(1)+1,RVpoint(2)+1,:);
    blood = blood + deltaSIcurves(RVpoint(1),RVpoint(2)+1,:);
    blood = blood + deltaSIcurves(RVpoint(1)+1,RVpoint(2),:);
    blood = blood + deltaSIcurves(RVpoint(1),RVpoint(2)-1,:);
    blood = blood + deltaSIcurves(RVpoint(1)-1,RVpoint(2),:);
    blood = blood + deltaSIcurves(RVpoint(1)-1,RVpoint(2)+1,:);
    blood = blood + deltaSIcurves(RVpoint(1)-1,RVpoint(2)-1,:);
    blood = blood + deltaSIcurves(RVpoint(1)+1,RVpoint(2)-1,:);
    blood = squeeze(blood / 9);
    temp = vertcat(blood,repmat(mean(blood((end-5):end)),[5,1]));
    temp = filter2(h,temp);
    blood = squeeze(temp(1:st));
    
    %find foot of RV bolus
    beginningFrame = -1;
    maxSlope = 0;
    for t=1:15
        values = blood(t:(t+5));
        b = robustfit(1:6,values);
        if(b(2) > maxSlope)
            beginningFrame = t;
            maxSlope = b(2);
        end
    end
    beginningFrame = beginningFrame + 3;
    %construct all shifts of AIF
    bloodvariations = zeros(st,last);
    integral = zeros(st,last);
    for t=TimeDelayWindow
        temp = blood(1:(st-t+1));
        bloodvariations((st-length(temp)+1):st,t) = temp;
        integral(:,t) = -cumtrapz(bloodvariations(:,t));
    end
    
    %for each pixel retrieve the signal and find the parameters that best
    %fit it
    smoothMovie = zeros(sx,sy,st);
    for y=1:sy
        %disp(num2str(y));
        for x=1:sx
            tissue = squeeze(deltaSIcurves(x,y,:));
            if(any(isnan(tissue)))
                keyboard
            end
%             for j=1:length(randomPoints)
%                 if(randomPoints(j,1) == x && y == randomPoints(j,2))
%                     stop = 0;
%                     debugg = 1;
%                     figure(2), clf
%                     plot(tissue), hold on;
%                 end
%             end
            A(:,2) = -cumtrapz(tissue);
            err = zeros(last,1) + 9999999;
            for t=TimeDelayWindow
                A(:,1) = integral(:,t);
                A(:,3) = bloodvariations(:,t);
                par = A\tissue;
                fit = A*par;
                %if(debugg)
                %    plot(fit,'Color',rand(3,1));
                %end
                err(t) = sum(abs(fit-tissue));
            end
            [bestError, besti] = min(err);
            
            %now that we have the best delay calculate our optimal fit
            A(:,1) = integral(:,besti);
            A(:,3) = bloodvariations(:,besti);
            %disp(['x: ' num2str(x)]);
            par = A\tissue;
            fit = A*par;
            
            %smooth the fit
%             if(debugg), subplot(2,1,1),cla, plot(fit), hold on, end
%             temp = vertcat(fit,repmat(mean(fit((end-5):end)),[5,1]));
%             temp = filter2(h,temp);
%             smoothedfit = temp(1:st);
%             if(debugg), subplot(2,1,1), plot(smoothedfit,'g');plot(tissue,'r'), end
%             fit = smoothedfit;
            %debugg = 1;
            smoothMovie(x,y,:) = fit;
            if(debugg), figure(1), plot(tissue), hold on, plot(fit,'g'); hold off, end
            pixelwisevariance(x,y) = bestError;
            if(debugg), pause(.1);  debugg = 0;end
        end
    end
%     %smooth the smooth movie
%     h = fspecial('gaussian',[sx,sy],1);
%     for t=1:st
%         temp = smoothMovie(:,:,t);
%         temp = filter2(h,temp);
%         smoothMovie(:,:,t) = temp;
%     end
    
    %register cinemri to smoothMovie
    if(any(isnan(cinemri1(:))))
            keyboard
    end
    %add in the offset
    for t=1:st
        cinemri1(:,:,t) = cinemri1(:,:,t) + InitialTissue(rangex,rangey);
        smoothMovie(:,:,t) = smoothMovie(:,:,t) + InitialTissue(rangex,rangey);
    end
    for t=beginningFrame:st
        
        c = normxcorr2(cinemri1(:,:,t),smoothMovie(:,:,t));
        %figure(2), subplot(2,1,1), imagesc(cinemri1(:,:,t)), colormap gray, colorbar, subplot(2,1,2),imagesc(smoothMovie(:,:,t)), colormap gray, colorbar,pause(.1)
        %subplot(1,2,1), surf(c)
        c = c .* diskFilter;
        %subplot(1,2,2), surf(c), pause();
        [~, imax] = max(abs(c(:)));
        [xpeak, ypeak] = ind2sub(size(c),imax(1));
        if(norm([(xpeak-sx) (ypeak-sy)])<= 20)
            shifts(itter+1,t,:) = [(xpeak-sx) (ypeak-sy)];
        else
            x = 0;
        end
        if(norm(squeeze(sum(shifts(:,t,:),1))) > 20)
            % we've shifted it way too far now.  We need to simply reset it
            % and bring it back home
            upsampled(:,:,t) = originalUpsampled(:,:,t);
            shifts(:,t,:) = 0;
        end
        if(any(isnan(shifts(:))))
            keyboard
        end
        
        %shift the whole frame and reextract
        upsampled(:,:,t) = circshift(upsampled(:,:,t), squeeze(shifts(itter,t,:)));
    end
    cinemri1 = upsampled(rangex, rangey, :);
    
    stop = 0;
end
shifts = squeeze(sum(shifts,1));

disp('------------');
disp(abs(old - shifts));