function shifts = TemporalSmoothing(divisor, Nitter,cinemri1, rangex, rangey)
global waith
if((exist('waith','var')>0))
    if(~isempty(ishandle(waith)) && ishandle(waith))
        close(waith);
    end
end
close all
clc
waith = waitbar(0,'upsampling');
saved = cinemri1;
[sx sy st] = size(cinemri1);
shifts = zeros(Nitter,st,2);
for itter = 1:Nitter
    disp(num2str(itter));
    [sx sy st] = size(cinemri1);
    waitbar(0,waith,'smoothing');
    smoothedmovie = zeros(size(cinemri1));
    window = round(st/divisor);   %must be even
    if(mod(window,2) ~= 0)
        window = window + 1;
    end
    padded = cat(3,mean(cinemri1(:,:,1:5),3), cinemri1, mean(cinemri1(:,:,(end-5):end),3));
    padded = padarray(padded,[0 0 (window-1)],'replicate');
    DC = mean(padded,3);
    for t=1:st
        padded(:,:,t) = padded(:,:,t) - DC;
    end
    h = fspecial('gaussian',[1 (window+1)],st/divisor/3);
    h(window/2+1) = 0;
    summ = sum(h);
    h = h/summ; %make it sum to 1
    h = reshape(h,1,1,window+1);
    h = repmat(h,[sx sy]);
    for t=1:st
        waitbar(t/st,waith);
        smoothedmovie(:,:,t) = smoothimg(padded, window, t, h);
        smoothedmovie(:,:,t) = smoothedmovie(:,:,t) + DC;
    end
    
    [sx sy st] = size(smoothedmovie(rangex,rangey,:));
    waitbar(0,waith,'fixing');
    myccs = zeros(1,st);
    for t=1:st
        waitbar(t/st,waith);
        cc = normxcorr2(smoothedmovie(rangex,rangey,t),cinemri1(rangex,rangey,t)); 
        [max_cc, imax] = max(abs(cc(:)));
        myccs(t) = max_cc;
        [xpeak, ypeak] = ind2sub(size(cc),imax(1));
        shiftAmount = -[ (xpeak-sx) (ypeak-sy) ];
        shifts(itter,t,:) = shiftAmount;
        disp(num2str(shifts(itter,t,:)));
        %cinemri1(:,:,t) = circshift(cinemri1(:,:,t),shiftAmount);

    end
    
    %fix frames that are completely wrong
    [fatnesses,b] = findWeakestLink(squeeze(sum(shifts,1)));
    threshold = (Nitter-itter-1)/Nitter*3.5*std(fatnesses);
    onesToFix = (abs(fatnesses) > threshold) ;
   
    for tt=1:length(onesToFix)
        if(norm(squeeze(sum(shifts(:,tt,:),1)))> 20)
            onesToFix(tt) = 1;
        end
    end
    for t=1:length(onesToFix)
        for mydivisor=linspace(4,2*st,4)
            figure(1), clf, hold on
            myshifts = squeeze(sum(shifts,1));
            for tt=1:st
                if(onesToFix(tt) == 0)
                    color = [0 0 1];
                else
                    color = [1 0 0];
                end
                plot(myshifts(tt,1),myshifts(tt,2),'o','Color',color);
            end
            myx = min(myshifts(:,1)):max(myshifts(:,1));
            plot(myx, b(2)*myx + b(1));
            hold off
            
            if(onesToFix(t) == 0) continue; end;
            mywindow = round(st/mydivisor);   %must be even
            if(mod(mywindow,2) ~= 0)
                mywindow = mywindow + 1;
            end
            h = fspecial('gaussian',[1 (mywindow+1)],st/mydivisor/3);
            h(mywindow/2+1) = 0;
            summ = sum(h);
            h = h/summ; %make it sum to 1
            h = reshape(h,1,1,mywindow+1);
            h = repmat(h,[size(smoothedmovie,1) size(smoothedmovie,2)]);
            figure(2),subplot(1,2,1), imagesc(smoothedmovie(:,:,t)), colormap gray
            smoothedmovie(:,:,t) = smoothimg(padded, window, t, h);
            figure(2),subplot(1,2,2), imagesc(smoothedmovie(:,:,t)), colormap gray
            cc = normxcorr2(smoothedmovie(rangex,rangey,t),cinemri1(rangex,rangey,t)); 
            [max_cc, imax] = max(abs(cc(:)));
            myccs(t) = max_cc;
            [xpeak, ypeak] = ind2sub(size(cc),imax(1));
            shiftAmount = -[ (xpeak-sx) (ypeak-sy) ];
            shifts(itter,t,:) = shiftAmount;
            [fatnesses,b] = findWeakestLink(squeeze(sum(shifts,1)));
            threshold = (Nitter-itter-1)/Nitter*3.5*std(fatnesses);
            onesToFix = (fatnesses > threshold);
            for tt=1:length(onesToFix)
                if(norm(squeeze(sum(shifts(:,tt,:),1)))> 20)
                    onesToFix(tt) = 1;
                end
            end
        end
    end
end
shifts = squeeze(sum(shifts,1));
close(waith);

